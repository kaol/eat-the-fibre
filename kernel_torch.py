"""kernel_torch.py: the kernel stage on the GPU (PyTorch: ROCm/HIP, CUDA or CPU). Peer of pcensus.kernel_folds_ref +
pcensus.torus_reps_ref. It returns the SAME lists in the SAME order, so checkpoints (fibres keyed by rep index) and
representative-dependent output fields don't change.

  subalgebras  Every subspace of K = F_p^kdim has exactly one RREF. The RREFs with k rows and pivot set π form a Schubert
               cell ≅ F_p^f, where f counts the free entries (row i, column j ∉ π, j > π_i). Index t ∈ [0, p^f) unranks to
               free values in itertools.product order (the first free entry is the most significant digit), so cell
               order × t order is the reference order. S is a subalgebra iff every product r_a·r_b of its rows lies in
               span S. With S in RREF, membership is linear: v ∈ span S ⇔ v − v[π]·S = 0. No rank computation is
               needed, and the whole cell is one batched bilinear form.
  pairs        Folds of K ↔ pairs (H, Q) of complementary subalgebras: P = Bt E Bt⁻¹ = Bt[:, :k] (Bt⁻¹)[:k] with
               Bt = [H; Q]ᵀ. Complementary ⇔ Bt is invertible, which the batched Gauss–Jordan reports (fibre_torch.rref).
               Order: k ascending, then H, then Q (row-major), as in the reference double loop.
  torus        The grading torus l ∈ F_p^× acts on folds of K by PK ↦ A_l PK A_l⁻¹. All conjugates go through one
               batched contraction. Each conjugate is mapped to its index in Kf (folds are coded base p, then binary
               search). The orbit label is the minimum index, so the reps (first member of each orbit in Kf order) and
               their l-lists (first l of each distinct conjugate) are those of the reference loop.
Arithmetic: float32 residues through Fp.mm (exactness asserted), as in fibre_torch.

Env: ETF_MEM (MiB per working set, shared with fibre_torch)."""
import itertools
import numpy as np
import torch
from fibre_torch import Fp, MEM, F32, pick_device, rref


class KPlan(Fp):
    """Static side of the kernel stage: K's multiplication in K coordinates, the Schubert cells."""

    def __init__(self, C, dev):
        p, kdim = C['p'], C['kdim']
        super().__init__(p, dev)
        Kmat, M = np.asarray(C['Kmat']) % p, np.asarray(C['M']) % p
        Kinv = np.asarray(C['Binv'])[:kdim] % p                              # old coords (on K) → K coords
        prod = np.einsum('xi,yj,xyc->ijc', Kmat, Kmat, M) % p               # K_i·K_j in old coords
        MK = np.einsum('ijc,kc->ijk', prod, Kinv) % p                       # (kdim, kdim, kdim)
        assert ((np.einsum('ijk,ck->ijc', MK, Kmat) - prod) % p == 0).all(), "K is not closed under multiplication"
        self.n = kdim
        i, j, c = np.nonzero(MK)                                            # MK is sparse (e.g. 10 of 216 for 1⇉2⇉3)
        ix = lambda a: torch.as_tensor(a, dtype=torch.long, device=dev)
        self.I, self.J = ix(i), ix(j)
        self.coef = torch.as_tensor(MK[i, j, c], dtype=torch.float32, device=dev)
        self.onehot = torch.as_tensor(np.eye(kdim)[c], dtype=torch.float32, device=dev)   # (nnz, kdim): term → c
        self.cells = []                                                     # (k, π, free entries (i, j))
        for k in range(kdim + 1):
            for piv in itertools.combinations(range(kdim), k):
                free = [(i, j) for i in range(k) for j in range(kdim) if j not in piv and j > piv[i]]
                self.cells.append((k, piv, free))

    def batch(self, k):
        """subspaces per closure-test batch: S, a gathered copy, nnz terms, product, residual (floats per item, ×3)."""
        n = self.n
        return max(1, MEM // (3 * F32 * (1 + 2 * k * n + len(self.coef) + 3 * n)))


def unrank(KP, k, piv, free, t):
    """t (N,) int64 cell indices → S (N, k, n) float residues, the RREFs in itertools.product order."""
    N, n, p = t.shape[0], KP.n, KP.p
    S = torch.zeros(N, k, n, device=KP.dev)
    for i, c in enumerate(piv): S[:, i, c] = 1
    for f, (i, j) in enumerate(free):
        S[:, i, j] = (t // p ** (len(free) - 1 - f) % p).to(torch.float32)
    return S


def closed(KP, S, piv):
    """S (N, k, n) RREF with pivots piv → (N,) bool: span S is a subalgebra of K. The k² membership conditions combine
    by ∧, so they run one at a time on the survivors of the previous ones: the first condition sees the whole batch,
    the rest see only what is left (most of a cell fails early)."""
    N, k, n = S.shape
    alive = torch.arange(N, device=KP.dev)
    for a, b in itertools.product(range(k), repeat=2):
        Sa = S[alive]                                                       # (N', k, n)
        t = KP.md(Sa[:, a, KP.I] * Sa[:, b, KP.J] * KP.coef)                # (N', nnz): MK's nonzero terms of r_a·r_b
        v = KP.mm(t, KP.onehot)                                             # (N', n) = r_a·r_b
        res = KP.md(v - KP.mm(v[:, None, list(piv)], Sa)[:, 0])             # v − v[π]·S
        alive = alive[(res == 0).all(1)]
        if not len(alive): break
    ok = torch.zeros(N, dtype=torch.bool, device=KP.dev); ok[alive] = True
    return ok


def subalgebras(KP):
    """→ list over k of (m_k, k, n) float tensors: all subalgebras of dimension k, in reference order."""
    out = [[] for _ in range(KP.n + 1)]
    for k, piv, free in KP.cells:
        total, CB = KP.p ** len(free), KP.batch(k)
        for s in range(0, total, CB):
            S = unrank(KP, k, piv, free, torch.arange(s, min(total, s + CB), device=KP.dev))
            out[k].append(S[closed(KP, S, piv)])
    return [torch.cat(o) for o in out]


def inverse(F, A):
    """A (N, n, n) float residues → (A⁻¹ (N, n, n), invertible (N,)) by batched Gauss–Jordan on [A | I]."""
    N, n, _ = A.shape
    I = torch.eye(n, device=F.dev).expand(N, n, n)
    R, prow, _ = rref(F, torch.cat([A, I, torch.zeros(N, n, 1, device=F.dev)], 2), ncols=n)
    ok = (prow[:, :n] >= 0).all(1)
    rows = prow[:, :n].clamp(min=0)                                         # row holding the pivot of column j
    inv = torch.gather(R[:, :, n:2 * n], 1, rows[:, :, None].expand(N, n, n))
    return inv, ok


def folds(KP, subs):
    """Complementary subalgebra pairs → Kf (nK, n, n) int8 numpy, in the reference order (k, H, Q)."""
    n, out = KP.n, []
    per = 3 * F32 * 4 * n * (2 * n + 1)
    for k in range(n + 1):
        H, Q = subs[k], subs[n - k]
        a, b = H.shape[0], Q.shape[0]
        if not a or not b: continue
        hb = max(1, MEM // per // b)                                        # H rows per batch (all Q each)
        for s in range(0, a, hb):
            Hs = H[s:s + hb]; m = Hs.shape[0]
            Bm = torch.cat([Hs[:, None].expand(m, b, k, n), Q[None].expand(m, b, n - k, n)], 2).reshape(m * b, n, n)
            Bt = Bm.transpose(1, 2)
            inv, ok = inverse(KP, Bt)
            P = KP.mm(Bt[ok][:, :, :k], inv[ok][:, :k, :])                  # Bt E Bt⁻¹
            out.append(P.to(torch.int8).cpu().numpy())
    return np.concatenate(out) if out else np.zeros((0, n, n), np.int8)


def torus_reps(KP, C, Kf, phis):
    """Grading-torus orbits of Kf → reps = [(PK int64 (n, n), [l, …])], the reference torus loop's result."""
    p, n, nK = KP.p, KP.n, len(Kf)
    if not nK: return []
    Binv, Kmat = np.asarray(C['Binv']) % p, np.asarray(C['Kmat']) % p
    ls = sorted(phis); L = len(ls)
    A = np.stack([(Binv @ phis[l][0] @ Kmat % p)[:n] for l in ls])          # (L, n, n): A_l, and A_l⁻¹ below
    Ai = np.stack([(Binv @ phis[l][1] @ Kmat % p)[:n] for l in ls])
    # vec(A PK A⁻¹) = vec(PK) · Kr with Kr[(j,k), (l,i,m)] = A_l[i,j] A_l⁻¹[k,m]: all conjugates in one GEMM
    Kr = torch.as_tensor(np.einsum('lij,lkm->jklim', A, Ai).reshape(n * n, L * n * n) % p,
                         dtype=torch.float32, device=KP.dev)
    conj = np.empty((nK, L, n * n), np.int8)
    CB = max(1, MEM // (3 * F32 * (L + 1) * n * n))
    for s in range(0, nK, CB):
        PK = torch.as_tensor(Kf[s:s + CB].reshape(-1, n * n), dtype=torch.float32, device=KP.dev)
        conj[s:s + CB] = KP.mm(PK, Kr).reshape(-1, L, n * n).to(torch.int8).cpu().numpy()
    # index of every conjugate in Kf: one np.unique over the byte rows of Kf and of the conjugates
    row = np.dtype((np.void, n * n))
    keys = np.concatenate([np.ascontiguousarray(Kf.reshape(nK, n * n)).view(row)[:, 0],
                           np.ascontiguousarray(conj.reshape(-1, n * n)).view(row)[:, 0]])
    uniq, inv = np.unique(keys, return_inverse=True)
    at = np.full(len(uniq), -1, np.int64); at[inv[:nK]] = np.arange(nK)
    assert len(uniq) == nK, "duplicate kernel folds, or a torus conjugate of a kernel fold outside Kf"
    idx = at[inv[nK:]].reshape(nK, L)                                       # idx[i, q] = index of A_q PK_i A_q⁻¹
    lab = idx.min(1)                                                        # orbit label = first member in Kf order
    reps = []
    for i in np.nonzero(lab == np.arange(nK))[0]:
        r = idx[i]
        reps.append((Kf[i].astype(np.int64), [ls[q] for q in range(L) if r[q] not in r[:q]]))   # torb.setdefault order
    return reps


def kernel_stage(C, phis, device=None, log=print):
    """→ (Kf (nK, kdim, kdim) int8, reps), equal to (pcensus.kernel_folds_ref, pcensus.torus_reps_ref)."""
    dev = pick_device(device)
    KP = KPlan(C, dev)
    subs = subalgebras(KP)
    log(f"  kernel stage on {dev}: {sum(KP.p ** len(f) for _, _, f in KP.cells)} subspaces, "
        f"subalgebras by dimension {[s.shape[0] for s in subs]}")
    Kf = folds(KP, subs)
    return Kf, torus_reps(KP, C, Kf, phis)
