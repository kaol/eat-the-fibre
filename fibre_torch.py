"""fibre_torch.py: batched, level-synchronous fibre engine for pcensus.py (PyTorch: ROCm/HIP, CUDA or CPU).

Peer of pcensus.extend. Observable contract: for each parent X (columns < c filled) the engine returns as a SET
    children(X, c) = { X[:, c] := L + k : k ∈ K, pcensus.check(·, c+1, c) passes }
where L = lift(X, c). ETF_XCHECK=1 asserts this per parent against pcensus.lift/check (numpy reference).

Per level c, for a batch of parents (all fibres mixed; the column schedule is data-independent):
  lift        Δ-lift L, batched: acc = (P⊗P)Δ(b_c), y = Linv·acc, keep parents with DC·y = acc.
  prefilter   Every check condition is h ↦ [P h = h] (or [P h = 0] for the complement variants), checked only where
              Binv h has no support ≥ c+1 (the mask). Write k = Kmat κ. A condition with h = h0 + H κ affine in κ and
              (Binv H)[c:] = 0 has a κ-independent mask and an affine residual, so it is a linear system A κ = r mod p.
              Admitted: the batch-constant conditions (H = 0), the new-column idempotency (H = Kmat), the HL/HR
              conditions U_i·W, W·U_j (H_t = ±U_i·K_t), each one per parent if it is eligible. The result: batched RREF, then
              the candidates are κ0 + span(nullspace) with p^f points, f = kdim − rank, a superset of the survivors.
  check       The full masked check (pcensus.check's condition set, masks and targets) on the candidates.
Soundness: prefilter equations ⊆ check conditions, encoded exactly ⇒ survivors are unchanged (only work drops).
Arithmetic: residues in [0, p) stored as float32. Every contraction goes through mm(), which asserts
len · (p−1)² < 2^24 (exact in float32). Plan() also runs a device self-test against int64 numpy.

Env: ETF_MEM (MiB per working set, default 1024), ETF_GROUP (fibres per checkpoint group, default 4096),
     ETF_OUT (MiB of children per chunk yielded by level(), default ETF_MEM),
     ETF_HOST_STACK=0 (keep pending batches on the device instead of in host memory),
     ETF_XCHECK=1 (per-parent comparison with pcensus.lift/check; slow, for validation)."""
import os, time, itertools, collections
import numpy as np
import torch

EXACT = 1 << 24
MEM = int(os.environ.get('ETF_MEM', '1024') or 1024) << 20
GROUP = int(os.environ.get('ETF_GROUP', '4096') or 4096)
XCHECK = int(os.environ.get('ETF_XCHECK', '0') or 0)
OUT = (int(os.environ.get('ETF_OUT_BYTES', '0') or 0)                 # children (int8 bytes) per yielded chunk;
       or int(os.environ.get('ETF_OUT', '0') or 0) << 20 or MEM)      # ETF_OUT_BYTES is for the tests
HOST_STACK = int(os.environ.get('ETF_HOST_STACK', '1') or 0)          # pending batches wait in host memory
QUOTIENT = int(os.environ.get('ETF_QUOTIENT', '1') or 0)             # 0: no frontier quotient (fibre_torch.quotient)
F32 = 4


def pick_device(spec=None):
    if spec: return torch.device(spec)
    if torch.cuda.is_available(): return torch.device('cuda')          # ROCm builds expose HIP as 'cuda'
    return torch.device('cpu')


def describe(dev):
    if dev.type == 'cuda':
        be = f"HIP {torch.version.hip}" if getattr(torch.version, 'hip', None) else f"CUDA {torch.version.cuda}"
        return f"{dev} ({torch.cuda.get_device_name(dev)}, {be}, torch {torch.__version__})"
    return f"{dev} (torch {torch.__version__}, {torch.get_num_threads()} threads)"


class Fp:
    """F_p on float32 residues in [0, p): exact as long as every contraction stays below 2^24 (asserted in mm)."""

    def __init__(self, p, dev):
        torch.backends.cuda.matmul.allow_tf32 = False
        torch.set_float32_matmul_precision('highest')
        self.p, self.dev = p, dev
        self.inv = torch.as_tensor([0] + [pow(a, p - 2, p) for a in range(1, p)], dtype=torch.float32, device=dev)

    def md(self, x): return torch.remainder(x, self.p)

    def mm(self, a, b):
        assert a.shape[-1] * (self.p - 1) ** 2 < EXACT, ("contraction too long for exact float32", a.shape)
        return torch.remainder(a @ b, self.p)


class Plan(Fp):
    """The static (Mu) side: everything that depends only on (quiver, p, column) and not on the frontier."""

    def __init__(self, C, dev):
        p, d, kdim = C['p'], C['d'], C['kdim']
        super().__init__(p, dev)
        assert p < 128, "frontier stored as int8"
        self.d, self.kdim, self.C = d, kdim, C
        t = lambda a: torch.as_tensor(np.asarray(a, dtype=np.int64) % p, dtype=torch.float32, device=dev)
        self.Binv, self.Kmat, self.B = t(C['Binv']), t(C['Kmat']), t(C['B'])                 # Kmat: (d, kdim)
        M = np.asarray(C['M']) % p
        self.Mx = t(M.reshape(d, d * d))                                     # x ↦ Lop(x)[y,k] = Σ_x x_x M[x,y,k]
        self.My = t(M.transpose(1, 0, 2).reshape(d, d * d))                  # y ↦ Rop(y)[x,k] = Σ_y y_y M[x,y,k]
        DC = np.asarray(C['DC']) % p; nC = DC.shape[1]
        R, piv = C['F'].rref(np.concatenate([DC, np.eye(DC.shape[0], dtype=np.int64)], axis=1))
        assert piv[:nC] == list(range(nC)), "DC not of full column rank"
        self.Linv, self.DC = t(R[:nC, nC:]), t(DC)                           # Linv·DC = I
        self.Cidx = torch.as_tensor(C['Cidx'], dtype=torch.long, device=dev)
        self.col = {}                                                        # c ↦ Δ(b_c) as (rows, cols, coefs)
        Binv = np.asarray(C['Binv']) % p
        for c in range(kdim, d):
            T = np.asarray(C['DL'][C['Cidx'][c - kdim]]) % p
            ar, br = np.nonzero(T)
            assert not Binv[c:, ar].any() and not Binv[c:, br].any(), "Δ-lift needs an unknown column"
            ix = lambda a: torch.as_tensor(a, dtype=torch.long, device=dev)
            self.col[c] = (ix(ar), ix(br), t(T[ar, br]))
        self.combos = {f: t(list(itertools.product(range(p), repeat=f))).reshape(p ** f, f) for f in range(kdim + 1)}
        self.RopK = self.rop(self.Kmat.T)                                    # (kdim, d_x, d_k): Rop(K_t)
        self.LopK = self.lop(self.Kmat.T)                                    # (kdim, d_y, d_k): Lop(K_t)
        assert d * d * (p - 1) ** 2 < EXACT, "float32 contractions would not be exact"
        self.selftest()
        # frontier quotient (pcensus.symmetry_plan): N = ker(H → H|_K), H = G × T, acting as X ↦ Lh X Rh
        N = C.get('N'); self.sym = N is not None and len(N[0]) > 1 and QUOTIENT
        if self.sym:
            Binv_, B_ = np.asarray(C['Binv']) % p, np.asarray(C['B']) % p
            Rh = np.einsum('ij,hjk,kl->hil', Binv_, np.asarray(N[1]) % p, B_) % p
            self.HL, self.HR = t(N[0]), t(Rh)                                    # (|N|, d, d)
            self.Nc = {c: torch.as_tensor(np.nonzero(~Rh[:, c:, :c].any((1, 2)))[0], dtype=torch.long, device=dev)
                       for c in range(kdim + 1, d)}                              # N_c = {n : Rn[c:, :c] = 0}
        self.packw = max(1, int(62 / np.log2(p)))                            # base-p digits per int64 key word
        self.pw = torch.as_tensor([p ** e for e in range(self.packw - 1, -1, -1)], dtype=torch.long, device=dev)

    def lop(self, x):  # (..., d) → (..., d_y, d_k): left multiplication by x, i.e. y ↦ x·y is y @ lop(x)
        return self.mm(x, self.Mx).reshape(*x.shape[:-1], self.d, self.d)

    def rop(self, y):  # (..., d) → (..., d_x, d_k): right multiplication by y, i.e. x ↦ x·y is x @ rop(y)
        return self.mm(y, self.My).reshape(*y.shape[:-1], self.d, self.d)

    def selftest(self):
        g = np.random.default_rng(1); L = self.d * self.d
        a, b = g.integers(0, self.p, (64, L)), g.integers(0, self.p, (L, 32))
        got = self.mm(torch.as_tensor(a, dtype=torch.float32, device=self.dev),
                      torch.as_tensor(b, dtype=torch.float32, device=self.dev)).cpu().numpy().astype(np.int64)
        assert (got == (a @ b) % self.p).all(), f"device {self.dev}: float32 matmul is not exact (TF32/xf32 enabled?)"

    # ---- budgets (bytes are estimates of the per-item working set; MEM bounds a single batch)
    def parent_batch(self, c):
        d, k, m = self.d, self.kdim, c
        Qc, Qv = m + 2 * m * m, 4 * m + 1
        per = F32 * (6 * Qc * d + 2 * m * d * d + 3 * Qv * d * k + (Qv + 1) * d * (k + 1) * 3 + 4 * d * d)
        return max(1, MEM // per)

    def cand_batch(self, c):
        d, m = self.d, c
        Qc, Qv = m + 2 * m * m, 4 * m + 3
        per = F32 * (2 * Qc * d + 4 * m * d + 8 * d * d + 4 * Qv * d + d * d)
        return max(1, MEM // per)


def lift(P, X, c):
    """X (n,d,d) float residues, columns < c filled. → (L (n,d), ok (n,)): unique Δ-preimage in the complement."""
    ar, br, coef = P.col[c]
    n, d = X.shape[0], P.d
    Pa = P.md(P.mm(X[:, :, :c], P.Binv[:c, ar]) * coef)                     # (n, d, nnz)
    Pb = P.mm(X[:, :, :c], P.Binv[:c, br])
    bv = P.mm(Pa, Pb.transpose(1, 2)).reshape(n, d * d)                    # (P⊗P)Δ(b_c)
    y = P.mm(bv, P.Linv.T)
    ok = (P.mm(y, P.DC.T) == bv).all(1)
    L = torch.zeros(n, d, device=P.dev); L[:, P.Cidx] = y
    return L, ok


class Parents:
    """Per-parent constants at column c, shared by prefilter and check."""

    def __init__(self, P, X, L, c):
        p, d, m = P.p, P.d, c
        self.X, self.L, self.c = X, L, c
        U0 = X[:, :, :m]; U1 = P.md(P.B[:, :m] - U0)
        self.U = torch.stack([U0, U1])                                       # (2, n, d, m): P b_i, (1−P) b_i
        LU = P.lop(self.U.transpose(2, 3))                                   # (2, n, m, d_y, d_k): Lop(U_a,i)
        prod = P.mm(LU.transpose(3, 4), self.U[:, :, None])                  # (2, n, i, d_k, j) = U_a,i·U_a,j
        prod = prod.permute(0, 1, 2, 4, 3).reshape(2, X.shape[0], m * m, d)
        Hc = torch.cat([U0.transpose(1, 2), prod[0], prod[1]], 1)            # (n, Qc, d)
        self.Tc = torch.cat([U0.transpose(1, 2), prod[0], torch.zeros_like(prod[1])], 1)
        self.cc = P.mm(Hc, P.Binv.T)                                         # new-basis coordinates
        act = (self.cc[:, :, c + 1:] == 0).all(2)
        ver = (self.cc[:, :, c:] == 0).all(2) if c > P.kdim else torch.zeros_like(act)
        self.sel = act & ~ver                                                # exactly pcensus.check's `sel`
        self.R0 = P.mm(self.cc[:, :, :c], X[:, :, :c].transpose(1, 2))       # (n, Qc, d): Σ_{j<c} cc_j X[:, j]
        self.bcol = P.B[:, c]


def prefilter(P, S):
    """→ (consistent (n,), κ0 (n,kdim), N (n,kdim,kdim) with rows = nullspace generators, free (n,kdim) bool)."""
    p, d, k, c = P.p, P.d, P.kdim, S.c
    n = S.X.shape[0]; Xc = S.X[:, :, :c]; L = S.L
    # batch-constant conditions: c0[c]·(L + Kκ) = T − R0. If c0[c] ≠ 0 this pins Kκ, otherwise it is a consistency bit.
    cz = S.cc[:, :, c]                                                       # (n, Qc)
    rhs = P.md(S.Tc - S.R0)                                                  # (n, Qc, d)
    zero = S.sel & (cz == 0)
    ok = ~(zero & (rhs != 0).any(2)).any(1)
    piv = S.sel & (cz != 0)
    v = P.md(rhs * P.inv[cz.long()][:, :, None] - L[:, None, :])             # Kκ = v_q for each pinning q
    has = piv.any(1)
    q0 = piv.to(torch.uint8).argmax(1)
    v0 = v[torch.arange(n, device=P.dev), q0]                                # (n, d)
    ok &= ~(piv & (v != v0[:, None, :]).any(2)).any(1)
    rows = [torch.where(has[:, None, None], torch.cat([P.Kmat.expand(n, d, k), v0[:, :, None]], 2),
                        torch.zeros(n, d, k + 1, device=P.dev))]
    # κ-affine variable conditions: (h0, H) with h = h0 + Hκ; target h (tg) or 0.
    hs, Hs, tg = [], [], []
    for a, sgn in ((0, 1), (1, p - 1)):
        Ua = S.U[a].transpose(1, 2)                                          # (n, m, d)
        w0 = L if a == 0 else P.md(S.bcol - L)
        Rw, Lw = P.rop(w0), P.lop(w0)                                        # (n, d, d)
        hs += [P.mm(Ua, Rw), P.mm(Ua, Lw)]                                   # U_i·w0, w0·U_j: (n, m, d)
        HL = torch.einsum('nmx,txk->nmkt', Ua, P.RopK)                       # U_i·K_t
        HR = torch.einsum('nmy,tyk->nmkt', Ua, P.LopK)                       # K_t·U_j
        Hs += [P.md(sgn * P.md(HL)), P.md(sgn * P.md(HR))]
        tg += [a == 0] * (2 * c)
    hs.append(L[:, None]); Hs.append(P.Kmat.expand(n, 1, d, k)); tg.append(True)
    h0 = torch.cat(hs, 1); H = torch.cat(Hs, 1)                              # (n, Qv, d), (n, Qv, d, k)
    tgt = torch.as_tensor(tg, device=P.dev)
    c0 = P.mm(h0, P.Binv.T)                                                  # (n, Qv, d)
    CH = P.mm(P.Binv, H)                                                     # (n, Qv, d, k)
    elig = (c0[:, :, c + 1:] == 0).all(2) & (CH[:, :, c:] == 0).all(3).all(2)
    A = P.md(P.mm(Xc[:, None], CH[:, :, :c]) + P.Kmat * c0[:, :, c, None, None])
    b = P.md(P.mm(c0[:, :, None, :c], Xc[:, None].transpose(2, 3))[:, :, 0] + L[:, None] * c0[:, :, c, None])
    A = P.md(A - H * tgt[None, :, None, None]); b = P.md(b - h0 * tgt[None, :, None])
    Ab = torch.cat([A, P.md(-b)[..., None]], 3) * elig[:, :, None, None]
    rows.append(Ab.reshape(n, -1, k + 1))
    cons, kap0, N, free = solve_affine(P, torch.cat(rows, 1))
    return ok & cons, kap0, N, free


def solve_affine(P, Ab):
    """Ab (n, R, k+1) = [E | r]. Solutions of E κ = r mod p are {κ0 + Σ_{f free} c_f N[f] : c ∈ F_p^free}, or ∅ if not
    consistent. → (consistent (n,), κ0 (n,k), N (n,k,k), free (n,k) bool); rows N[f] for pivot f are meaningless."""
    n, k = Ab.shape[0], Ab.shape[2] - 1
    if Ab.shape[1] == 0: Ab = torch.zeros(n, 1, k + 1, device=P.dev)         # no equations ≡ the zero equation
    A, prow, cons = rref(P, Ab)
    isp = prow >= 0                                                          # pivot columns
    R = A[torch.arange(n, device=P.dev)[:, None], prow.clamp(min=0)]         # (n, k, k+1): pivot row of column j
    kap0 = torch.where(isp, R[:, :, k], torch.zeros_like(R[:, :, k]))
    eye = torch.eye(k, device=P.dev)
    N = torch.where(isp[:, None, :], P.md(-R[:, :, :k].transpose(1, 2)), eye)  # N[f] = e_f − Σ_j R[j,f] e_j
    return cons, kap0, N, ~isp


def rref(P, A, ncols=None):
    """Batched reduced row echelon form mod p of A (n, R, k+1) = [E | r]. → (A, pivot row per column or −1, consistent).
    ncols < k: pivot only in the first ncols columns (the rest are carried along, e.g. [B | I | 0] for B⁻¹)."""
    n, Rn, K1 = A.shape; k = K1 - 1
    ar = torch.arange(n, device=P.dev)
    used = torch.zeros(n, Rn, dtype=torch.bool, device=P.dev)
    prow = torch.full((n, k), -1, dtype=torch.long, device=P.dev)
    for j in range(min(k, k if ncols is None else ncols) if Rn else 0):
        cand = (A[:, :, j] != 0) & ~used
        has = cand.any(1)
        pr_i = cand.to(torch.uint8).argmax(1)
        pv = A[ar, pr_i]                                                     # (n, k+1)
        pr = P.md(pv * P.inv[pv[:, j].long()][:, None])
        An = P.md(A - A[:, :, j:j + 1] * pr[:, None, :])
        An[ar, pr_i] = pr
        A = torch.where(has[:, None, None], An, A)
        used[ar[has], pr_i[has]] = True
        prow[:, j] = torch.where(has, pr_i, torch.full_like(pr_i, -1))
    cons = ~(((A[:, :, :k] == 0).all(2)) & (A[:, :, k] != 0)).any(1)
    return A, prow, cons


def check(P, S, pi, w):
    """pcensus.check on candidates: parent index pi (nc,), new column w (nc, d). → survivor mask (nc,)."""
    p, d, c = P.p, P.d, S.c
    # batch-constant conditions: Σ_{j<c} cc_j X_j + cc_c w = T on `sel`
    R = P.md(S.R0[pi] + S.cc[:, :, c][pi][:, :, None] * w[:, None, :])
    bad = (S.sel[pi] & (R != S.Tc[pi]).any(2)).any(1)
    # conditions involving column c (same list, targets and masks as pcensus.check)
    W = torch.stack([w, P.md(S.bcol - w)])                                   # (2, nc, d)
    U = S.U[:, pi].transpose(2, 3)                                           # (2, nc, m, d)
    Rw, Lw = P.rop(W), P.lop(W)
    HL, HR = P.mm(U, Rw), P.mm(U, Lw)                                        # (2, nc, m, d)
    Hnn = P.mm(W[:, :, None], Lw)[:, :, 0]                                   # (2, nc, d)
    H = torch.cat([HL[0], HR[0], HL[1], HR[1], Hnn[0, :, None], Hnn[1, :, None], w[:, None]], 1)
    fl = torch.as_tensor([1.] * (2 * c) + [0.] * (2 * c) + [1., 0., 1.], device=P.dev)
    cb = P.mm(H, P.Binv.T)                                                   # (nc, Qv, d)
    unk = (cb[:, :, c + 1:] != 0).any(2)
    Xp = S.X[pi][:, :, :c]
    Rb = P.md(P.mm(cb[:, :, :c], Xp.transpose(1, 2)) + cb[:, :, c, None] * w[:, None, :])
    bad |= ((Rb != P.md(H * fl[None, :, None])).any(2) & ~unk).any(1)
    return ~bad


STATS = collections.defaultdict(lambda: np.zeros(4, dtype=np.int64))   # c ↦ [parents, lifted, candidates, children]
QSTATS = collections.defaultdict(lambda: np.zeros(2, dtype=np.int64))  # c ↦ [frontier, kept] of the N-quotient


def pack(P, v):
    """(..., L) residues → (..., ⌈L/w⌉) int64 words of w base-p digits: lexicographic order on rows is preserved."""
    w = P.packw; L = v.shape[-1]
    v = torch.nn.functional.pad(v, (0, (-L) % w)).long().reshape(*v.shape[:-1], -1, w)
    return (v * P.pw).sum(-1)


def lexmin(keys):
    """(n, m, q) int64 → (n,) index along m of the lexicographically smallest row."""
    cand = torch.ones(keys.shape[:2], dtype=torch.bool, device=keys.device)
    for j in range(keys.shape[2]):
        col = torch.where(cand, keys[:, :, j], torch.iinfo(torch.long).max)
        cand &= col == col.min(1, keepdim=True).values
    return cand.to(torch.uint8).argmax(1)


def first_unique(rows):
    """(n, q) int64 → indices of one row per distinct value (the first), ascending."""
    _, inv = torch.unique(rows, dim=0, return_inverse=True)
    first = torch.full((int(inv.max()) + 1,), rows.shape[0], dtype=torch.long, device=rows.device)
    first.scatter_reduce_(0, inv, torch.arange(rows.shape[0], device=rows.device), 'amin')
    return first.sort().values


def orbit_keys(P, X, g, c):
    """X (n, d, d) with columns < c filled, g: indices into N with Rn[c:, :c] = 0 → (n, words) int64: the
    lexicographically smallest image over n ∈ g of (n·X)[:, kdim:c] = Ln (X[:, :c] Rn[:c, kdim:c]), packed. The K columns
    are left out: n fixes PK, so (n·X)[:, :kdim] = X[:, :kdim]."""
    k = P.kdim
    L_, R_ = P.HL[g], P.HR[g][:, :c, k:c]; m = len(g); keys = []
    CB = max(1, MEM // (3 * F32 * m * P.d * (c - k + 1)))
    for s in range(0, X.shape[0], CB):
        img = P.mm(L_[None], P.mm(X[s:s + CB, None, :, :c].to(torch.float32), R_[None]))    # (n, m, d, c − kdim)
        kw = pack(P, img.transpose(2, 3).reshape(*img.shape[:2], -1))                      # (n, m, words)
        keys.append(kw[torch.arange(kw.shape[0], device=P.dev), lexmin(kw)])
    return torch.cat(keys)


def quotient(P, rid, X, c):
    """Frontier at column c → one member per (fibre, N_c-orbit). N = ker(H → H|_K) fixes every PK and maps folds to
    folds; for n ∈ N_c = {n : Rn[c:, :c] = 0}, (n·X)[:, :c] depends on X[:, :c] alone, so completions(n·X) =
    n·completions(X). Of the frontier members in one orbit only the first is kept. Kept members are genuine frontier
    nodes, never images, so the engine's skipping of already-verified conditions stays exact. The orbit key is the
    lexicographically smallest image. Invariant: fibre = N · completions(frontier); expand() restores it."""
    g = P.Nc[c]
    if len(g) <= 1: return rid, X
    keep = first_unique(torch.cat([rid[:, None].long(), orbit_keys(P, X, g, c)], 1))
    QSTATS[c] += [len(rid), len(keep)]
    return rid[keep], X[keep]


def expand(P, rid, X):
    """Complete folds (c = d) → their N-orbits, deduplicated per fibre (inverse of the quotients on the way)."""
    if not P.sym or not len(rid): return rid, X
    m = P.HL.shape[0]; out_r, out_X = [], []
    CB = max(1, MEM // (3 * F32 * m * P.d * P.d))
    for s in range(0, len(rid), CB):
        img = P.mm(P.mm(P.HL[None], X[s:s + CB, None].to(torch.float32)), P.HR[None])   # (n, m, d, d)
        out_r.append(rid[s:s + CB, None].expand(-1, m).reshape(-1)); out_X.append(img.reshape(-1, P.d, P.d).to(torch.int8))
    rid2, X2 = torch.cat(out_r), torch.cat(out_X)
    keep = first_unique(torch.cat([rid2[:, None].long(), pack(P, X2.transpose(1, 2).reshape(len(rid2), -1))], 1))
    return rid2[keep], X2[keep]


def level(P, rid, X, c):
    """One column for a batch of parents: (rid (n,), X (n,d,d) int8) → children (rid', X') at column c+1, yielded in
    chunks of at most about OUT bytes. The children of one parent batch are unbounded (up to p^kdim per parent; run_batch job 14:
    44 M children = 13.4 GiB on 1⇉2⇇3 + 1→3 + 2→4 over 𝔽₅), so they are never concatenated as a whole."""
    Xf = X.to(torch.float32)
    L, ok = lift(P, Xf, c)
    st = STATS[c]; st[0] += len(rid); st[1] += int(ok.sum())
    if XCHECK: _xcheck_lift(P, Xf, L, ok, c)
    keep = ok.nonzero()[:, 0]
    rid, Xf, L = rid[keep], Xf[keep], L[keep]
    out_r, out_X, acc = [], [], 0
    seen_r, seen_X = [], []                                              # only for XCHECK
    def flush():
        nonlocal out_r, out_X, acc
        r2, x2 = torch.cat(out_r), torch.cat(out_X)
        out_r, out_X, acc = [], [], 0
        st[3] += len(r2)
        if XCHECK: seen_r.append(r2); seen_X.append(x2)
        return r2, x2
    if len(rid):
        S = Parents(P, Xf, L, c)
        cons, kap0, N, free = prefilter(P, S)
        fn = free.sum(1)
        CB = P.cand_batch(c)
        for f in range(P.kdim + 1):
            idx = (cons & (fn == f)).nonzero()[:, 0]
            if not len(idx): continue
            fc = torch.argsort((~free[idx]).to(torch.uint8), dim=1, stable=True)[:, :f]   # free columns, ascending
            Nsel = torch.gather(N[idx], 1, fc[:, :, None].expand(-1, -1, P.kdim))         # (n_f, f, k)
            p_pow = P.p ** f; total = len(idx) * p_pow
            for s in range(0, total, CB):
                flat = torch.arange(s, min(total, s + CB), device=P.dev)
                li, ci = flat // p_pow, flat % p_pow
                kap = P.md(kap0[idx[li]] + P.mm(P.combos[f][ci][:, None, :], Nsel[li])[:, 0])
                w = P.md(L[idx[li]] + P.mm(kap, P.Kmat.T))
                sv = check(P, S, idx[li], w); st[2] += len(flat)
                pi = idx[li][sv]
                if not len(pi): continue
                Xn = Xf[pi].clone(); Xn[:, :, c] = w[sv]
                out_r.append(rid[pi]); out_X.append(Xn.to(torch.int8)); acc += Xn.shape[0] * P.d * P.d
                if acc >= OUT: yield flush()
    if out_r: yield flush()
    if XCHECK:
        rid2 = torch.cat(seen_r) if seen_r else rid[:0]
        X2 = torch.cat(seen_X) if seen_X else X[:0]
        _xcheck_end(P, rid, Xf, L, c, rid2, X2)


# ---- validation against the numpy reference (pcensus.lift / pcensus.check over all of K)
def _ref(P):
    import pcensus as pc                                   # under `python3 pcensus.py` this is a second module copy
    if not pc.C: pc.C.update(P.C)
    return pc


def _xcheck_lift(P, Xf, L, ok, c):
    pc = _ref(P)
    for Xn, Ln, o in zip(Xf.cpu().numpy().astype(np.int64), L.cpu().numpy().astype(np.int64), ok.cpu().numpy()):
        Lr = pc.lift(Xn, c)
        assert (Lr is None) == (not o), f"column {c}: lift solvability disagrees"
        assert Lr is None or ((Lr % P.p) == Ln).all(), f"column {c}: lift disagrees"


def _xcheck_end(P, rid, Xf, L, c, rid2, X2):
    pc = _ref(P)
    p = P.p; allK = P.C['allK']
    got = collections.defaultdict(set)
    for r, x in zip(rid2.cpu().numpy(), X2.cpu().numpy()): got[int(r)].add(x.astype(np.int8).tobytes())
    want = collections.defaultdict(set)
    for r, Xn, Ln in zip(rid.cpu().numpy(), Xf.cpu().numpy().astype(np.int64), L.cpu().numpy().astype(np.int64)):
        Xb = np.repeat(Xn[None], len(allK), axis=0); Xb[:, :, c] = (Ln + allK) % p
        for x in Xb[pc.check(Xb, c + 1, c)]: want[int(r)].add(x.astype(np.int8).tobytes())
    assert got == want, f"column {c}: engine children differ from pcensus.check"


def run(C, reps, fibres, save=lambda force=False: None, log=print, device=None, only=None):
    """Fill fibres[r] (list of int16 (d,d) arrays, columns in the new basis) for every torus rep r (r ∈ only, if given)
    not yet in fibres.
    Same checkpoint content as the process-pool engine. Work order: groups of GROUP fibres, each driven to completion
    by a stack of (column, batch) frontiers, so memory is bounded by the batch budgets, not by the fibre sizes."""
    dev = pick_device(device)
    P = Plan(C, dev)
    log(f"  torch engine on {describe(dev)}; mem budget {MEM >> 20} MiB, group {GROUP}")
    d, kdim, p = P.d, P.kdim, P.p
    Kmat = np.asarray(C['Kmat'])
    todo = [r for r in (range(len(reps)) if only is None else only) if r not in fibres]
    t0 = last = time.time(); ndone = 0
    for g0 in range(0, len(todo), GROUP):
        grp = todo[g0:g0 + GROUP]
        X0 = np.zeros((len(grp), d, d), dtype=np.int8)
        for i, r in enumerate(grp): X0[i, :, :kdim] = Kmat @ reps[r][0] % p
        home = torch.device('cpu') if HOST_STACK else dev               # where pending batches wait
        stack = [(kdim, torch.as_tensor(grp, device=home), torch.as_tensor(X0, device=home))]
        out = collections.defaultdict(list)
        while stack:
            c, rid, X = stack.pop()
            rid, X = rid.to(dev), X.to(dev)
            if c == d:                   # expand multiplies by |N|: slice so each image set fits the budget
                ne = max(1, MEM // (P.HL.shape[0] * d * d)) if P.sym else len(rid)
                for s in range(0, len(rid), ne):
                    r_, X_ = expand(P, rid[s:s + ne], X[s:s + ne])
                    for r, x in zip(r_.cpu().numpy(), X_.cpu().numpy()): out[int(r)].append(x.astype(np.int16))
                continue
            nb = P.parent_batch(c)
            if len(rid) > nb:
                stack += [(c, rid[s:s + nb].to(home), X[s:s + nb].to(home)) for s in range(0, len(rid), nb)][::-1]
                continue
            for rid2, X2 in level(P, rid, X, c):
                if P.sym and c + 1 < d and len(rid2): rid2, X2 = quotient(P, rid2, X2, c + 1)
                if len(rid2): stack.append((c + 1, rid2.to(home), X2.to(home)))
            del rid, X
            if time.time() - last > 60:
                log(f"  {ndone}/{len(todo)} fibres computed, group {g0 // GROUP + 1}/{-(-len(todo) // GROUP)}, "
                    f"{len(stack)} stacked batches at column {c}, t={time.time() - t0:.0f}s"); last = time.time()
        for r in grp:                    # orbits split across batches can expand to the same folds: dedupe per fibre
            fibres[r] = list({x.tobytes(): x for x in out.get(r, [])}.values())
        ndone += len(grp); save()
    log(f"  done in {time.time() - t0:.0f}s; per column [parents, lifted, candidates, children] (brute force would "
        f"check lifted·{p}^{kdim}): " + ", ".join(f"{c}: {list(map(int, STATS[c]))}" for c in sorted(STATS))
        + (f"; N-quotient (|N| = {P.HL.shape[0]}) per column [frontier, kept]: "
           + ", ".join(f"{c}: {list(map(int, QSTATS[c]))}" for c in sorted(QSTATS)) if QSTATS else ""))
    return fibres
