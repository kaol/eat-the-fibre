"""Parallel, resumable census of Δ-compatible folds of kQ (Q acyclic, any parallel classes) over F_p.

Same mathematics, same output JSON, same checkpoint format as the original pcensus.py. Optimizations (validated
against the reference implementation on K2 and 1⇉2→3 over F_3: identical fold sets, and K2 reproduces N=64):

  1. Field.rref: row elimination vectorized (one numpy op per pivot column instead of a Python loop over rows).
  2. lift: the per-task Gauss solve of DC·x = acc is replaced by a precomputed left inverse of DC (one matmul);
     the Δ-contraction is vectorized over the nonzero entries of Δ(b_col).
  3. check: fully vectorized (all masked conditions computed in a handful of batched einsums/matmuls instead of
     ~200 per-pair numpy calls), and it SKIPS conditions that were already fully determined and passed when the
     parent column was checked (a condition's value never changes once it becomes determined, so this is exact).
  4. extend: a task extends STRIDE columns before returning its survivors as one int16 array (children shared
     columns are re-verified nowhere), cutting task count and IPC by the survival branching factor per stride.

Env: PCENSUS_STRIDE (default 4), PCENSUS_MAXB (default 1024; a task stops extending and returns early once it
holds more than MAXB survivors), PCENSUS_XCHECK=1 to cross-check lift/check against the reference implementations
at every step (slow, for validation runs on small primes).

Engines (--engine): 'pool' = the process pool above (default); 'torch' = fibre_torch.py (GPU via ROCm/CUDA, or CPU):
same fibres as sets, same checkpoint, so a run may be resumed with either engine. Fibres are sorted canonically before
materialization, so the output JSON (including representative-dependent fields) is engine- and schedule-independent.

Usage:
  python3 pcensus.py QUIVER p [--engine pool|torch] [--device DEV] [--jobs J] [--checkpoint FILE] [--out results.jsonl]"""
import sys, os, json, time, pickle, argparse, itertools, collections
import numpy as np
from concurrent.futures import ProcessPoolExecutor, FIRST_COMPLETED, wait

exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'gcensus.py')).read().split("def census(")[0])
QUIVERS.update({
    'K3':        (2, [(0, 1)] * 3),                                     # 1⇛2  (v2.18: 34 orbits)
    'K3_after':  (3, [(0, 1)] * 3 + [(1, 2)]),                          # 1⇛2→3
    'K3_srcout': (3, [(0, 1)] * 3 + [(0, 2)]),                          # 1⇛2, 1→3
    # v2.66: multiplicity variants (vertex-ideal comparison) and the one-arrow control
    'K4': (2, [(0, 1)] * 4),
    'A2': (2, [(0, 1)]),
    'A3': (3, [(0, 1), (1, 2)]),
    'K3_before': (3, [(1, 2)] * 3 + [(0, 1)]),
    'K3_tgtin': (3, [(0, 1)] * 3 + [(2, 1)]),
    'K3_after2': (4, [(0, 1)] * 3 + [(1, 2), (2, 3)]),
    'K3_mid': (4, [(1, 2)] * 3 + [(0, 1), (2, 3)]),
    'K3_tgtfork': (4, [(0, 1)] * 3 + [(1, 2), (1, 3)]),
    'K3_srcboth': (4, [(0, 1)] * 3 + [(0, 2), (3, 1)]),
    # v2.66: simplifications (every parallel class reduced to one arrow)
    'S_srcout': (3, [(0, 1), (0, 2)]),
    'S_in': (3, [(0, 1), (2, 1)]),
    'S_short': (3, [(0, 1), (2, 1), (0, 2)]),
    'A4': (4, [(0, 1), (1, 2), (2, 3)]),
    'S_tgtfork': (4, [(0, 1), (1, 2), (1, 3)]),
    'S_srcboth': (4, [(0, 1), (0, 2), (3, 1)]),
    'S_short_c4': (4, [(0, 1), (2, 1), (0, 2), (2, 3)]),
    'two_pairs': (3, [(0, 1), (0, 1), (1, 2), (1, 2)]),                 # 1⇉2⇉3
    'two_pairs_sep': (4, [(0, 1), (0, 1), (1, 2), (2, 3), (2, 3)]),     # 1⇉2→3⇉4
    'pair_after3': (5, [(0, 1), (0, 1), (1, 2), (2, 3), (3, 4)]),       # 1⇉2→3→4→5
    # v2.27: consecutive parallel classes (Theorem 2) and controls
    'pair_triple':  (3, [(0, 1), (0, 1), (1, 2), (1, 2), (1, 2)]),      # 1⇉2⇛3
    'two_pairs_pre':  (4, [(0, 1), (1, 2), (1, 2), (2, 3), (2, 3)]),    # 1→2⇉3⇉4
    'two_pairs_post': (4, [(0, 1), (0, 1), (1, 2), (1, 2), (2, 3)]),    # 1⇉2⇉3→4
    'two_pairs_short': (3, [(0, 1), (0, 1), (1, 2), (1, 2), (0, 2)]),   # 1⇉2⇉3 with 1→3 (Theorem 2's hypothesis fails)
    'pairs_in':     (3, [(0, 1), (0, 1), (2, 1), (2, 1)]),              # 1⇉2⇇3 (two classes, not consecutive)
    # v2.27: converging / diverging parallel classes (Theorem 3) and controls
    'pairs_out':    (3, [(1, 0), (1, 0), (1, 2), (1, 2)]),              # 1⇇2⇉3 (opposite of pairs_in)
    'pairs_in_tail': (4, [(0, 1), (0, 1), (2, 1), (2, 1), (1, 3)]),     # 1⇉2⇇3, 2→4
    'pairs_in_extra': (4, [(0, 1), (0, 1), (2, 1), (2, 1), (3, 1)]),    # 1⇉2⇇3, 4→2
    'pairs_in_pre': (4, [(3, 0), (0, 1), (0, 1), (2, 1), (2, 1)]),      # 4→1⇉2⇇3
    'pairs_in_side': (4, [(0, 3), (0, 1), (0, 1), (2, 1), (2, 1)]),     # 4←1⇉2⇇3
    'pairs_in_short': (3, [(0, 1), (0, 1), (2, 1), (2, 1), (0, 2)]),    # 1⇉2⇇3, 1→3 (hypothesis fails)
    'triple_pair_in': (3, [(0, 1), (0, 1), (0, 1), (2, 1), (2, 1)]),    # 1⇛2⇇3
    # v2.28: tests of the post hoc reading (converging classes with adjacent outer ends)
    'triple_pair_in_short': (3, [(0, 1), (0, 1), (0, 1), (2, 1), (2, 1), (0, 2)]),   # 1⇛2⇇3, 1→3
    'pairs_in_short_tail': (4, [(0, 1), (0, 1), (2, 1), (2, 1), (0, 2), (1, 3)]),   # 1⇉2⇇3, 1→3, 2→4
    # v2.64: a second leaf extension of pairs_in_short (sink leaf at 3)
    'pairs_in_short_c4': (4, [(0, 1), (0, 1), (2, 1), (2, 1), (0, 2), (2, 3)]),     # 1⇉2⇇3, 1→3, 3→4
    'two_pairs_dense': (4, [(0, 1), (0, 1), (1, 2), (1, 2), (0, 2), (0, 3), (3, 2), (1, 3)]),  # 1⇉2⇉3 + 1→3, 1→4→3, 2→4
})

def _rref_fast(self, A):
    """Same result as the original rref; the row-elimination step is one vectorized numpy op."""
    p = self.p; A = A.copy() % p; r = 0; piv = []
    for c in range(A.shape[1]):
        pr = next((i for i in range(r, A.shape[0]) if A[i, c]), None)
        if pr is None: continue
        if pr != r: A[[r, pr]] = A[[pr, r]]
        A[r] = A[r] * pow(int(A[r, c]), p - 2, p) % p
        nz = np.nonzero(A[:, c])[0]; nz = nz[nz != r]
        if len(nz): A[nz] = (A[nz] - A[nz, c][:, None] * A[r]) % p
        piv.append(c); r += 1
    return A[:r], piv
Field.rref = _rref_fast

C = {}   # worker context
STRIDE = int(os.environ.get('PCENSUS_STRIDE', '4') or 4)
MAXB = int(os.environ.get('PCENSUS_MAXB', '1024') or 1024)
XCHECK = int(os.environ.get('PCENSUS_XCHECK', '0') or 0)

def corders():
    """Column orders from the environment, each a list of basis indices to put first. ETF_CORDERS = o1;o2;... (an order
    portfolio, v2.51: fibres are tried with o1, then o2, ..., see fibre_torch.portfolio) or ETF_CORDER = o1 (one order)."""
    spec = os.environ.get('ETF_CORDERS') or os.environ.get('ETF_CORDER') or ''
    return [[int(x) for x in o.split(',') if x.strip()] for o in spec.split(';') if o.strip()]

def column_order(basis, d, pref=None):
    """Order in which complement basis vectors are tried (and hence the engine's column schedule). Default: basis order,
    or the first order of corders(). pref (basis indices) goes first, in that order, then the rest in basis order. The
    Δ-lift needs every column's Δ-support earlier (fibre_torch.Plan asserts it)."""
    if pref is None: pref = (corders() or [[]])[0]
    return list(pref) + [i for i in range(d) if i not in pref]

def setup(name, p, order=None):
    nv, arrows = QUIVERS[name]; F = Field(p)
    basis, M, DL = build(nv, arrows); d = len(basis); I = np.eye(d, dtype=np.int64)
    Dmat = DL.reshape(d, d * d).T % p
    Kb = F.nullspace(Dmat); kdim = len(Kb)
    Cidx = []; cur = list(Kb)
    for i in column_order(basis, d, order):
        if F.rank(np.array(cur + [I[i]])) > len(cur): cur.append(I[i]); Cidx.append(i)
    B = np.array(cur).T % p; Binv = F.inv(B)
    DC = Dmat[:, Cidx]; assert F.rank(DC) == len(Cidx)
    Kmat = np.array(Kb).T
    allK = np.array(list(itertools.product(range(p), repeat=kdim)), dtype=np.int64) @ np.array(Kb) % p
    return dict(name=name, p=p, F=F, basis=basis, arrows=arrows, M=M, DL=DL, d=d, I=I, Kb=Kb, kdim=kdim, Cidx=Cidx,
                B=B, Binv=Binv, DC=DC, Kmat=Kmat, allK=allK)

def init_worker(name, p):
    C.update(setup(name, p))

def kernel_folds_staged(ctx):
    """(Slower alternative, kept for cross-checks.) Algebra folds of K = ker Δ, enumerated column by column (K coordinates), with the same masked checks as the fibres."""
    p, F, kdim, Kmat, M, d = ctx['p'], ctx['F'], ctx['kdim'], ctx['Kmat'], ctx['M'], ctx['d']
    Kinv = ctx['Binv'][:kdim]                                            # old coords (on K) -> K coords
    MK = np.zeros((kdim, kdim, kdim), dtype=np.int64)
    for i in range(kdim):
        for j in range(kdim):
            prod = np.einsum('i,j,ijk->k', Kmat[:, i], Kmat[:, j], M) % p
            MK[i, j] = Kinv @ prod % p
            assert ((Kmat @ MK[i, j] - prod) % p == 0).all(), "K is not closed under multiplication"
    I = np.eye(kdim, dtype=np.int64)
    allv = np.array(list(itertools.product(range(p), repeat=kdim)), dtype=np.int64)
    def kcheck(Xb, known, new):
        st = {'Xb': Xb, 'idx': np.arange(len(Xb))}
        def Pap(v):
            Xb = st['Xb']; unk = v[..., known:].any(axis=-1)
            cc = v[:, :known] if v.ndim == 2 else np.broadcast_to(v[:known], (Xb.shape[0], known))
            return np.einsum('nij,nj->ni', Xb[:, :, :known], cc) % p, unk
        def keep(ok): st['Xb'] = st['Xb'][ok]; st['idx'] = st['idx'][ok]; return len(st['idx'])
        col = lambda i: st['Xb'][:, :, i] if i == new else st['Xb'][0, :, i]
        for i in range(known):
            v = col(i); r, unk = Pap(v)
            if not keep(~(((r - v) % p).any(axis=1) & ~unk)): return st['idx']
        for i in range(known):
            for j in range(known):
                for comp in (False, True):
                    xi, yj = col(i), col(j)
                    if comp: xi, yj = (I[:, i] - xi) % p, (I[:, j] - yj) % p
                    if xi.ndim == 1 and yj.ndim == 1: h = np.einsum('i,j,ijk->k', xi, yj, MK) % p
                    elif xi.ndim == 1: h = np.einsum('i,nj,ijk->nk', xi, yj, MK) % p
                    elif yj.ndim == 1: h = np.einsum('ni,j,ijk->nk', xi, yj, MK) % p
                    else: h = np.einsum('ni,nj,ijk->nk', xi, yj, MK) % p
                    r, unk = Pap(h)
                    if not keep(~((r if comp else (r - h) % p).any(axis=1) & ~unk)): return st['idx']
        return st['idx']
    layer = [np.zeros((kdim, kdim), dtype=np.int64)]
    for c in range(kdim):
        nxt = []
        for X in layer:
            Xb = np.repeat(X[None], len(allv), axis=0); Xb[:, :, c] = allv
            nxt += list(Xb[kcheck(Xb, c + 1, c)])
        layer = nxt
    return layer

def kernel_folds(ctx):
    """Algebra folds of K = ker Δ as complementary pairs of subalgebras (subspace enumeration)."""
    p, F, kdim, Kmat, M = ctx['p'], ctx['F'], ctx['kdim'], ctx['Kmat'], ctx['M']
    mul = lambda x, y: np.einsum('i,j,ijk->k', x, y, M) % p
    def subspaces(n, k):
        for piv in itertools.combinations(range(n), k):
            free = [(i, j) for i in range(k) for j in range(n) if j not in piv and j > piv[i]]
            for vals in itertools.product(range(p), repeat=len(free)):
                S = np.zeros((k, n), dtype=np.int64)
                for i, c in enumerate(piv): S[i, c] = 1
                for (i, j), v in zip(free, vals): S[i, j] = v
                yield S
    def is_subalg(S):
        vecs = [Kmat @ r % p for r in S]
        return all(F.rank(np.array(vecs + [mul(x, y)])) == len(vecs) for x in vecs for y in vecs) if vecs else True
    subalgs = [S for k in range(kdim + 1) for S in subspaces(kdim, k) if is_subalg(S)]
    out = []
    for H in subalgs:
        for Q in subalgs:
            if H.shape[0] + Q.shape[0] != kdim: continue
            Bm = np.vstack([H, Q]).reshape(kdim, kdim)
            if F.rank(Bm) != kdim: continue
            Bt = Bm.T % p; E = np.diag([1] * H.shape[0] + [0] * Q.shape[0])
            out.append(Bt @ E @ F.inv(Bt) % p)
    return out

def grading(ctx):
    """φ_l = diag(l^len(b)) for l ∈ F_p^×, with inverses: the grading torus."""
    p, basis = ctx['p'], ctx['basis']; lens = [0 if b[0] == 'v' else len(b[1]) for b in basis]
    return {l: (np.diag([pow(l, k, p) for k in lens]), np.diag([pow(pow(l, p - 2, p), k, p) for k in lens])) for l in range(1, p)}

def torus_reps(ctx, Kf, phis):
    """Reference: one rep per grading-torus orbit of Kf (first in Kf order), with the first l of each distinct conjugate."""
    p, kdim, Binv, Kmat = ctx['p'], ctx['kdim'], ctx['Binv'], ctx['Kmat']
    kkey = lambda A: (A % p).astype(np.int8).tobytes(); done_k = set(); reps = []      # (PK, {kkey: l})
    for PK in Kf:
        PK = np.asarray(PK, dtype=np.int64)
        if kkey(PK) in done_k: continue
        torb = {}
        for l, (ph, phi_) in phis.items():
            PKl = (Binv @ ph @ (Kmat @ PK % p) % p)[:kdim] % p
            PKl = PKl @ (Binv @ phi_ @ Kmat % p)[:kdim] % p
            torb.setdefault(kkey(PKl), l)
        done_k |= set(torb); reps.append((PK, list(torb.values())))
    return reps

KERNEL_CACHE_VERSION = 1

def kernel_digest(ctx):
    """Version manifest of the kernel stage's seed: (Kf, reps) is a pure function of these, so a cache keyed by it is sound."""
    import hashlib
    h = hashlib.sha256(f"kernel-v{KERNEL_CACHE_VERSION} p={ctx['p']}".encode())
    for k in ('M', 'DL', 'Kmat', 'Binv'): h.update(np.ascontiguousarray(np.asarray(ctx[k], dtype=np.int64) % ctx['p']).tobytes())
    return h.hexdigest()

def kernel_stage(ctx, phis, engine='auto', device=None, cache=None):
    """(Kf, reps): algebra folds of K and their grading-torus reps. Engines 'ref' (kernel_folds + torus_reps, numpy) and
    'torch' (kernel_torch, GPU) give equal lists in equal order. Memoized in `cache` (npz) under kernel_digest."""
    p, kdim = ctx['p'], ctx['kdim']; dig = kernel_digest(ctx)
    if cache and os.path.exists(cache):
        z = np.load(cache)
        if str(z['digest']) == dig:
            Kf = z['Kf']; reps = [(Kf[i].astype(np.int64), [int(l) for l in ls[:n]]) for i, ls, n in zip(z['ri'], z['rl'], z['rn'])]
            print(f"  kernel stage: cached in {cache}", flush=True)
            return Kf, reps
        print(f"  kernel cache {cache} is for another setup (digest mismatch), recomputing", flush=True)
    if engine == 'auto':
        try: import torch; engine = 'torch'
        except ImportError: engine = 'ref'
    if engine == 'torch':
        import kernel_torch
        Kf, reps = kernel_torch.kernel_stage(ctx, phis, device=device, log=lambda s: print(s, flush=True))
    else:
        Kf = np.array(kernel_folds(ctx), dtype=np.int8).reshape(-1, kdim, kdim); reps = torus_reps(ctx, Kf, phis)
    if cache:
        pos = {(Kf[i] % p).astype(np.int8).tobytes(): i for i in range(len(Kf))}
        ri = np.array([pos[(PK % p).astype(np.int8).tobytes()] for PK, _ in reps], dtype=np.int64)
        rl = np.zeros((len(reps), max(p - 1, 1)), np.int64); rn = np.array([len(ls) for _, ls in reps], dtype=np.int64)
        for i, (_, ls) in enumerate(reps): rl[i, :len(ls)] = ls
        tmp = cache + '.tmp.npz'
        np.savez(tmp, Kf=Kf, ri=ri, rl=rl, rn=rn, digest=np.array(dig)); os.replace(tmp, cache)
    return Kf, reps

def symmetry_plan(ctx, reps, phis):
    """H = G × T (substitutions × grading torus). Every h ∈ H is an algebra automorphism with Δh = (h⊗h)Δ (column sums
    1), so h maps folds to folds and fibre(h·PK) = h·fibre(PK)·h⁻¹; K is H-stable, with h|_K = hK.
      transport  Each H-orbit of torus reps is computed once, at its first rep (the source); the other reps get
                 X' = Lh X Rh with Lh = h, Rh = Binv h⁻¹ B (X = P B in the new basis).
      kernel     N = ker(H → H|_K) fixes every PK, so it acts on every fibre; the torch engine quotients its frontier by N.
    → (sources, {r: (source, Lh, Rh)}, N = (h (m,d,d), h⁻¹ (m,d,d)) for h ∈ N).
    Streaming version (v2.28): H is never stored. Only one h per restriction hK and the elements of N are kept, and the
    rep-to-rep table is recomputed in chunks (two passes), so memory is O(|H_K| + |N| + chunk). Same output as
    symmetry_plan_ref: same sources, and for each transported rep the same source and the same choice of h (the first h,
    in G × T enumeration order, with the needed restriction)."""
    p, d, kdim = ctx['p'], ctx['d'], ctx['kdim']
    B, Binv, Kmat = ctx['B'] % p, ctx['Binv'] % p, ctx['Kmat'] % p; Kinv = Binv[:kdim]
    HKd = {}; firsth = []; Nd = {}; seen_h = set()
    Id8 = np.eye(d, dtype=np.int8).tobytes()
    for g, gi in build_G(ctx):
        for l, (ph, phi_) in phis.items():
            h, hi = ph @ g % p, gi @ phi_ % p
            hb = h.astype(np.int8).tobytes()
            if hb in seen_h: continue                                   # (the reference dedups h the same way)
            seen_h.add(hb)
            hK = Kinv @ h @ Kmat % p
            assert ((h @ Kmat - Kmat @ hK) % p == 0).all(), "K is not H-stable"
            k = hK.astype(np.int8).tobytes()
            if k not in HKd: HKd[k] = len(HKd); firsth.append((h, hi, hK))
            if k == np.eye(kdim, dtype=np.int8).tobytes(): Nd.setdefault(hb, (h, hi))
    seen_h = None
    Nl = list(Nd.values()); N = (np.array([x[0] for x in Nl]), np.array([x[1] for x in Nl]))
    R = len(reps)
    if R == 0: return [], {}, N
    # generators of H_K: greedily add an element outside the subgroup generated so far (closure by BFS on restrictions)
    key = lambda A: (A % p).astype(np.int8).tobytes()
    gens = []; sub = {key(np.eye(kdim, dtype=np.int64))}
    for h, hi, hK in firsth:
        if key(hK) in sub: continue
        gens.append((h, hi, hK)); I_ = np.eye(kdim, dtype=np.int64); frontier = [I_]; sub = {key(I_)}
        while frontier:
            nxt = []
            for A in frontier:
                for _, _, gK in gens:
                    Bk = gK @ A % p; kb = key(Bk)
                    if kb not in sub: sub.add(kb); nxt.append(Bk)
            frontier = nxt
    assert len(sub) == len(firsth), "generators do not generate H_K"
    PK = np.array([PK_ for PK_, _ in reps], dtype=np.int64)
    rk = fold_keys(PK % p); order = np.argsort(rk); srk = rk[order]
    # the image of a torus rep is a kernel fold in some torus orbit; φ_l moves it to that orbit's rep. So each edge is
    # labelled by (generator, l) and carries the element φ_l·g.
    tor = [(ph, phi_, Kinv @ ph @ Kmat % p, Kinv @ phi_ @ Kmat % p) for ph, phi_ in phis.values()]
    def image(gK, gKi):                                                  # (rep index, torus index) of the image of every r
        out = np.full(R, -1, np.int64); tl = np.full(R, -1, np.int64); CH = 1 << 16
        for s0 in range(0, R, CH):
            Cg = (gK[None] @ PK[s0:s0 + CH] @ gKi[None]) % p
            for li, (_, _, tK, tKi) in enumerate(tor):
                ck = fold_keys((tK[None] @ Cg @ tKi[None]) % p)
                jj = np.minimum(np.searchsorted(srk, ck), R - 1)
                hit = (srk[jj] == ck) & (out[s0:s0 + len(ck)] < 0)
                out[s0:s0 + len(ck)][hit] = order[jj][hit]; tl[s0:s0 + len(ck)][hit] = li
        assert (out >= 0).all(), "an image lies in no torus orbit of a rep"
        return out, tl
    imgs_tl = [image(gK, Kinv @ hi @ Kmat % p) for _, hi, gK in gens]
    imgs = [x[0] for x in imgs_tl]
    # orbits = connected components of the generator graph; scanning r in index order makes each source the least
    # index of its orbit (as in symmetry_plan_ref). A BFS tree from the source gives, for every other r, an h ∈ H with
    # h·PK_source = PK_r (a product of generators). The choice of h may differ from the reference; any such h transports
    # the fibre, which is the same set.
    par = np.full(R, -1, np.int64); pgen = np.full(R, -1, np.int64); root = np.full(R, -1, np.int64)
    ptor = np.full(R, -1, np.int64)
    sources = []; bfs = []
    for r0 in range(R):
        if root[r0] >= 0: continue
        sources.append(r0); root[r0] = r0; frontier = np.array([r0], np.int64)
        while len(frontier):
            nxt = []
            for gi_, im in enumerate(imgs):
                tg = im[frontier]; new_ = root[tg] < 0
                tg_new, par_new = tg[new_], frontier[new_]
                tg_new, first = np.unique(tg_new, return_index=True); par_new = par_new[first]
                root[tg_new] = r0; par[tg_new] = par_new; pgen[tg_new] = gi_; ptor[tg_new] = imgs_tl[gi_][1][par_new]
                bfs.append(tg_new); nxt.append(tg_new)
            frontier = np.concatenate(nxt)
    Lh = {}; trans = {}
    for blk in bfs:                                                      # parents precede children in BFS order
        for r in blk:
            r = int(r); g0, gi0, _ = gens[int(pgen[r])]; pr = int(par[r]); ph, phi_ = tor[int(ptor[r])][:2]
            g, gi_ = ph @ g0 % p, gi0 @ phi_ % p                               # edge element φ_l·g
            h0, hi0 = Lh.get(pr, (None, None))
            h, hi = (g, gi_) if h0 is None else (g @ h0.astype(np.int64) % p, hi0.astype(np.int64) @ gi_ % p)
            h8 = h.astype(np.int8); Lh[r] = (h8, hi.astype(np.int8))           # entries < p ≤ 127: int8 is exact
            trans[r] = (int(root[r]), h8, (Binv @ hi % p @ B % p).astype(np.int8))
    return sources, trans, N

def symmetry_plan_ref(ctx, reps, phis):
    """(Reference implementation, kept for cross-checks: stores all of H and an R × |H_K| table.)"""
    """H = G × T (substitutions × grading torus). Every h ∈ H is an algebra automorphism with Δh = (h⊗h)Δ (column sums
    1), so h maps folds to folds and fibre(h·PK) = h·fibre(PK)·h⁻¹; K is H-stable, with h|_K = hK.
      transport  Each H-orbit of torus reps is computed once, at its first rep (the source); the other reps get
                 X' = Lh X Rh with Lh = h, Rh = Binv h⁻¹ B (X = P B in the new basis).
      kernel     N = ker(H → H|_K) fixes every PK, so it acts on every fibre; the torch engine quotients its frontier by N.
    → (sources, {r: (source, Lh, Rh)}, N = (h (m,d,d), h⁻¹ (m,d,d)) for h ∈ N)."""
    p, d, kdim = ctx['p'], ctx['d'], ctx['kdim']
    B, Binv, Kmat = ctx['B'] % p, ctx['Binv'] % p, ctx['Kmat'] % p; Kinv = Binv[:kdim]
    Hd, HKd = {}, {}                                                    # h bytes → (h, h⁻¹, q); hK bytes → q
    for g, gi in build_G(ctx):
        for l, (ph, phi_) in phis.items():
            h, hi = ph @ g % p, gi @ phi_ % p
            hK = Kinv @ h @ Kmat % p
            assert ((h @ Kmat - Kmat @ hK) % p == 0).all(), "K is not H-stable"
            q = HKd.setdefault(hK.astype(np.int8).tobytes(), len(HKd))
            Hd.setdefault(h.astype(np.int8).tobytes(), (h, hi, q))
    Hl = list(Hd.values())
    H, Hi, q = np.array([x[0] for x in Hl]), np.array([x[1] for x in Hl]), np.array([x[2] for x in Hl])
    first = {}                                                          # q → one h with that restriction
    for i, qi in enumerate(q): first.setdefault(int(qi), i)
    rep_h = [first[k] for k in range(len(HKd))]
    HK = np.array([Kinv @ H[i] @ Kmat % p for i in rep_h]); HKi = np.array([Kinv @ Hi[i] @ Kmat % p for i in rep_h])
    qid = q[(H == np.eye(d, dtype=np.int64)).all((1, 2))][0]
    N = (H[q == qid], Hi[q == qid])
    R = len(reps)
    if R == 0: return [], {}, N
    PK = np.array([PK_ for PK_, _ in reps], dtype=np.int64)
    rk = fold_keys(PK % p); order = np.argsort(rk); srk = rk[order]
    Mx = np.empty((R, len(HK)), np.int64)                              # Mx[r, q] = rep index of hK PK_r hK⁻¹, or −1
    CH = max(1, (1 << 22) // (len(HK) * kdim * kdim))
    for s0 in range(0, R, CH):
        conj = (HK[None] @ PK[s0:s0 + CH, None] @ HKi[None]) % p               # (c, nq, kdim, kdim)
        ck = fold_keys(conj.reshape(-1, kdim, kdim))
        j = np.minimum(np.searchsorted(srk, ck), R - 1)
        Mx[s0:s0 + CH] = np.where(srk[j] == ck, order[j], -1).reshape(-1, len(HK))
    src = np.where(Mx >= 0, Mx, R).min(1)                               # first rep of the H-orbit (the orbit is H·PK_r)
    assert (src <= np.arange(R)).all()
    sources = [int(r) for r in np.nonzero(src == np.arange(R))[0]]
    trans = {}
    for r in np.nonzero(src != np.arange(R))[0]:
        s_ = int(src[r]); i = rep_h[int(np.nonzero(Mx[s_] == r)[0][0])]       # h·PK_source = PK_r
        trans[int(r)] = (s_, H[i], Binv @ Hi[i] @ B % p)
    return sources, trans, N

def transport(ctx, reps, fibres, trans):
    """Fill fibres[r] for transported reps whose source is done (see symmetry_plan)."""
    p, kdim, Kmat = ctx['p'], ctx['kdim'], ctx['Kmat'] % ctx['p']
    for r, (s_, Lh, Rh) in trans.items():
        if r in fibres or s_ not in fibres: continue
        Xs = np.array(fibres[s_], dtype=np.int64).reshape(-1, ctx['d'], ctx['d'])
        Xt = (Lh.astype(np.int64) @ Xs @ Rh.astype(np.int64) % p).astype(np.int16)
        assert (Xt[:, :, :kdim] == (Kmat @ reps[r][0] % p)).all(), "transported fibre is not over its kernel fold"
        fibres[r] = list(Xt)

def lift_ref(X, col):
    """Reference lift (per-task Gauss solve), kept for PCENSUS_XCHECK."""
    p, d, I, DL, Binv, DC, Cidx, kdim, F = (C[k] for k in ('p', 'd', 'I', 'DL', 'Binv', 'DC', 'Cidx', 'kdim', 'F'))
    T = DL[Cidx[col - kdim]]; acc = np.zeros((d, d), dtype=np.int64)
    for a, b in zip(*np.nonzero(T)):
        ca, cb = Binv @ I[a] % p, Binv @ I[b] % p
        assert not ca[col:].any() and not cb[col:].any(), "Δ-lift needs an unknown column"
        acc = (acc + T[a, b] * np.outer(X[:, :col] @ ca[:col] % p, X[:, :col] @ cb[:col] % p)) % p
    sol = F.solve(DC, acc.ravel())
    if sol is None: return None
    return sum((sol[k] * I[Cidx[k]] for k in range(len(Cidx))), np.zeros(d, dtype=np.int64)) % p

def lift(X, col):
    """Fast lift: precomputed left inverse of DC (one matmul instead of a Gauss solve per task)."""
    p, d, DL, Binv, Cidx, kdim, DC = (C[k] for k in ('p', 'd', 'DL', 'Binv', 'Cidx', 'kdim', 'DC'))
    if 'Linv' not in C:
        nC = DC.shape[1]
        R, piv = C['F'].rref(np.concatenate([DC % p, np.eye(DC.shape[0], dtype=np.int64)], axis=1))
        assert piv[:nC] == list(range(nC)), "DC not of full column rank"
        C['Linv'] = R[:nC, nC:] % p                    # Linv @ DC = I  (unique Δ-preimage in the complement)
    T = DL[Cidx[col - kdim]]
    ar, br = np.nonzero(T); coefs = T[ar, br]
    ca, cb = Binv[:, ar], Binv[:, br]
    assert not ca[col:].any() and not cb[col:].any(), "Δ-lift needs an unknown column"
    Pa = (X[:, :col] @ Binv[:col, ar]) % p             # P(prefix) in old coords, per Δ entry
    Pb = (X[:, :col] @ Binv[:col, br]) % p
    acc = ((Pa * coefs[None, :]) @ Pb.T) % p           # (P⊗P)Δ(b_col) accumulated
    bv = acc.ravel()
    y = (C['Linv'] @ bv) % p
    if ((DC @ y) % p != bv).any(): return None         # no Δ-preimage in the complement
    out = np.zeros(d, dtype=np.int64); out[Cidx] = y
    return out

def check_ref(Xb, known, new):
    """Reference masked check (original implementation), kept for PCENSUS_XCHECK."""
    p, M, Binv, B = C['p'], C['M'], C['Binv'], C['B']
    st = {'Xb': Xb, 'idx': np.arange(Xb.shape[0])}
    def Papply(v):
        Xb = st['Xb']; c = np.einsum('ij,...j->...i', Binv, v) % p; unk = c[..., known:].any(axis=-1)
        cc = c[:, :known] if c.ndim == 2 else np.broadcast_to(c[:known], (Xb.shape[0], known))
        return np.einsum('nij,nj->ni', Xb[:, :, :known], cc) % p, unk
    def keep(ok):
        st['Xb'] = st['Xb'][ok]; st['idx'] = st['idx'][ok]; return len(st['idx'])
    col = lambda i: st['Xb'][:, :, i] if i == new else st['Xb'][0, :, i]
    for i in range(known):
        v = col(i); r, unk = Papply(v)
        if not keep(~(((r - v) % p).any(axis=1) & ~unk)): return st['idx']
    for i in range(known):
        for j in range(known):
            for comp in (False, True):
                xi, yj = col(i), col(j)
                if comp: xi, yj = (B[:, i] - xi) % p, (B[:, j] - yj) % p
                if xi.ndim == 1 and yj.ndim == 1: h = np.einsum('i,j,ijk->k', xi, yj, M) % p
                elif xi.ndim == 1: h = np.einsum('i,nj,ijk->nk', xi, yj, M) % p
                elif yj.ndim == 1: h = np.einsum('ni,j,ijk->nk', xi, yj, M) % p
                else: h = np.einsum('ni,nj,ijk->nk', xi, yj, M) % p
                r, unk = Papply(h)
                if not keep(~((r if comp else (r - h) % p).any(axis=1) & ~unk)): return st['idx']
    return st['idx']

def check(Xb, known, new):
    """Vectorized masked check; returns exactly the survivor set of check_ref. Conditions whose operands are all
    in columns < col (constant across the batch) and that were fully determined (mask c[col:] = 0) were already
    verified when column col was checked, and never change afterwards — so they are skipped. Validated against
    check_ref on random batches and in full pipeline runs (identical fold sets)."""
    p, d, M, B, Binv, kdim = C['p'], C['d'], C['M'], C['B'], C['Binv'], C['kdim']
    n = Xb.shape[0]
    if n == 0: return np.zeros(0, dtype=np.int64)
    col = new; parent_ran = col > kdim; m = col
    X0 = Xb[0]
    bad = np.zeros(n, dtype=bool)
    sel = None; cs = None; ts = None; U = None
    if m:
        U = np.stack([X0[:, :m], (B[:, :m] - X0[:, :m]) % p])            # P(b_i), (1-P)(b_i) for i < col
        Hoo = np.einsum('axi,ayj,xyk->aijk', U, U, M) % p                # b_i·b_j and (1-P)b_i·(1-P)b_j
        Hc = np.concatenate([U[0].T, Hoo[0].reshape(m * m, d), Hoo[1].reshape(m * m, d)])
        Tc = np.concatenate([U[0].T, Hoo[0].reshape(m * m, d), np.zeros((m * m, d), np.int64)])
        cc = (Hc @ Binv.T) % p                                           # new-basis coords (batch-constant)
        act = ~cc[:, known:].any(axis=1)                                 # determined now
        ver = (~cc[:, col:].any(axis=1)) if parent_ran else np.zeros(len(Hc), dtype=bool)
        sel = np.nonzero(act & ~ver)[0]                                  # to (re)check
        if len(sel): cs, ts = cc[sel][:, :known], Tc[sel]
    CH = 4096
    for s0 in range(0, n, CH):
        Xs = Xb[s0:s0 + CH]; ch = Xs.shape[0]
        b = np.zeros(ch, dtype=bool)
        if m and len(sel):                                               # constant conditions
            R = np.einsum('nxk,qk->qnx', Xs[:, :, :known], cs) % p
            b |= (((R - ts[:, None, :]) % p) != 0).any(axis=2).any(axis=0)
        W0 = Xs[:, :, new]                                               # conditions involving column col
        W = np.stack([W0, (B[None, :, new] - W0) % p])                    # (2, ch, d)
        Hbl, fl = [], []
        if m:
            WL = np.einsum('axi,xyk->aiyk', U, M) % p
            WR = np.einsum('ayj,xyk->ajxk', U, M) % p
            HL = np.einsum('any,aiyk->anik', W, WL) % p                   # b_i·b_col and (1-P)b_i·(1-P)b_col
            HR = np.einsum('anx,ajxk->anjk', W, WR) % p
            for a in (0, 1):
                Hbl += [HL[a].transpose(1, 0, 2), HR[a].transpose(1, 0, 2)]
                fl += [a == 0] * (2 * m)                                 # one flag per condition, not per block
        Hnn = np.einsum('anx,any,xyk->ank', W, W, M) % p
        for a in (0, 1):
            Hbl.append(Hnn[a][None]); fl.append(a == 0)                  # P b_col·P b_col, (1-P) variant
        Hbl.append(W0[None]); fl.append(True)                            # P(P b_col) = P b_col
        H = np.concatenate(Hbl, axis=0)                                  # (Qb, ch, d)
        Tg = H * np.asarray(fl, np.int64)[:, None, None]                 # target: h (comp=False) or 0
        cb = (H @ Binv.T) % p
        unk = cb[:, :, known:].any(axis=2)
        Rb = np.einsum('nxk,qnk->qnx', Xs[:, :, :known], cb[:, :, :known]) % p
        b |= ((((Rb - Tg) % p) != 0).any(axis=2) & ~unk).any(axis=0)
        bad[s0:s0 + ch] = b
    return np.nonzero(~bad)[0]

def extend(task):
    """task = (rep, X, col) with columns 0..col-1 filled: extend by up to STRIDE columns; return
    (rep, next column to fill, survivors as one int16 (m, d, d) array)."""
    rep, X, col = task
    d, p, allK = C['d'], C['p'], C['allK']
    Xs = np.asarray(X, dtype=np.int64)[None]
    steps = STRIDE
    while col < d and steps > 0 and Xs.shape[0] <= MAXB:
        nxt = []
        for Xc in Xs:
            L = lift(Xc, col)
            if XCHECK:
                L2 = lift_ref(Xc, col)
                assert (L is None) == (L2 is None) and (L is None or (L == L2).all()), "fast lift disagrees with reference"
            if L is None: continue
            Xb = np.repeat(Xc[None], len(allK), axis=0)
            Xb[:, :, col] = (L + allK) % p
            keep = check(Xb, col + 1, col)
            if XCHECK:
                keep_ref = check_ref(Xb, col + 1, col)
                assert np.array_equal(keep, keep_ref), "fast check disagrees with reference check"
            if len(keep): nxt.append(Xb[keep])
        if not nxt: return rep, col, np.zeros((0, d, d), np.int16)
        Xs = np.concatenate(nxt); col += 1; steps -= 1
    return rep, col, Xs.astype(np.int16)

def fold_keys(A):
    """(N, d, d) residues → (N,) byte-row keys (the int8 bytes used as fold keys everywhere), for np.unique / np.isin."""
    A = np.ascontiguousarray(A, dtype=np.int8); n = A.shape[1] * A.shape[2]
    return A.reshape(len(A), n).view(np.dtype((np.void, n)))[:, 0]

def orbit_label(folds, mats, p, chunk=4096):
    """G-orbits: dict fold-key -> orbit id, list of representatives, list of orbit sizes."""
    g, gi = mats; lab = {}; reps = []; sizes = []; keys = fold_keys(folds)
    for i, P in enumerate(folds):
        if keys[i].tobytes() in lab: continue
        o = set()
        for s_ in range(0, len(g), chunk):
            Q = (g[s_:s_ + chunk] @ P.astype(np.int64)) @ gi[s_:s_ + chunk] % p        # g P g⁻¹ (matmul, not 5-index einsum)
            o.update(k.tobytes() for k in np.unique(fold_keys(Q)))
        for k in o: lab[k] = len(reps)
        reps.append(P); sizes.append(len(o))
    return lab, reps, sizes

def merge(lab, reps, sizes, extra, p):
    """orbits of ⟨G, extra⟩ (extra normalizes G): union of G-orbits under the extra elements. Returns count, size histogram,
    and one representative fold per merged orbit."""
    par = list(range(len(reps)))
    def find(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    key = lambda A: (A % p).astype(np.int8).tobytes()
    for i, P in enumerate(reps):
        for t, ti in extra:
            j = lab[key(t @ P @ ti)]; a, b = find(i), find(j)
            if a != b: par[a] = b
    comp = collections.defaultdict(int); rep = {}
    for i in range(len(reps)):
        r = find(i); comp[r] += sizes[i]; rep.setdefault(r, reps[i])
    return len(comp), dict(sorted(collections.Counter(comp.values()).items())), list(rep.values())

def orbit_stats(folds, mats, p, chunk=4096):
    """orbit count and size histogram of the group given as stacked (g, g⁻¹) arrays; vectorized over the group."""
    g, gi = mats; seen = set(); sizes = collections.Counter(); key = lambda A: (A % p).astype(np.int8).tobytes()
    for P in folds:
        if key(P) in seen: continue
        o = set()
        for s in range(0, len(g), chunk):
            Q = (g[s:s + chunk] @ P.astype(np.int64)) @ gi[s:s + chunk] % p
            o.update(q.astype(np.int8).tobytes() for q in Q)
        seen |= o; sizes[len(o)] += 1
    return sum(sizes.values()), dict(sorted(sizes.items()))

def lie_G(ctx):
    """Lie(G) on A: for each parallel class (a_1..a_m), derivations a_j ↦ a_i − a_m (i < m), extended to paths; plus E (grading)."""
    basis, arrows, d = ctx['basis'], ctx['arrows'], ctx['d']; idxb = {b: i for i, b in enumerate(basis)}
    classes = collections.defaultdict(list)
    for ai, (s_, t_) in enumerate(arrows): classes[(s_, t_)].append(ai)
    gens = []
    for c in (c for c in classes.values() if len(c) > 1):
        for i in range(len(c) - 1):
            for j in range(len(c)):
                D = np.zeros((d, d), dtype=np.int64)          # derivation with D(a_j) = a_i − a_m
                for col, b in enumerate(basis):
                    if b[0] == 'v': continue
                    for pos, x in enumerate(b[1]):
                        if x != c[j]: continue
                        for y, coef in ((c[i], 1), (c[-1], -1)):
                            D[idxb[('p', b[1][:pos] + (y,) + b[1][pos + 1:])], col] += coef
                gens.append(D)
    E = np.diag([0 if b[0] == 'v' else len(b[1]) for b in basis]).astype(np.int64)
    return gens, E

def rank_Q(rows):
    """Exact rank over Q of an integer matrix (rows of ints): fraction-free elimination in Python ints, each new row
    divided by its content (gcd) so the entries stay small. Replaces sympy.Matrix.rank (same value, ~10³× faster)."""
    import math
    A = [[int(x) for x in r] for r in rows]; A = [r for r in A if any(r)]; r = 0
    for c in range(len(A[0]) if A else 0):
        piv = next((i for i in range(r, len(A)) if A[i][c]), None)
        if piv is None: continue
        A[r], A[piv] = A[piv], A[r]; a = A[r][c]
        for i in range(r + 1, len(A)):
            b = A[i][c]
            if b:
                row = [a * x - b * y for x, y in zip(A[i], A[r])]
                g = math.gcd(*row) or 1; A[i] = [x // g for x in row]
        r += 1
    return r

def q_fold(P, ctx):
    """The representative read over Z with p−1 ↦ −1, if that is an exact fold over Q (idempotent, image and kernel
    subalgebras, Δ-compatible); else None. All conditions as batched integer contractions."""
    p, M, DL, d = ctx['p'], ctx['M'], ctx['DL'], ctx['d']
    Pq = np.where(P == p - 1, -1, P).astype(np.int64); I = np.eye(d, dtype=np.int64)
    if (Pq @ Pq != Pq).any(): return None
    h = np.einsum('ia,jb,ijk->abk', Pq, Pq, M); q = np.einsum('ia,jb,ijk->abk', I - Pq, I - Pq, M)   # P b_a·P b_b, ...
    if (np.einsum('xk,abk->abx', Pq, h) != h).any() or np.einsum('xk,abk->abx', Pq, q).any(): return None
    if (np.einsum('ml,mab->lab', Pq, DL) != np.einsum('ai,lij,bj->lab', Pq, DL, Pq)).any(): return None
    return Pq

def transversal(reps, ctx):
    """For one representative per (G × grading)-orbit: r_G = rank{[D,P]}, r_H = rank with E added (mod p, and over Q
    when the representative lifts to an exact Q-fold with entries p−1 ↦ −1). Counts orbits with r_H > r_G."""
    p, F = ctx['p'], ctx['F']; gens, E = lie_G(ctx)
    cnt = collections.Counter()
    for P in reps:
        P = np.asarray(P, dtype=np.int64)
        rows = [((D @ P - P @ D) % p).ravel() for D in gens]
        rG = F.rank(np.array(rows)) if rows else 0; rH = F.rank(np.array(rows + [((E @ P - P @ E) % p).ravel()]))
        Pq = q_fold(P, ctx)
        if Pq is not None:
            rowsq = [(D @ Pq - Pq @ D).ravel() for D in gens]
            rGq = rank_Q(rowsq); rHq = rank_Q(rowsq + [(E @ Pq - Pq @ E).ravel()])
            cnt['lifted'] += 1; cnt['lifted_transversal'] += int(rHq > rGq)
        cnt['transversal_modp'] += int(rH > rG)
    return dict(cnt)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('quiver'); ap.add_argument('p', type=int)
    ap.add_argument('--jobs', type=int, default=os.cpu_count()); ap.add_argument('--checkpoint'); ap.add_argument('--out')
    ap.add_argument('--no-orbits', action='store_true'); ap.add_argument('--save-folds')
    ap.add_argument('--save-reps')                                       # pickle G-orbit reps + local-torus components
    ap.add_argument('--engine', choices=('pool', 'torch'), default='pool'); ap.add_argument('--device')
    ap.add_argument('--kernel', choices=('auto', 'torch', 'ref'), default='auto')   # auto: torch if importable
    ap.add_argument('--kernel-cache')                                    # default: CHECKPOINT minus .pkl, + .kernel.npz
    ap.add_argument('--no-symmetry', action='store_true')                # compute every fibre (no G × T transport)
    a = ap.parse_args(); name, p = a.quiver, a.p
    t0 = time.time(); ctx = setup(name, p); C.update(ctx)
    d, kdim = ctx['d'], ctx['kdim']
    phis = grading(ctx)
    cache = a.kernel_cache or (a.checkpoint[:-4] if a.checkpoint and a.checkpoint.endswith('.pkl') else a.checkpoint)
    Kf, reps = kernel_stage(ctx, phis, a.kernel, a.device, cache and cache + '.kernel.npz')
    print(f"  kernel stage: {time.time() - t0:.1f}s", flush=True)
    eng = f"jobs={a.jobs}" if a.engine == 'pool' else "engine=torch"
    print(f"{name} F_{p}: d={d}, dim K={kdim}, {len(Kf)} kernel folds in {len(reps)} torus orbits; {eng}", flush=True)
    fibres = {}
    if a.checkpoint and os.path.exists(a.checkpoint):
        ck = pickle.load(open(a.checkpoint, 'rb'))
        assert (ck['name'], ck['p'], ck['nreps']) == (name, p, len(reps)), "checkpoint is for another run"
        # fibres are stored in the column basis B, which depends on the column order (ETF_CORDER); checkpoints from
        # before v2.51 carry no Cidx and were made with the default order
        assert ck.get('Cidx', None if corders() else list(ctx['Cidx'])) == list(ctx['Cidx']), \
            "checkpoint was made with another column order (ETF_CORDER / first of ETF_CORDERS)"
        fibres = ck['fibres']; print(f"  resumed: {len(fibres)}/{len(reps)} fibres done", flush=True)
    last_save = [time.time()]
    def save(force=False):
        if not a.checkpoint or (not force and time.time() - last_save[0] < 30): return
        tmp = a.checkpoint + '.tmp'
        pickle.dump(dict(name=name, p=p, nreps=len(reps), fibres=fibres, Cidx=list(ctx['Cidx'])), open(tmp, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp, a.checkpoint); last_save[0] = time.time()
    if a.no_symmetry: src, trans = list(range(len(reps))), {}
    else:
        src, trans, ctx['N'] = symmetry_plan(ctx, reps, phis)
        print(f"  symmetry: {len(reps)} torus orbits in {len(src)} H-orbits (H = G × T); {len(trans)} fibres by transport; "
              f"|N| = |ker(H → H|_K)| = {len(ctx['N'][0])}", flush=True)
    if a.engine == 'torch':
        import fibre_torch
        orders = corders()
        if len(orders) > 1:                  # order portfolio: per-fibre breadth-first search, each order in turn
            alts = []
            for o in orders[1:]:
                c2 = setup(name, p, order=o); c2['N'] = ctx.get('N'); alts.append(c2)
            fibre_torch.portfolio(ctx, alts, reps, fibres, save, log=lambda s: print(s, flush=True), device=a.device,
                                  only=src)
        fibre_torch.run(ctx, reps, fibres, save, log=lambda s: print(s, flush=True), device=a.device, only=src)
    else:
        pool_fibres(name, p, a.jobs, ctx, reps, fibres, save, t0, only=src)
    transport(ctx, reps, fibres, trans)
    assert len(fibres) == len(reps)
    save(force=True)
    finish(a, name, p, ctx, Kf, reps, phis, fibres, t0)

def pool_fibres(name, p, jobs, ctx, reps, fibres, save, t0, only=None):
    d, kdim, Kmat = ctx['d'], ctx['kdim'], ctx['Kmat']
    stack = []; outstanding = collections.Counter(); partial = collections.defaultdict(list)
    todo = [r for r in (range(len(reps)) if only is None else only) if r not in fibres]
    for r in reversed(todo):
        X = np.zeros((d, d), dtype=np.int64); X[:, :kdim] = Kmat @ reps[r][0] % p
        if kdim == d: fibres[r] = [X.astype(np.int16)]; continue
        stack.append((r, X.astype(np.int16), kdim)); outstanding[r] += 1
    with ProcessPoolExecutor(jobs, initializer=init_worker, initargs=(name, p)) as ex:
        inflight = set(); lastlog = time.time(); ntasks = 0
        while stack or inflight:
            while stack and len(inflight) < 4 * jobs:
                inflight.add(ex.submit(extend, stack.pop()))
            fin, inflight = wait(inflight, return_when=FIRST_COMPLETED)
            for fu in fin:
                r, col, kids = fu.result(); outstanding[r] -= 1; ntasks += 1
                if col == d: partial[r] += list(kids)
                else:
                    for k in kids: stack.append((r, k, col)); outstanding[r] += 1
                if outstanding[r] == 0:
                    fibres[r] = partial.pop(r, []); del outstanding[r]; save()
            if time.time() - lastlog > 60:
                print(f"  {len(fibres)}/{len(reps)} fibres done, {len(outstanding)} open, {len(stack)} queued tasks, "
                      f"{ntasks} tasks run, t={time.time()-t0:.0f}s", flush=True); lastlog = time.time()

def build_G(ctx):
    """G = ∏ over parallel classes of the column-sum-1 substitution group, as (Φ, Φ⁻¹) pairs acting on A = kQ.
    Memoized in ctx (the symmetry plan and finish both use it)."""
    if 'G' not in ctx: ctx['G'] = _build_G(ctx)
    return ctx['G']

def _build_G(ctx):
    p, d, basis = ctx['p'], ctx['d'], ctx['basis']
    # column-sum-1 substitutions on every parallel class (as in gcensus.py), stacked
    F = ctx['F']; arrows = ctx['arrows']
    classes = collections.defaultdict(list)
    for ai, (s_, t_) in enumerate(arrows): classes[(s_, t_)].append(ai)
    classes = [c for c in classes.values() if len(c) > 1]
    def mats(m):      # column sums 1: the last row is determined by the first m − 1
        out = []
        for vals in itertools.product(range(p), repeat=(m - 1) * m):
            top = np.array(vals, dtype=np.int64).reshape(m - 1, m)
            g = np.vstack([top, (1 - top.sum(axis=0)) % p])
            if F.rank(g) == m: out.append(g)
        return out
    idxb = {b: i for i, b in enumerate(basis)}
    def phi(subs):
        Phi = np.zeros((d, d), dtype=np.int64)
        for i, b in enumerate(basis):
            if b[0] == 'v': Phi[i, i] = 1; continue
            for choice in itertools.product(*[list(subs.get(x, {x: 1}).items()) for x in b[1]]):
                j = idxb[('p', tuple(x for x, _ in choice))]
                Phi[j, i] = (Phi[j, i] + int(np.prod([c for _, c in choice]))) % p
        return Phi
    G = []
    for gs in itertools.product(*[mats(len(c)) for c in classes]):
        subs = {}
        for c, g in zip(classes, gs):
            for jc, x in enumerate(c): subs[x] = {c[ir]: int(g[ir, jc]) for ir in range(len(c)) if g[ir, jc]}
        Phi = phi(subs); G.append((Phi, F.inv(Phi)))
    return G

def finish(a, name, p, ctx, Kf, reps, phis, fibres, t0):
    d, kdim, B, Binv, Kmat, basis = ctx['d'], ctx['kdim'], ctx['B'], ctx['Binv'], ctx['Kmat'], ctx['basis']
    # materialize all folds (old coordinates) via the grading on each torus orbit; canonical order within each fibre
    # (one int8 (N, d, d) array; φ_l is diagonal, so φ_l P φ_l⁻¹ = P ∘ (diag φ_l ⊗ diag φ_l⁻¹), batched per fibre)
    fib = collections.Counter(); n_out = sum(len(fibres[r]) * len(ls) for r, (_, ls) in enumerate(reps))
    folds = np.empty((n_out, d, d), np.int8); at = 0
    scale = {l: np.outer(np.diag(ph), np.diag(phi_)) for l, (ph, phi_) in phis.items()}
    for r, (PK, ls) in enumerate(reps):
        Xs = sorted(fibres[r], key=lambda X: (X % p).astype(np.int8).tobytes())
        base = np.array(Xs, dtype=np.int64).reshape(-1, d, d) @ Binv % p
        for l in ls:
            fib[len(base)] += 1
            folds[at:at + len(base)] = base * scale[l] % p; at += len(base)
    assert sum(fib.values()) == len(Kf) and at == n_out
    assert len(np.unique(fold_keys(folds))) == len(folds)
    res = dict(name=name, p=p, d=d, kdim=kdim, nK=len(Kf), reps=len(reps), N=len(folds), fib=dict(sorted(fib.items())))
    if not a.no_orbits:
        F = ctx['F']; arrows = ctx['arrows']
        G = build_G(ctx)
        Gg = (np.array([g for g, _ in G]), np.array([gi for _, gi in G]))
        lab, greps, gsizes = orbit_label(folds, Gg, p)
        res['G'] = len(G); res['orbits_G'] = len(greps); res['sizes_G'] = dict(sorted(collections.Counter(gsizes).items()))
        res['orbits_GxGr'], res['sizes_GxGr'], hreps = merge(lab, greps, gsizes, list(phis.values()), p)
        res['transversal'] = transversal(hreps, ctx)
        # class-weight torus: φ_w scales each parallel class (single arrows included) by w_c; keep the w preserving Folds
        cls = collections.defaultdict(list)
        for ai, st in enumerate(arrows): cls[st].append(ai)
        # Folds is a union of G-orbits (FS = set(lab)) and φ_w commutes with G (it scales each parallel class by one
        # scalar), so φ_w preserves Folds iff it maps every G-orbit representative into FS.
        cl = list(cls.values()); Tpres = []; FS = set(lab)
        for w in itertools.product(range(1, p), repeat=len(cl)):
            wa = {a: w[k] for k, c in enumerate(cl) for a in c}
            diag = [1 if b[0] == 'v' else int(np.prod([wa[x] for x in b[1]])) % p for b in basis]
            ph, phi_ = np.diag(diag), np.diag([pow(x, p - 2, p) for x in diag])
            if all(((ph @ P_ @ phi_) % p).astype(np.int8).tobytes() in FS for P_ in greps): Tpres.append((ph, phi_))
        res['classes'] = len(cl); res['torus_preserving'] = len(Tpres); res['torus_all'] = (p - 1) ** len(cl)
        res['orbits_GxT'], res['sizes_GxT'], _ = merge(lab, greps, gsizes, Tpres, p)
        # local class torus (v2.24): for each G-orbit representative P, every weighting w with P compatible with
        # Δ_w = Σ_c w_c Δ_c moves P to a fold φ_w P φ_w⁻¹ (v2.23 Lemma, applied to Δ_w directly). W(P) is constant on
        # G-orbits (G substitutes within classes, so Δ_w g = (g⊗g)Δ_w). Merge G-orbits along these moves.
        idxb = {b: i for i, b in enumerate(basis)}
        DLc = []
        for c in cl:
            D_ = np.zeros((d, d, d), dtype=np.int64)
            for l, b in enumerate(basis):
                if b[0] == 'v': continue
                for k, x in enumerate(b[1]):
                    if x not in c: continue
                    pre = idxb[('p', b[1][:k])] if k > 0 else idxb[('v', arrows[x][0])]
                    suf = idxb[('p', b[1][k + 1:])] if k + 1 < len(b[1]) else idxb[('v', arrows[x][1])]
                    D_[l, pre, suf] += 1
            DLc.append(D_)
        DLc = np.array(DLc)                                              # (classes, d, d, d)
        par = list(range(len(greps)))
        def find(x):
            while par[x] != x: par[x] = par[par[x]]; x = par[x]
            return x
        wcount = collections.Counter()
        Ws = list(itertools.product(range(1, p), repeat=len(cl)))
        Warr = np.array(Ws, dtype=np.int64)                                                # (W, classes)
        scs = []
        for w in Ws:
            wa = {a: w[k] for k, c in enumerate(cl) for a in c}
            diag = [1 if b[0] == 'v' else int(np.prod([wa[x] for x in b[1]])) % p for b in basis]
            scs.append(np.outer(diag, [pow(x, p - 2, p) for x in diag]))
        for i, P_ in enumerate(greps):
            P_ = P_.astype(np.int64)
            # P compatible with Δ_w: Σ_m P[m,l] Δ_w[m] = P Δ_w[l] Pᵀ for every l. The defect is linear in w:
            # Σ_c w_c E_c with E_c = Σ_m P[m,l] Δ_c[m] − P Δ_c[l] Pᵀ, so all weightings are tested by one product
            # (v2.56; before, two einsums over (W, d, d, d), the second one naive, cost W·d⁵ per representative).
            Ecl = (np.einsum('ml,cmab->clab', P_, DLc) - np.einsum('ai,clij,bj->clab', P_, DLc, P_, optimize=True)) % p
            ok = ((Warr @ Ecl.reshape(len(cl), -1)) % p == 0).all(1)
            for q in np.nonzero(ok)[0]:
                kq = (P_ * scs[q] % p).astype(np.int8).tobytes(); assert kq in lab, "local torus move left the fold set"
                a_, b_ = find(i), find(lab[kq])
                if a_ != b_: par[a_] = b_
            wcount[int(ok.sum())] += 1
        comp = collections.Counter(find(i) for i in range(len(greps)))
        res['orbits_GxTloc'] = len(comp)
        if getattr(a, 'save_reps', None):      # G-orbit representatives, sizes, and their local-torus components
            pickle.dump(dict(name=name, p=p, basis=basis, greps=np.array(greps, dtype=np.int8), gsizes=gsizes,
                             loc=[find(i) for i in range(len(greps))]), open(a.save_reps, 'wb'))
        res['local_torus_sizes'] = dict(sorted(wcount.items()))         # how many weightings each G-orbit admits
    res['t'] = round(time.time() - t0)
    line = json.dumps(res); print(line, flush=True)
    if a.out: open(a.out, 'a').write(line + "\n")
    if a.save_folds: pickle.dump(dict(name=name, p=p, folds=list(folds)), open(a.save_folds, 'wb'))

if __name__ == '__main__':
    main()