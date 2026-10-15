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
    'two_pairs': (3, [(0, 1), (0, 1), (1, 2), (1, 2)]),                 # 1⇉2⇉3
    'two_pairs_sep': (4, [(0, 1), (0, 1), (1, 2), (2, 3), (2, 3)]),     # 1⇉2→3⇉4
    'pair_after3': (5, [(0, 1), (0, 1), (1, 2), (2, 3), (3, 4)]),       # 1⇉2→3→4→5
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

def setup(name, p):
    nv, arrows = QUIVERS[name]; F = Field(p)
    basis, M, DL = build(nv, arrows); d = len(basis); I = np.eye(d, dtype=np.int64)
    Dmat = DL.reshape(d, d * d).T % p
    Kb = F.nullspace(Dmat); kdim = len(Kb)
    Cidx = []; cur = list(Kb)
    for i in range(d):
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

def orbit_label(folds, mats, p, chunk=4096):
    """G-orbits: dict fold-key -> orbit id, list of representatives, list of orbit sizes."""
    g, gi = mats; lab = {}; reps = []; sizes = []; key = lambda A: (A % p).astype(np.int8).tobytes()
    for P in folds:
        if key(P) in lab: continue
        o = set()
        for s_ in range(0, len(g), chunk):
            Q = np.einsum('gij,jk,gkl->gil', g[s_:s_ + chunk], P, gi[s_:s_ + chunk]) % p
            o.update(q.astype(np.int8).tobytes() for q in Q)
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
            Q = np.einsum('gij,jk,gkl->gil', g[s:s + chunk], P, gi[s:s + chunk]) % p
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

def transversal(reps, ctx):
    """For one representative per (G × grading)-orbit: r_G = rank{[D,P]}, r_H = rank with E added (mod p, and over Q
    when the representative lifts to an exact Q-fold with entries p−1 ↦ −1). Counts orbits with r_H > r_G."""
    import sympy as sp
    p, F = ctx['p'], ctx['F']; gens, E = lie_G(ctx)
    cnt = collections.Counter()
    M, DL, d = ctx['M'], ctx['DL'], ctx['d']
    def q_fold(P):
        Pq = np.where(P == p - 1, -1, P).astype(np.int64)
        if (Pq @ Pq != Pq).any(): return None
        I = np.eye(d, dtype=np.int64)
        for i in range(d):
            for j in range(d):
                h = np.einsum('i,j,ijk->k', Pq[:, i], Pq[:, j], M); q = np.einsum('i,j,ijk->k', (I - Pq)[:, i], (I - Pq)[:, j], M)
                if (Pq @ h != h).any() or (Pq @ q).any(): return None
        for l in range(d):
            lhs = np.einsum('m,mab->ab', Pq[:, l], DL)
            if (lhs != Pq @ DL[l] @ Pq.T).any(): return None
        return Pq
    for P in reps:
        rows = [((D @ P - P @ D) % p).ravel() for D in gens]
        rG = F.rank(np.array(rows)) if rows else 0; rH = F.rank(np.array(rows + [((E @ P - P @ E) % p).ravel()]))
        Pq = q_fold(P); tag = 'modp'
        if Pq is not None:
            rowsq = [(D @ Pq - Pq @ D).ravel().tolist() for D in gens]
            rGq = sp.Matrix(rowsq).rank() if rowsq else 0; rHq = sp.Matrix(rowsq + [(E @ Pq - Pq @ E).ravel().tolist()]).rank()
            cnt['lifted'] += 1; cnt['lifted_transversal'] += int(rHq > rGq)
        cnt['transversal_modp'] += int(rH > rG)
    return dict(cnt)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('quiver'); ap.add_argument('p', type=int)
    ap.add_argument('--jobs', type=int, default=os.cpu_count()); ap.add_argument('--checkpoint'); ap.add_argument('--out')
    ap.add_argument('--no-orbits', action='store_true'); ap.add_argument('--save-folds')
    ap.add_argument('--engine', choices=('pool', 'torch'), default='pool'); ap.add_argument('--device')
    a = ap.parse_args(); name, p = a.quiver, a.p
    t0 = time.time(); ctx = setup(name, p); C.update(ctx)
    d, kdim, B, Binv, Kmat, basis = ctx['d'], ctx['kdim'], ctx['B'], ctx['Binv'], ctx['Kmat'], ctx['basis']
    Kf = kernel_folds(ctx)
    lens = [0 if b[0] == 'v' else len(b[1]) for b in basis]
    phis = {l: (np.diag([pow(l, k, p) for k in lens]), np.diag([pow(pow(l, p - 2, p), k, p) for k in lens])) for l in range(1, p)}
    kkey = lambda A: (A % p).tobytes(); done_k = set(); reps = []      # (PK, {kkey: l})
    for PK in Kf:
        if kkey(PK) in done_k: continue
        torb = {}
        for l, (ph, phi_) in phis.items():
            PKl = (Binv @ ph @ (Kmat @ PK % p) % p)[:kdim] % p
            PKl = PKl @ (Binv @ phi_ @ Kmat % p)[:kdim] % p
            torb.setdefault(kkey(PKl), l)
        done_k |= set(torb); reps.append((PK, list(torb.values())))
    eng = f"jobs={a.jobs}" if a.engine == 'pool' else "engine=torch"
    print(f"{name} F_{p}: d={d}, dim K={kdim}, {len(Kf)} kernel folds in {len(reps)} torus orbits; {eng}", flush=True)
    fibres = {}
    if a.checkpoint and os.path.exists(a.checkpoint):
        ck = pickle.load(open(a.checkpoint, 'rb'))
        assert (ck['name'], ck['p'], ck['nreps']) == (name, p, len(reps)), "checkpoint is for another run"
        fibres = ck['fibres']; print(f"  resumed: {len(fibres)}/{len(reps)} fibres done", flush=True)
    last_save = [time.time()]
    def save(force=False):
        if not a.checkpoint or (not force and time.time() - last_save[0] < 30): return
        tmp = a.checkpoint + '.tmp'
        pickle.dump(dict(name=name, p=p, nreps=len(reps), fibres=fibres), open(tmp, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp, a.checkpoint); last_save[0] = time.time()
    if a.engine == 'torch':
        import fibre_torch
        fibre_torch.run(ctx, reps, fibres, save, log=lambda s: print(s, flush=True), device=a.device)
    else:
        pool_fibres(name, p, a.jobs, ctx, reps, fibres, save, t0)
    save(force=True)
    finish(a, name, p, ctx, Kf, reps, phis, fibres, t0)

def pool_fibres(name, p, jobs, ctx, reps, fibres, save, t0):
    d, kdim, Kmat = ctx['d'], ctx['kdim'], ctx['Kmat']
    stack = []; outstanding = collections.Counter(); partial = collections.defaultdict(list)
    todo = [r for r in range(len(reps)) if r not in fibres]
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
    """G = ∏ over parallel classes of the column-sum-1 substitution group, as (Φ, Φ⁻¹) pairs acting on A = kQ."""
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
    folds = []; fib = collections.Counter()
    for r, (PK, ls) in enumerate(reps):
        base = [X.astype(np.int64) @ Binv % p for X in sorted(fibres[r], key=lambda X: (X % p).astype(np.int8).tobytes())]
        for l in ls:
            ph, phi_ = phis[l]; fib[len(base)] += 1
            folds += [ph @ P_ @ phi_ % p for P_ in base]
    assert sum(fib.values()) == len(Kf)
    assert len({P_.astype(np.int8).tobytes() for P_ in folds}) == len(folds)
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
        cl = list(cls.values()); FS = set(lab); Tpres = []
        for w in itertools.product(range(1, p), repeat=len(cl)):
            wa = {a: w[k] for k, c in enumerate(cl) for a in c}
            diag = [1 if b[0] == 'v' else int(np.prod([wa[x] for x in b[1]])) % p for b in basis]
            ph, phi_ = np.diag(diag), np.diag([pow(x, p - 2, p) for x in diag])
            if all(((ph @ P_ @ phi_) % p).astype(np.int8).tobytes() in FS for P_ in folds): Tpres.append((ph, phi_))
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
        for i, P_ in enumerate(greps):
            ok_w = 0
            for w in itertools.product(range(1, p), repeat=len(cl)):
                Dw = np.einsum('c,cmab->mab', np.array(w), DLc) % p
                if not all(((np.einsum('m,mab->ab', P_[:, l], Dw) - P_ @ Dw[l] @ P_.T) % p == 0).all() for l in range(d)): continue
                ok_w += 1
                wa = {a: w[k] for k, c in enumerate(cl) for a in c}
                diag = [1 if b[0] == 'v' else int(np.prod([wa[x] for x in b[1]])) % p for b in basis]
                Q_ = (np.diag(diag) @ P_ @ np.diag([pow(x, p - 2, p) for x in diag])) % p
                kq = Q_.astype(np.int8).tobytes(); assert kq in lab, "local torus move left the fold set"
                a_, b_ = find(i), find(lab[kq])
                if a_ != b_: par[a_] = b_
            wcount[ok_w] += 1
        comp = collections.Counter(find(i) for i in range(len(greps)))
        res['orbits_GxTloc'] = len(comp)
        res['local_torus_sizes'] = dict(sorted(wcount.items()))         # how many weightings each G-orbit admits
    res['t'] = round(time.time() - t0)
    line = json.dumps(res); print(line, flush=True)
    if a.out: open(a.out, 'a').write(line + "\n")
    if a.save_folds: pickle.dump(dict(name=name, p=p, folds=[P_.astype(np.int8) for P_ in folds]), open(a.save_folds, 'wb'))

if __name__ == '__main__':
    main()