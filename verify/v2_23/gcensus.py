"""General census of Δ-compatible folds of the path algebra kQ of an acyclic quiver Q over F_p (kernel fibration).

A fold is an idempotent P on A = kQ whose image H and kernel Q are subalgebras, with ΔP = (P⊗P)Δ, Δ = arrow deletion.
Method (generalizing q123s2.py):
  1. K = ker Δ (computed mod p). Every fold restricts to an algebra fold of K (P(K) = H∩K, ker = Q∩K).
  2. New basis: a basis of K, then complement basis paths in length order. Δ is injective on the complement, so for each
     complement path x, P(x) = (unique Δ-preimage of (P⊗P)Δx in the complement) + (unknown K-component).
  3. K-components are enumerated one column at a time (p^dimK candidates); masked necessary conditions on the known
     columns are applied after each column, progressively filtering; the last column gets the full checks.
Orbits: G = ∏ over parallel classes of {g ∈ GL_m : column sums 1}, acting by substituting parallel arrows inside paths.
Usage: python3 gcensus.py QUIVER p     (QUIVER is a key of QUIVERS below)"""
import sys, itertools, time, collections
import numpy as np

QUIVERS = {   # name: (number of vertices, arrows as (src, tgt))
    'K2':        (2, [(0, 1), (0, 1)]),                      # 1⇉2
    'pair_after':  (3, [(0, 1), (0, 1), (1, 2)]),            # 1⇉2→3
    'pair_before': (3, [(1, 2), (1, 2), (0, 1)]),            # 1→2⇉3  (opposite of pair_after)
    'pair_srcout': (3, [(0, 1), (0, 1), (0, 2)]),            # 1⇉2, 1→3
    'pair_tgtin':  (3, [(0, 1), (0, 1), (2, 1)]),            # 1⇉2, 3→2  (opposite of pair_srcout)
    'pair_after2': (4, [(0, 1), (0, 1), (1, 2), (2, 3)]),    # 1⇉2→3→4
    'pair_mid':    (4, [(1, 2), (1, 2), (0, 1), (2, 3)]),    # 1→2⇉3→4
    'pair_tgtfork': (4, [(0, 1), (0, 1), (1, 2), (1, 3)]),   # 1⇉2, 2→3, 2→4  (fork after the pair's target)
    'pair_srcboth': (4, [(0, 1), (0, 1), (0, 2), (3, 1)]),   # 1⇉2, 1→3, 4→2  (forks at both ends)
}

def build(nv, arrows):
    """Path basis (vertices first, then paths by length), multiplication tensor M, arrow-deletion tensor DL."""
    paths = [()]; frontier = [(a,) for a in range(len(arrows))]; allp = []
    while frontier:
        allp += frontier
        frontier = [p_ + (a,) for p_ in frontier for a in range(len(arrows)) if arrows[p_[-1]][1] == arrows[a][0]]
    basis = [('v', v) for v in range(nv)] + [('p', p_) for p_ in allp]
    idx = {b: i for i, b in enumerate(basis)}; d = len(basis)
    src = lambda b: b[1] if b[0] == 'v' else arrows[b[1][0]][0]
    tgt = lambda b: b[1] if b[0] == 'v' else arrows[b[1][-1]][1]
    M = np.zeros((d, d, d), dtype=np.int64)
    for i, x in enumerate(basis):
        for j, y in enumerate(basis):
            if tgt(x) != src(y): continue
            if x[0] == 'v': M[i, j, j] = 1
            elif y[0] == 'v': M[i, j, i] = 1
            else: M[i, j, idx[('p', x[1] + y[1])]] = 1
    DL = np.zeros((d, d, d), dtype=np.int64)
    for l, b in enumerate(basis):
        if b[0] == 'v': continue
        s = b[1]
        for k, a in enumerate(s):
            pre = idx[('p', s[:k])] if k > 0 else idx[('v', arrows[a][0])]
            suf = idx[('p', s[k + 1:])] if k + 1 < len(s) else idx[('v', arrows[a][1])]
            DL[l, pre, suf] += 1
    return basis, M, DL

class Field:
    def __init__(self, p): self.p = p
    def rref(self, A):
        p = self.p; A = A.copy() % p; r = 0; piv = []
        for c in range(A.shape[1]):
            pr = next((i for i in range(r, A.shape[0]) if A[i, c]), None)
            if pr is None: continue
            A[[r, pr]] = A[[pr, r]]; A[r] = A[r] * pow(int(A[r, c]), p - 2, p) % p
            for i in range(A.shape[0]):
                if i != r and A[i, c]: A[i] = (A[i] - A[i, c] * A[r]) % p
            piv.append(c); r += 1
        return A[:r], piv
    def rank(self, A): return len(self.rref(A)[1]) if A.size else 0
    def nullspace(self, A):
        R, piv = self.rref(A); n = A.shape[1]; free = [c for c in range(n) if c not in piv]; out = []
        for f in free:
            v = np.zeros(n, dtype=np.int64); v[f] = 1
            for i, c in enumerate(piv): v[c] = (-R[i, f]) % self.p
            out.append(v)
        return out
    def inv(self, A):
        n = A.shape[0]; R, piv = self.rref(np.concatenate([A % self.p, np.eye(n, dtype=np.int64)], axis=1))
        assert piv[:n] == list(range(n)); return R[:, n:]
    def solve(self, A, b):
        R, piv = self.rref(np.concatenate([A % self.p, (b % self.p).reshape(-1, 1)], axis=1))
        n = A.shape[1]
        if n in piv: return None
        x = np.zeros(n, dtype=np.int64)
        for i, c in enumerate(piv): x[c] = R[i, n]
        return x

def census(name, p, verbose=True):
    nv, arrows = QUIVERS[name]; F = Field(p)
    basis, M, DL = build(nv, arrows); d = len(basis); I = np.eye(d, dtype=np.int64)
    mul = lambda x, y: np.einsum('i,j,ijk->k', x, y, M) % p
    Dmat = DL.reshape(d, d * d).T % p                                   # Δ as (d²×d)
    Kb = F.nullspace(Dmat); kdim = len(Kb)
    # complement: greedily add basis paths (length order) independent of K
    Cidx = []; cur = list(Kb)
    for i in range(d):
        if F.rank(np.array(cur + [I[i]])) > len(cur): cur.append(I[i]); Cidx.append(i)
    B = np.array(cur).T % p; Binv = F.inv(B)                            # columns: new basis in old coords
    DC = Dmat[:, Cidx]                                                  # Δ on complement (injective)
    assert F.rank(DC) == len(Cidx)
    Kmat = np.array(Kb).T
    # ---- 1. algebra folds of K (by complementary subalgebra pairs; coordinates in Kb)
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
    Kfolds = []
    for H in subalgs:
        for Q in subalgs:
            if H.shape[0] + Q.shape[0] != kdim: continue
            Bm = np.vstack([H, Q]).reshape(kdim, kdim)
            if F.rank(Bm) != kdim: continue
            Bt = Bm.T % p; E = np.diag([1] * H.shape[0] + [0] * Q.shape[0])
            Kfolds.append(Bt @ E @ F.inv(Bt) % p)
    # ---- 2. staged fibres. X[:, i] = P(new basis vector i) in old coords; P v = X @ (Binv v)
    newc = lambda v: (Binv @ v) % p                                      # new-basis coordinates
    order = list(range(kdim, d))                                         # complement columns, length order
    def lift(X, known, col):
        x = I[Cidx[col - kdim]]; T = DL[Cidx[col - kdim]]
        acc = np.zeros((d, d), dtype=np.int64)
        for a, b in zip(*np.nonzero(T)):
            ca, cb = newc(I[a]), newc(I[b])
            assert not ca[known:].any() and not cb[known:].any(), "Δ-lift needs an unknown column"
            Pa, Pb_ = X[:, :known] @ ca[:known] % p, X[:, :known] @ cb[:known] % p
            acc = (acc + T[a, b] * np.outer(Pa, Pb_)) % p
        sol = F.solve(DC, acc.ravel())
        if sol is None: return None
        return sum((sol[k] * I[Cidx[k]] for k in range(len(Cidx))), np.zeros(d, dtype=np.int64)) % p
    allK = np.array(list(itertools.product(range(p), repeat=kdim)), dtype=np.int64) @ np.array(Kb) % p
    NB = [newc(I[i]) for i in range(d)]
    def check(Xb, known, new):
        """masked conditions on columns < known; only column `new` varies across the batch. Returns surviving indices."""
        idx = np.arange(Xb.shape[0]); X0 = Xb[0]
        state = {'Xb': Xb, 'idx': idx}
        def Papply(v):   # batched or single vector v (old coords) -> P v, plus mask of "uses unknown columns"
            Xb = state['Xb']
            c = (np.einsum('ij,...j->...i', Binv, v)) % p
            unk = c[..., known:].any(axis=-1)
            r = np.einsum('nij,nj->ni', Xb[:, :, :known], c[:, :known] if c.ndim == 2 else np.broadcast_to(c[:known], (Xb.shape[0], known))) % p
            return r, unk
        def keep(ok):
            state['Xb'] = state['Xb'][ok]; state['idx'] = state['idx'][ok]
            return len(state['idx'])
        col = lambda i: state['Xb'][:, :, i] if i == new else state['Xb'][0, :, i]
        for i in range(known):                                           # idempotency P(P b_i) = P b_i
            v = col(i); r, unk = Papply(v)
            bad = ((r - v) % p).any(axis=1) & ~unk
            if not keep(~bad): return state['idx']
        for i in range(known):
            for j in range(known):
                for comp in (False, True):
                    xi, yj = col(i), col(j)
                    if comp: xi, yj = (B[:, i] - xi) % p, (B[:, j] - yj) % p
                    if xi.ndim == 1 and yj.ndim == 1: h = mul(xi, yj)
                    elif xi.ndim == 1: h = np.einsum('i,nj,ijk->nk', xi, yj, M) % p
                    elif yj.ndim == 1: h = np.einsum('ni,j,ijk->nk', xi, yj, M) % p
                    else: h = np.einsum('ni,nj,ijk->nk', xi, yj, M) % p
                    r, unk = Papply(h)
                    bad = (r if comp else (r - h) % p).any(axis=1) & ~unk
                    if not keep(~bad): return state['idx']
        return state['idx']
    folds = []; fib = collections.Counter(); t0 = time.time()
    lens = np.array([0 if b[0] == 'v' else len(b[1]) for b in basis])
    phis = {l: (np.diag([pow(l, int(k), p) for k in lens]), np.diag([pow(pow(l, p - 2, p), int(k), p) for k in lens])) for l in range(1, p)}
    Kpinv = np.array([newc(v)[:kdim] for v in np.eye(d, dtype=np.int64)]).T      # old coords -> K coords (on K)
    kkey = lambda A: (A % p).tobytes()
    done = {}
    for ki, PK in enumerate(Kfolds):
        if kkey(PK) in done: continue
        torb = {}
        for l, (ph, phi_) in phis.items():
            PKl = (Binv @ ph @ (Kmat @ PK % p) % p)[:kdim] % p                      # columns: φ P(k_i), K coords
            PKl = PKl @ (Binv @ phi_ @ Kmat % p)[:kdim] % p                          # ∘ φ⁻¹ on K
            torb.setdefault(kkey(PKl), l)
        for kk in torb: done[kk] = True
        X = np.zeros((d, d), dtype=np.int64); X[:, :kdim] = Kmat @ PK % p
        layer = [X]
        for col in order:
            nxt = []
            for Xc in layer:
                L = lift(Xc, col, col)
                if L is None: continue
                Xb = np.repeat(Xc[None], len(allK), axis=0); Xb[:, :, col] = (L + allK) % p
                nxt += list(Xb[check(Xb, col + 1, col)])
            layer = nxt
            if not layer: break
        base = [Xf @ Binv % p for Xf in layer]                          # P in old coordinates
        for kk, l in torb.items():
            ph, phi_ = phis[l]; fib[len(layer)] += 1
            folds += [ph @ P_ @ phi_ % p for P_ in base]
    # ---- 3. orbits under G (column-sum-1 substitutions on each parallel class)
    classes = collections.defaultdict(list)
    for a, (s, t) in enumerate(arrows): classes[(s, t)].append(a)
    classes = [c for c in classes.values() if len(c) > 1]
    def mats(m):
        out = []
        for vals in itertools.product(range(p), repeat=m * m):
            g = np.array(vals, dtype=np.int64).reshape(m, m)
            if (g.sum(axis=0) % p == 1).all() and F.rank(g) == m: out.append(g)
        return out
    idxb = {b: i for i, b in enumerate(basis)}
    def phi(subs):        # subs: arrow -> {arrow: coeff}; algebra automorphism on paths
        Phi = np.zeros((d, d), dtype=np.int64)
        for i, b in enumerate(basis):
            if b[0] == 'v': Phi[i, i] = 1; continue
            for choice in itertools.product(*[list(subs.get(a, {a: 1}).items()) for a in b[1]]):
                coef = int(np.prod([c for _, c in choice])) % p
                Phi[idxb[('p', tuple(a for a, _ in choice))], i] = (Phi[idxb[('p', tuple(a for a, _ in choice))], i] + coef) % p
        return Phi
    G = []
    for gs in itertools.product(*[mats(len(c)) for c in classes]):
        subs = {}
        for c, g in zip(classes, gs):
            for jcol, a in enumerate(c): subs[a] = {c[irow]: int(g[irow, jcol]) for irow in range(len(c)) if g[irow, jcol]}
        Phi = phi(subs); G.append((Phi, F.inv(Phi)))
    seen = set(); sizes = collections.Counter(); key = lambda A: (A % p).tobytes()
    for P in folds:
        if key(P) in seen: continue
        o = {key(g @ P @ gi) for g, gi in G}; seen |= o; sizes[len(o)] += 1
    assert len({key(P) for P in folds}) == len(folds)
    assert sum(fib.values()) == len(Kfolds), (sum(fib.values()), len(Kfolds))
    res = dict(name=name, p=p, d=d, kdim=kdim, nK=len(Kfolds), N=len(folds), G=len(G), orbits=sum(sizes.values()),
               sizes=dict(sorted(sizes.items())), fib=dict(sorted(fib.items())), t=round(time.time() - t0))
    if verbose: print(res, flush=True)
    return res, folds

if __name__ == '__main__':
    census(sys.argv[1], int(sys.argv[2]))
