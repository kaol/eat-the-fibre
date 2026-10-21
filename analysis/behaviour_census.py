"""v2.33: the behaviour census. Folds up to conjugation by all algebra automorphisms of kQ.

Behaviour of a fold P := its Aut(kQ)-orbit (equivalently, the splitting (im P, ker P) up to ambient isomorphism).
For acyclic Q without a path of length ≥ 2 parallel to an arrow, Aut(kQ) = Inn ⋊ (Sub ⋊ Aut(Q)), where
- Inn  = conjugations by u = 1 + r, r ∈ rad kQ (normal in Aut(kQ));
- Sub  = arrow substitutions within each parallel class, any invertible matrix (single arrows: scalars);
- Aut(Q) = quiver automorphisms (vertex permutations preserving arrow multiplicities), reported separately.
The census group G is the column-sum-1 part of Sub (exactly the Δ-compatible substitutions). Right cosets G·s of G in
Sub are labelled by the row vector 1ᵀs, so with S a set of right-coset representatives,
    Aut⁰(kQ) · P ∩ Folds = ⋃_{s ∈ S, u ∈ Inn} G · (Inn_u s P s⁻¹) ∩ Folds,
because g Inn_u s = Inn_{g(u)} g s and Inn ⊴ Aut. So behaviours are the components of: G-orbits, joined whenever
Inn_u s P s⁻¹ (P a G-orbit member, s ∈ S, u ∈ Inn) is a fold.

Usage: python3 analysis/behaviour_census.py QUIVER p FOLDS.pkl [CLASSES.jsonl]
       (FOLDS from `pcensus.py QUIVER p --checkpoint COPY --save-folds FOLDS.pkl`)
"""
import os, sys, pickle, itertools, collections, json
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_argv = sys.argv; sys.argv = ['x']; _cwd = os.getcwd(); os.chdir(ROOT); __file__ = os.path.join(ROOT, 'pcensus.py')
exec(open(__file__).read().split("if __name__ == '__main" + "__':")[0])
sys.argv = _argv; os.chdir(_cwd)


class DSU:
    def __init__(s, n): s.p = np.arange(n)
    def find(s, x):
        r = x
        while s.p[r] != r: r = s.p[r]
        while s.p[x] != r: s.p[x], x = r, s.p[x]
        return r
    def union(s, a, b):
        a, b = s.find(a), s.find(b)
        if a != b: s.p[max(a, b)] = min(a, b)


def main():
    name, p, ffile = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    nv, arrows = QUIVERS[name]
    basis, M, DL = build(nv, arrows)
    M = np.array(M, dtype=np.int64); d = len(basis)
    idxb = {b: i for i, b in enumerate(basis)}
    folds = np.array(pickle.load(open(ffile, 'rb'))['folds'], dtype=np.int64) % p
    N = len(folds)
    # ---- hashing folds
    rng = np.random.default_rng(1); W = rng.integers(1, 2**62, size=d * d, dtype=np.int64)
    hsh = lambda X: (X.reshape(len(X), -1) * W).sum(1)                    # wraps mod 2^64
    H = hsh(folds); order = np.argsort(H); Hs = H[order]
    assert len(np.unique(Hs)) == N, "hash collision among folds"
    def lookup(X):
        h = hsh(X); pos = np.searchsorted(Hs, h).clip(0, N - 1)
        idx = order[pos]; hit = (Hs[pos] == h) & (folds[idx] == X).all((1, 2))   # exact confirmation
        return np.where(hit, idx, -1)
    # ---- substitution matrices
    classes = collections.defaultdict(list)
    for ai, st in enumerate(arrows): classes[st].append(ai)
    cl = list(classes.values())
    def phi(Ms, perm=None):
        """vertex-fixing (or vertex-permuting, perm = (σ_v, σ_a)) automorphism: arrow x ↦ Σ_y Ms[class][y, x] y"""
        sv, sa = perm if perm else (list(range(nv)), list(range(len(arrows))))
        img_ar = {}
        for c, Mc in zip(cl, Ms):
            for jx, x in enumerate(c):
                img_ar[x] = {sa[c[iy]]: int(Mc[iy, jx]) % p for iy in range(len(c)) if Mc[iy, jx] % p}
        Phi = np.zeros((d, d), dtype=np.int64)
        for i, b in enumerate(basis):
            if b[0] == 'v': Phi[idxb[('v', sv[b[1]])], i] = 1; continue
            for choice in itertools.product(*[list(img_ar[x].items()) for x in b[1]]):
                j = idxb[('p', tuple(x for x, _ in choice))]
                Phi[j, i] = (Phi[j, i] + int(np.prod([v for _, v in choice]))) % p
        return Phi
    def inv(A):
        import sympy as sp
        return np.array(sp.Matrix(A.tolist()).inv_mod(p).tolist(), dtype=np.int64)
    def gl(m):
        out = []
        for vals in itertools.product(range(p), repeat=m * m):
            A = np.array(vals, dtype=np.int64).reshape(m, m)
            if round(np.linalg.det(A)) % p: out.append(A)
        return out
    GL = [gl(len(c)) for c in cl]
    Gc = [[A for A in g if ((A.sum(0) - 1) % p == 0).all()] for g in GL]  # column sums 1
    cos = []                                                               # right-coset reps: one per 1ᵀs
    for g in GL:
        seen = {}
        for A in g: seen.setdefault(tuple(A.sum(0) % p), A)
        cos.append(list(seen.values()))
    # quiver automorphisms (vertex permutations with matching arrow multiplicities; arrows mapped in class order)
    qaut = []
    for sv in itertools.permutations(range(nv)):
        ok = all(len(classes.get((sv[s], sv[t]), [])) == len(c) for (s, t), c in classes.items())
        if not ok: continue
        sa = [0] * len(arrows)
        for (s, t), c in classes.items():
            for x, y in zip(c, classes[(sv[s], sv[t])]): sa[x] = y
        qaut.append((list(sv), sa))
    # ---- G-orbits by union over generators (random elements generating each class group)
    dsu = DSU(N)
    gens = []
    for k, c in enumerate(cl):
        if len(Gc[k]) <= 1: continue
        for A in [Gc[k][i] for i in rng.choice(len(Gc[k]), size=min(6, len(Gc[k])), replace=False)]:
            Ms = [np.eye(len(cc), dtype=np.int64) for cc in cl]; Ms[k] = A; gens.append(phi(Ms))
    for g in gens:
        gi = inv(g); img = lookup(np.einsum('ij,njk,kl->nil', g, folds, gi) % p)
        assert (img >= 0).all(), "G does not preserve folds"
        for a, b in zip(range(N), img): dsu.union(a, b)
    Gorb = np.array([dsu.find(i) for i in range(N)])
    nG = len(set(Gorb))
    # ---- inner automorphisms
    vert = [i for i, b in enumerate(basis) if b[0] == 'v']; rad = [i for i, b in enumerate(basis) if b[0] == 'p']
    one = np.zeros(d, np.int64); one[vert] = 1
    mulv = lambda x, y: np.einsum('...i,...j,ijm->...m', x, y, M) % p
    R = np.array(list(itertools.product(range(p), repeat=len(rad))), dtype=np.int64)
    r = np.zeros((len(R), d), np.int64); r[:, rad] = R
    u = (one + r) % p
    pw = np.broadcast_to(one, r.shape).copy(); uinv = pw.copy()
    for _ in range(d): pw = mulv(pw, (-r) % p); uinv = (uinv + pw) % p        # r is nilpotent
    assert (mulv(u, uinv) == one).all()
    def innmat(a, b):                                                      # x ↦ a x b
        La = np.einsum('ni,ijm->nmj', a, M) % p; Rb = np.einsum('nj,ijm->nmi', b, M) % p
        return np.einsum('nmk,nkj->nmj', Rb, La) % p
    Inn, Innv = innmat(u, uinv), innmat(uinv, u)
    # ---- behaviour edges from one member per G-orbit
    reps = {}
    for i, g in enumerate(Gorb): reps.setdefault(g, i)
    def components(cosets):
        dsu2 = DSU(N)
        for i in range(N): dsu2.p[i] = Gorb[i]
        S = [(phi(list(Ms), perm), None) for perm in cosets for Ms in itertools.product(*cos)]
        S = [(s, inv(s)) for s, _ in S]
        Sm = np.array([a for a, _ in S], dtype=np.float64); Si = np.array([b for _, b in S], dtype=np.float64)
        InnF, InnvF = Inn.astype(np.float64), Innv.astype(np.float64)
        CH = max(1, (1 << 27) // (len(Inn) * d * d * 8))                    # coset reps per batch (~128 MiB)
        for gi, i in reps.items():
            P = folds[i].astype(np.float64)
            for c0 in range(0, len(S), CH):
                X = np.mod(Sm[c0:c0 + CH] @ P @ Si[c0:c0 + CH], p)                       # (s, d, d)
                Q = np.mod(InnF[None] @ X[:, None] @ InnvF[None], p).astype(np.int64)    # exact: entries < 2^53
                hit = lookup(Q.reshape(-1, d, d)); hit = hit[hit >= 0]
                for j in np.unique(hit): dsu2.union(gi, Gorb[j])
        comp = np.array([dsu2.find(Gorb[i]) for i in range(N)])
        return comp
    comp0 = components([None])
    compQ = components(qaut) if len(qaut) > 1 else comp0
    res = dict(name=name, p=p, N=N, orbits_G=nG, behaviours_Aut0=len(set(comp0)), behaviours_Aut=len(set(compQ)),
               quiver_automorphisms=len(qaut), inner=len(R), coset_reps=int(np.prod([len(c) for c in cos])))
    # fibre widths: G-orbits per behaviour, as a histogram
    gper = collections.Counter(); seen = set()
    for i in range(N):
        if (comp0[i], Gorb[i]) not in seen: seen.add((comp0[i], Gorb[i])); gper[comp0[i]] += 1
    res['G_orbits_per_behaviour'] = dict(sorted(collections.Counter(gper.values()).items()))
    out = os.environ.get('BEH_OUT')
    if out: pickle.dump(dict(res=res, Gorb=Gorb, comp0=comp0, compQ=compQ), open(out, 'wb'))
    print(json.dumps(res), flush=True)


if __name__ == '__main__':
    main()
