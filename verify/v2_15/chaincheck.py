"""Chain theorem, intermediate claim: for linearly oriented A_n and every arrow-rule vertex fold, every element of a
diagonal component 1_X Q_m 1_X (m ≥ 1) is supported on paths having an interior vertex v with P₀(e_v) = 0 (a Podge leaf).
(For m = 1 there are no interior vertices, so the claim says the diagonal components vanish.)"""
import sys, itertools
from sage.all__sagemath_singular import QQ, matrix, vector, identity_matrix
from quiver import build
from extend import vertex_folds
n = int(sys.argv[1]); arrows = [(i, i + 1) for i in range(n - 1)]; A = set(arrows)
d, mult, Dl, allp = build(n, arrows)
plen = lambda p: 0 if p[0] != 'a' else len(p) - 1
vseq = lambda p: [p[0]] if p[0] != 'a' else [arrows[p[1]][0]] + [arrows[k][1] for k in p[1:]]
bys = {tuple(vseq(p)): j for j, p in enumerate(allp)}
ok = True; nf = 0; maxdim = 0
for P0 in vertex_folds(n):
    if not all(not [a for a in range(n) if P0[a, s]] or not [b for b in range(n) if P0[b, t]] or
               all((a, b) in A for a in range(n) if P0[a, s] for b in range(n) if P0[b, t]) for s, t in arrows): continue
    nf += 1
    T = matrix(QQ, d, d)
    for j, p in enumerate(allp):
        vs = vseq(p); sup = [[a for a in range(n) if P0[a, v]] for v in vs]
        for seq in itertools.product(*sup):
            c = 1
            for a, v in zip(seq, vs): c *= P0[a, v]
            T[bys[seq], j] += c
    zero = {v for v in range(n) if not any(P0.column(v))}
    Q0 = (identity_matrix(QQ, n) - matrix(QQ, P0)).column_space()
    ind = []
    for S in itertools.chain.from_iterable(itertools.combinations(range(n), r) for r in range(1, n + 1)):
        v = vector(QQ, [1 if i in S else 0 for i in range(n)])
        if v in Q0 and not any(set(S) > {i for i in range(n) if w[i]} for w in ind): ind.append(v)
    iota = sum(ind, vector(QQ, [0] * n)); idem = ind + [vector(QQ, [1] * n) - iota]
    for m in range(1, n):
        Am = [j for j, p in enumerate(allp) if plen(p) == m]
        Qsp = (identity_matrix(QQ, d) - T).matrix_from_columns(Am).column_space()
        for X in idem:
            D = matrix.diagonal(QQ, [X[vseq(p)[0]] * X[vseq(p)[-1]] for p in allp])
            comp = D.column_space().intersection(Qsp); maxdim = max(maxdim, comp.dimension())
            for b in comp.basis():
                for j in range(d):
                    if b[j] and not (set(vseq(allp[j])[1:-1]) & zero): ok = False
print(f"A_{n}: {nf} vertex folds (= a({n})); every diagonal Podge element lies on paths through an interior Podge leaf: {ok}; max diagonal dim {maxdim}")
