import sys, itertools
from sage.all__sagemath_singular import QQ, matrix, vector, identity_matrix
from quiver import build
from extend import vertex_folds
from arrowrule import QS
name = sys.argv[1]; nv, arrows = QS[name]; A = set(arrows)
d, mult, Dl, allp = build(nv, arrows)
plen = lambda p: 0 if p[0] != 'a' else len(p) - 1
vseq = lambda p: [p[0]] if p[0] != 'a' else [arrows[p[1]][0]] + [arrows[k][1] for k in p[1:]]
shown = 0
for P0 in vertex_folds(nv):
    if not all(not [a for a in range(nv) if P0[a, s]] or not [b for b in range(nv) if P0[b, t]] or
               all((a, b) in A for a in range(nv) if P0[a, s] for b in range(nv) if P0[b, t]) for s, t in arrows): continue
    bys = {tuple(vseq(p)): j for j, p in enumerate(allp)}
    T = matrix(QQ, d, d)
    for j, p in enumerate(allp):
        vs = vseq(p); sup = [[a for a in range(nv) if P0[a, v]] for v in vs]
        for seq in itertools.product(*sup):
            c = 1
            for a, v in zip(seq, vs): c *= P0[a, v]
            T[bys[seq], j] += c
    Q0 = (identity_matrix(QQ, nv) - matrix(QQ, P0)).column_space()
    ind = []
    for S in itertools.chain.from_iterable(itertools.combinations(range(nv), r) for r in range(1, nv + 1)):
        v = vector(QQ, [1 if i in S else 0 for i in range(nv)])
        if v in Q0 and not any(set(S) > {i for i in range(nv) if w[i]} for w in ind): ind.append(v)
    iota = sum(ind, vector(QQ, [0] * nv)); idem = ind + [vector(QQ, [1] * nv) - iota]
    for m in (2, 3):
        Am = [j for j, p in enumerate(allp) if plen(p) == m]
        if not Am: continue
        Qsp = (identity_matrix(QQ, d) - T).matrix_from_columns(Am).column_space()
        for X in idem:
            D = matrix.diagonal(QQ, [X[vseq(p)[0]] * X[vseq(p)[-1]] for p in allp])
            comp = D.column_space().intersection(Qsp)
            for b in comp.basis():
                if shown < 6:
                    print(f"P0 cols={[list(P0.column(v)) for v in range(nv)]}  F2∪∅ block {list(X)}  m={m}:",
                          {tuple(vseq(allp[j])): b[j] for j in range(d) if b[j]}, " T(paths):",
                          {tuple(vseq(allp[j])): {tuple(vseq(allp[i])): T[i, j] for i in range(d) if T[i, j]} for j in range(d) if b[j]})
                    shown += 1
