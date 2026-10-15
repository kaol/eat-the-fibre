"""Podge half, structure.  ψ: Q_m → H₀ is a Q₀-bimodule map (q⋆h = π_H(qh)).  Decompose Q_m = ⊕ 1_X Q_m 1_Y over
X, Y ∈ F2 ∪ {∅} (1_∅ = 1 − ι).  Mixed components are killed by orthogonality.  Test, for every arrow-rule vertex
fold and every length:  (D1) the diagonal components 1_C Q_m 1_C and (1−ι)Q_m(1−ι) vanish for m = 1;
(D2) for m ≥ 2 they lie in the span of glued products Q_i·Q_{m−i} (i, m−i ≥ 1)."""
import sys, itertools
from sage.all__sagemath_singular import QQ, matrix, vector, identity_matrix
from quiver import build
from extend import vertex_folds
from arrowrule import QS
def analyse(nv, arrows, P0):
    d, mult, Dl, allp = build(nv, arrows)
    plen = lambda p: 0 if p[0] != 'a' else len(p) - 1
    vseq = lambda p: [p[0]] if p[0] != 'a' else [arrows[p[1]][0]] + [arrows[k][1] for k in p[1:]]
    bys = {tuple(vseq(p)): j for j, p in enumerate(allp)}
    T = matrix(QQ, d, d)
    for j, p in enumerate(allp):
        vs = vseq(p); sup = [[a for a in range(nv) if P0[a, v]] for v in vs]
        for seq in itertools.product(*sup):
            c = 1
            for a, v in zip(seq, vs): c *= P0[a, v]
            T[bys[seq], j] += c
    I = identity_matrix(QQ, d); Qm = I - T
    Pm = matrix(QQ, P0); Q0 = (identity_matrix(QQ, nv) - Pm).column_space()
    # F2 block indicators: Q0 is spanned by indicators of disjoint blocks; recover them as the primitive idempotents
    basisQ0 = [vector(QQ, b) for b in Q0.basis()]
    # primitive idempotents of the subalgebra Q0 of k^n: indicator vectors in Q0 with minimal support
    ind = []
    for S in itertools.chain.from_iterable(itertools.combinations(range(nv), r) for r in range(1, nv + 1)):
        v = vector(QQ, [1 if i in S else 0 for i in range(nv)])
        if v in Q0 and not any(set(S) > set(T_) for T_ in [tuple(i for i in range(nv) if w[i]) for w in ind]): ind.append(v)
    iota = sum(ind, vector(QQ, [0] * nv))
    idem = ind + [vector(QQ, [1] * nv) - iota]
    def Lmul(f):   # left multiplication by a vertex function on A
        M = matrix(QQ, d, d)
        for j, p in enumerate(allp): M[j, j] = f[vseq(p)[0]]
        return M
    def Rmul(f):
        M = matrix(QQ, d, d)
        for j, p in enumerate(allp): M[j, j] = f[vseq(p)[-1]]
        return M
    out = {}
    for m in range(1, max(plen(p) for p in allp) + 1):
        Am = [j for j, p in enumerate(allp) if plen(p) == m]
        Qsp = Qm.matrix_from_columns(Am).column_space()
        diag = []
        for X in idem:
            comp = (Lmul(X) * Rmul(X)).column_space().intersection(Qsp)        # 1_X Q_m 1_X (Q_m is an A0-bimodule)
            diag.append(comp)
        tot = sum(c.dimension() for c in diag)
        if m == 1: out[m] = ('D1', tot == 0, tot)
        else:
            prods = []
            for i in range(1, m):
                Ai = [j for j, p in enumerate(allp) if plen(p) == i]; Aj = [j for j, p in enumerate(allp) if plen(p) == m - i]
                for a in Ai:
                    ua = Qm.column(a)
                    for b in Aj:
                        ub = Qm.column(b)
                        w = vector(QQ, d)
                        for x in range(d):
                            if ua[x] == 0: continue
                            for y in range(d):
                                if ub[y] == 0: continue
                                for z in range(d):
                                    if mult[x, y, z]: w[z] += ua[x] * ub[y]
                        prods.append(w)
            span = matrix(QQ, prods).row_space() if prods else None
            ok = all(span is not None and c.is_subspace(span) for c in diag if c.dimension() > 0)
            out[m] = ('D2', ok, tot)
    return out
if __name__ == '__main__':
    for name in sys.argv[1].split(';'):
        nv, arrows = QS[name]; A = set(arrows); allok = True; cnt = 0; dims = set()
        for P0 in vertex_folds(nv):
            if not all(not [a for a in range(nv) if P0[a, s]] or not [b for b in range(nv) if P0[b, t]] or
                       all((a, b) in A for a in range(nv) if P0[a, s] for b in range(nv) if P0[b, t]) for s, t in arrows): continue
            cnt += 1
            for m, (kind, ok, tot) in analyse(nv, arrows, P0).items():
                allok &= ok; dims.add((m, tot))
        print(f"{name}: {cnt} vertex folds; D1 (m=1 diagonal components vanish) and D2 (m≥2 diagonal ⊂ Q·Q): {allok}; (m, total diagonal dim) seen: {sorted(dims)}", flush=True)
