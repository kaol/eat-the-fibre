"""Δ-compatible derivations of path algebras with loops, determined on generators (vertices, arrows), filtered
(D(arrow) of length ≤ 1).  A derivation of a path algebra is determined by D(e_v), D(α); Δ-compatibility on generators
implies it everywhere (both ΔD and (D⊗1+1⊗D)Δ are D-twisted ε-derivations).  Orbit dimension at a fold P:
rank of P ↦ [D, P] on the truncated algebra, compared with the fibre dimension."""
import sys, itertools
from sage.all__sagemath_singular import QQ, matrix, vector, PolynomialRing
from loops_build import Q, tbuild
from cfibre import fibre
from extend import vertex_folds
def derivations(nv, arrows, L):
    d, mult, Dl, allp = tbuild(nv, arrows, L)
    plen = lambda p: 0 if p[0] != 'a' else len(p) - 1
    low = [i for i in range(d) if plen(allp[i]) <= 1]        # D(generator) ∈ A_{≤1}
    gens = list(range(nv)) + [nv + k for k in range(len(arrows))]   # vertices then arrows (build order)
    nvar = len(gens) * len(low); S = PolynomialRing(QQ, nvar, 'x'); X = S.gens()
    Dg = {g: [0] * d for g in gens}
    for a, g in enumerate(gens):
        for b, i in enumerate(low): Dg[g][i] = X[a * len(low) + b]
    MT = lambda i, j: [(k, int(mult[i, j, k])) for k in range(d) if mult[i, j, k]]
    def mulv(u, v):
        out = [0] * d
        for i in range(d):
            if u[i] == 0: continue
            for j in range(d):
                if v[j] == 0: continue
                for k, c in MT(i, j): out[k] += c * u[i] * v[j]
        return out
    E = lambda i: [1 if k == i else 0 for k in range(d)]
    eqs = []
    # derivation on generator relations: e_v e_w = δ e_v ; e_s α = α = α e_t ; e_v α = 0 (v≠s), α e_v = 0 (v≠t)
    for v in range(nv):
        for w in range(nv):
            lhs = Dg[v] if v == w else [0] * d
            r = [x + y for x, y in zip(mulv(Dg[v], E(w)), mulv(E(v), Dg[w]))]
            eqs += [l - x for l, x in zip(lhs, r)]
    for k, (s, t) in enumerate(arrows):
        a = nv + k
        for v in range(nv):
            lhs = Dg[a] if v == s else [0] * d
            r = [x + y for x, y in zip(mulv(Dg[v], E(a)), mulv(E(v), Dg[a]))]; eqs += [l - x for l, x in zip(lhs, r)]
            lhs = Dg[a] if v == t else [0] * d
            r = [x + y for x, y in zip(mulv(Dg[a], E(v)), mulv(E(a), Dg[v]))]; eqs += [l - x for l, x in zip(lhs, r)]
    # Δ-compatibility on generators: Δ D g = (D⊗1 + 1⊗D) Δ g
    for g in gens:
        lhs = [[0] * d for _ in range(d)]
        for l in range(d):
            if Dg[g][l] == 0: continue
            for a in range(d):
                for b in range(d):
                    if Dl[l, a, b]: lhs[a][b] += int(Dl[l, a, b]) * Dg[g][l]
        rhs = [[0] * d for _ in range(d)]
        for a in range(d):
            for b in range(d):
                c = int(Dl[g, a, b])
                if not c: continue
                for k in range(d):
                    rhs[k][b] += c * Dg[a][k]; rhs[a][k] += c * Dg[b][k]
        eqs += [lhs[a][b] - rhs[a][b] for a in range(d) for b in range(d)]
    eqs = [e for e in eqs if e != 0]
    A = matrix(QQ, [[S(e).monomial_coefficient(x) for x in X] for e in eqs]) if eqs else matrix(QQ, 0, nvar)
    K = A.right_kernel().basis()
    return d, gens, low, K
def extend_der(nv, arrows, L, gens, low, kvec, d, mult, allp):
    """matrix of the derivation on the truncated algebra (paths via Leibniz)"""
    D = matrix(QQ, d, d)
    for a, g in enumerate(gens):
        for b, i in enumerate(low): D[i, g] = kvec[a * len(low) + b]
    for i, p in enumerate(allp):
        if p[0] != 'a' or len(p) <= 2: continue
        head, tail = ('a', p[1]), ('a',) + p[2:]
        idx = {q: j for j, q in enumerate(allp)}; h, t = idx[head], idx[tail]
        col = vector(QQ, d)
        for k in range(d):                                  # D(h·t) = D(h)·t + h·D(t)
            if D[k, h]:
                for m in range(d):
                    if mult[k, t, m]: col[m] += D[k, h] * int(mult[k, t, m])
            if D[k, t]:
                for m in range(d):
                    if mult[h, k, m]: col[m] += D[k, t] * int(mult[h, k, m])
        D.set_column(i, col)
    return D
if __name__ == '__main__':
    name, L = sys.argv[1], int(sys.argv[2]); nv, arrows = Q[name]; Aset = set(arrows)
    d, gens, low, K = derivations(nv, arrows, L)
    _, mult, Dl, allp = tbuild(nv, arrows, L)
    plen = lambda p: 0 if p[0] != 'a' else len(p) - 1
    Ds = [extend_der(nv, arrows, L, gens, low, k, d, mult, allp) for k in K]
    print(f"{name}: dim of filtered Δ-compatible derivations (on generators) = {len(K)}")
    for P0 in vertex_folds(nv):
        if not all(not [a for a in range(nv) if P0[a, s]] or not [b for b in range(nv) if P0[b, t]] or
                   all((a, b) in Aset for a in range(nv) if P0[a, s] for b in range(nv) if P0[b, t]) for s, t in arrows): continue
        R, J, n = fibre(nv, arrows, P0, L - 1)
        fdim = J.dimension()
        # canonical fold T on the truncated algebra
        import itertools as it
        T = matrix(QQ, d, d)
        vseq = lambda p: [p[0]] if p[0] != 'a' else [arrows[p[1]][0]] + [arrows[k][1] for k in p[1:]]
        bys = {tuple(vseq(p)): j for j, p in enumerate(allp)}
        for j, p in enumerate(allp):
            vs = vseq(p); sup = [[a for a in range(nv) if P0[a, v]] for v in vs]
            for seq in it.product(*sup):
                c = 1
                for a, v in zip(seq, vs): c *= P0[a, v]
                T[bys[seq], j] += c
        keep = [i for i in range(d) if plen(allp[i]) <= L - 1]
        orb = matrix(QQ, [list((Dm * T - T * Dm).matrix_from_rows_and_columns(keep, keep).list()) for Dm in Ds]).rank() if Ds else 0
        print(f"   P0={list(P0)}: fibre dim {fdim}, orbit dim of the canonical fold {orb}")
