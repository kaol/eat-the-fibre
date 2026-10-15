"""Exact fibre over a vertex fold for quivers WITH cycles (no parallel arrows), truncated honestly:
unknowns a_x (x of length ≤ L); equations = Δ-existence for all paths of length ≤ L+1, idempotency and both
subalgebra conditions for every product of total length ≤ L.  P is length-filtered, so these are exactly the fold
equations that involve only a_{≤L}: every genuine fold restricts to a point of this fibre."""
import sys, itertools
from sage.all__sagemath_singular import PolynomialRing, QQ, matrix, vector
import loops_build
from extend import vertex_folds
def fibre(nv, arrows, P0, L):
    d, mult, Dl, allp = loops_build.tbuild(nv, arrows, L + 1)
    plen = lambda p: 0 if p[0] != 'a' else len(p) - 1
    n = nv; paths = [i for i in range(n, d)]
    Dm = matrix(QQ, d * d, d, lambda ab, l: int(Dl[l, ab // d, ab % d]))
    Dp = Dm.matrix_from_columns(paths); Lp = (Dp.transpose() * Dp).inverse() * Dp.transpose()
    unk = [x for x in paths if plen(allp[x]) <= L]
    names = [f'a{x}_{v}' for x in unk for v in range(n)]
    R = PolynomialRing(QQ, len(names), names, order='degrevlex'); g = dict(zip(names, R.gens()))
    MT = {(i, j): [(k, int(mult[i, j, k])) for k in range(d) if mult[i, j, k]] for i in range(d) for j in range(d)}
    def mul(u, v):
        out = [R(0)] * d
        for i in range(d):
            if u[i] == 0: continue
            for j in range(d):
                if v[j] == 0: continue
                for k, c in MT[(i, j)]: out[k] += c * u[i] * v[j]
        return vector(R, out)
    P = matrix(R, d, d); eqs = []
    for i in range(n):
        for j in range(n): P[i, j] = P0[i, j]
    for x in paths:
        T = [R(0)] * (d * d)
        for a in range(d):
            for b in range(d):
                c = int(Dl[x, a, b])
                if c:
                    for i in range(d):
                        if P[i, a] == 0: continue
                        for j in range(d):
                            if P[j, b] != 0: T[i * d + j] += c * P[i, a] * P[j, b]
        Tv = vector(R, T); y = Lp * Tv
        eqs += [e for e in (Dp * y - Tv) if e != 0]
        for k, p in enumerate(paths): P[p, x] = y[k]
        if plen(allp[x]) <= L:
            for v in range(n): P[v, x] = g[f'a{x}_{v}']
    low = [i for i in range(d) if plen(allp[i]) <= L]
    I_ = matrix.identity(R, d)
    for i in low:
        col = P.column(i); eqs += [e for e in (P * col - col) if e != 0]
    for i in low:
        for j in low:
            if plen(allp[i]) + plen(allp[j]) > L: continue
            h = mul(P.column(i), P.column(j)); eqs += [e for e in (P * h - h) if e != 0]
            q = mul(I_.column(i) - P.column(i), I_.column(j) - P.column(j)); eqs += [e for e in P * q if e != 0]
    J = R.ideal(eqs) if eqs else R.ideal([R(0)])
    return R, J, len(names)
if __name__ == '__main__':
    from loops_build import Q
    name, L = sys.argv[1], int(sys.argv[2]); nv, arrows = Q[name]; A = set(arrows)
    for P0 in vertex_folds(nv):
        if not all(not [a for a in range(nv) if P0[a, s]] or not [b for b in range(nv) if P0[b, t]] or
                   all((a, b) in A for a in range(nv) if P0[a, s] for b in range(nv) if P0[b, t]) for s, t in arrows): continue
        R, J, nvars = fibre(nv, arrows, P0, L)
        dim = J.dimension()
        extra = ''
        if dim == 0: extra = f", points = {J.radical().vector_space_dimension()}"
        print(f"{name} L={L} P0={list(P0)}: {nvars} unknowns, fibre dimension {dim}{extra}", flush=True)
