"""Exact fold enumeration on path algebras WITHOUT parallel arrows, via the vertex algebra.
Lemma (proved): with no parallel arrows, ker Δ = A_0 (vertex span); since ΔP = (P⊗P)Δ and Δe_v = 0, P(A_0) ⊂ A_0.
So P|A_0 is a fold of k^n = A_0 (a Cayley forest), and for a path x of length ≥ 1, ΔP(x) = (P⊗P)Δx involves only
shorter paths: P(x) = L((P⊗P)Δx) + a_x with a_x ∈ A_0 unknown (L = left inverse of the injective Δ on A_≥1).
All remaining fold equations are polynomial in the a's; each vertex splitting's fibre is solved exactly (Singular)."""
import sys, itertools
from fractions import Fraction
from sage.all__sagemath_singular import PolynomialRing, QQ, matrix, vector, ideal
from quiver import build
from collections import Counter
def forests(n):
    for par in itertools.product(list(range(n)) + ['R1', 'R2'], repeat=n):
        ok = True
        for v in range(n):
            seen = set(); x = v
            while isinstance(x, int):
                if x in seen: ok = False; break
                seen.add(x); x = par[x]
            if not ok: break
        if ok: yield par
def splitting(par):
    n = len(par); F = {1: [], 2: []}
    for w in range(n):
        d = 0; x = w
        while isinstance(x, int): d += 1; x = par[x]
        rs = 1 if x == 'R1' else 2
        F[rs if d % 2 == 1 else 3 - rs].append(frozenset([w] + [v for v in range(n) if par[v] == w]))
    return F[1], F[2]
def vertex_folds(n):
    seen = set()
    for par in forests(n):
        F1, F2 = splitting(par)
        B = matrix(QQ, [[1 if v in s else 0 for s in F1 + F2] for v in range(n)])
        P0 = B * matrix.diagonal(QQ, [1] * len(F1) + [0] * len(F2)) * B.inverse()
        key = tuple(P0.list())
        if key not in seen: seen.add(key); yield P0
def fibres(nv, arrows, verbose=False):
    assert len(set(arrows)) == len(arrows), "parallel arrows: ker Δ ≠ A_0"
    d, mult, Dl, allp = build(nv, arrows)
    n = nv; paths = list(range(n, d))
    Dm = matrix(QQ, d * d, d, lambda ab, l: int(Dl[l, ab // d, ab % d]))
    Dp = Dm.matrix_from_columns(paths); Lp = (Dp.transpose() * Dp).inverse() * Dp.transpose()
    names = [f'a{x}_{v}' for x in paths for v in range(n)] or ['dummy']
    R = PolynomialRing(QQ, len(names), names, order='degrevlex'); gens = dict(zip(names, R.gens()))
    MT = {(i, j): [(k, int(mult[i, j, k])) for k in range(d) if mult[i, j, k]] for i in range(d) for j in range(d)}
    def mul(u, v):
        out = [R(0)] * d
        for i in range(d):
            if u[i] == 0: continue
            for j in range(d):
                if v[j] == 0: continue
                for k, c in MT[(i, j)]: out[k] += c * u[i] * v[j]
        return vector(R, out)
    total = 0; stats = Counter(); ranks = Counter(); points = []
    for P0 in vertex_folds(n):
        P = matrix(R, d, d); eqs = []
        for i in range(n):
            for j in range(n): P[i, j] = P0[i, j]
        for x in paths:                                     # increasing length (build order)
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
            eqs += list(Dp * y - Tv)                         # (P⊗P)Δx must lie in Δ(A_≥1)
            for k, p in enumerate(paths): P[p, x] = y[k]
            for v in range(n): P[v, x] = gens[f'a{x}_{v}']
        eqs += (P * P - P).list()
        I_ = matrix.identity(R, d)
        cols = [P.column(i) for i in range(d)]; qcols = [I_.column(i) - cols[i] for i in range(d)]
        for i in range(d):
            for j in range(d):
                h = mul(cols[i], cols[j]); eqs += list(P * h - h)
                eqs += list(P * mul(qcols[i], qcols[j]))
        nz = [e for e in eqs if e != 0]
        if not paths:
            if not nz: total += 1; stats[1] += 1; stats['the point is a = 0 (canonical tensor extension)'] += 1
            continue
        J = R.ideal(nz) if nz else R.ideal([R(0)])
        dim = J.dimension()
        if dim < 0: continue
        if dim > 0: stats['positive-dim fibre'] += 1; continue
        cnt = J.radical().vector_space_dimension()
        total += cnt; stats[cnt] += 1
        if cnt == 1 and J.radical() == R.ideal(list(R.gens())): stats['the point is a = 0 (canonical tensor extension)'] += 1
    return d, total, stats
if __name__ == '__main__':
    QS = {'A2': (2, [(0, 1)]), 'A3': (3, [(0, 1), (1, 2)]), '1→2←3': (3, [(0, 1), (2, 1)]), '1←2→3': (3, [(1, 0), (1, 2)]),
          'A4': (4, [(0, 1), (1, 2), (2, 3)]), 'A4 alt': (4, [(0, 1), (2, 1), (2, 3)]), 'D4': (4, [(0, 3), (1, 3), (2, 3)]),
          '1→2→3, 1→3': (3, [(0, 1), (1, 2), (0, 2)]), 'diamond': (4, [(0, 1), (0, 2), (1, 3), (2, 3)]),
          'A5': (5, [(0, 1), (1, 2), (2, 3), (3, 4)]), 'A1': (1, []), 'A4 →→←': (4, [(0, 1), (1, 2), (3, 2)]), 'A4 ←←→': (4, [(1, 0), (2, 1), (2, 3)]), 'A4 →←←': (4, [(0, 1), (2, 1), (3, 2)])}
    for name in sys.argv[1:]:
        d, total, stats = fibres(*QS[name])
        print(f"{name}: dim A = {d}; folds = {total}; fibre sizes over vertex splittings: {dict(stats)}", flush=True)
