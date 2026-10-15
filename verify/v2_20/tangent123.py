"""Lift the F_3 folds of k(1⇉2→3) to Q (2 ↦ −1), keep those that are exact folds over Q, and compare at each the dimension of
the Zariski tangent space of the fold variety (linearized equations, exact) with the dimension of the orbit of
G = {column-sum-1 substitutions of α1, α2} (Lie algebra: D(α1) = aδ, D(α2) = bδ, δ = α1 − α2, extended to α_iβ)."""
import sys, numpy as np
sys.argv = ['x', '3']
src = open('q123v.py').read().replace('exec("G = []" + open(\'q123.py\').read().split("G = []")[1])', '')
exec(src)
from sage.all__sagemath_singular import QQ, matrix, vector
Mi = M.astype(int); DLi = DL.astype(int)
def qmul(x, y):
    out = [QQ(0)] * d
    for i in range(d):
        if x[i] == 0: continue
        for j in range(d):
            if y[j] == 0: continue
            for k in range(d):
                if Mi[i, j, k]: out[k] += x[i] * y[j]
    return vector(QQ, out)
def is_fold(P):
    Id = matrix.identity(QQ, d)
    if P * P != P: return False
    for i in range(d):
        for j in range(d):
            h = qmul(P.column(i), P.column(j))
            if P * h != h: return False
            q = qmul(Id.column(i) - P.column(i), Id.column(j) - P.column(j))
            if P * q != 0: return False
    for j in range(d):
        lhs = matrix(QQ, d, d); 
        for l in range(d):
            if P[l, j]:
                for a in range(d):
                    for b in range(d):
                        if DLi[l, a, b]: lhs[a, b] += P[l, j] * DLi[l, a, b]
        T = matrix(QQ, d, d, lambda a, b: DLi[j, a, b])
        if lhs != P * T * P.transpose(): return False
    return True
def residual_jacobian(P):
    """Jacobian of the fold equations at P w.r.t. the 64 entries (exact, by bilinearity: forward difference of degree-2 maps)."""
    from sage.all__sagemath_singular import PolynomialRing
    R = PolynomialRing(QQ, d * d, 'x'); X = matrix(R, d, d, R.gens())
    Pp = P.change_ring(R) + X
    Id = matrix.identity(R, d); eqs = (Pp * Pp - Pp).list()
    def rmul(x, y):
        out = [R(0)] * d
        for i in range(d):
            for j in range(d):
                for k in range(d):
                    if Mi[i, j, k]: out[k] += x[i] * y[j]
        return vector(R, out)
    for i in range(d):
        for j in range(d):
            h = rmul(Pp.column(i), Pp.column(j)); eqs += list(Pp * h - h)
            q = rmul(Id.column(i) - Pp.column(i), Id.column(j) - Pp.column(j)); eqs += list(Pp * q)
    for j in range(d):
        lhs = matrix(R, d, d)
        for l in range(d):
            for a in range(d):
                for b in range(d):
                    if DLi[l, a, b]: lhs[a, b] += Pp[l, j] * DLi[l, a, b]
        T = matrix(R, d, d, lambda a, b: DLi[j, a, b]); eqs += (lhs - Pp * T * Pp.transpose()).list()
    lin = [[e.monomial_coefficient(g) if e != 0 else 0 for g in R.gens()] for e in eqs]
    return matrix(QQ, lin)
def derivs():
    out = []
    for a, b in ((1, 0), (0, 1)):
        Dm = matrix(QQ, d, d)
        Dm[3, 3] += a; Dm[4, 3] -= a; Dm[3, 4] += b; Dm[4, 4] -= b      # D(α1) = aδ, D(α2) = bδ
        Dm[6, 6] += a; Dm[7, 6] -= a; Dm[6, 7] += b; Dm[7, 7] -= b      # on α1β, α2β
        out.append(Dm)
    return out
Ds = derivs()
lifted = 0; rows = []
import os
for f in folds[:int(os.environ.get("LIM", "100000"))]:
    P = matrix(QQ, [[(-1 if x == 2 else x) for x in row] for row in f.tolist()])
    if not is_fold(P): continue
    lifted += 1
    J = residual_jacobian(P); tdim = d * d - J.rank()
    orb = matrix(QQ, [(Dm * P - P * Dm).list() for Dm in Ds]).rank()
    rows.append((tdim, orb))
from collections import Counter
print(f"{len(folds)} F_3 folds; {lifted} lift to exact folds over Q; (tangent dim, orbit dim) distribution: {sorted(Counter(rows).items())}")
print("tangent dim = orbit dim at every lifted fold:", all(t == o for t, o in rows))
