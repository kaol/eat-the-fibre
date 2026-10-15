"""Closed-form candidate for the 2-cycle+loop fold family, checked exactly against ALL fold equations on the length-
truncated algebra (products of total length ≤ L, Δ-compatibility on paths of length ≤ L):
   P(p) = Π_{vertices v of p} σ(v) · Π_{arrows a of p} w(a) · e_0,   σ(0)=1, σ(1)=−1,  w(u)=u, w(α)=u−c₁, w(β)=u−c₂."""
import sys
from fractions import Fraction as Fr
from sympy import symbols, expand, Poly
from loops_build import Q, tbuild
def check(name, w_of, sigma, L):
    nv, arrows = Q[name]; d, mult, Dl, allp = tbuild(nv, arrows, L)
    plen = lambda p: 0 if p[0] != 'a' else len(p) - 1
    vseq = lambda p: [p[0]] if p[0] != 'a' else [arrows[p[1]][0]] + [arrows[k][1] for k in p[1:]]
    loopidx = {('a',) + (0,) * k if k else (0,): None for k in range(L + 1)}
    upow = {}                                          # index of u^k e_0 (loop is arrow with (0,0))
    lk = [k for k, a in enumerate(arrows) if a == (0, 0)][0]
    for i, p in enumerate(allp):
        if p == (0,): upow[0] = i
        elif p[0] == 'a' and all(x == lk for x in p[1:]): upow[len(p) - 1] = i
    t = symbols('t')
    P = [[0] * d for _ in range(d)]                    # P[i][j]: coefficient of basis i in P(basis j)
    for j, p in enumerate(allp):
        s = 1
        for v in vseq(p): s *= sigma[v]
        poly = s
        if p[0] == 'a':
            for k in p[1:]: poly = poly * w_of[k](t)
        pl = Poly(expand(poly), t)
        for (e,), c in zip(pl.monoms(), pl.coeffs()):
            if e > L: return None
            P[upow[e]][j] += c
    def apply(vec):
        out = [0] * d
        for j, c in enumerate(vec):
            if c:
                for i in range(d):
                    if P[i][j]: out[i] += P[i][j] * c
        return [expand(x) for x in out]
    def mulv(u, v):
        out = [0] * d
        for i in range(d):
            if u[i] == 0: continue
            for j in range(d):
                if v[j] == 0 or plen(allp[i]) + plen(allp[j]) > L: continue
                for k in range(d):
                    if mult[i, j, k]: out[k] += u[i] * v[j]
        return [expand(x) for x in out]
    col = lambda j: [P[i][j] for i in range(d)]
    ok = True
    for j in range(d): ok &= apply(col(j)) == col(j)                                  # idempotent
    E = lambda j: [1 if i == j else 0 for i in range(d)]
    for i in range(d):
        for j in range(d):
            if plen(allp[i]) + plen(allp[j]) > L: continue
            h = mulv(col(i), col(j)); ok &= apply(h) == h                              # H subalgebra
            q = mulv([a - b for a, b in zip(E(i), col(i))], [a - b for a, b in zip(E(j), col(j))])
            ok &= all(x == 0 for x in apply(q))                                        # Q subalgebra
    for j in range(d):                                                                 # Δ P = (P⊗P) Δ
        lhs = [[0] * d for _ in range(d)]
        for l in range(d):
            if P[l][j]:
                for a in range(d):
                    for b in range(d):
                        if Dl[l, a, b]: lhs[a][b] += P[l][j] * int(Dl[l, a, b])
        rhs = [[0] * d for _ in range(d)]
        for a in range(d):
            for b in range(d):
                if Dl[j, a, b]:
                    for x in range(d):
                        for y in range(d):
                            if P[x][a] and P[y][b]: rhs[x][y] += int(Dl[j, a, b]) * P[x][a] * P[y][b]
        ok &= all(expand(lhs[a][b] - rhs[a][b]) == 0 for a in range(d) for b in range(d))
    return ok
if __name__ == '__main__':
    c1, c2 = symbols('c1 c2')
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    r = check('2-cycle+loop', {0: lambda t: t - c1, 1: lambda t: t - c2, 2: lambda t: t}, {0: 1, 1: -1}, L)
    print(f"2-cycle+loop closed form, symbolic c1, c2, truncation L={L}: all fold equations hold: {r}")
    r2 = check('loop+arrow', {0: lambda t: t, 1: lambda t: t - c1}, {0: 1, 1: -1}, L)
    print(f"loop+arrow closed form α ↦ c1 − u: {r2}")
