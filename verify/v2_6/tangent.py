import sys, time
from fold import *
import sympy as sp
def spec(kind, Lmax):
    if kind == 'ident':   # character splitting, 1 in H: H = k, Q = (u)
        return Splitting([(1, 0)], [(u**k, k) for k in range(1, Lmax + 1)], Lmax)
    if kind == 'framed':  # roles swapped: H = (u), Q = k
        return Splitting([(u**k, k) for k in range(1, Lmax + 1)], [(1, 0)], Lmax)
    if kind == 'dihedral':  # H = k[(u-a)^2], Q = ((u-b)^2 - s) k[(u-b)^2]; a=0, b=1, s=2
        import os; a, b, s = map(sp.Rational, os.environ.get("DIH", "0,1,2").split(",")); w = (u - b)**2
        return Splitting([((u - a)**(2 * k), k) for k in range(Lmax + 1)], [((w - s) * w**k, k + 1) for k in range(Lmax)], Lmax)
def gvec(F, h):
    S = F.S; vec = {}
    for key, j in F.var.items():
        if key[0] != 'c': continue
        _, (x, y), o = key
        val = {}
        if F.isH(x) and F.isH(y):
            for k, m in S.mul(x, h).items():
                for k2, m2 in S.mul(k, y).items(): val[k2] = (val.get(k2, 0) + m * m2) % MOD
        xh = {}
        for k, m in S.mul(x, h).items():
            for k2, m2 in S.mul(k, y).items(): xh[k2] = (xh.get(k2, 0) + m * m2) % MOD
        for k, m in xh.items():
            if F.isH(k): val[k] = (val.get(k, 0) - m) % MOD
        if val.get(o, 0): vec[j] = val[o]
    return vec
def run(kind, Lw, delta, Llow):
    t0 = time.time()
    S = spec(kind, Lw + 2 * delta + 2)
    F = Fold(S, Lw, delta)
    eqs = equations(F, Lw)
    J = linear_part(eqs)
    K, E = kernel(J, len(F.names))
    low = [j for j, key in enumerate(F.names) if key[0] == 'c' and F.wlev(key[1]) <= Llow]
    lowc2 = [j for j, key in enumerate(F.names) if key[0] == 'c2' and F.wlev(key[1]) <= Llow]
    proj = lambda vecs, idx: [{j: v[j] for j in idx if j in v} for v in vecs]
    G = [gvec(F, h) for h in range(len(S.B)) if S.lev[h] <= delta]
    inker = all(not E.reduce({}) and all(sum(r.get(j, 0) * g.get(j, 0) for j in r) % MOD == 0 for r in J if r) for g in G)
    rK = rank(proj(K, low)); rG = rank(proj(G, low)); rK2 = rank(proj(K, lowc2))
    # c2-only tangent directions (zero on all low c coords)
    print(f"{kind:9s} Lw={Lw} δ={delta} low≤{Llow}: vars={len(F.names)} eqs={len(eqs)} tangent(c,low)={rK} G-orbit(c,low)={rG} extra={rK-rG} | G⊂ker:{inker} | {time.time()-t0:.0f}s", flush=True)
    return F, K, G, low
if __name__ == '__main__':
    kind = sys.argv[1]; Lw, delta, Llow = map(int, sys.argv[2:5])
    run(kind, Lw, delta, Llow)
