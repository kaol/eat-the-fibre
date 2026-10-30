"""v2.76: which vertex parts c can the path family carry? As v2_76_pathfamily, but N(u·a0·π_s) = u·(c_s + m_s ℓ₀) with
c_s = Σ_v x_{s,v} e_v fully symbolic (one vector per tail path); prints a Gröbner basis of the fold conditions in the
x's and m's (the empty ideal means every choice is a fold).
Usage: python3 analysis/v2_76_vertexpart.py QUIVER A0,A1,A2 PI1;PI2 [OUT.json]"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE)); sys.path.insert(0, HERE)
import numpy as np, sympy as sp
from sympy import QQ
from sympy.polys.rings import ring
import pcensus as pc, fold_ideal as fi
quiver = sys.argv[1]; cls = tuple(map(int, sys.argv[2].split(',')))
paths = [tuple(map(int, p.split('.'))) for p in sys.argv[3].split(';')]; out = sys.argv[4] if len(sys.argv) > 4 else None
ctx = pc.setup(quiver, 7); d = ctx['d']; basis = ctx['basis']; idx = {b: k for k, b in enumerate(basis)}
arrows = ctx['arrows']; a0, a1, a2 = cls; i, j = arrows[a0]
nv = sum(1 for b in basis if b[0] == 'v')
names = [f'x{s}_{v}' for s in range(len(paths)) for v in range(nv)] + [f'm{s}' for s in range(len(paths))]
R, *gens = ring(','.join(names), QQ); G = dict(zip(names, gens))
M = np.asarray(ctx['M'], dtype=np.int64); DL = np.asarray(ctx['DL'], dtype=np.int64)
nzM = [(a, b, k, int(M[a, b, k])) for a in range(d) for b in range(d) for k in range(d) if M[a, b, k]]
def mul(x, y):
    o = [R(0)] * d
    for a, b, k, c in nzM:
        if x[a] and y[b]: o[k] += c * x[a] * y[b]
    return o
unit = lambda k: [R(1) if r == k else R(0) for r in range(d)]
P = [[R(1) if r == c else R(0) for c in range(d)] for r in range(d)]
for col, bb in enumerate(basis):
    if bb[0] != 'p' or a0 not in bb[1]: continue
    k = bb[1].index(a0); u, r = bb[1][:k], bb[1][k + 1:]
    v = unit(idx[('p', u + (a1,) + r)])
    if r in paths:
        s = paths.index(r); corr = [R(0)] * d
        for w in range(nv): corr[idx[('v', w)]] += G[f'x{s}_{w}']
        corr[idx[('p', (a1,))]] += G[f'm{s}']; corr[idx[('p', (a2,))]] -= G[f'm{s}']
        if u: corr = mul(unit(idx[('p', u)]), corr)
        v = [v[q] + corr[q] for q in range(d)]
    for rr in range(d): P[rr][col] = v[rr]
eqs = fi._equations(P, R, d, nzM, DL)
syms = [sp.Symbol(n) for n in names]
gb = sp.groebner([fi.to_expr(e) for e in eqs], *syms, order='grevlex').exprs if eqs else []
res = dict(quiver=quiver, paths=paths, i=i, j=j, nvert=nv, groebner=[str(g) for g in gb])
print(json.dumps(res))
if out: json.dump(res, open(out, 'w'), indent=1)
