"""v2.79: infinitesimal moduli of the path family with a general vertex part c = e_X (the same X for every path):
N(u·a0·π_s) = u·(λ_s e_X + m_s ℓ₀) (u = e_i only matters when i is a source). As v2_77_tangent: at a random rational
member, moduli = rank(T + F) − rank(T), T the tangent of the H-orbit, F the four family directions. Also checks that
the member is a fold (all fold equations vanish at the point; exact over ℚ).
Usage: python3 analysis/v2_79_tangent.py QUIVER A0,A1,A2 PI1;PI2 X1,X2,... [SEED] [OUT.json]"""
import os, sys, json, itertools, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import numpy as np, sympy as sp
from sympy import QQ
from sympy.polys.rings import ring
import pcensus as pc, fold_ideal as fi
quiver = sys.argv[1]; cls = tuple(map(int, sys.argv[2].split(',')))
paths = [tuple(map(int, p.split('.'))) for p in sys.argv[3].split(';')]; X = [int(v) for v in sys.argv[4].split(',')]
seed = int(sys.argv[5]) if len(sys.argv) > 5 else 1; out = sys.argv[6] if len(sys.argv) > 6 else None
ctx = pc.setup(quiver, 7); d = ctx['d']; basis = ctx['basis']; idx = {b: k for k, b in enumerate(basis)}
arrows = ctx['arrows']; a0, a1, a2 = cls; i, j = arrows[a0]
R, *gens = ring(','.join(f'l{s},m{s}' for s in range(len(paths))), QQ)
M_ = np.asarray(ctx['M'], dtype=np.int64); DL = np.asarray(ctx['DL'], dtype=np.int64)
nzM = [(a, b, k, int(M_[a, b, k])) for a in range(d) for b in range(d) for k in range(d) if M_[a, b, k]]
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
        s = paths.index(r); l, m = gens[2 * s], gens[2 * s + 1]; corr = [R(0)] * d
        for w in X: corr[idx[('v', w)]] += l
        corr[idx[('p', (a1,))]] += m; corr[idx[('p', (a2,))]] -= m
        if u: corr = mul(unit(idx[('p', u)]), corr)
        v = [v[q] + corr[q] for q in range(d)]
    for rr in range(d): P[rr][col] = v[rr]
eqs = fi._equations(P, R, d, nzM, DL)
rng = random.Random(seed); syms = [sp.Symbol(str(g)) for g in R.gens]
val = {g: sp.Rational(rng.choice([-1, 1]) * rng.randint(1, 9), rng.randint(1, 7)) for g in syms}
fold_at_point = all(fi.to_expr(e).subs(val) == 0 for e in eqs)
Pexpr = sp.Matrix(d, d, lambda r, c: P[r][c].as_expr()); P0 = Pexpr.subs(val)
F = [sp.diff(Pexpr, g).subs(val) for g in syms]
classes = {}
for a, st in enumerate(arrows): classes.setdefault(st, []).append(a)
cl = list(classes.values()); cls_of = {a: k for k, arr in enumerate(cl) for a in arr}
th = []; G = []; s_ = []
for k, arr in enumerate(cl):
    m = len(arr); Mk = sp.zeros(m, m)
    for a in range(m - 1):
        for b in range(m):
            v = sp.Symbol(f'g{k}_{a}{b}'); th.append((v, 1 if a == b else 0)); Mk[a, b] = v
    for b in range(m): Mk[m - 1, b] = 1 - sum(Mk[a, b] for a in range(m - 1))
    G.append(Mk); v = sp.Symbol(f's{k}'); th.append((v, 1)); s_.append(v)
S = sp.zeros(d, d)
for col, b in enumerate(basis):
    if b[0] == 'v': S[col, col] = 1; continue
    choices = [[(cl[cls_of[a]][q], s_[cls_of[a]] * G[cls_of[a]][q, cl[cls_of[a]].index(a)]) for q in range(len(cl[cls_of[a]]))] for a in b[1]]
    for combo in itertools.product(*choices):
        cf = sp.Integer(1)
        for _, x in combo: cf *= x
        S[idx[('p', tuple(a for a, _ in combo))], col] += cf
at1 = dict(th); D = [sp.diff(S, v).subs(at1) for v, _ in th]
T = [(Dk * P0 - P0 * Dk) for Dk in D]
rT = sp.Matrix([list(M) for M in T]).rank(); rTF = sp.Matrix([list(M) for M in T + F]).rank()
res = dict(quiver=quiver, paths=paths, X=X, seed=seed, fold_identically=len(eqs) == 0, fold_at_point=fold_at_point,
           dim_h=len(D), rank_T=rT, moduli=rTF - rT)
print(json.dumps(res))
if out: json.dump(res, open(out, 'w'), indent=1)
