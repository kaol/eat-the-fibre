"""v2.77: infinitesimal moduli of the path family (v2_76_pathfamily): at a rational member P (parameters from SEED),
the number of family directions ∂P/∂λ_s, ∂P/∂m_s outside the tangent of the H-orbit, T = ⟨[D, P] : D ∈ 𝔥⟩ with
𝔥 = the derivatives at 1 of the census symmetries (class substitutions with column sums 1, class weightings).
Exact over ℚ. moduli = rank(T + F) − rank(T).
Usage: python3 analysis/v2_77_tangent.py QUIVER A0,A1,A2 PI1;PI2 [SEED] [OUT.json]"""
import os, sys, json, itertools, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sympy as sp
from v2_76_pathfamily import build
quiver = sys.argv[1]; cls = tuple(map(int, sys.argv[2].split(',')))
paths = [tuple(map(int, p.split('.'))) for p in sys.argv[3].split(';')]
seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1; out = sys.argv[5] if len(sys.argv) > 5 else None
ctx, R, Pr = build(quiver, cls, paths)
d = ctx['d']; basis = ctx['basis']; arrows = ctx['arrows']; idx = {b: k for k, b in enumerate(basis)}
rng = random.Random(seed)
gens = [sp.Symbol(str(g)) for g in R.gens]
val = {g: sp.Rational(rng.choice([-1, 1]) * rng.randint(1, 9), rng.randint(1, 7)) for g in gens}
Pexpr = sp.Matrix(d, d, lambda r, c: Pr[r][c].as_expr())
P0 = Pexpr.subs(val)
F = [sp.diff(Pexpr, g).subs(val) for g in gens]
classes = {}
for a, st in enumerate(arrows): classes.setdefault(st, []).append(a)
cl = list(classes.values()); cls_of = {a: k for k, arr in enumerate(cl) for a in arr}
th = []; G = []; s = []
for k, arr in enumerate(cl):
    m = len(arr)
    Mk = sp.zeros(m, m)
    for a in range(m - 1):
        for b in range(m):
            v = sp.Symbol(f'g{k}_{a}{b}'); th.append((v, 1 if a == b else 0)); Mk[a, b] = v
    for b in range(m): Mk[m - 1, b] = 1 - sum(Mk[a, b] for a in range(m - 1))
    G.append(Mk); v = sp.Symbol(f's{k}'); th.append((v, 1)); s.append(v)
S = sp.zeros(d, d)
for col, b in enumerate(basis):
    if b[0] == 'v': S[col, col] = 1; continue
    choices = [[(cl[cls_of[a]][q], s[cls_of[a]] * G[cls_of[a]][q, cl[cls_of[a]].index(a)]) for q in range(len(cl[cls_of[a]]))] for a in b[1]]
    for combo in itertools.product(*choices):
        cf = sp.Integer(1)
        for _, x in combo: cf *= x
        S[idx[('p', tuple(a for a, _ in combo))], col] += cf
at1 = dict(th)
assert S.subs(at1) == sp.eye(d)
D = [sp.diff(S, v).subs(at1) for v, _ in th]
T = [(Dk * P0 - P0 * Dk) for Dk in D]
flat = lambda M: list(M)
rT = sp.Matrix([flat(M) for M in T]).rank()
rTF = sp.Matrix([flat(M) for M in T + F]).rank()
res = dict(quiver=quiver, paths=paths, seed=seed, dim_h=len(D), rank_T=rT, family_dirs=len(F), moduli=rTF - rT, point={str(k): str(v) for k, v in val.items()})
print(json.dumps(res))
if out: json.dump(res, open(out, 'w'), indent=1)
