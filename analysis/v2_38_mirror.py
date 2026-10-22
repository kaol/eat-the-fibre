"""v2.38: on 1⇉2⇇3, identify the inner-curve families of v2.36 with F1–F4 (Lean tables thm3F{1..4}P in
lean/EatTheFibre/Gen/pairs_in.lean) and test whether the source swap σ (vertex 1 ↔ 3, a_i ↔ b_i) maps families to
families. σ ∈ Aut(kQ) commutes with Δ, so σPσ⁻¹ is a fold; it is not in the census group G.
Usage: python3 analysis/v2_38_mirror.py PAIRS_IN_REPS_q.pkl
"""
import os, sys, io, re, ast, contextlib
import numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
pin_file = sys.argv[1]
sys.argv = ['family_curves.py', pin_file]
g = {'__name__': 'fc', '__file__': os.path.join(HERE, 'family_curves.py')}
with contextlib.redirect_stdout(io.StringIO()):
    exec(open(os.path.join(HERE, 'family_curves.py')).read(), g)
p, basis, key, key2loc = g['p'], g['basis'], g['key'], g['key2loc']
d = len(basis); idx = {b: i for i, b in enumerate(basis)}
arrows = g['arrows']
# σ: vertex 0 ↔ 2, arrow class (0,1) ↔ class (2,1), i-th arrow ↦ i-th arrow
cls = {}
for ai, st in enumerate(arrows): cls.setdefault(st, []).append(ai)
amap = {}
for x, y in zip(cls[(0, 1)], cls[(2, 1)]): amap[x] = y; amap[y] = x
vmap = {0: 2, 1: 1, 2: 0}
S = np.zeros((d, d), np.int64)
for i, b in enumerate(basis):
    img = ('v', vmap[b[1]]) if b[0] == 'v' else ('p', tuple(amap[a] for a in b[1]))
    S[idx[img], i] = 1
fams = [f for f in g['fams'] if len(f[0]) == p]
hits = {i: frozenset(h) for h, (i, vec) in fams}
def curve_pts(i, vec):
    Pz = g['lifted'](i); Inn, Innv = g['INN'][vec]
    Ps = (Inn * sp.Matrix(Pz.tolist()) * Innv).applyfunc(sp.expand)
    return [np.array(Ps.subs(g['t'], s).tolist(), np.int64) % p for s in range(p)]
print("σ-images of the families:")
for h, (i, vec) in fams:
    img = frozenset(key2loc.get(key(S @ P @ S.T % p)) for P in curve_pts(i, vec))
    match = [j for j, hj in hits.items() if hj == img]
    print(f"  family from rep {i}: σ-image is family {match if match else img}")
# Lean tables
src = open(os.path.join(HERE, '..', 'lean', 'EatTheFibre', 'Gen', 'pairs_in.lean')).read()
def table(name):
    m = re.search(rf"def {name} : Fin 7 → Fin 7 → Poly :=\s*(.*?)\n\s*\n", src, re.S)
    return ast.literal_eval(m.group(1).replace("![", "["))
def ev(poly, s): return sum(c * s ** e for e, c in enumerate(poly))
for F in range(1, 5):
    name = f"thm3F{F}P"
    if name not in src: continue
    T = table(name)
    locs = {key2loc.get(key(np.array([[ev(T[i][j], s) for j in range(d)] for i in range(d)], np.int64) % p))
            for s in range(p)}
    best = max(hits, key=lambda i: len(hits[i] & locs))
    print(f"F{F}: Lean table meets {len(locs)} local orbits; shares {len(hits[best] & locs)} with family from rep {best}")
