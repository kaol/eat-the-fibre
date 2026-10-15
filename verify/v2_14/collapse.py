"""Theorem (loop-collapse folds) checked exactly on truncations, with symbolic parameters, plus negative controls."""
import sys
from sympy import symbols
from closedform import check
L = int(sys.argv[1]) if len(sys.argv) > 1 else 4
c = symbols('c0:6')
res = {}
# 2-cycle+loop with a tail through vertex 2 (σ(2)=0): arrows α:0→1, β:1→0, u:0→0, γ:1→2, δ:2→0
w = {0: lambda t: t - c[0], 1: lambda t: t - c[1], 2: lambda t: t, 3: lambda t: t - c[3], 4: lambda t: t - c[4]}
res['2-cycle+loop+tail, σ=(1,−1,0)'] = check('2-cycle+loop+tail', w, {0: 1, 1: -1, 2: 0}, L)
res['NEG three-vertex support σ=(1,−1,1)'] = check('2-cycle+loop+tail', w, {0: 1, 1: -1, 2: 1}, L)
res['NEG equal signs σ=(1,1,0)'] = check('2-cycle+loop+tail', w, {0: 1, 1: 1, 2: 0}, L)
w2 = dict(w); w2[0] = lambda t: 2 * t - c[0]
res['NEG non-monic weight 2u − c'] = check('2-cycle+loop+tail', w2, {0: 1, 1: -1, 2: 0}, L)
w3 = dict(w); w3[2] = lambda t: t - c[2]
res['NEG translated loop at s, w(u)=u−c'] = check('2-cycle+loop+tail', w3, {0: 1, 1: -1, 2: 0}, L)
# triangle with loop at 0 and back-arrow: support {0,1}, σ(2)=0; arrows (0,1),(1,2),(2,0),(0,0),(1,0)
wt = {0: lambda t: t - c[0], 1: lambda t: t - c[1], 2: lambda t: t - c[2], 3: lambda t: t, 4: lambda t: t - c[4]}
res['triangle+loop, σ=(1,−1,0)'] = check('triangle+loop', wt, {0: 1, 1: -1, 2: 0}, L)
for k, v in res.items(): print(f"{k}: fold equations hold: {v}")
ok = res['2-cycle+loop+tail, σ=(1,−1,0)'] and res['triangle+loop, σ=(1,−1,0)'] and not any(v for k, v in res.items() if k.startswith('NEG'))
print("collapse theorem and all negative controls behave as predicted:", ok)
