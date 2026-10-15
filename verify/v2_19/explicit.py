"""Explicit cyc A family: solve the exact fibre with the free coordinate fixed to c = 1 and list P − T (the deviation
from the canonical fold) on all paths of length ≤ L."""
import sys, itertools
from sage.all__sagemath_singular import QQ
import cfibre
from loops_build import Q, tbuild
from extend import vertex_folds
name, L = sys.argv[1], int(sys.argv[2]); nv, arrows = Q[name]
d, mult, Dl, allp = tbuild(nv, arrows, L + 1)
vseq = lambda p: (p[0],) if p[0] != 'a' else (arrows[p[1]][0],) + tuple(arrows[k][1] for k in p[1:])
want = eval(sys.argv[3])
for P0 in vertex_folds(nv):
    if [list(P0.column(v)) for v in range(nv)] != want: continue
    R, J, n = cfibre.fibre(nv, arrows, P0, L)
    free = [v for v in R.gens() if v not in J.radical()]
    pts = (J + R.ideal([free[0] - 1] if free else [])).variety(QQ) if free else J.variety(QQ)
    print("free coordinate:", free, " solutions with it set to 1:", len(pts))
    sol = pts[0]
    for k, v in sol.items():
        if v != 0:
            nm = str(k); x = int(nm[1:].split('_')[0]); w = int(nm.split('_')[1])
            print(f"   a({'→'.join(map(str, vseq(allp[x])))}) has e_{w}-component {v}")
