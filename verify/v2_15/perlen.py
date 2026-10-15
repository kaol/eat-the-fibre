import sys
from lin import system
from extend import vertex_folds
from arrowrule import QS
name, use = sys.argv[1], tuple(sys.argv[2].split(','))
nv, arrows = QS[name]; A = set(arrows); worst = {}
for P0 in vertex_folds(nv):
    if not all(not [a for a in range(nv) if P0[a, s]] or not [b for b in range(nv) if P0[b, t]] or
               all((a, b) in A for a in range(nv) if P0[a, s] for b in range(nv) if P0[b, t]) for s, t in arrows): continue
    for m in range(1, nv):
        k = system(nv, arrows, P0, m, use)
        if k is not None: worst[m] = max(worst.get(m, 0), k)
print(name, use, "max nullity per length:", worst)
