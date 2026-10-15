import sys
from cfibre import fibre
from loops_build import Q, tbuild
from extend import vertex_folds
name, L = sys.argv[1], int(sys.argv[2])
nv, arrows = Q[name]; A = set(arrows)
d, mult, Dl, allp = tbuild(nv, arrows, L + 1)
lab = {f"a{i}": allp[i] for i in range(nv, d)}
for P0 in vertex_folds(nv):
    if not all(not [a for a in range(nv) if P0[a, s]] or not [b for b in range(nv) if P0[b, t]] or
               all((a, b) in A for a in range(nv) if P0[a, s] for b in range(nv) if P0[b, t]) for s, t in arrows): continue
    R, J, n = fibre(nv, arrows, P0, L)
    if J.dimension() <= 0: continue
    print(f"{name} L={L} P0={list(P0)} dim {J.dimension()}")
    for P in J.minimal_associated_primes():
        gb = [g for g in P.groebner_basis() if not (len(g.monomials()) == 1 and g.degree() == 1)]
        free = sorted(str(v) for v in R.gens() if v not in P)
        print(f"   component dim {P.dimension()}: nonlinear/other generators {gb}; nonzero coordinates {[(v, lab[v.split('_')[0]], v.split('_')[1]) for v in free]}")
