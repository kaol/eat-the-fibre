import sys
from classify_loops import predicted
from cfibre import fibre
from extend import vertex_folds
from loops_build import tbuild
nv = 3; arrows = eval(sys.argv[1]); L = int(sys.argv[2]); A = set(arrows)
import loops_build; loops_build.Q['__q'] = (nv, arrows)
d, mult, Dl, allp = tbuild(nv, arrows, L + 1)
lab = {f"a{i}": allp[i] for i in range(nv, d)}
for P0 in vertex_folds(nv):
    if not all(not [a for a in range(nv) if P0[a, s]] or not [b for b in range(nv) if P0[b, t]] or
               all((a, b) in A for a in range(nv) if P0[a, s] for b in range(nv) if P0[b, t]) for s, t in arrows): continue
    R, J, n = fibre(nv, arrows, P0, L)
    dim = J.dimension()
    if dim > 0 and predicted(nv, arrows, P0) != dim:
        print("P0 columns:", [list(P0.column(v)) for v in range(nv)], "fibre dim", dim, flush=True)
        for P in J.minimal_associated_primes():
            free = [str(v) for v in R.gens() if v not in P]
            gb = [g for g in P.groebner_basis() if not (len(g.monomials()) == 1 and g.degree() == 1)]
            print("  component dim", P.dimension(), "nonzero:", [(lab[v.split('_')[0]], 'vertex ' + v.split('_')[1]) for v in free], "relations:", gb[:6], flush=True)
print("arrow indices:", list(enumerate(arrows)))
