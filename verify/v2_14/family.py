import sys
from cfibre import fibre
from loops_build import Q, tbuild
from extend import vertex_folds
name, L, which = sys.argv[1], int(sys.argv[2]), sys.argv[3]
nv, arrows = Q[name]
d, mult, Dl, allp = tbuild(nv, arrows, L + 1)
for P0 in vertex_folds(nv):
    if str(list(P0)) != which: continue
    R, J, nvars = fibre(nv, arrows, P0, L)
    print("dimension", J.dimension())
    rad = J.radical()
    for P in rad.minimal_associated_primes():
        gb = P.groebner_basis()
        nz = sorted({str(v) for v in R.gens()} - {str(g) for g in gb if g.is_monomial() or len(g.monomials()) == 1})
        print("component dim", P.dimension(), " generators:", [str(g) for g in gb if not (len(g.monomials()) == 1 and g.degree() == 1)])
        print("   variables not forced to 0:", nz)
    pathname = {i: p for i, p in enumerate(allp)}
    print("path labels:", {f"a{i}": pathname[i] for i in range(nv, d)})
