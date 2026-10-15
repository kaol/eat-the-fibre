import sys
import cfibre, lin
from loops_build import Q, tbuild
from extend import vertex_folds
nv, arrows = Q['obstructed']
want = [[1, 0, 0, 0], [-1, 0, 0, 0], [0, 0, 1, 0], [0, 0, -1, 0]]
for P0 in vertex_folds(nv):
    if [list(P0.column(v)) for v in range(nv)] != want: continue
    for L in map(int, sys.argv[1:]):
        lin.build = lambda a, b, L=L: tbuild(a, b, L)
        nul = [lin.system(nv, arrows, P0, m) for m in range(1, L + 1)]
        R, J, n = cfibre.fibre(nv, arrows, P0, L)
        dim = J.dimension()
        extra = f", points {J.radical().vector_space_dimension()}" if dim == 0 else ""
        print(f"L={L}: linear nullities {nul}; exact fibre dimension {dim}{extra}", flush=True)
