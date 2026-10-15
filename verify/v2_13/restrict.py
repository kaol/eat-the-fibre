"""Do folds on kQ preserve the vertex algebra A_0, and are they determined by their restriction to A_0?"""
import sys, sympy as sp
from foldcount import run
from exactfold import QUIVERS
for name in sys.argv[1:] or ['A2', 'A3', '1→2←3']:
    nv = QUIVERS[name][0]; folds = run(name)
    pres = all(all(P[i, j] == 0 for i in range(nv, P.rows) for j in range(nv)) for P in folds)
    res = [tuple(P[:nv, :nv]) for P in folds]
    print(f"{name}: P(A_0) ⊂ A_0 for every fold: {pres}; distinct restrictions to A_0: {len(set(res))} of {len(folds)} folds", flush=True)
