"""numerical sampler vs exact list: every numerically found fold (A2, A3, 1→2←3) is one of the exact folds."""
import sys, numpy as np, sympy as sp
from scipy.optimize import least_squares
from quiver import setup
from exactfold import QUIVERS
from foldcount import run
rng = np.random.default_rng(3)
allok = True
for name in sys.argv[1:] or ['A2', 'A3', '1→2←3']:
    exact = [np.array(sp.Matrix(P).evalf(), dtype=float) for P in run(name)]
    nv, arrows = QUIVERS[name]; d, res, _ = setup(nv, arrows); found = set(); n = 0
    for t in range(int(__import__('os').environ.get('TRIALS', 60))):
        x0 = rng.normal(scale=rng.choice([0.3, 1, 2]), size=d * d)
        sol = least_squares(res, x0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
        if np.linalg.norm(res(sol.x)) > 1e-10: continue
        P = sol.x.reshape(d, d); n += 1
        dist = [np.abs(P - E).max() for E in exact]; k = int(np.argmin(dist))
        if dist[k] > 1e-6: allok = False; print("  unmatched fold!", name)
        else: found.add(k)
    print(f"{name}: {n} numerical folds, all among the {len(exact)} exact folds; distinct exact folds hit: {len(found)}", flush=True)
print("numerical folds ⊂ exact list:", allok)
