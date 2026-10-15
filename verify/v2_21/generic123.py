"""At each lifted fold where local dim > (Aut_Δ + grading) orbit dim, step to a generic nearby fold on the local family
(converged solution at displacement 1e-2) and compare there: numerical local dim vs numerical orbit rank."""
import numpy as np
from scipy.optimize import least_squares
exec(open('localdim123.py').read().split("Edeg = matrix.diagonal")[0])
Edeg = np.diag([0, 0, 0, 1, 1, 1, 2, 2]).astype(float)
Dn = [np.array(Dm, dtype=float) for Dm in Ds] + [Edeg]
def orbit_rank(Pn):
    A = np.array([(Dm @ Pn - Pn @ Dm).ravel() for Dm in Dn]); s = np.linalg.svd(A, compute_uv=False)
    return int((s > 1e-6 * max(1, s[0])).sum())
rows = []
for f in folds:
    P = matrix(QQ, [[(-1 if x == 2 else x) for x in row] for row in f.tolist()])
    if not is_fold(P): continue
    J = residual_jacobian(P); tdim = d * d - J.rank()
    orbE = matrix(QQ, [(Dm * P - P * Dm).list() for Dm in Ds + [matrix.diagonal(QQ, [0, 0, 0, 1, 1, 1, 2, 2])]]).rank()
    if tdim == orbE: continue
    Pn = np.array(P, dtype=float); ld, _ = local_dim(Pn)
    if ld <= orbE: continue
    x = Pn.ravel(); Jn = jac(res, x); U, s, Vt = np.linalg.svd(Jn); rk = int((s > 1e-6 * s[0]).sum()); N = Vt[rk:].T
    best = None
    for _ in range(12):
        sol = least_squares(res, x + 1e-2 * (N @ rng.normal(size=N.shape[1])), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=5000)
        if np.linalg.norm(res(sol.x)) < 1e-11 and np.linalg.norm(sol.x - x) > 1e-4: best = sol.x; break
    if best is None: rows.append((ld, orbE, None, None)); continue
    Pg = best.reshape(d, d); ldg, _ = local_dim(Pg); og = orbit_rank(Pg)
    rows.append((ld, orbE, ldg, og)); print(f"  at P: local {ld}, orbit {orbE};  at generic nearby fold: local {ldg}, orbit {og}", flush=True)
from collections import Counter
print("flagged folds:", len(rows), " (local, orbit) at P → (local, orbit) nearby:", sorted(Counter(rows).items(), key=str))
print("generic local dim = generic orbit dim at all flagged folds:", all(r[2] is not None and r[2] == r[3] for r in rows))
