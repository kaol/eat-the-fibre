import numpy as np
from scipy.optimize import least_squares
from quiver import setup, jac
d, res, derivs = setup(2, [(0, 1), (0, 1)]); rng = np.random.default_rng(1)
orbit = lambda P: np.linalg.matrix_rank(np.array([(D @ P - P @ D).ravel() for D in derivs]), tol=1e-7)
cnt = 0
for t in range(400):
    x0 = rng.normal(scale=rng.choice([0.3, 1, 2]), size=d * d)
    sol = least_squares(res, x0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
    if np.linalg.norm(res(sol.x)) > 1e-11: continue
    x = sol.x; J = jac(res, x); U, s, Vt = np.linalg.svd(J); rk = int((s > 1e-6 * s[0]).sum())
    if d * d - rk != 3: continue
    N = Vt[rk:].T; out = []
    for eps in (1e-4, 1e-5, 1e-6):
        diffs = []
        for _ in range(60):
            sol2 = least_squares(res, x + eps * (N @ rng.normal(size=3)), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
            if np.linalg.norm(res(sol2.x)) < 1e-12: diffs.append(sol2.x - x)
        sv = np.linalg.svd(np.array(diffs), compute_uv=False) / eps
        out.append(np.round(sv[:4], 6))
    P = x.reshape(d, d)
    print(f"tangent 3, orbit {orbit(P)}: normalized sv of re-solved displacements for eps=1e-4,1e-5,1e-6:", [list(o) for o in out])
    cnt += 1
    if cnt >= 6: break
