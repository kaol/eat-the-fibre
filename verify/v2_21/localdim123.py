"""At every lifted Q-fold of k(1⇉2→3) whose tangent dimension exceeds its orbit dimension, measure the numerical local
dimension (v2.12 ε-scaling test: solutions of the fold equations near P + εN that move ∝ ε)."""
import sys, os, numpy as np
from scipy.optimize import least_squares
from quiver import setup, jac
exec(open('tangent123.py').read().split("lifted = 0; rows = []")[0])
dq, res, derivs = setup(3, [(0, 1), (0, 1), (1, 2)])
rng = np.random.default_rng(1)
def local_dim(Pn, samples=24):
    x = Pn.ravel(); J = jac(res, x); U, s, Vt = np.linalg.svd(J); rk = int((s > 1e-6 * s[0]).sum()); N = Vt[rk:].T
    if N.shape[1] == 0: return 0, 0
    svs = []
    for eps in (1e-4, 1e-5):
        diffs = []
        for _ in range(samples):
            sol = least_squares(res, x + eps * (N @ rng.normal(size=N.shape[1])), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
            if np.linalg.norm(res(sol.x)) < 1e-11: diffs.append(sol.x - x)
        sv = np.zeros(N.shape[1])
        if diffs:
            v = np.linalg.svd(np.array(diffs), compute_uv=False)[:N.shape[1]] / eps; sv[:len(v)] = v
        svs.append(sv)
    genuine = [(a > 0.2) and (b > 0.2) and (0.25 < b / a < 4) for a, b in zip(svs[0], svs[1])]
    return int(sum(genuine)), N.shape[1]
Edeg = matrix.diagonal(QQ, [0, 0, 0, 1, 1, 1, 2, 2])     # grading generator: length
out = []
for f in folds:
    P = matrix(QQ, [[(-1 if x == 2 else x) for x in row] for row in f.tolist()])
    if not is_fold(P): continue
    J = residual_jacobian(P); tdim = d * d - J.rank()
    orb = matrix(QQ, [(Dm * P - P * Dm).list() for Dm in Ds]).rank()
    orbE = matrix(QQ, [(Dm * P - P * Dm).list() for Dm in Ds + [Edeg]]).rank()
    if tdim == orb: out.append((tdim, orb, orbE, orb)); continue
    ld, nt = local_dim(np.array(P, dtype=float))
    out.append((tdim, orb, orbE, ld)); print(f"  tangent {tdim}, Aut orbit {orb}, Aut+grading orbit {orbE}, numerical local dim {ld}", flush=True)
from collections import Counter
print("lifted folds:", len(out), " (tangent, Aut orbit, Aut+grading orbit, local) distribution:", sorted(Counter(out).items()))
print("local dim ≤ Aut+grading orbit dim at every lifted fold:", all(l <= oe for t, o, oe, l in out))
print("folds with local dim > Aut+grading orbit dim:", sum(1 for t, o, oe, l in out if l > oe))
