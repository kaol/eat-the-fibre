"""Numerical local dimension of the fold variety at transversal Q-folds where tangent dim > r_H, and at a generic nearby fold.
Residuals: P² − P; P(P_i·P_j) − P_i·P_j; P((1−P)_i·(1−P)_j); ΔP − (P⊗P)Δ.  Local dim: v2.12 ε-scaling test (solutions near
P + εN, N = tangent directions, whose displacement scales ∝ ε). Generic point: a converged solution at displacement ~1e-2.
Verdict per fold: generic local dim vs generic numerical rank of Lie(G) + E (r_H). local > r_H at a generic point = moduli
beyond Aut_Δ ⋊ grading (Conjecture 3″ fails there)."""
import os as _os, sys as _sys
if 'ROOT' not in globals():
    AN = _os.path.dirname(_os.path.abspath(__file__)); ROOT = _os.path.dirname(AN); _os.chdir(ROOT)
    __file__ = _os.path.join(ROOT, 'pcensus.py')   # pcensus.py locates gcensus.py relative to its own __file__
import sys, collections, numpy as np
from scipy.optimize import least_squares
name, p, maxn = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 30
sys.argv = ['x']
exec(open(_os.path.join(AN, 'tangent.py')).read().split("tab = collections.Counter(); seen = set()")[0].replace("name, p = sys.argv[1], int(sys.argv[2])\nsys.argv = ['x']\n", ""))
Mf, DLf = M.astype(float), DL.astype(float); If = np.eye(d)
def res(x):
    P = x.reshape(d, d)
    H = np.einsum('ai,bj,abk->ijk', P, P, Mf); Q = If - P; Qp = np.einsum('ai,bj,abk->ijk', Q, Q, Mf)
    r1 = P @ P - P
    r2 = np.einsum('kl,ijl->ijk', P, H) - H
    r3 = np.einsum('kl,ijl->ijk', P, Qp)
    r4 = np.einsum('ml,mab->lab', P, DLf) - np.einsum('xa,lab,yb->lxy', P, DLf, P)
    return np.concatenate([r1.ravel(), r2.ravel(), r3.ravel(), r4.ravel()])
def jacn(x, h=1e-6):
    cols = []
    for k in range(len(x)):
        e = np.zeros_like(x); e[k] = h; cols.append((res(x + e) - res(x - e)) / (2 * h))
    return np.array(cols).T
rng = np.random.default_rng(7)
def tangent_basis(x):
    J = jacn(x); U, s, Vt = np.linalg.svd(J); rk = int((s > 1e-7 * s[0]).sum()); return Vt[rk:].T
def local_dim(x, samples=24):
    N = tangent_basis(x)
    if N.shape[1] == 0: return 0
    svs = []
    for eps in (1e-4, 1e-5):
        diffs = []
        for _ in range(samples):
            sol = least_squares(res, x + eps * (N @ rng.normal(size=N.shape[1])), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=2000)
            if np.linalg.norm(res(sol.x)) < 1e-11: diffs.append(sol.x - x)
        sv = np.zeros(N.shape[1])
        if diffs:
            v = np.linalg.svd(np.array(diffs), compute_uv=False)[:N.shape[1]] / eps; sv[:len(v)] = v
        svs.append(sv)
    return int(sum((a > 0.2) and (b > 0.2) and (0.25 < b / a < 4) for a, b in zip(svs[0], svs[1])))
Gf = [g.astype(float) for g in gens] + [E.astype(float)]
def rankH(x, tol=1e-6):
    P = x.reshape(d, d); A = np.array([(D @ P - P @ D).ravel() for D in Gf]); s = np.linalg.svd(A, compute_uv=False)
    return int((s > tol * max(1, s[0])).sum())
def generic_near(x):
    N = tangent_basis(x)
    for _ in range(20):
        sol = least_squares(res, x + 1e-2 * (N @ rng.normal(size=N.shape[1])), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=5000)
        if np.linalg.norm(res(sol.x)) < 1e-11 and np.linalg.norm(sol.x - x) > 1e-4: return sol.x
    return None
cands = collections.defaultdict(list)
for P in folds:
    Pq = qfold(P)
    if Pq is None: continue
    rows = [(D @ Pq - Pq @ D).ravel().tolist() for D in gens]
    rG = sp.Matrix(rows).rank(); rH = sp.Matrix(rows + [(E @ Pq - Pq @ E).ravel().tolist()]).rank()
    if rH <= rG: continue
    t = jac(Pq)
    if t > rH: cands[(rG, rH, t)].append(Pq)
out = collections.Counter(); bad = 0
for key_, lst in sorted(cands.items()):
    for Pq in lst[:max(1, maxn // len(cands))]:
        x = Pq.astype(float).ravel(); ld = local_dim(x)
        if ld <= key_[1]: out[(key_, 'local ≤ r_H at P')] += 1; continue
        y = generic_near(x)
        if y is None: out[(key_, 'no nearby solution')] += 1; continue
        ldy, rHy = local_dim(y), rankH(y)
        verdict = f'generic: local {ldy}, r_H {rHy}'
        out[(key_, verdict)] += 1; bad += int(ldy > rHy)
        print(f"  {key_}: local at P {ld}; {verdict}", flush=True)
print(f"{name} F_{p}: (r_G, r_H, tangent) -> local-dimension verdicts: {sorted(out.items())}")
print(f"generic points with local dim > r_H (moduli beyond Aut_Δ ⋊ grading): {bad}")
