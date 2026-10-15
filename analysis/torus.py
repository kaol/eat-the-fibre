"""Class-weight torus T_Q: φ_w scales every arrow of parallel class c by w_c (single arrows are classes of size 1).
Exact at lifted Q-folds: compatibility with each Δ_c, and ranks r_G (Lie Aut_Δ), r_H (+ grading E), r_T (+ every class
generator E_c), against the Zariski tangent dimension. Numerical at generic nearby folds: Δ_c-residuals, local dimension,
and rank r_T. Verdict: local dim ≤ r_T everywhere means no moduli beyond Aut_Δ ⋊ T_Q at the points tested."""
import os as _os, sys as _sys
if 'ROOT' not in globals():
    AN = _os.path.dirname(_os.path.abspath(__file__)); ROOT = _os.path.dirname(AN); _os.chdir(ROOT)
    __file__ = _os.path.join(ROOT, 'pcensus.py')   # pcensus.py locates gcensus.py relative to its own __file__
import sys, collections, numpy as np
name, p, maxn = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 24
sys.argv = ['localdim.py', name, str(p), str(maxn)]
exec(open(_os.path.join(AN, 'localdim.py')).read().split("cands = collections.defaultdict(list)")[0])
idx = {b: i for i, b in enumerate(basis)}
cls = collections.defaultdict(list)
for a, st in enumerate(arrows): cls[st].append(a)
Ec, DLc = [], []
for st, c in cls.items():
    Ec.append(np.diag([0 if b[0] == 'v' else sum(1 for x in b[1] if x in c) for b in basis]).astype(np.int64))
    D_ = np.zeros((d, d, d), dtype=np.int64)
    for l, b in enumerate(basis):
        if b[0] == 'v': continue
        s = b[1]
        for k, a in enumerate(s):
            if a not in c: continue
            pre = idx[('p', s[:k])] if k > 0 else idx[('v', arrows[a][0])]
            suf = idx[('p', s[k + 1:])] if k + 1 < len(s) else idx[('v', arrows[a][1])]
            D_[l, pre, suf] += 1
    DLc.append(D_)
def split_ok(P): return all((np.einsum('m,mab->ab', P[:, l], D_) == P @ D_[l] @ P.T).all() for D_ in DLc for l in range(d))
def split_res(x):
    P = x.reshape(d, d); return max(np.abs(np.einsum('ml,mab->lab', P, D_.astype(float)) - np.einsum('xa,lab,yb->lxy', P, D_.astype(float), P)).max() for D_ in DLc)
Tf = [g.astype(float) for g in gens] + [e.astype(float) for e in Ec]
def rankT(x, tol=1e-6):
    P = x.reshape(d, d); A = np.array([(D @ P - P @ D).ravel() for D in Tf]); s = np.linalg.svd(A, compute_uv=False)
    return int((s > tol * max(1, s[0])).sum())
exact = collections.Counter(); numer = collections.Counter(); cands = collections.defaultdict(list); nsplit = [0, 0]
for P in folds:
    Pq = qfold(P)
    if Pq is None: continue
    rows = [(D @ Pq - Pq @ D).ravel().tolist() for D in gens]
    rG = sp.Matrix(rows).rank() if rows else 0
    rH = sp.Matrix(rows + [(E @ Pq - Pq @ E).ravel().tolist()]).rank()
    if rH <= rG: continue
    rT = sp.Matrix(rows + [(e @ Pq - Pq @ e).ravel().tolist() for e in Ec]).rank()
    t = jac(Pq); s_ok = split_ok(Pq); nsplit[s_ok] += 1
    exact[(rG, rH, rT, t, s_ok)] += 1
    if t > rT: cands[(rG, rH, rT, t)].append(Pq)
print(f"{name}: transversal Q-folds by (r_G, r_H, r_T, tangent, Δ_c-split compatible): {sorted(exact.items())}", flush=True)
for key_, lst in sorted(cands.items()):
    for Pq in lst[:max(1, maxn // max(1, len(cands)))]:
        y = generic_near(Pq.astype(float).ravel())
        if y is None: numer[(key_, 'no nearby solution')] += 1; continue
        ldy, rTy, sr = local_dim(y), rankT(y), split_res(y)
        numer[(key_, f'generic: local {ldy}, r_T {rTy}, split residual {"<1e-9" if sr < 1e-9 else f"{sr:.1e}"}')] += 1
print(f"{name}: at folds with tangent > r_T, generic nearby verdicts: {sorted(numer.items())}", flush=True)
