"""Local test of Conjecture 3* (finitely many orbits of Aut_Δ ⋊ T_Q) at every G-orbit of folds over F_p.

For each G-orbit: find a member that lifts to an exact Q-fold (entries p−1 ↦ −1); this is an orbit property ("some member
lifts"), unlike pcensus's transversal.lifted. At that Q-fold compute exactly:
  r_G = rank{[D, P] : D ∈ Lie G},  r_T = rank with the torus generators added — every class generator E_c if P is
  compatible with each class deletion Δ_c (Lemma of v2.23: then the whole class torus maps P to folds), otherwise the
  grading E only — and the Zariski tangent dimension t (product rule, rank mod 1000003).
t = r_T: the fold variety near P is one orbit of Aut_Δ ⋊ T (no moduli there).
t > r_T: walk (12 projected tangent steps of 0.05) to a generic fold and compare numerical local dimension with numerical r_T there (torus generators
  again chosen by the split residual). Generic local dim > generic r_T = moduli beyond Aut_Δ ⋊ T_Q (Conjecture 3* fails).
Usage (from anywhere): python3 analysis/orbitcheck.py QUIVER p [max generic tests per class]
Needs logs/folds_QUIVER_p.pkl (pcensus.py ... --save-folds logs/folds_QUIVER_p.pkl)."""
import os as _os, sys as _sys
if 'ROOT' not in globals():
    AN = _os.path.dirname(_os.path.abspath(__file__)); ROOT = _os.path.dirname(AN); _os.chdir(ROOT)
    __file__ = _os.path.join(ROOT, 'pcensus.py')
import sys, collections, numpy as np
name, p, per = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 2
sys.argv = ['torus.py', name, str(p), '0']
exec(open(_os.path.join(AN, 'torus.py')).read().split("exact = collections.Counter()")[0])
G = build_G(ctx)
Gg = (np.array([g for g, _ in G]), np.array([gi for _, gi in G]))
lab, greps, gsizes = orbit_label(folds, Gg, p)
members = collections.defaultdict(list)
for P in folds: members[lab[(P % p).astype(np.int8).tobytes()]].append(P)
Tgen = lambda split: Ec if split else [E]
def generic_far(x, steps=12, h=0.05):
    """Random walk on the fold variety: tangent step of size h, project back (least squares), repeat. Moves O(steps·h) away
    from special strata, so orbit ranks at the end point are well conditioned (a single 1e-2 step is not enough there)."""
    for _ in range(steps):
        N = tangent_basis(x)
        if N.shape[1] == 0: return x
        for _try in range(10):
            sol = least_squares(res, x + h * (N @ rng.normal(size=N.shape[1])) / np.sqrt(N.shape[1]), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=5000)
            if np.linalg.norm(res(sol.x)) < 1e-11: x = sol.x; break
        else: return None
    return x
tab = collections.Counter(); cands = collections.defaultdict(list); nolift = 0
for oid in range(len(greps)):
    Pq = next((q for q in (qfold(P) for P in members[oid]) if q is not None), None)
    if Pq is None: nolift += 1; continue
    s = split_ok(Pq)
    rows = [(D @ Pq - Pq @ D).ravel().tolist() for D in gens]
    rG = sp.Matrix(rows).rank() if rows else 0
    rT = sp.Matrix(rows + [(e @ Pq - Pq @ e).ravel().tolist() for e in Tgen(s)]).rank()
    t = jac(Pq); tab[(rG, rT, t, s)] += 1
    if t > rT: cands[(rG, rT, t, s)].append(Pq)
print(f"{name} F_{p}: {len(greps)} G-orbits; {nolift} with no member lifting to Q", flush=True)
print(f"  lifted orbits by (r_G, r_T, tangent, split compatible): {sorted(tab.items())}", flush=True)
verdicts = collections.Counter(); bad = 0
for key_, lst in sorted(cands.items()):
    for Pq in lst[:per]:
        y = generic_far(Pq.astype(float).ravel())
        if y is None: verdicts[(key_, 'no nearby solution')] += 1; continue
        split_y = split_res(y) < 1e-9
        Tf_y = [g.astype(float) for g in gens] + [e.astype(float) for e in Tgen(split_y)]
        A = np.array([(D @ y.reshape(d, d) - y.reshape(d, d) @ D).ravel() for D in Tf_y]); sv = np.linalg.svd(A, compute_uv=False)
        rTy = int((sv > 1e-4 * max(1, sv[0])).sum()); ldy = local_dim(y)
        verdicts[(key_, f'generic: local {ldy}, r_T {rTy}, split {split_y}, min orbit sv {sv[rTy-1] if rTy else 0:.1e}')] += 1; bad += int(ldy > rTy)
print(f"  generic-point verdicts where tangent > r_T: {sorted(verdicts.items())}", flush=True)
print(f"  generic points with local dim > r_T (moduli beyond Aut_Δ ⋊ T_Q): {bad}", flush=True)
