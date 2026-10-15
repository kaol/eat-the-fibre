"""Deconstructing folds on path algebras of acyclic quivers: sample, then compare local dimension with the
orbit dimension of Delta-compatible automorphisms (Lie algebra = derivations D with Delta D = (D(x)1 + 1(x)D) Delta)."""
import numpy as np, sys, itertools
from scipy.optimize import least_squares
from collections import Counter
def build(nv, arrows):
    paths = [(v,) for v in range(nv)]                       # vertex paths: (v,)
    def extend(p):
        for k, (s, t) in enumerate(arrows):
            if (len(p) == 1 and p[0] == s) or (len(p) > 1 and arrows[p[-1]][1] == s):
                yield (p + (k,)) if len(p) > 1 else ('a', k)
    # represent arrow paths as tuples ('a', k1, k2, ...)
    allp = [(v,) for v in range(nv)]; frontier = [('a', k) for k in range(len(arrows))]
    while frontier:
        allp += frontier; nf = []
        for p in frontier:
            last = arrows[p[-1]][1]
            for k, (s, t) in enumerate(arrows):
                if s == last: nf.append(p + (k,))
        frontier = nf
    idx = {p: i for i, p in enumerate(allp)}; d = len(allp)
    start = lambda p: p[0] if p[0] != 'a' else arrows[p[1]][0]
    end = lambda p: p[0] if p[0] != 'a' else arrows[p[-1]][1]
    mult = np.zeros((d, d, d))
    for p in allp:
        for q in allp:
            if end(p) != start(q): continue
            if p[0] != 'a': r = q
            elif q[0] != 'a': r = p
            else: r = p + q[1:]
            mult[idx[p], idx[q], idx[r]] = 1
    Dl = np.zeros((d, d, d))                                 # Delta(path) = sum over its arrows prefix (x) suffix
    for p in allp:
        if p[0] != 'a': continue
        ks = p[1:]
        for i, k in enumerate(ks):
            pre = ('a',) + ks[:i] if i > 0 else (arrows[k][0],)
            suf = ('a',) + ks[i + 1:] if i + 1 < len(ks) else (arrows[k][1],)
            Dl[idx[p], idx[pre], idx[suf]] += 1
    return d, mult, Dl, allp
def setup(nv, arrows):
    d, mult, Dl, allp = build(nv, arrows); I = np.eye(d)
    mul = lambda x, y: np.einsum('i,j,ijk->k', x, y, mult)
    def residual(pv):
        P = pv.reshape(d, d); r = [(P @ P - P).ravel()]
        for i in range(d):
            for j in range(d):
                x, y = P[:, i], P[:, j]; xy = mul(x, y); r.append(P @ xy - xy)
                r.append(P @ mul(I[:, i] - x, I[:, j] - y))
        for j in range(d):
            lhs = np.einsum('i,ijk->jk', P[:, j], Dl)
            rhs = np.einsum('ab,ja,kb->jk', Dl[j], P, P)
            r.append((lhs - rhs).ravel())
        return np.concatenate(r)
    # Lie algebra of Delta-compatible automorphisms: derivations D (D(xy) = D(x)y + xD(y)) with Delta D = (D(x)1 + 1(x)D) Delta
    rows = []
    for i in range(d):
        for j in range(d):
            # derivation: D(b_i b_j) - D(b_i) b_j - b_i D(b_j) = 0, linear in D entries D[k, l] (D[:, l] = D(b_l))
            for out in range(d):
                row = np.zeros((d, d))
                for k in range(d): row[out, k] += mult[i, j, k]          # D(b_i b_j)_out = sum_k mult[i,j,k] D[out,k]
                row[:, i] -= mult[:, j, out]                             # D(b_i)_a * (b_a b_j)_out
                row[:, j] -= mult[i, :, out]
                rows.append(row.ravel())
    for l in range(d):
        for (a_, b_) in itertools.product(range(d), repeat=2):
            row = np.zeros((d, d))
            for k in range(d): row[k, l] += Dl[k, a_, b_]                # Delta(D b_l)[a,b]
            for x in range(d):
                for y in range(d):
                    c = Dl[l, x, y]
                    if c: row[a_, x] -= c * (1 if y == b_ else 0); row[b_, y] -= c * (1 if x == a_ else 0)
            rows.append(row.ravel())
    M = np.array(rows); _, s, Vt = np.linalg.svd(M); rk = int((s > 1e-9 * s[0]).sum())
    derivs = [v.reshape(d, d) for v in Vt[rk:]]
    return d, residual, derivs
def jac(res, x, h=1e-7):
    f0 = res(x); J = np.zeros((len(f0), len(x)))
    for k in range(len(x)):
        e = np.zeros(len(x)); e[k] = h; J[:, k] = (res(x + e) - res(x - e)) / (2 * h)
    return J
def analyse(name, nv, arrows, trials, seed=0):
    d, res, derivs = setup(nv, arrows); rng = np.random.default_rng(seed)
    orbit = lambda P: np.linalg.matrix_rank(np.array([(D @ P - P @ D).ravel() for D in derivs]), tol=1e-7) if derivs else 0
    def local_dim(P, samples=30):
        """count displacement directions that scale with eps (genuine), discarding a constant noise floor."""
        x = P.ravel(); J = jac(res, x); U, s, Vt = np.linalg.svd(J); rk = int((s > 1e-6 * s[0]).sum()); N = Vt[rk:].T
        if N.shape[1] == 0: return 0, 0
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
        genuine = [(a > 0.2) and (b > 0.2) and (0.25 < b / a < 4) for a, b in zip(svs[0], svs[1])]
        return int(sum(genuine)), d * d - rk
    stats = Counter(); seen = 0
    for t in range(trials):
        x0 = rng.normal(scale=rng.choice([0.3, 1, 2]), size=d * d)
        sol = least_squares(res, x0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
        if np.linalg.norm(res(sol.x)) > 1e-10: continue
        P = sol.x.reshape(d, d); ld, td = local_dim(P); od = orbit(P); seen += 1
        stats[(int(np.linalg.matrix_rank(P, tol=1e-6)), td, ld, od)] += 1
    mx = max((ld - od for (_, _, ld, od) in stats), default=0)
    print(f"{name}: dim A = {d}, dim G = {len(derivs)}, {seen} folds sampled; max (local dim − orbit dim) = {mx}")
    for k, v in sorted(stats.items()): print(f"     (rank P, tangent, local, orbit) = {k}: {v}")
    return mx
if __name__ == '__main__':
    Q = {'Kronecker 3': (2, [(0, 1), (0, 1), (0, 1)]), 'A2 (1→2)': (2, [(0, 1)]), 'Kronecker 2': (2, [(0, 1), (0, 1)]), 'A3 (1→2→3)': (3, [(0, 1), (1, 2)]),
         '1⇉2→3': (3, [(0, 1), (0, 1), (1, 2)]), '1→2←3': (3, [(0, 1), (2, 1)])}
    which = sys.argv[1:] or list(Q)
    import os
    for nm in which: analyse(nm, *Q[nm], trials=int(os.environ.get('TRIALS', 60)))
