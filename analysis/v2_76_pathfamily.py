"""v2.76: the path version of Theorem 5. Class arrows a0, a1, a2 : i→j (census indices), tail paths π₁, …, π_r from j
(tuples of arrow indices). P = Π₀ + N, Π₀ = substitution a0 ↦ a1 on every path, N(u·a0·π_s) = u·(λ_s c + m_s ℓ₀)
(an empty u is not multiplied in), N = 0 elsewhere; c = e_i ('src') and ℓ₀ = a1 − a2.
Checks the fold equations identically in ℚ[λ, m] (Gröbner basis of what fails), and, for INV=1, the H-invariance of
the cross-ratio φ for the first two tail paths, with S generic in H (symbolic) or random (INV_RANDOM=seed).
Usage: [INV=1] [INV_RANDOM=SEED] python3 analysis/v2_76_pathfamily.py QUIVER A0,A1,A2 PI1;PI2[;...] [OUT.json]
       (a path as arrow indices joined by '.', e.g. '3.5;4')"""
import os, sys, json, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE)); sys.path.insert(0, HERE)
import numpy as np, sympy as sp
from sympy import QQ
from sympy.polys.rings import ring
import pcensus as pc, fold_ideal as fi
from v2_75_phi import phi

def build(quiver, cls, paths):
    ctx = pc.setup(quiver, 7); d = ctx['d']; basis = ctx['basis']; idx = {b: i for i, b in enumerate(basis)}
    arrows = ctx['arrows']; a0, a1, a2 = cls; i, j = arrows[a0]
    assert arrows[a1] == arrows[a2] == (i, j)
    for pth in paths:
        assert arrows[pth[0]][0] == j and all(arrows[pth[k]][1] == arrows[pth[k + 1]][0] for k in range(len(pth) - 1))
    R, *gens = ring(','.join(f'l{s},m{s}' for s in range(len(paths))), QQ)
    par = {tuple(pth): (gens[2 * s], gens[2 * s + 1]) for s, pth in enumerate(paths)}
    M = np.asarray(ctx['M'], dtype=np.int64)
    nzM = [(a, b, k, int(M[a, b, k])) for a in range(d) for b in range(d) for k in range(d) if M[a, b, k]]
    def mul(x, y):
        out = [R(0)] * d
        for a, b, k, c in nzM:
            if x[a] and y[b]: out[k] += c * x[a] * y[b]
        return out
    unit = lambda k: [R(1) if r == k else R(0) for r in range(d)]
    P = [[R(1) if r == c else R(0) for c in range(d)] for r in range(d)]
    for col, bb in enumerate(basis):
        if bb[0] != 'p' or a0 not in bb[1]: continue
        k = bb[1].index(a0); u, r = bb[1][:k], bb[1][k + 1:]
        v = unit(idx[('p', u + (a1,) + r)])
        if r in par:
            l, m = par[r]; corr = [R(0)] * d
            corr[idx[('v', i)]] += l; corr[idx[('p', (a1,))]] += m; corr[idx[('p', (a2,))]] -= m
            if u: corr = mul(unit(idx[('p', u)]), corr)
            v = [v[q] + corr[q] for q in range(d)]
        for rr in range(d): P[rr][col] = v[rr]
    return ctx, R, P

if __name__ == '__main__':
    quiver = sys.argv[1]; cls = tuple(map(int, sys.argv[2].split(',')))
    paths = [tuple(map(int, p.split('.'))) for p in sys.argv[3].split(';')]; out = sys.argv[4] if len(sys.argv) > 4 else None
    ctx, R, P = build(quiver, cls, paths)
    d = ctx['d']; M = np.asarray(ctx['M'], dtype=np.int64); DL = np.asarray(ctx['DL'], dtype=np.int64)
    nzM = [(a, b, k, int(M[a, b, k])) for a in range(d) for b in range(d) for k in range(d) if M[a, b, k]]
    res = dict(quiver=quiver, cls=cls, paths=paths, d=d, params=[str(g) for g in R.gens])
    for part in ('idem', 'delta', 'closure'):
        res[part] = len(fi._equations(P, R, d, nzM, DL, (part,)))
    eqs = fi._equations(P, R, d, nzM, DL)
    if eqs: res['groebner'] = [str(g) for g in sp.groebner([fi.to_expr(e) for e in eqs], *[sp.Symbol(str(g)) for g in R.gens], order='grevlex').exprs]
    if os.environ.get('INV') == '1':
        basis = ctx['basis']; arrows = ctx['arrows']; idx = {b: k for k, b in enumerate(basis)}
        Pm = sp.Matrix(d, d, lambda r, c: P[r][c].as_expr())
        classes = {}
        for a, st in enumerate(arrows): classes.setdefault(st, []).append(a)
        cl = list(classes.values()); cls_of = {a: k for k, arr in enumerate(cl) for a in arr}
        s = list(sp.symbols(f's0:{len(cl)}')); G = []
        for k, arr in enumerate(cl):
            m = len(arr)
            Mk = sp.Matrix(m, m, lambda a, b: sp.Symbol(f'g{k}_{a}{b}') if a < m - 1 else 0) if m > 1 else sp.Matrix([[1]])
            if m > 1:
                for b in range(m): Mk[m - 1, b] = 1 - sum(Mk[a, b] for a in range(m - 1))
            G.append(Mk)
        RND = os.environ.get('INV_RANDOM')
        if RND:
            import random; rng = random.Random(int(RND))
            sub = {v: sp.Rational(rng.randint(-9, 9), rng.randint(1, 5)) for Mk in G for v in Mk.free_symbols}
            sub.update({v: sp.Rational(rng.choice([-1, 1]) * rng.randint(1, 9), rng.randint(1, 5)) for v in s})
            G = [Mk.subs(sub) for Mk in G]; s = [sub[v] for v in s]
        Gi = [Mk.inv() for Mk in G]
        def sigma(Gs, sc):
            S = sp.zeros(d, d)
            for col, b in enumerate(basis):
                if b[0] == 'v': S[col, col] = 1; continue
                choices = [[(cl[cls_of[a]][q], sc[cls_of[a]] * Gs[cls_of[a]][q, cl[cls_of[a]].index(a)]) for q in range(len(cl[cls_of[a]]))] for a in b[1]]
                for combo in itertools.product(*choices):
                    cf = sp.Integer(1)
                    for _, x in combo: cf *= x
                    S[idx[('p', tuple(a for a, _ in combo))], col] += cf
            return S
        S = sigma(G, s); Si = sigma(Gi, [1 / v for v in s])
        P2 = (S * Pm * Si).applyfunc(sp.cancel)
        carr = cl[cls_of[cls[0]]]
        f0, e0 = phi(Pm, basis, carr, paths[0], paths[1]); f1, e1 = phi(P2, basis, carr, paths[0], paths[1])
        res.update(phi_P=str(f0), err_P=e0, phi_SPSinv=str(f1), err_SPSinv=e1, inv_random=RND,
                   invariant=(f0 is not None and f1 is not None and sp.cancel(f0 - f1) == 0),
                   formula=(f0 is not None and sp.cancel(f0 - sp.Symbol('l0') * sp.Symbol('m1') / (sp.Symbol('l1') * sp.Symbol('m0'))) == 0))
    print(json.dumps(res))
    if out: json.dump(res, open(out, 'w'), indent=1)
