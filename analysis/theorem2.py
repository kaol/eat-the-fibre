"""Theorem 2 (v2.27): consecutive parallel classes give a one-parameter family of folds with an intrinsic invariant.

Let A: u → v and B: v → w be parallel classes with |A|, |B| ≥ 2, u ≠ w, and no arrow u → w, no arrow w → u, no loop at
u or w. Take δ_a = a₁ − a₂ ∈ A, δ_b = b₁ − b₂ ∈ B. For t ∈ k:
    P_t(e_u) = e_u − δ_a − (1 + t)·δ_aδ_b,   P_t(e_w) = e_w − δ_b + t·δ_aδ_b,   P_t = 0 on every other basis path.
Then P_t is a Δ-compatible fold. Writing P(e_u) = e_u − x − y, P(e_w) = e_w − z + w' (x, z of length 1; y, w' of
length 2), t = w'/(x·z) is preserved by every graded automorphism fixing the vertices (G, class weightings), so the P_t
are pairwise inequivalent.
This script checks the fold equations symbolically (over Q(t)) on several quivers, including one that violates the
hypothesis (an arrow u → w), where it must fail. Usage: python3 analysis/theorem2.py"""
import os as _os, sys as _sys
if 'ROOT' not in globals():
    AN = _os.path.dirname(_os.path.abspath(__file__)); ROOT = _os.path.dirname(AN); _os.chdir(ROOT)
    __file__ = _os.path.join(ROOT, 'pcensus.py')
import sys, collections, numpy as np, sympy as sp
sys.argv = ['x']
exec(open(_os.path.join(ROOT, 'pcensus.py')).read().split("if __name__ == '__main__':")[0])
t = sp.Symbol('t')
def check(name):
    nv, arrows = QUIVERS[name]; basis, M, DL = build(nv, arrows); d = len(basis); idx = {b: i for i, b in enumerate(basis)}
    cls = collections.defaultdict(list)
    for a, st in enumerate(arrows): cls[st].append(a)
    out = []
    for (u, v), A in cls.items():
        for (v2, w), B in cls.items():
            if v2 != v or len(A) < 2 or len(B) < 2: continue
            hyp = u != w and not any((s_, t_) in ((u, w), (w, u), (u, u), (w, w)) for s_, t_ in arrows)
            P = sp.zeros(d, d); eu, ew = idx[('v', u)], idx[('v', w)]
            a1, a2, b1, b2 = A[0], A[1], B[0], B[1]
            P[eu, eu] = 1; P[idx[('p', (a1,))], eu] = -1; P[idx[('p', (a2,))], eu] = 1
            P[ew, ew] = 1; P[idx[('p', (b1,))], ew] = -1; P[idx[('p', (b2,))], ew] = 1
            for (x, y), sgn in (((a1, b1), 1), ((a1, b2), -1), ((a2, b1), -1), ((a2, b2), 1)):
                P[idx[('p', (x, y))], eu] = -(1 + t) * sgn; P[idx[('p', (x, y))], ew] = t * sgn
            Id = sp.eye(d)
            def mul(xv, yv):
                o = [0] * d
                for i in range(d):
                    if xv[i] == 0: continue
                    for j in range(d):
                        if yv[j] == 0: continue
                        for kk in np.nonzero(M[i, j])[0]: o[kk] += xv[i] * yv[j] * int(M[i, j, kk])
                return sp.Matrix(o)
            ok = sp.expand(P * P - P) == sp.zeros(d, d)
            for i in range(d):
                for j in range(d):
                    h = mul(P[:, i], P[:, j]); q = mul((Id - P)[:, i], (Id - P)[:, j])
                    ok &= sp.expand(P * h - h) == sp.zeros(d, 1) and sp.expand(P * q) == sp.zeros(d, 1)
            for l in range(d):
                lhs = sum((P[m, l] * sp.Matrix(DL[m].tolist()) for m in range(d)), sp.zeros(d, d))
                ok &= sp.expand(lhs - P * sp.Matrix(DL[l].tolist()) * P.T) == sp.zeros(d, d)
            out.append(((u + 1, v + 1, w + 1), hyp, ok))
    return out
for q in ('two_pairs', 'pair_triple', 'two_pairs_pre', 'two_pairs_post', 'two_pairs_short', 'pairs_in', 'two_pairs_sep'):
    r = check(q)
    print(f"{q:16s}: " + ("; ".join(f"classes {u}⇉{v}⇉{w}: hypothesis {h}, P_t is a fold for all t: {ok}" for (u, v, w), h, ok in r)
                           if r else "no consecutive parallel classes"))
