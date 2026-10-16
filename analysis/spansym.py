"""Exact Lie algebra of derivations D of A = kQ that preserve the span of the class deletions:
Δ_c D − (D⊗1 + 1⊗D) Δ_c = Σ_d m_cd Δ_d for every class c (unknowns: D and the matrix m).
Its flows map folds compatible with every Δ_c (split-compatible folds) to folds: they preserve the span, which contains Δ.
Contains Lie(Aut_Δ) (m = 0) and the class generators E_c (m diagonal). Reports its dimension, to compare with
dim Lie G + #classes. Usage: python3 analysis/spansym.py QUIVER"""
import os as _os, sys as _sys
if 'ROOT' not in globals():
    AN = _os.path.dirname(_os.path.abspath(__file__)); ROOT = _os.path.dirname(AN); _os.chdir(ROOT)
    __file__ = _os.path.join(ROOT, 'pcensus.py')
import sys, collections, sympy as sp, numpy as np
name = sys.argv[1]; sys.argv = ['x']
exec(open(_os.path.join(ROOT, 'pcensus.py')).read().split("if __name__ == '__main__':")[0])
nv, arrows = QUIVERS[name]; basis, M, DL = build(nv, arrows); d = len(basis); idx = {b: i for i, b in enumerate(basis)}
cls = collections.defaultdict(list)
for a, st in enumerate(arrows): cls[st].append(a)
cl = list(cls.values()); nc = len(cl)
DLc = []
for c in cl:
    D_ = np.zeros((d, d, d), dtype=np.int64)
    for l, b in enumerate(basis):
        if b[0] == 'v': continue
        for k, x in enumerate(b[1]):
            if x not in c: continue
            pre = idx[('p', b[1][:k])] if k > 0 else idx[('v', arrows[x][0])]
            suf = idx[('p', b[1][k + 1:])] if k + 1 < len(b[1]) else idx[('v', arrows[x][1])]
            D_[l, pre, suf] += 1
    DLc.append(D_)
Ds = sp.symbols(f'D0:{d*d}'); ms = sp.symbols(f'm0:{nc*nc}'); D = sp.Matrix(d, d, Ds)
eqs = []
for i in range(d):
    for j in range(d):
        lhs = sum((M[i, j, k] * D[:, k] for k in range(d) if M[i, j, k]), sp.zeros(d, 1)); rhs = sp.zeros(d, 1)
        for a in range(d):
            for k in range(d):
                if M[a, j, k]: rhs[k] += D[a, i]
                if M[i, a, k]: rhs[k] += D[a, j]
        eqs += list(lhs - rhs)
for c in range(nc):
    for l in range(d):
        T = sp.zeros(d, d)
        for m_ in range(d): T += D[m_, l] * sp.Matrix(DLc[c][m_])
        S = sp.Matrix(DLc[c][l]); T -= D * S + S * D.T
        for e in range(nc): T -= ms[c * nc + e] * sp.Matrix(DLc[e][l])
        eqs += list(T)
eqs = [e for e in eqs if e != 0]
A, _ = sp.linear_eq_to_matrix(eqs, list(Ds) + list(ms))
ns = A.nullspace()
mrank = sp.Matrix([[v[d*d + t] for t in range(nc * nc)] for v in ns]).rank() if ns else 0
print(f"{name}: dim(span-preserving derivations) = {len(ns)}; classes = {nc}; rank of the induced action on span{{Δ_c}} = {mrank}")
