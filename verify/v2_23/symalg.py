"""Exact symmetry Lie algebra of the fold problem for each quiver: derivations D of A = kQ with
ΔD − (D⊗1 + 1⊗D)Δ = c·Δ for a scalar c (c = 0: Aut_Δ; c ≠ 0 allowed: conformal, e.g. the path-length grading E).
Reports dim{c = 0} (Lie Aut_Δ), dim of the full conformal algebra, and the dimension of Lie(G) for the column-sum-1 group."""
import sys, sympy as sp, numpy as np
exec(open('gcensus.py').read().split("class Field")[0])
exec(open('pcensus.py').read().split('C = {}')[0].split('exec(open(os.path.join')[0])
QUIVERS.update(dict(K3=(2, [(0, 1)] * 3), K3_after=(3, [(0, 1)] * 3 + [(1, 2)]), K3_srcout=(3, [(0, 1)] * 3 + [(0, 2)]), two_pairs=(3, [(0, 1), (0, 1), (1, 2), (1, 2)])))
def symalg(name):
    nv, arrows = QUIVERS[name]; basis, M, DL = build(nv, arrows); d = len(basis)
    Ds = sp.symbols(f'D0:{d*d}'); c = sp.Symbol('c'); D = sp.Matrix(d, d, Ds)
    eqs = []
    for i in range(d):                       # derivation: D(b_i b_j) = D(b_i) b_j + b_i D(b_j)
        for j in range(d):
            lhs = sum((M[i, j, k] * D[:, k] for k in range(d) if M[i, j, k]), sp.zeros(d, 1))
            rhs = sp.zeros(d, 1)
            for a in range(d):
                if D[a, i] != 0:
                    for k in range(d):
                        if M[a, j, k]: rhs[k] += D[a, i]
                for k in range(d):
                    if M[i, a, k]: rhs[k] += D[a, j]
            eqs += list(lhs - rhs)
    for l in range(d):                       # Δ D b_l − (D⊗1 + 1⊗D) Δ b_l − c Δ b_l = 0
        T = sp.zeros(d, d)
        for m in range(d):
            T += D[m, l] * sp.Matrix(DL[m])
        S = sp.Matrix(DL[l])
        T -= D * S + S * D.T + c * S
        eqs += list(T)
    eqs = [e for e in eqs if e != 0]
    A, _ = sp.linear_eq_to_matrix(eqs, list(Ds) + [c])
    ns = A.nullspace()
    full = len(ns); delta_only = sum(1 for v in ns if v[-1] == 0)
    # Aut_Δ part = kernel of the c-coordinate restricted to the solution space
    cvals = sp.Matrix([[v[-1] for v in ns]])
    aut = full - cvals.rank()
    return full, aut
for name in sys.argv[1:]:
    full, aut = symalg(name)
    print(f"{name}: dim(conformal symmetry algebra) = {full}, dim Lie Aut_Δ (c = 0) = {aut}")
