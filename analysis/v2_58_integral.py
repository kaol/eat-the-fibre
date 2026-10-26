"""v2.58 (Program 2, step 0–1): the integral structure behind the stratum coordinates, and the census group over ℤ.

For QUIVER (default pairs_in_short):
1. Δ as an integer matrix Dmat (d² × d). Its Smith normal form: if every elementary divisor is 1, then for every
   prime p, ker(Δ ⊗ 𝔽_p) = (ker Δ over ℤ) ⊗ 𝔽_p, so K_p is spanned by the integer basis Kmat mod p.
2. Kmat (as used by fold_ideal.stratum): integral, of rank kdim mod every p (gcd of its maximal minors), and spans
   the saturated lattice ker_ℤ Δ.
3. B = [Kmat | e_C] for the complement used by fold_ideal: det B. For p ∤ det B the stratum coordinates (P_K, P(C))
   are an isomorphism over 𝔽_p.
4. The census group G (column-sum-1 substitutions on each parallel class, pcensus.build_G) with generic parameters:
   Φ(g) commutes with Δ, Σ_m Φ[m,l] DL[m] = Φ DL[l] Φᵀ, as polynomial identities over ℤ[g]. Its Lie algebra (column
   sums 0) contains the Der_Δ basis used by the pipeline.
5. Each class weighting φ_w (diagonal, w ∈ 𝔾_m^classes) is an algebra automorphism over ℤ[w^±1] and satisfies
   Δ ∘ φ_w = (φ_w ⊗ φ_w) ∘ Δ_w with Δ_w = Σ_c w_c Δ_c (checked as identities). Hence (char-free) if P is a fold
   compatible with Δ_w, φ_w P φ_w⁻¹ is a fold: ΔP' = (φ⊗φ)Δ_w P φ⁻¹ = (φ⊗φ)(P⊗P)Δ_w φ⁻¹ = (P'⊗P')Δ.
Usage: python3 analysis/v2_58_integral.py [QUIVER]"""
import os, sys, itertools, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import numpy as np, sympy as sp
from sympy.matrices.normalforms import smith_normal_form
import pcensus as pc, fold_ideal as fi
q = sys.argv[1] if len(sys.argv) > 1 else 'pairs_in_short'
ctx = pc.setup(q, 7); d, kdim, basis, arrows = ctx['d'], ctx['kdim'], ctx['basis'], ctx['arrows']
M = np.asarray(ctx['M']); DL = np.asarray(ctx['DL'])
Dmat = sp.Matrix(DL.reshape(d, d * d).T.tolist())
snf = smith_normal_form(Dmat, domain=sp.ZZ)
ed = [snf[i, i] for i in range(min(snf.shape)) if snf[i, i] != 0]
print(f"{q}: d = {d}, kdim = {kdim}; Δ: rank {len(ed)} = d − kdim: {len(ed) == d - kdim}; elementary divisors {sorted(set(abs(x) for x in ed))}")
Kmat = sp.Matrix(((np.asarray(ctx['Kmat']) + 3) % 7 - 3).tolist())
assert (Dmat * Kmat).is_zero_matrix, "Kmat not in ker Δ over ℤ"
mins = [Kmat.extract(list(r), list(range(kdim))).det() for r in itertools.combinations(range(d), kdim)]
g = 0
for m_ in mins: g = sp.gcd(g, m_)
print(f"Kmat: integral, in ker Δ over ℤ; gcd of maximal minors {g} (1: saturated basis, rank kdim mod every p)")
comp = [list(basis).index(b) for b in fi.PAIRS_IN_SHORT_C] if q == 'pairs_in_short' else list(ctx['Cidx'])
B = sp.Matrix.hstack(Kmat, sp.Matrix([[1 if i == j else 0 for j in comp] for i in range(d)]))
print(f"B = [Kmat | e_C], C = {[basis[i] for i in comp]}: det B = {B.det()}")
# 4. the census group with generic parameters
cls = collections.defaultdict(list)
for ai, st in enumerate(arrows): cls[st].append(ai)
big = [c for c in cls.values() if len(c) > 1]
idxb = {b: i for i, b in enumerate(basis)}
params = []; subs = {}
for ci, c in enumerate(big):
    m = len(c); top = sp.Matrix(m - 1, m, lambda i, j: sp.Symbol(f'g{ci}_{i}{j}')); params += list(top)
    gm = top.col_join(sp.Matrix(1, m, lambda _, j: 1 - sum(top[:, j])))
    for j, x in enumerate(c): subs[x] = {c[i]: gm[i, j] for i in range(m)}
Phi = sp.zeros(d, d)
for i, b in enumerate(basis):
    if b[0] == 'v': Phi[i, i] = 1; continue
    for choice in itertools.product(*[list(subs.get(x, {x: 1}).items()) for x in b[1]]):
        Phi[idxb[('p', tuple(x for x, _ in choice))], i] += sp.Mul(*[c_ for _, c_ in choice])
DLs = [sp.Matrix(DL[l].tolist()) for l in range(d)]
ok = all(sp.expand(sum((Phi[m, l] * DLs[m] for m in range(d)), sp.zeros(d, d)) - Phi * DLs[l] * Phi.T) == sp.zeros(d, d) for l in range(d))
mul = lambda x, y: sp.Matrix([sum(int(M[i, j, k]) * x[i] * y[j] for i in range(d) for j in range(d) if M[i, j, k]) for k in range(d)])
alg = all(sp.expand(Phi * mul(sp.eye(d)[:, i], sp.eye(d)[:, j]) - mul(Phi[:, i], Phi[:, j])) == sp.zeros(d, 1) for i in range(d) for j in range(d))
print(f"census G ({len(params)} parameters over ℤ): Φ is multiplicative: {alg}; Δ-compatible as an identity over ℤ[g]: {ok}")
# its Lie algebra at the identity (g = identity substitution) vs the pipeline's Der_Δ basis
idpt = {}
for ci, c in enumerate(big):
    for i in range(len(c) - 1):
        for j in range(len(c)): idpt[sp.Symbol(f'g{ci}_{i}{j}')] = 1 if i == j else 0
tang = [sp.Matrix(d, d, lambda a, b: sp.diff(Phi[a, b], x).subs(idpt)) for x in params]
from sympy.polys.matrices import DomainMatrix
def nullsp(rows, n):
    Nm = DomainMatrix.from_list_sympy(len(rows), n, rows).convert_to(sp.QQ).nullspace().to_Matrix()
    return [Nm.row(k).T for k in range(Nm.rows)]
rowsD = []
for i in range(d):
    for j in range(d):
        for m in range(d):
            row = [0] * (d * d)
            for k in range(d):
                if M[i, j, k]: row[m * d + k] += int(M[i, j, k])
            for r_ in range(d):
                if M[r_, j, m]: row[r_ * d + i] -= int(M[r_, j, m])
                if M[i, r_, m]: row[r_ * d + j] -= int(M[i, r_, m])
            if any(row): rowsD.append(row)
Der = [sp.Matrix(d, d, list(v)) for v in nullsp(rowsD, d * d)]
cond = [sum((list(sum((D[m, l] * DLs[m] for m in range(d)), sp.zeros(d, d)) - D * DLs[l] - DLs[l] * D.T) for l in range(d)), []) for D in Der]
cx = nullsp([list(r) for r in sp.Matrix(cond).T.tolist()], len(Der))
DerD = [sum((c[k] * Der[k] for k in range(len(Der))), sp.zeros(d, d)) for c in cx]
T = sp.Matrix([list(t) for t in tang]).T
sol = [T.gauss_jordan_solve(sp.Matrix(list(D)))[0] for D in DerD]
integral = all(all(x.is_integer for x in s_) for s_ in sol)
print(f"Der_Δ basis (dim {len(DerD)}) = integral combinations of the tangents of G at 1: {integral} "
      f"(so mod every p the pipeline's 𝔥 basis lies in Lie G_p)")
# 5. class weightings are algebra automorphisms (diagonal, multiplicative on paths) over ℤ[w^±]
ws = sp.symbols(f'w0:{len(cls)}'); wa = {a: ws[k] for k, c in enumerate(cls.values()) for a in c}
Ph = sp.diag(*[1 if b[0] == 'v' else sp.Mul(*[wa[x] for x in b[1]]) for b in basis])
algw = all(sp.expand(Ph * mul(sp.eye(d)[:, i], sp.eye(d)[:, j]) - mul(Ph[:, i], Ph[:, j])) == sp.zeros(d, 1) for i in range(d) for j in range(d))
print(f"class weightings φ_w: algebra automorphisms over ℤ[w^±1]: {algw}")
DLc = []
for c in cls.values():
    D_ = [sp.zeros(d, d) for _ in range(d)]
    for l, b in enumerate(basis):
        if b[0] == 'v': continue
        for k, x in enumerate(b[1]):
            if x not in c: continue
            pre = idxb[('p', b[1][:k])] if k > 0 else idxb[('v', arrows[x][0])]
            suf = idxb[('p', b[1][k + 1:])] if k + 1 < len(b[1]) else idxb[('v', arrows[x][1])]
            D_[l][pre, suf] += 1
    DLc.append(D_)
assert all(sum((DLc[c][l] for c in range(len(DLc))), sp.zeros(d, d)) == DLs[l] for l in range(d)), "Δ = Σ_c Δ_c"
Dw = [sum((ws[c] * DLc[c][l] for c in range(len(DLc))), sp.zeros(d, d)) for l in range(d)]
# Δ φ_w (column l) = Σ_m Ph[m,l] DL[m];  (φ⊗φ) Δ_w (column l) = Ph Dw[l] Phᵀ
okw = all(sp.expand(sum((Ph[m, l] * DLs[m] for m in range(d)), sp.zeros(d, d)) - Ph * Dw[l] * Ph.T) == sp.zeros(d, d) for l in range(d))
print(f"Δ ∘ φ_w = (φ_w ⊗ φ_w) ∘ Δ_w over ℤ[w]: {okw}")
