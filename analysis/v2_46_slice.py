"""v2.46: the dominant-slice lemma for V_A (rep 5107 line), exactly over ℚ.

(21) P_t = P₀ + t·N is compatible with each class-split coproduct Δ_c over ℚ[t] (so class weightings act on it).
(22) 𝔥 = {D ∈ Der_Δ ⊕ span(E_c) : D(A) ⊆ A}; rank of {[D, P_t] : D ∈ 𝔥} at t = 1, and with N added.
Der(kQ) and Der_Δ are computed exactly over ℚ from their linear defining equations:
  D(e_i e_j) = D(e_i) e_j + e_i D(e_j);   Σ_m D[m, l] DL[m] = D·DL[l] + DL[l]·Dᵀ  (Δ∘D = (D⊗1 + 1⊗D)∘Δ).
E_c = diag(number of arrows of class c in each basis path).
Usage: python3 analysis/v2_46_slice.py REPS.pkl REP_INDEX LINES_LOG [which-N]"""
import os, sys, io, re, ast, contextlib, collections
import numpy as np, sympy as sp
from sympy.polys.matrices import DomainMatrix
HERE = os.path.dirname(os.path.abspath(__file__))
rf, ri, logf = sys.argv[1], int(sys.argv[2]), sys.argv[3]; which = int(sys.argv[4]) if len(sys.argv) > 4 else 0
sys.argv = ['family_curves.py', rf]
g = {'__name__': 'fc', '__file__': os.path.join(HERE, 'family_curves.py')}
with contextlib.redirect_stdout(io.StringIO()):
    exec(open(os.path.join(HERE, 'family_curves.py')).read().split("fams = []; nlift")[0], g)
d, M, DL, basis, arrows = g['d'], g['M'], g['DL'], g['basis'], g['arrows']
idx = {b: i for i, b in enumerate(basis)}
QQ = sp.QQ
def nullsp(rows, n):
    if not rows: return [sp.Matrix([1 if i == j else 0 for i in range(n)]) for j in range(n)]
    A = DomainMatrix.from_list_sympy(len(rows), n, rows).convert_to(QQ)
    Nm = A.nullspace().to_Matrix()
    return [Nm.row(k).T for k in range(Nm.rows)]
def rank(vecs):
    if not vecs: return 0
    return DomainMatrix.from_list_sympy(len(vecs), len(vecs[0]), [list(v) for v in vecs]).convert_to(QQ).rank()
# class-split coproduct
cls = collections.defaultdict(list)
for a_, st in enumerate(arrows): cls[st].append(a_)
cl = list(cls.values())
DLc = []
for c in cl:
    D_ = [[[0] * d for _ in range(d)] for _ in range(d)]
    for l, b in enumerate(basis):
        if b[0] == 'v': continue
        for k, x in enumerate(b[1]):
            if x not in c: continue
            pre = idx[('p', b[1][:k])] if k > 0 else idx[('v', arrows[x][0])]
            suf = idx[('p', b[1][k + 1:])] if k + 1 < len(b[1]) else idx[('v', arrows[x][1])]
            D_[l][pre][suf] += 1
    DLc.append([sp.Matrix(D_[l]) for l in range(d)])
DLs = [sp.Matrix(DL[l].tolist()) for l in range(d)]
assert all(sum((DLc[c][l] for c in range(len(cl))), sp.zeros(d, d)) == DLs[l] for l in range(d)), "class split ≠ Δ"
# the line
Nd = [ast.literal_eval(m.group(1)) for m in re.finditer(r"N \(nonzero entries\): (\{.*\})", open(logf).read())][which]
N = sp.zeros(d, d)
for (i, j), v in Nd.items(): N[i, j] = sp.Rational(v)
P0 = sp.Matrix(g['lifted'](ri).tolist()); t = sp.Symbol('t')
Pt = (P0 + t * N).applyfunc(sp.expand)
# (21)
ok21 = []
for c in range(len(cl)):
    good = True
    for l in range(d):
        lhs = sum((Pt[m, l] * DLc[c][m] for m in range(d)), sp.zeros(d, d))
        if (lhs - Pt * DLc[c][l] * Pt.T).applyfunc(sp.expand) != sp.zeros(d, d): good = False; break
    ok21.append(good)
print(f"(21) classes {[[arrows[a_] for a_ in c][0] for c in cl]}: P_t compatible with Δ_c over Q[t]: {ok21}", flush=True)
# Der over Q
Mi = [[[int(M[i, j, m]) for m in range(d)] for j in range(d)] for i in range(d)]
rows = []
for i in range(d):
    for j in range(d):
        for m in range(d):
            row = [0] * (d * d)
            for k in range(d):
                if Mi[i][j][k]: row[m * d + k] += Mi[i][j][k]
            for r_ in range(d):
                if Mi[r_][j][m]: row[r_ * d + i] -= Mi[r_][j][m]
                if Mi[i][r_][m]: row[r_ * d + j] -= Mi[i][r_][m]
            if any(row): rows.append(row)
Der = [sp.Matrix(d, d, list(v)) for v in nullsp(rows, d * d)]
# Der_Δ: coefficients x with Σ x_k (Δ-condition of Der[k]) = 0
cond = []
for D in Der:
    v = []
    for l in range(d):
        E = sum((D[m, l] * DLs[m] for m in range(d)), sp.zeros(d, d)) - D * DLs[l] - DLs[l] * D.T
        v += list(E)
    cond.append(v)
cx = nullsp([list(r) for r in sp.Matrix(cond).T.tolist()], len(Der))
DerD = [sum((c[k] * Der[k] for k in range(len(Der))), sp.zeros(d, d)) for c in cx]
Ec = [sp.diag(*[0 if b[0] == 'v' else sum(1 for x in b[1] if x in c) for b in basis]) for c in cl]
print(f"dim Der = {len(Der)}, dim Der_Δ = {len(DerD)}, torus derivations {len(Ec)}", flush=True)
# admissible weights W = {w : P_t compatible with Δ_w for all t}: linear in w; coefficients of t must vanish
wv = sp.symbols(f'w0:{len(cl)}')
lin = []
for l in range(d):
    DLw = [sum((wv[c] * DLc[c][m] for c in range(len(cl))), sp.zeros(d, d)) for m in range(d)]
    E = (sum((Pt[m, l] * DLw[m] for m in range(d)), sp.zeros(d, d)) - Pt * DLw[l] * Pt.T).applyfunc(sp.expand)
    for e in E:
        if e != 0:
            for coeff in sp.Poly(e, t).all_coeffs(): lin.append([sp.expand(coeff).coeff(w) for w in wv])
W = nullsp(lin, len(cl))
print(f"admissible weights W (all t): {[list(w) for w in W]}", flush=True)
Ew = [sum((w[c] * Ec[c] for c in range(len(cl))), sp.zeros(d, d)) for w in W]
# 𝔥′ = stabiliser of A in Der_Δ ⊕ {Σ w_c E_c : w ∈ W}
L = DerD + Ew
I = sp.eye(d); Q0 = I - P0
hc = nullsp([list(r) for r in sp.Matrix([list(Q0 * D * P0) for D in L]).T.tolist()], len(L))
h = [sum((c[k] * L[k] for k in range(len(L))), sp.zeros(d, d)) for c in hc]
print(f"dim 𝔥 (stabiliser of A in the census Lie algebra) = {len(h)}", flush=True)
for tv in (1, 2, 3):
    P1 = Pt.subs(t, tv)
    brs = [list((D * P1 - P1 * D)) for D in h]
    r1 = rank(brs); r2 = rank(brs + [list(N)])
    print(f"(22) t = {tv}: rank [𝔥, P_t] = {r1}; with N: {r2}", flush=True)
