"""v2.44: exact checks for a line family P(t) = P₀ + t·N found by analysis/v2_44_lines.py.
(a) P(t) is a fold over ℚ[t]: idempotent, image and kernel closed under multiplication (basis products), Δ-compatible.
(b) N ∉ [Der(kQ), P₀] over ℚ, with Der(kQ) computed exactly over ℚ from D(xy) = D(x)y + xD(y) on basis products.
Usage: python3 analysis/v2_44_exact.py REPS.pkl REP_INDEX LINES_LOG [which-N (0-based)]"""
import os, sys, io, re, ast, contextlib
import numpy as np, sympy as sp
from sympy.polys.matrices import DomainMatrix
HERE = os.path.dirname(os.path.abspath(__file__))
rf, ri, logf = sys.argv[1], int(sys.argv[2]), sys.argv[3]; which = int(sys.argv[4]) if len(sys.argv) > 4 else 0
sys.argv = ['family_curves.py', rf]
g = {'__name__': 'fc', '__file__': os.path.join(HERE, 'family_curves.py')}
with contextlib.redirect_stdout(io.StringIO()):
    exec(open(os.path.join(HERE, 'family_curves.py')).read().split("fams = []; nlift")[0], g)
d, M, smul = g['d'], g['M'], g['smul']
Ns = [ast.literal_eval(m.group(1)) for m in re.finditer(r"N \(nonzero entries\): (\{.*\})", open(logf).read())]
Nd = Ns[which]
N = sp.zeros(d, d)
for (i, j), v in Nd.items(): N[i, j] = sp.Rational(v)
P0 = sp.Matrix(g['lifted'](ri).tolist())
t = sp.Symbol('t'); I = sp.eye(d)
Pt = (P0 + t * N).applyfunc(sp.expand); Qt = (I - Pt).applyfunc(sp.expand)
idem = (Pt * Pt - Pt).applyfunc(sp.expand) == sp.zeros(d, d)
img = ker = True
for a in range(d):
    for b in range(d):
        u = smul(list(Pt[:, a]), list(Pt[:, b]))
        if any(sp.expand(x) for x in (Pt * sp.Matrix(u) - sp.Matrix(u))): img = False
        w = smul(list(Qt[:, a]), list(Qt[:, b]))
        if any(sp.expand(x) for x in (Pt * sp.Matrix(w))): ker = False
delta = g['is_fold_family'](Pt)
print(f"(a) over Q[t]: idempotent {idem}, image closed {img}, kernel closed {ker}, Δ-compatible {delta}")
# (b) derivations over Q: D (d×d, column j = D e_j); D(e_i e_j) = D(e_i) e_j + e_i D(e_j)
Ms = [[[int(M[i, j, m]) for m in range(d)] for j in range(d)] for i in range(d)]
rows = []
for i in range(d):
    for j in range(d):
        for m in range(d):
            row = [0] * (d * d)                                  # unknown D[r, c] at index r*d + c
            for r in range(d):
                row[r * d + m] += 0
            # (D(e_i e_j))_m = Σ_k M[i,j,k] D[m,k]
            for k in range(d):
                if Ms[i][j][k]: row[m * d + k] += Ms[i][j][k]
            # − (D e_i · e_j)_m = − Σ_r D[r,i] M[r,j,m];  − (e_i · D e_j)_m = − Σ_r D[r,j] M[i,r,m]
            for r in range(d):
                if Ms[r][j][m]: row[r * d + i] -= Ms[r][j][m]
                if Ms[i][r][m]: row[r * d + j] -= Ms[i][r][m]
            if any(row): rows.append(row)
A = DomainMatrix.from_list_sympy(len(rows), d * d, rows).convert_to(sp.QQ)
Dbasis = A.nullspace().to_Matrix()
print(f"(b) dim Der(kQ) over Q = {Dbasis.rows}")
brs = []
for k in range(Dbasis.rows):
    D = sp.Matrix(d, d, list(Dbasis.row(k)))
    brs.append(list((D * P0 - P0 * D).reshape(d * d, 1)))
B = DomainMatrix.from_list_sympy(len(brs), d * d, brs).convert_to(sp.QQ)
rB = B.rank()
B2 = DomainMatrix.from_list_sympy(len(brs) + 1, d * d, brs + [list(N.reshape(d * d, 1))]).convert_to(sp.QQ)
print(f"    rank [Der, P0] = {rB}; with N: {B2.rank()} -> N in [Der, P0]: {B2.rank() == rB}")
