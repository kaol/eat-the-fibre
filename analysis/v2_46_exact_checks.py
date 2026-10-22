"""v2.46: exact checks over ℚ.
(25) At P₂ = P₀ + 2N: dim T_{P₂}V_A (Jacobian of the V_A equations at Z = 2·Z_N) and whether N ∈ span [Der(kQ), P₂].
(24) Images of the line families in the given logs: im P₀ for each line's base representative, grouped by equality.
Usage: python3 analysis/v2_46_exact_checks.py REPS.pkl REP_INDEX LINES_LOG LOG..."""
import os, sys, io, re, ast, contextlib, collections
import numpy as np, sympy as sp
from sympy.polys.matrices import DomainMatrix
HERE = os.path.dirname(os.path.abspath(__file__))
rf, ri, logf = sys.argv[1], int(sys.argv[2]), sys.argv[3]; more = sys.argv[4:]
sys.argv = ['family_curves.py', rf]
g = {'__name__': 'fc', '__file__': os.path.join(HERE, 'family_curves.py')}
with contextlib.redirect_stdout(io.StringIO()):
    exec(open(os.path.join(HERE, 'family_curves.py')).read().split("fams = []; nlift")[0], g)
d, M, DL = g['d'], g['M'], g['DL']
QQ = sp.QQ
def rank_rows(rows):
    return DomainMatrix.from_list_sympy(len(rows), len(rows[0]), rows).convert_to(QQ).rank() if rows else 0
I = sp.eye(d)
P0 = sp.Matrix(g['lifted'](ri).tolist())
Nd = [ast.literal_eval(m.group(1)) for m in re.finditer(r"N(?: \(nonzero entries\):| =) (\{.*\})", open(logf).read())][0]
N = sp.zeros(d, d)
for (i, j), v in Nd.items(): N[i, j] = sp.Rational(v)
# --- (25a) tangent of V_A at P₂
BA = P0.columnspace(); BK = (I - P0).columnspace(); nA, nK = len(BA), len(BK)
B = sp.Matrix.hstack(*BA, *BK); Binv = B.inv(); RK = Binv[nA:, :]; BAm = sp.Matrix.hstack(*BA)
Z = sp.Matrix(nA, nK, lambda i, j: sp.Symbol(f'z{i}_{j}'))
ZN = Binv[:nA, :] * N * sp.Matrix.hstack(*BK)
assert (BAm * ZN * RK - N).applyfunc(sp.simplify) == sp.zeros(d, d), "N is not of the form B_A Z R_K"
P = (P0 + BAm * Z * RK).applyfunc(sp.expand); Qm = (I - P).applyfunc(sp.expand)
def mul(x, y): return [sp.expand(sum(int(M[i, j, m]) * x[i] * y[j] for i in range(d) for j in range(d) if M[i, j, m])) for m in range(d)]
eqs = set()
for a in range(d):
    for b in range(d):
        for e in (P * sp.Matrix(mul(list(Qm[:, a]), list(Qm[:, b])))):
            e = sp.expand(e)
            if e != 0: eqs.add(e)
DLs = [sp.Matrix(DL[m].tolist()) for m in range(d)]
for l in range(d):
    for e in (sum((P[m, l] * DLs[m] for m in range(d)), sp.zeros(d, d)) - P * DLs[l] * P.T):
        e = sp.expand(e)
        if e != 0: eqs.add(e)
eqs = list(eqs); zs = list(Z)
pt = {z: 2 * ZN[i, j] for (i, j), z in zip([(i, j) for i in range(nA) for j in range(nK)], zs)}
assert all(sp.expand(e.subs(pt)) == 0 for e in eqs), "P₂ is not on V_A"
J = [[sp.diff(e, z).subs(pt) for z in zs] for e in eqs]
print(f"(25a) V_A: {len(eqs)} equations in {len(zs)} unknowns; dim T_{{P2}} V_A = {len(zs) - rank_rows(J)}", flush=True)
# --- (25b) N ∉ [Der, P₂]
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
Dn = DomainMatrix.from_list_sympy(len(rows), d * d, rows).convert_to(QQ).nullspace().to_Matrix()
Der = [sp.Matrix(d, d, list(Dn.row(k))) for k in range(Dn.rows)]
P2 = P0 + 2 * N
brs = [list(D * P2 - P2 * D) for D in Der]
r1 = rank_rows(brs); r2 = rank_rows(brs + [list(N)])
print(f"(25b) dim Der = {len(Der)}; rank [Der, P2] = {r1}; with N: {r2}; N in [Der, P2]: {r1 == r2}", flush=True)
# --- (24) images of the lines
items = []
for lf in [logf] + more:
    txt = open(lf).read()
    reps = [int(m.group(1)) for m in re.finditer(r"from orbit \d+ rep (\d+)", txt)] or [ri]
    for r_ in sorted(set(reps)): items.append((os.path.basename(lf), r_))
classes = []
for src, r_ in items:
    Pz = g['lifted'](r_)
    if Pz is None: print(f"  rep {r_}: not liftable"); continue
    A = sp.Matrix.hstack(*sp.Matrix(Pz.tolist()).columnspace())
    k = next((k for k, (A2, _) in enumerate(classes) if sp.Matrix.hstack(A, A2).rank() == A.rank() == A2.rank()), None)
    if k is None: classes.append((A, [(src, r_)]))
    else: classes[k][1].append((src, r_))
print(f"(24) {len(items)} base representatives, {len(classes)} distinct images (exact over Q):")
for k, (A, mem) in enumerate(classes):
    print(f"  image #{k} (dim {A.rank()}): {mem}")
