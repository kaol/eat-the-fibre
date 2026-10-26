"""v2.57 (Program 2, step 1): an explicit basis of Der_Δ for a quiver, whether D² = 0, and whether D(x)D(y) = 0 on the
basis (so that 1 + tD is multiplicative). Usage: python3 analysis/v2_57_derd.py QUIVER"""
import sys, io, contextlib, os
sys.path.insert(0, 'analysis'); sys.path.insert(0, '.')
q = sys.argv[1]
# only the group part: build ctx without the stratum ideal
import numpy as np, sympy as sp, pcensus as pc
from sympy.polys.matrices import DomainMatrix
ctx = pc.setup(q, 7); d = ctx['d']; basis = ctx['basis']; M = np.asarray(ctx['M']); DL = np.asarray(ctx['DL'])
def nullsp(rows, n):
    if not rows: return [sp.Matrix([1 if i == j else 0 for i in range(n)]) for j in range(n)]
    Nm = DomainMatrix.from_list_sympy(len(rows), n, rows).convert_to(sp.QQ).nullspace().to_Matrix()
    return [Nm.row(k).T for k in range(Nm.rows)]
rows = []
for i in range(d):
    for j in range(d):
        for m in range(d):
            row = [0] * (d * d)
            for k in range(d):
                if M[i, j, k]: row[m * d + k] += int(M[i, j, k])
            for r_ in range(d):
                if M[r_, j, m]: row[r_ * d + i] -= int(M[r_, j, m])
                if M[i, r_, m]: row[r_ * d + j] -= int(M[i, r_, m])
            if any(row): rows.append(row)
Der = [sp.Matrix(d, d, list(v)) for v in nullsp(rows, d * d)]
DLs = [sp.Matrix(DL[l].tolist()) for l in range(d)]
cond = [sum((list(sum((D[m, l] * DLs[m] for m in range(d)), sp.zeros(d, d)) - D * DLs[l] - DLs[l] * D.T) for l in range(d)), []) for D in Der]
cx = nullsp([list(r) for r in sp.Matrix(cond).T.tolist()], len(Der))
DerD = [sum((c[k] * Der[k] for k in range(len(Der))), sp.zeros(d, d)) for c in cx]
print(q, "d", d, "dim Der", len(Der), "dim Der_Δ", len(DerD))
lab = lambda i: ('e%d' % basis[i][1]) if basis[i][0] == 'v' else 'p' + ''.join(map(str, basis[i][1]))
for D in DerD:
    terms = [f"{D[a,b]}·{lab(a)}←{lab(b)}" for a in range(d) for b in range(d) if D[a, b]]
    print("  D:", ", ".join(terms), "| D²=0:", (D * D).is_zero_matrix, "| [D,D']=0 all:", all((D*E - E*D).is_zero_matrix for E in DerD))
# 1 + tD automorphism? check multiplicativity of (1+tD) symbolically: (1+tD)(xy) = (1+tD)x (1+tD)y  <=> t^2 D(x)D(y) = 0
for D in DerD:
    ok = True
    for i in range(d):
        for j in range(d):
            Dx = D[:, i]; Dy = D[:, j]
            prod = [sum(int(M[a, b, k]) * Dx[a] * Dy[b] for a in range(d) for b in range(d) if M[a, b, k]) for k in range(d)]
            if any(prod): ok = False
    print("  D(x)D(y) = 0 for all basis x, y:", ok)
