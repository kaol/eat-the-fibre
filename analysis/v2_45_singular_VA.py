"""v2.45: the variety V_A of folds with image A = im P₀ (P₀ a lifted census representative), as an ideal for Singular.

N = B_A · Z · R_K with B_A a basis of A (columns), R_K the rows of B⁻¹ dual to a basis of ker P₀, Z a 7×5 (generally
dim A × dim ker) matrix of unknowns. P = P₀ + N is idempotent with image A for every Z. The remaining fold conditions:
  kernel closure: P((1 − P)e_a · (1 − P)e_b) = 0 for all basis a, b;
  Δ-compatibility: Σ_m P[m, l] DL[m] = P DL[l] Pᵀ for all l.
Writes a Singular script, runs it, and reports dim V_A and the dimension of the tangent space at Z = 0.
Usage: python3 analysis/v2_45_singular_VA.py REPS.pkl REP_INDEX OUT_PREFIX"""
import os, sys, io, contextlib, subprocess
import numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
rf, ri, outp = sys.argv[1], int(sys.argv[2]), sys.argv[3]
sys.argv = ['family_curves.py', rf]
g = {'__name__': 'fc', '__file__': os.path.join(HERE, 'family_curves.py')}
with contextlib.redirect_stdout(io.StringIO()):
    exec(open(os.path.join(HERE, 'family_curves.py')).read().split("fams = []; nlift")[0], g)
d, M, DL = g['d'], g['M'], g['DL']
P0 = sp.Matrix(g['lifted'](ri).tolist()); I = sp.eye(d)
BA = P0.columnspace(); BK = (I - P0).columnspace()
nA, nK = len(BA), len(BK)
B = sp.Matrix.hstack(*BA, *BK); Binv = B.inv()
RK = Binv[nA:, :]
Z = sp.Matrix(nA, nK, lambda i, j: sp.Symbol(f'z{i}_{j}'))
N = sp.Matrix.hstack(*BA) * Z * RK
P = (P0 + N).applyfunc(sp.expand); Qm = (I - P).applyfunc(sp.expand)
def mul(x, y): return [sp.expand(sum(int(M[i, j, m]) * x[i] * y[j] for i in range(d) for j in range(d) if M[i, j, m])) for m in range(d)]
eqs = set()
assert (P * P - P).applyfunc(sp.expand) == sp.zeros(d, d)
for a in range(d):
    for b in range(a, d) if False else range(d):
        w = sp.Matrix(mul(list(Qm[:, a]), list(Qm[:, b])))
        for e in (P * w): e = sp.expand(e); eqs.add(e) if e != 0 else None
DLs = [sp.Matrix(DL[m].tolist()) for m in range(d)]
for l in range(d):
    lhs = sum((P[m, l] * DLs[m] for m in range(d)), sp.zeros(d, d))
    for e in (lhs - P * DLs[l] * P.T): e = sp.expand(e); eqs.add(e) if e != 0 else None
eqs = sorted(eqs, key=str)
vars_ = [str(v) for v in Z]
print(f"dim A = {nA}, dim ker = {nK}, unknowns {len(vars_)}, equations {len(eqs)}, max degree "
      f"{max(sp.Poly(e, *Z).total_degree() for e in eqs)}", flush=True)
# tangent at Z = 0: linear parts
lin = sp.Matrix([[sp.diff(e, v).subs({w: 0 for w in Z}) for v in Z] for e in eqs])
print(f"tangent dimension at Z = 0: {len(vars_) - lin.rank()}", flush=True)
src = f"""ring r = 0,({','.join(vars_)}),dp;
ideal I = {',\n  '.join(str(e).replace('**', '^') for e in eqs)};
ideal G = std(I);
print("dim " + string(dim(G)));
ideal L = minAssGTZ(I);
""" if False else f"""LIB "primdec.lib";
ring r = 0,({','.join(vars_)}),dp;
ideal I = {','.join(str(e).replace('**', '^') for e in eqs)};
ideal G = std(I);
print("dim " + string(dim(G)));
list L = minAssGTZ(I);
print("components " + string(size(L)));
int i;
for (i = 1; i <= size(L); i++) {{ print("component " + string(i) + " dim " + string(dim(std(L[i]))) + " contains 0: " + string(reduce(0, std(L[i])) == 0 && vdim(std(L[i] + maxideal(1))) == 1)); }}
quit;
"""
open(outp + '.sing', 'w').write(src)
res = subprocess.run(["Singular", "-q", "--no-rc", outp + '.sing'], capture_output=True, text=True, timeout=3000)
print(res.stdout[-3000:], res.stderr[-1000:])
