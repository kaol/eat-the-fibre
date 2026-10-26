"""v2.58 (Program 2): the vertex-block theorem of v2.52 in every characteristic.

The v2.52 system (P_K idempotent, image and kernel closed under products on K, and the Δ-conditions of the given
arrows), written over ℤ: the K-coordinates of a product use the integral rows of B⁻¹ (B = [Kmat | e_C], det ±1)
instead of v2.52's rational left inverse, so the equations are the correct conditions over every 𝔽_p.
MODE int: eliminate to the 9 vertex-block entries over ℤ (Singular, ring integer) and print the strong Gröbner basis
of the elimination ideal with its leading coefficients, and its vdim over ℚ.
MODE p=P: the same elimination over 𝔽_P, with the points of the vertex-block variety over 𝔽_P (minAssGTZ).
Usage: python3 analysis/v2_58_vertex_integral.py QUIVER OUT_PREFIX MODE [--arrow 1,3,4]"""
import os, sys, subprocess, argparse
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import numpy as np, sympy as sp
import pcensus as pc, fold_ideal as fi
ap = argparse.ArgumentParser(); ap.add_argument('quiver'); ap.add_argument('out'); ap.add_argument('mode'); ap.add_argument('--arrow', default='1,3,4')
a = ap.parse_args()
ctx = pc.setup(a.quiver, 7)
d, kdim, basis, arrows = ctx['d'], ctx['kdim'], ctx['basis'], ctx['arrows']
M = np.asarray(ctx['M']); DL = np.asarray(ctx['DL'])
Kmat = sp.Matrix(((np.asarray(ctx['Kmat']) + 3) % 7 - 3).tolist())
comp = [list(basis).index(b) for b in fi.PAIRS_IN_SHORT_C] if a.quiver == 'pairs_in_short' else list(ctx['Cidx'])
B = sp.Matrix.hstack(Kmat, sp.Matrix([[1 if i == j else 0 for j in comp] for i in range(d)]))
assert abs(B.det()) == 1
Kl = B.inv()[:kdim, :]                                   # integral left inverse on K
vid = [i for i, b in enumerate(basis) if b[0] == 'v']
ks = sp.Matrix(kdim, kdim, lambda i, j: sp.Symbol(f'k{i}{j}'))
def mul(x, y): return sp.Matrix([sp.expand(sum(int(M[i, j, m]) * x[i] * y[j] for i in range(d) for j in range(d) if M[i, j, m])) for m in range(d)])
eqs = list((ks * ks - ks).applyfunc(sp.expand)); I5 = sp.eye(kdim)
for i in range(kdim):
    for j in range(kdim):
        for Pm, target in ((ks, 'img'), (I5 - ks, 'ker')):
            m = mul(Kmat * Pm[:, i], Kmat * Pm[:, j])
            kap = (Kl * m).applyfunc(sp.expand)
            assert (Kmat * kap - m).applyfunc(sp.expand) == sp.zeros(d, 1)
            eqs += list(((ks * kap - kap) if target == 'img' else ks * kap).applyfunc(sp.expand))
ys = []
for ar in [int(x) for x in a.arrow.split(',') if x.strip()]:
    s_, t_ = arrows[ar]; yv = sp.symbols(f'y{ar}_0:{d}'); ys += list(yv)
    Pe = lambda v: Kmat * ks[:, vid.index(v)]
    lhs = sum((yv[l] * sp.Matrix(DL[l].tolist()) for l in range(d)), sp.zeros(d, d))
    eqs += list((lhs - Pe(s_) * Pe(t_).T).applyfunc(sp.expand))
eqs = sorted({e for e in eqs if e != 0}, key=str)
assert all(c.is_integer for e in eqs for c in sp.Poly(e, *(list(ks) + ys)).coeffs())
allv = list(ks) + ys; keep = [ks[i, j] for i in range(len(vid)) for j in range(len(vid))]
elim = [v for v in allv if v not in keep]
ch = 'integer' if a.mode == 'int' else a.mode.split('=')[1]
src = [f'LIB "primdec.lib"; ring r = {ch},({",".join(map(str, allv))}),dp;',
       f"ideal I = {','.join(str(e).replace('**', '^') for e in eqs)};",
       f"ideal E = eliminate(I, {'*'.join(map(str, elim))});",
       f"ring rv = {ch},({','.join(map(str, keep))}),dp; option(redSB); ideal E2 = std(imap(r, E));"]
if a.mode == 'int':
    src += ['int ii; for (ii = 1; ii <= size(E2); ii++) { print("LC " + string(leadcoef(E2[ii])) + " LM " + string(leadmonom(E2[ii]))); }',
            f"ring rq = 0,({','.join(map(str, keep))}),dp; ideal E3 = std(imap(rv, E2));",
            'print("QDIM " + string(dim(E3)) + " VDIM " + string(vdim(E3)));',
            # radical certificate: the ideal Rq of the rational points of E3 (its associated primes), integral generators
            'list LL = minAssGTZ(E3); print("POINTS " + string(size(LL))); ideal Rq = LL[1]; int jj;',
            'for (jj = 2; jj <= size(LL); jj++) { Rq = intersect(Rq, LL[jj]); } Rq = std(Rq);',
            'for (jj = 1; jj <= nvars(rq); jj++) { Rq = Rq + ideal(var(jj)^2 - var(jj)); }',
            'setring rv; ideal Rz = imap(rq, Rq); int kk; int ok;',
            # every generator ρ of Rz: the least k ≤ 6 with ρ^k ∈ E (reduction by the strong GB E2 over ℤ, all LC 1)
            'for (jj = 1; jj <= size(Rz); jj++) { ok = 0; for (kk = 1; kk <= 6 && ok == 0; kk++) { if (reduce(Rz[jj]^kk, E2) == 0) { ok = kk; } } '
            'print("RAD " + string(jj) + " k " + string(ok) + " : " + string(Rz[jj])); }']
else:
    src += ['print("PDIM " + string(dim(E2)));', 'if (dim(E2) == 0) { print("VDIM " + string(vdim(E2))); list L = minAssGTZ(E2); print("POINTS " + string(size(L))); '
            'int i; for (i = 1; i <= size(L); i++) { print(std(L[i])); print("--"); } }']
src += ['print("DONE"); quit;']
open(a.out + '.sing', 'w').write("\n".join(src))
print(f"{a.quiver}: {len(eqs)} equations in {len(allv)} unknowns, mode {a.mode}", flush=True)
res = subprocess.run(["Singular", "-q", "--no-rc", a.out + '.sing'], capture_output=True, text=True, timeout=20000)
print(res.stdout, res.stderr[-2000:])
