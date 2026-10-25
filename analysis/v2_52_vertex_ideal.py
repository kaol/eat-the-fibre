"""v2.52: necessary conditions on the vertex block of a fold, as an ideal over ℚ (Singular).

Unknowns:
- P_K, the 5×5 matrix of P on K = ker Δ in an integer basis of K (P(K) ⊆ K, since ΔP = (P⊗P)Δ);
- optionally y = P(c) for the arrow c = 1→3 (pairs_in_short).
Conditions:
- P_K² = P_K;
- for basis vectors x, y of K: P(P x · P y) = P x · P y and P((1−P)x · (1−P)y) = 0. K is a subalgebra, so these are
  conditions on P_K;
- for each arrow x: s → t given with --arrow: Δ(y_x) = P(e_s) ⊗ P(e_t) for y_x = P(x), i.e. Σ_l y_l DL[l] = (P e_s)(P e_t)ᵀ.
The ideal is eliminated down to the 9 vertex-block unknowns. Reported: the dimension of the elimination ideal and, if
it is zero-dimensional, its rational points (minAssGTZ, each component a point).
Usage: python3 analysis/v2_52_vertex_ideal.py QUIVER OUT_PREFIX [--arrow INDEX]"""
import os, sys, subprocess, argparse
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
import numpy as np, sympy as sp
import pcensus as pc
ap = argparse.ArgumentParser(); ap.add_argument('quiver'); ap.add_argument('out'); ap.add_argument('--arrow', default='')   # comma-separated arrow indices
a = ap.parse_args()
p = 7                     # only to read the basis, M, DL and K (integer lifts); setup enumerates p^dim K points
ctx = pc.setup(a.quiver, p)
d, kdim, basis, arrows = ctx['d'], ctx['kdim'], ctx['basis'], ctx['arrows']
M = np.asarray(ctx['M']); DL = np.asarray(ctx['DL'])
Kmat = sp.Matrix(((np.asarray(ctx['Kmat']) + p // 2) % p - p // 2).tolist())        # integer basis of K (columns)
vid = [i for i, b in enumerate(basis) if b[0] == 'v']
assert [int(np.nonzero(np.asarray(Kmat[:, j]).astype(int).ravel())[0][0]) for j in range(len(vid))] == vid, "vertices first in K"
Kpinv = (Kmat.T * Kmat).inv() * Kmat.T                                                  # left inverse
ks = sp.Matrix(kdim, kdim, lambda i, j: sp.Symbol(f'k{i}{j}'))
def mul(x, y): return sp.Matrix([sp.expand(sum(int(M[i, j, m]) * x[i] * y[j] for i in range(d) for j in range(d) if M[i, j, m])) for m in range(d)])
eqs = list((ks * ks - ks).applyfunc(sp.expand))
I5 = sp.eye(kdim)
for i in range(kdim):
    for j in range(kdim):
        for Pm, target in ((ks, 'img'), (I5 - ks, 'ker')):
            m = mul(Kmat * Pm[:, i], Kmat * Pm[:, j])
            kap = (Kpinv * m).applyfunc(sp.expand)
            assert (Kmat * kap - m).applyfunc(sp.expand) == sp.zeros(d, 1), "K is not closed under products"
            eqs += list(((ks * kap - kap) if target == 'img' else ks * kap).applyfunc(sp.expand))
ys = []
for ar in [int(x) for x in a.arrow.split(',') if x.strip()]:
    s_, t_ = arrows[ar]
    yv = sp.symbols(f'y{ar}_0:{d}'); ys += list(yv)                  # y = P(arrow)
    Pe = lambda v: Kmat * ks[:, vid.index(v)]
    lhs = sum((yv[l] * sp.Matrix(DL[l].tolist()) for l in range(d)), sp.zeros(d, d))
    eqs += list((lhs - Pe(s_) * Pe(t_).T).applyfunc(sp.expand))
eqs = sorted({e for e in eqs if e != 0}, key=str)
allv = list(ks) + list(ys); keep = [ks[i, j] for i in range(len(vid)) for j in range(len(vid))]
elim = [v for v in allv if v not in keep]
src = f"""LIB "primdec.lib";
ring r = 0,({','.join(map(str, allv))}),dp;
ideal I = {','.join(str(e).replace('**', '^') for e in eqs)};
ideal E = eliminate(I, {'*'.join(map(str, elim))});
ring rv = 0,({','.join(map(str, keep))}),dp;
ideal E2 = std(imap(r, E));
print("elimination ideal: dim " + string(dim(E2)));
if (dim(E2) == 0) {{
  list L = minAssGTZ(E2);
  print("points " + string(size(L)));
  int i;
  for (i = 1; i <= size(L); i++) {{ print(std(L[i])); print("--"); }}
}}
quit;
"""
open(a.out + '.sing', 'w').write(src)
print(f"{a.quiver}: {len(eqs)} equations in {len(allv)} unknowns (arrow conditions: {a.arrow or 'none'})", flush=True)
res = subprocess.run(["Singular", "-q", "--no-rc", a.out + '.sing'], capture_output=True, text=True, timeout=6000)
print(res.stdout, res.stderr[-2000:])
