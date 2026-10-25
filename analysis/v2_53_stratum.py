"""v2.53: the stratum S_V of folds of Q₀ = 1⇉2⇇3 + 1→3 with a given vertex block V, as an ideal over ℚ (Singular).

Parametrisation: P(K) ⊆ K (K = ker Δ), so P is given by
- P_K (5×5 in the integer K basis, with its vertex block fixed to V: 16 unknowns), and
- the images of a complement C = span(a₂, b₂, c, cb₁, cb₂) (5 columns of 10: 50 unknowns).
P in the standard basis is [Kmat·P_K | P(C)]·[Kmat | C]⁻¹. Equations: all fold conditions:
- P² = P;
- P(P e_a · P e_b) = P e_a · P e_b and P((1−P)e_a · (1−P)e_b) = 0 for all basis a, b;
- Σ_m P[m,l] DL[m] = P DL[l] Pᵀ for all l.
Reports dim of the ideal and, with --components, the components' dimensions (minAssGTZ).
Usage: python3 analysis/v2_53_stratum.py BLOCK OUT_PREFIX [--components]   (BLOCK like 100,001,001: the rows of V)"""
import os, sys, subprocess, argparse
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
import numpy as np, sympy as sp
import pcensus as pc
ap = argparse.ArgumentParser(); ap.add_argument('block'); ap.add_argument('out'); ap.add_argument('--components', action='store_true')
a = ap.parse_args()
V = [[int(ch) for ch in r] for r in a.block.split(',')]
QUIVER = os.environ.get('QUIVER', 'pairs_in_short')     # v2.55: any quiver
ctx = pc.setup(QUIVER, 7)
d, kdim, basis = ctx['d'], ctx['kdim'], ctx['basis']
M = np.asarray(ctx['M']); DL = np.asarray(ctx['DL'])
p = 7
Kmat = sp.Matrix(((np.asarray(ctx['Kmat']) + p // 2) % p - p // 2).tolist())
idx = {b: i for i, b in enumerate(basis)}
Cb = ([idx[('p', (1,))], idx[('p', (3,))], idx[('p', (4,))], idx[('p', (4, 2))], idx[('p', (4, 3))]]
      if QUIVER == 'pairs_in_short' else list(ctx['Cidx']))      # v2.53's complement for pairs_in_short; else pcensus's
Cmat = sp.Matrix([[1 if i == j else 0 for j in Cb] for i in range(d)])
Bm = sp.Matrix.hstack(Kmat, Cmat); assert Bm.rank() == d
Bi = Bm.inv()
nv = sum(1 for b in basis if b[0] == "v")
NV = sum(1 for b in basis if b[0] == 'v')
ks = sp.Matrix(kdim, kdim, lambda i, j: sp.Integer(V[i][j]) if i < nv and j < nv else sp.Symbol(f'k{i}{j}'))
ys = sp.Matrix(d, len(Cb), lambda i, j: sp.Symbol(f'y{i}_{j}'))
P = (sp.Matrix.hstack(Kmat * ks, ys) * Bi).applyfunc(sp.expand)
I = sp.eye(d); Q = (I - P).applyfunc(sp.expand)
def mul(x, y): return sp.Matrix([sp.expand(sum(int(M[i, j, m]) * x[i] * y[j] for i in range(d) for j in range(d) if M[i, j, m])) for m in range(d)])
eqs = set()
for e in (P * P - P): eqs.add(sp.expand(e))
for i in range(d):
    for j in range(d):
        u = mul(P[:, i], P[:, j]); w = mul(Q[:, i], Q[:, j])
        for e in (P * u - u): eqs.add(sp.expand(e))
        for e in (P * w): eqs.add(sp.expand(e))
DLs = [sp.Matrix(DL[l].tolist()) for l in range(d)]
for l in range(d):
    for e in (sum((P[m, l] * DLs[m] for m in range(d)), sp.zeros(d, d)) - P * DLs[l] * P.T): eqs.add(sp.expand(e))
eqs = sorted((e for e in eqs if e != 0), key=str)
# presolve: repeatedly solve an equation that is linear in some variable with a constant coefficient, and substitute
nsub = 0
while True:
    pick = None
    for e in sorted(eqs, key=lambda e: len(e.free_symbols)):
        for v in sorted(e.free_symbols, key=str):
            pe = sp.Poly(e, v)
            if pe.degree() == 1 and pe.coeff_monomial(v).is_number:
                pick = (v, sp.solve(e, v)[0]); break
        if pick: break
    if not pick: break
    v, val = pick; nsub += 1
    eqs = sorted({sp.expand(e.subs(v, val)) for e in eqs} - {0}, key=str)
print(f"presolve: {nsub} variables eliminated by linear equations", flush=True)
allv = set(ys) | {x for x in ks if isinstance(x, sp.Symbol)}
vars_ = sorted(set().union(*[e.free_symbols for e in eqs]), key=str)
nfree = len(allv) - nsub - len(vars_)          # unknowns that occur in no remaining equation: free coordinates
print(f"remaining: {len(eqs)} equations in {len(vars_)} unknowns; {nfree} unknowns free (dimension offset)", flush=True)
comp = """list L = minAssGTZ(I);
print("components " + string(size(L)));
int i;
for (i = 1; i <= size(L); i++) { print("component " + string(i) + " dim " + string(dim(std(L[i])))); }
int dm = dim(G); int na; int j;
for (i = 1; i <= size(L); i++) { if (dim(std(L[i])) == dm) { def A = absPrimdecGTZ(L[i]); setring A; na = size(absolute_primes); print("top component " + string(i) + ": absolute primes " + string(na) + ", extension degrees " + string(absolute_primes[1][2])); setring r; kill A; } }""" if a.components else ""
src = f"""LIB "primdec.lib";
ring r = 0,({','.join(map(str, vars_))}),dp;
ideal I = {','.join(str(e).replace('**', '^') for e in eqs)};
ideal G = std(I);
print("dim " + string(dim(G)) + " (add the offset)");
{comp}
quit;
"""
open(a.out + '.sing', 'w').write(src)
print(f"V = {a.block}: {len(eqs)} equations in {len(vars_)} unknowns", flush=True)
res = subprocess.run(["Singular", "-q", "--no-rc", a.out + '.sing'], capture_output=True, text=True, timeout=int(os.environ.get('TMO', '3000')))
print(res.stdout.strip(), res.stderr[-1500:], flush=True)
