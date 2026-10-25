"""v2.53: census-orbit dimension on each irreducible component of a vertex-block stratum S_V of Q₀ = 1⇉2⇇3 + 1→3.

For each component Z (a ℚ-prime from minAssGTZ, after the same linear presolve as v2_53_stratum.py):
- W_Z: the class weightings w with P compatible with Δ_w = Σ_c w_c Δ_c on all of Z. The class defects of P are reduced
  modulo a Gröbner basis of the prime (Singular), and W_Z is the kernel of w ↦ Σ_c w_c NF(defect_c).
- 𝔥_Z = Der_Δ ⊕ {Σ_c w_c E_c : w ∈ W_Z}. Its connected group (Δ-compatible automorphisms times the torus of W_Z)
  preserves Z: it preserves folds with W ⊆ W_P and fixes vertices, so it preserves S_V and permutes components.
- A rational point of Z: the variables of a maximal independent set get random small integers, the rest is solved.
  Retried until a rational solution exists. P is rebuilt and checked to be a fold exactly.
- rank [𝔥_Z, P] (the census-orbit dimension at the point) and rank [Der(kQ), P] (the Aut-orbit dimension), against
  dim Z.
rank [𝔥_Z, P] = dim Z means the census orbit of P is open in Z. A random point gives a lower bound for the generic
rank, so equality at one point already proves that Z has an open orbit.
Usage: python3 analysis/v2_53_orbits.py BLOCK OUT_PREFIX [--seed S] [--tries N]"""
import os, sys, subprocess, argparse, random, re, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
import numpy as np, sympy as sp
from sympy.polys.matrices import DomainMatrix
import pcensus as pc
ap = argparse.ArgumentParser(); ap.add_argument('block'); ap.add_argument('out')
ap.add_argument('--seed', type=int, default=53); ap.add_argument('--tries', type=int, default=40)
a = ap.parse_args(); random.seed(a.seed)
V = [[int(ch) for ch in r] for r in a.block.split(',')]
p = 7
QUIVER = os.environ.get('QUIVER', 'pairs_in_short')     # v2.55: any quiver
ctx = pc.setup(QUIVER, p)
d, kdim, basis, arrows = ctx['d'], ctx['kdim'], ctx['basis'], ctx['arrows']
M = np.asarray(ctx['M']); DL = np.asarray(ctx['DL'])
Kmat = sp.Matrix(((np.asarray(ctx['Kmat']) + p // 2) % p - p // 2).tolist())
idx = {b: i for i, b in enumerate(basis)}
Cb = ([idx[('p', (1,))], idx[('p', (3,))], idx[('p', (4,))], idx[('p', (4, 2))], idx[('p', (4, 3))]]
      if QUIVER == 'pairs_in_short' else list(ctx['Cidx']))      # v2.53's complement for pairs_in_short; else pcensus's
Cmat = sp.Matrix([[1 if i == j else 0 for j in Cb] for i in range(d)])
Bi = sp.Matrix.hstack(Kmat, Cmat).inv()
NV = sum(1 for b in basis if b[0] == 'v')
ks = sp.Matrix(kdim, kdim, lambda i, j: sp.Integer(V[i][j]) if i < NV and j < NV else sp.Symbol(f'k{i}{j}'))
ys = sp.Matrix(d, len(Cb), lambda i, j: sp.Symbol(f'y{i}_{j}'))
P = (sp.Matrix.hstack(Kmat * ks, ys) * Bi).applyfunc(sp.expand)
I = sp.eye(d); Q = (I - P).applyfunc(sp.expand)
def mul(x, y): return sp.Matrix([sp.expand(sum(int(M[i, j, m]) * x[i] * y[j] for i in range(d) for j in range(d) if M[i, j, m])) for m in range(d)])
def fold_eqs(P, Q):
    eqs = set(sp.expand(e) for e in (P * P - P))
    for i in range(d):
        for j in range(d):
            u = mul(P[:, i], P[:, j]); w = mul(Q[:, i], Q[:, j])
            eqs |= {sp.expand(e) for e in (P * u - u)} | {sp.expand(e) for e in (P * w)}
    DLs = [sp.Matrix(DL[l].tolist()) for l in range(d)]
    for l in range(d):
        eqs |= {sp.expand(e) for e in (sum((P[m, l] * DLs[m] for m in range(d)), sp.zeros(d, d)) - P * DLs[l] * P.T)}
    return sorted((e for e in eqs if e != 0), key=str)
eqs = fold_eqs(P, Q)
subs = []                                             # presolve, as in v2_53_stratum.py; the substitutions are kept
while True:
    pick = None
    for e in sorted(eqs, key=lambda e: len(e.free_symbols)):
        for v in sorted(e.free_symbols, key=str):
            pe = sp.Poly(e, v)
            if pe.degree() == 1 and pe.coeff_monomial(v).is_number:
                pick = (v, sp.solve(e, v)[0]); break
        if pick: break
    if not pick: break
    v, val = pick; subs.append(pick)
    eqs = sorted({sp.expand(e.subs(v, val)) for e in eqs} - {0}, key=str)
allv = sorted(set(ys) | {x for x in ks if isinstance(x, sp.Symbol)}, key=str)
elimd = {v for v, _ in subs}
rest = [v for v in allv if v not in elimd]            # the coordinates left after presolve (some may be free)
def full_point(vals):
    """values of `rest` → values of all unknowns (back-substitution in reverse order)"""
    env = dict(vals)
    for v, val in reversed(subs): env[v] = sp.nsimplify(val.subs(env))
    return env
# class-weighted Δ: defect_c(P) = Σ_m P[m,l] DL_c[m] − P DL_c[l] Pᵀ for each class c, with P in terms of `rest`
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
Pr = P
for v, val in subs: Pr = Pr.subs(v, val)
Pr = Pr.applyfunc(sp.expand)
defects = []
for c in range(len(cl)):
    ent = []
    for l in range(d):
        E = (sum((Pr[m, l] * DLc[c][m] for m in range(d)), sp.zeros(d, d)) - Pr * DLc[c][l] * Pr.T).applyfunc(sp.expand)
        ent += list(E)
    defects.append(ent)
nz = [k for k in range(len(defects[0])) if any(defects[c][k] != 0 for c in range(len(cl)))]
# Singular: primes, dims, independent sets, and normal forms of the class defects modulo each prime
vs = [str(v) for v in rest]
dfl = ",".join(str(defects[c][k]).replace('**', '^') for c in range(len(cl)) for k in nz)
src = f"""LIB "primdec.lib";
ring r = 0,({','.join(vs)}),dp;
ideal I = {','.join(str(e).replace('**', '^') for e in eqs) if eqs else '0'};
list L = minAssGTZ(I);
ideal DF = {dfl if dfl else '0'};
int i; int j;
for (i = 1; i <= size(L); i++) {{
  ideal Gi = std(L[i]);
  print("PRIME " + string(i) + " DIM " + string(dim(Gi)));
  print("GENS " + string(Gi));
  intvec iv = indepSet(Gi);
  print("INDEP " + string(iv));
  ideal NFd = reduce(DF, Gi);
  print("NFD " + string(NFd));
  kill iv; kill NFd; kill Gi;
}}
quit;
"""
open(a.out + '.sing', 'w').write(src)
res = subprocess.run(["Singular", "-q", "--no-rc", a.out + '.sing'], capture_output=True, text=True, timeout=6000)
out = res.stdout
loc = {v: sp.Symbol(v) for v in vs}
comps = []
for blk in out.split("PRIME ")[1:]:
    m = re.match(r"(\d+) DIM (-?\d+)", blk); i_, dim_ = int(m.group(1)), int(m.group(2))
    gens = [sp.sympify(g.replace('^', '**'), locals=loc) for g in re.search(r"GENS (.*)", blk).group(1).split(',') if g.strip()]
    indep = [int(x) for x in re.search(r"INDEP (.*)", blk).group(1).split(',')]
    nf = [sp.sympify(g.replace('^', '**'), locals=loc) for g in re.search(r"NFD (.*)", blk).group(1).split(',')]
    comps.append((i_, dim_, gens, indep, nf))
# Der(kQ) and Der_Δ over ℚ (as in v2_46_slice.py)
QQ = sp.QQ
def nullsp(rows, n):
    if not rows: return [sp.Matrix([1 if i == j else 0 for i in range(n)]) for j in range(n)]
    Nm = DomainMatrix.from_list_sympy(len(rows), n, rows).convert_to(QQ).nullspace().to_Matrix()
    return [Nm.row(k).T for k in range(Nm.rows)]
def rank(vecs): return DomainMatrix.from_list_sympy(len(vecs), len(vecs[0]), [list(v) for v in vecs]).convert_to(QQ).rank() if vecs else 0
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
cond = []
for D in Der:
    vv = []
    for l in range(d):
        vv += list(sum((D[m, l] * DLs[m] for m in range(d)), sp.zeros(d, d)) - D * DLs[l] - DLs[l] * D.T)
    cond.append(vv)
cx = nullsp([list(r) for r in sp.Matrix(cond).T.tolist()], len(Der))
DerD = [sum((c[k] * Der[k] for k in range(len(Der))), sp.zeros(d, d)) for c in cx]
Ec = [sp.diag(*[0 if b[0] == 'v' else sum(1 for x in b[1] if x in c) for b in basis]) for c in cl]
print(f"V = {a.block}: {len(comps)} components; dim Der = {len(Der)}, dim Der_Δ = {len(DerD)}, classes {len(cl)}", flush=True)
for i_, dim_, gens, indep, nf in comps:
    # W_Z: w with Σ_c w_c NF_c = 0 identically (as polynomials)
    n_ = len(nz); wv = sp.symbols(f'w0:{len(cl)}')
    lin = []
    for k in range(n_):
        e = sp.expand(sum(wv[c] * nf[c * n_ + k] for c in range(len(cl))))
        if e != 0:
            for co in (sp.Poly(e, *[sp.Symbol(v) for v in vs]).coeffs() if vs else [e]): lin.append([sp.expand(co).coeff(w) for w in wv])
    W = nullsp(lin, len(cl))
    h = DerD + [sum((w[c] * Ec[c] for c in range(len(cl))), sp.zeros(d, d)) for w in W]
    found = None; best = (-1, -1); npts = 0
    for t in range(a.tries):
        vals = {sp.Symbol(vs[k]): sp.Integer(random.randint(-4, 4)) for k in range(len(vs)) if indep[k] == 1}
        g2 = [sp.expand(g.subs(vals)) for g in gens]; g2 = [g for g in g2 if g != 0]
        unk = sorted(set().union(*[g.free_symbols for g in g2]), key=str) if g2 else []
        sols = sp.solve(g2, unk, dict=True) if g2 else [{}]
        for so in sols:
            if not all(x.is_rational for x in so.values()) or len(so) < len(unk): continue
            pt = dict(vals); pt.update(so)
            for v in vs:
                sv = sp.Symbol(v)
                if sv not in pt: pt[sv] = sp.Integer(random.randint(-4, 4))   # coordinates absent from the prime are free
            env = full_point({sp.Symbol(v): pt[sp.Symbol(v)] for v in vs})
            Pn = P.subs(env)
            if fold_eqs(Pn, I - Pn): continue
            npts += 1
            rh_ = rank([list(D * Pn - Pn * D) for D in h])
            if rh_ > best[0]: best = (rh_, rank([list(D * Pn - Pn * D) for D in Der])); found = Pn
            break
        if found is not None and best[0] == dim_: break
    if found is None:
        print(f"  component {i_}: dim {dim_}, W_Z = {[list(w) for w in W]}; no rational point found in {a.tries} tries", flush=True); continue
    rh, rD = best
    print(f"  component {i_}: dim {dim_}; W_Z = {[list(w) for w in W]}; best of {npts} rational points: rank [𝔥_Z, P] = {rh}, "
          f"rank [Der, P] = {rD}; open census orbit: {rh == dim_}", flush=True)
