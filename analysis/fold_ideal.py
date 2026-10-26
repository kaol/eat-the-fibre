"""Fast construction of a vertex-block stratum's fold ideal (v2.55), with sparse polynomials (sympy PolyRing over ℚ)
instead of sympy expressions. v2_53_stratum.py and v2_53_orbits.py build the same ideal with expressions, too slowly
for d ≥ 12.

stratum(quiver, V, complement=None) → dict with
  R, xs     the polynomial ring ℚ[unknowns] and its generators (k_ij: P_K outside the vertex block; y_i_j: P on the
            complement, column j);
  P         d×d list of ring elements: P in the standard basis, [Kmat·P_K | P(C)]·[Kmat | C]⁻¹;
  eqs       the fold equations (P² = P; image and kernel closed under products; Δ-compatibility), nonzero, deduplicated;
  ctx       pcensus.setup context (p = 7, used only for the basis, M, DL and K).
presolve(eqs, R) → (remaining eqs, substitutions [(x, value)] in order): repeatedly solve an equation that is linear in
  some variable with a constant coefficient, and substitute. (Order: fewest variables, fewest terms, leading monomial;
  v2_53 used fewest variables then the string order, so the surviving coordinates can differ, not the variety.)
to_expr(f) → sympy expression."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
import numpy as np, sympy as sp
from sympy import QQ
from sympy.polys.rings import ring
import pcensus as pc

PAIRS_IN_SHORT_C = [('p', (1,)), ('p', (3,)), ('p', (4,)), ('p', (4, 2)), ('p', (4, 3))]   # v2.53's complement

def stratum(quiver, V, complement=None):
    p = 7
    ctx = pc.setup(quiver, p)
    d, kdim, basis = ctx['d'], ctx['kdim'], ctx['basis']
    M = np.asarray(ctx['M'], dtype=np.int64); DL = np.asarray(ctx['DL'], dtype=np.int64)
    idx = {b: i for i, b in enumerate(basis)}
    nv = sum(1 for b in basis if b[0] == 'v')
    Kmat = sp.Matrix(((np.asarray(ctx['Kmat']) + p // 2) % p - p // 2).tolist())
    if complement is None:
        complement = [idx[b] for b in PAIRS_IN_SHORT_C] if quiver == 'pairs_in_short' else list(ctx['Cidx'])
    Cmat = sp.Matrix([[1 if i == j else 0 for j in complement] for i in range(d)])
    Bi = sp.Matrix.hstack(Kmat, Cmat).inv()
    names = [f'k{i}{j}' for i in range(kdim) for j in range(kdim) if not (i < nv and j < nv)] + \
            [f'y{i}_{j}' for i in range(d) for j in range(len(complement))]
    R, *xs = ring(','.join(names), QQ)
    X = dict(zip(names, xs))
    PK = [[R(V[i][j]) if i < nv and j < nv else X[f'k{i}{j}'] for j in range(kdim)] for i in range(kdim)]
    cols = [[sum((R(int(Kmat[r, t])) * PK[t][j] for t in range(kdim)), R(0)) for r in range(d)] for j in range(kdim)]
    cols += [[X[f'y{r}_{j}'] for r in range(d)] for j in range(len(complement))]
    BiQ = [[QQ(int(sp.fraction(Bi[i, j])[0]), int(sp.fraction(Bi[i, j])[1])) for j in range(d)] for i in range(d)]
    P = [[sum((cols[k][r] * BiQ[k][c] for k in range(d) if BiQ[k][c]), R(0)) for c in range(d)] for r in range(d)]
    Q = [[(R(1) if r == c else R(0)) - P[r][c] for c in range(d)] for r in range(d)]
    nzM = [(i, j, k, int(M[i, j, k])) for i in range(d) for j in range(d) for k in range(d) if M[i, j, k]]
    def mul(x, y):
        out = [R(0)] * d
        for i, j, k, c in nzM:
            if x[i] and y[j]: out[k] += c * x[i] * y[j]
        return out
    def matvec(A, v): return [sum((A[r][c] * v[c] for c in range(d) if A[r][c] and v[c]), R(0)) for r in range(d)]
    col = lambda A, j: [A[r][j] for r in range(d)]
    eqs = set()
    for c in range(d):
        Pc = matvec(P, col(P, c))
        for r in range(d): eqs.add(Pc[r] - P[r][c])
    for i in range(d):
        for j in range(d):
            u = mul(col(P, i), col(P, j)); w = mul(col(Q, i), col(Q, j))
            Pu = matvec(P, u); Pw = matvec(P, w)
            for r in range(d): eqs.add(Pu[r] - u[r]); eqs.add(Pw[r])
    for l in range(d):                                            # Σ_m P[m,l] DL[m] = P DL[l] Pᵀ
        DLl = DL[l]
        PD = [[sum((P[a][x] * int(DLl[x, y]) for x in range(d) if DLl[x, y]), R(0)) for y in range(d)] for a in range(d)]
        for a in range(d):
            for b in range(d):
                lhs = sum((P[m][l] * int(DL[m][a, b]) for m in range(d) if DL[m][a, b]), R(0))
                rhs = sum((PD[a][y] * P[b][y] for y in range(d) if PD[a][y]), R(0))
                eqs.add(lhs - rhs)
    eqs.discard(R(0))
    return dict(R=R, xs=xs, names=names, P=P, eqs=sorted(eqs, key=str), ctx=ctx, complement=complement)

def presolve(eqs, R):
    eqs = set(eqs); subs = []
    while True:
        pick = None
        for e in sorted(eqs, key=lambda e: (sum(1 for g in e.degrees() if g), len(e), e.LM)):
            for x in sorted((R.gens[i] for i, g in enumerate(e.degrees()) if g == 1), key=str):
                c = e.coeff_wrt(x, 1)
                if c.is_ground and c != 0:
                    rest = e - c * x
                    if rest.degree(x) <= 0:
                        pick = (x, -rest * (1 / c.LC)); break
            if pick: break
        if not pick: break
        x, val = pick; subs.append(pick)
        eqs = {e.compose(x, val) for e in eqs}; eqs.discard(R(0))
    return sorted(eqs, key=str), subs

def to_expr(f): return f.as_expr()

def presolve_fast(eqs, R):
    """Same substitution rule as presolve (an equation linear in some variable with a constant coefficient is solved
    for it), but with an index variable → equations, so a substitution touches only the equations that contain the
    variable. Equations are visited in passes, each in (variables, terms, leading monomial) order of the pass start,
    until a pass makes no substitution. The variety is the same as presolve's; the surviving coordinates can differ.
    (v2.57: presolve composed every equation at every step, 1.4 M compositions and 856 s on a 236-unknown stratum.)"""
    from collections import defaultdict
    E = {i: e for i, e in enumerate(set(eqs)) if e}
    gens = R.gens
    def vars_of(e): return [i for i, g in enumerate(e.degrees()) if g]
    occ = defaultdict(set)
    for i, e in E.items():
        for v in vars_of(e): occ[v].add(i)
    subs = []; changed = True
    while changed:
        changed = False
        order = sorted(E, key=lambda i: (len(vars_of(E[i])), len(E[i]), E[i].LM))
        for i in order:
            e = E.get(i)
            if e is None: continue
            pick = None
            for v in sorted(vars_of(e), key=lambda v: str(gens[v])):
                if e.degree(gens[v]) != 1: continue
                c = e.coeff_wrt(gens[v], 1)
                if c.is_ground and c != 0:
                    rest = e - c * gens[v]
                    if rest.degree(gens[v]) <= 0: pick = (v, -rest * (1 / c.LC)); break
            if not pick: continue
            v, val = pick; x = gens[v]; subs.append((x, val)); changed = True
            for j in list(occ[v]):
                old = E.pop(j, None)
                if old is None: continue
                for w in vars_of(old): occ[w].discard(j)
                new = old.compose(x, val)
                if new:
                    E[j] = new
                    for w in vars_of(new): occ[w].add(j)
            occ.pop(v, None)
    out = sorted(set(E.values()), key=str)
    return out, subs
