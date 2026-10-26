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
    if os.environ.get('FOLD_STAGED', '1') != '0':
        return dict(R=R, xs=xs, names=names, P=P, eqs=None, ctx=ctx, complement=complement,
                    gen=lambda P_: _equations(P_, R, d, nzM, DL))
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

def _equations(P, R, d, nzM, DL, parts=('idem', 'delta', 'closure')):
    """The fold equations of P (d×d ring elements): 'idem' P² = P; 'delta' Δ-compatibility; 'closure' image and
    kernel closed under products (cubic, the expensive part)."""
    Q = [[(R(1) if r == c else R(0)) - P[r][c] for c in range(d)] for r in range(d)]
    def mul(x, y):
        out = [R(0)] * d
        for i, j, k, c in nzM:
            if x[i] and y[j]: out[k] += c * x[i] * y[j]
        return out
    def matvec(A, v): return [sum((A[r][c] * v[c] for c in range(d) if A[r][c] and v[c]), R(0)) for r in range(d)]
    col = lambda A, j: [A[r][j] for r in range(d)]
    eqs = set()
    if 'idem' in parts:
        for c in range(d):
            Pc = matvec(P, col(P, c))
            for r in range(d): eqs.add(Pc[r] - P[r][c])
    if 'delta' in parts:
        for l in range(d):
            DLl = DL[l]
            PD = [[sum((P[a][x] * int(DLl[x, y]) for x in range(d) if DLl[x, y]), R(0)) for y in range(d)] for a in range(d)]
            for a in range(d):
                for b in range(d):
                    lhs = sum((P[m][l] * int(DL[m][a, b]) for m in range(d) if DL[m][a, b]), R(0))
                    eqs.add(lhs - sum((PD[a][y] * P[b][y] for y in range(d) if PD[a][y]), R(0)))
    if 'closure' in parts:
        for i in range(d):
            for j in range(d):
                u = mul(col(P, i), col(P, j)); w = mul(col(Q, i), col(Q, j))
                Pu = matvec(P, u); Pw = matvec(P, w)
                for r in range(d): eqs.add(Pu[r] - u[r]); eqs.add(Pw[r])
    eqs.discard(R(0))
    return eqs

def stage1_presolve(S):
    """Only the first stage of staged_presolve: P² = P and Δ-compatibility, presolved. Returns (remaining
    equations, substitutions). The closure equations are then generated elsewhere (v2_55_pipeline with
    CLOSURE_IN_SINGULAR=1 builds them inside Singular from the reduced P)."""
    R, P, ctx = S['R'], S['P'], S['ctx']
    d = ctx['d']; M = np.asarray(ctx['M'], dtype=np.int64); DL = np.asarray(ctx['DL'], dtype=np.int64)
    nzM = [(i, j, k, int(M[i, j, k])) for i in range(d) for j in range(d) for k in range(d) if M[i, j, k]]
    return presolve_fast(_equations(P, R, d, nzM, DL, ('idem', 'delta')), R)

def staged_presolve(S, batch=None):
    """v2.57: build and presolve in stages, to keep the polynomials small on d = 18.
    1. P² = P and Δ-compatibility (quadratic), presolved; the substitutions are applied to P.
    2. The closure equations (cubic), generated for one batch of column pairs (i, j) at a time from the current,
       already reduced P. Each batch is presolved together with the equations still pending, and the new
       substitutions are applied to P before the next batch.
    The result generates the same ideal as presolving all equations at once: substitution commutes with generating
    equations, and every equation is generated, so no condition is dropped. Only the order of substitutions differs.
    Returns (remaining equations, substitutions in order)."""
    R, P, ctx = S['R'], S['P'], S['ctx']
    d = ctx['d']; M = np.asarray(ctx['M'], dtype=np.int64); DL = np.asarray(ctx['DL'], dtype=np.int64)
    nzM = [(i, j, k, int(M[i, j, k])) for i in range(d) for j in range(d) for k in range(d) if M[i, j, k]]
    batch = batch or int(os.environ.get('FOLD_BATCH', str(d)))
    PS = presolve_markowitz if os.environ.get('PRESOLVE_RULE') == 'markowitz' else presolve_fast
    pend, subs = PS(_equations(P, R, d, nzM, DL, ('idem', 'delta')), R)
    def apply(P_, new):
        for x, val in new: P_ = [[e.compose(x, val) for e in row] for row in P_]
        return P_
    P = apply(P, subs)
    pairs = [(i, j) for i in range(d) for j in range(d)]
    Q = lambda P_: [[(R(1) if r == c else R(0)) - P_[r][c] for c in range(d)] for r in range(d)]
    def mul(x, y):
        out = [R(0)] * d
        for i, j, k, c in nzM:
            if x[i] and y[j]: out[k] += c * x[i] * y[j]
        return out
    def matvec(A, v): return [sum((A[r][c] * v[c] for c in range(d) if A[r][c] and v[c]), R(0)) for r in range(d)]
    for b0 in range(0, len(pairs), batch):
        Qm = Q(P); eqs = set(pend)
        for i, j in pairs[b0:b0 + batch]:
            u = mul([P[r][i] for r in range(d)], [P[r][j] for r in range(d)])
            w = mul([Qm[r][i] for r in range(d)], [Qm[r][j] for r in range(d)])
            Pu = matvec(P, u); Pw = matvec(P, w)
            for r in range(d): eqs.add(Pu[r] - u[r]); eqs.add(Pw[r])
        eqs.discard(R(0))
        pend, new = PS(eqs, R)
        subs += new; P = apply(P, new)
    return pend, subs

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

def presolve_markowitz(eqs, R):
    """Like presolve_fast, but each step picks, among all equations, the substitution x = val with the smallest fill
    estimate (terms of val) × (equations containing x), the Markowitz rule of sparse elimination. Fewer, smaller
    compositions on strata where the first-come order makes the polynomials grow (v2.57, Q₁ heavy blocks)."""
    from collections import defaultdict
    E = {i: e for i, e in enumerate(set(eqs)) if e}
    gens = R.gens
    def vars_of(e): return [i for i, g in enumerate(e.degrees()) if g]
    occ = defaultdict(set)
    for i, e in E.items():
        for v in vars_of(e): occ[v].add(i)
    subs = []
    while True:
        best = None
        for i, e in E.items():
            for v in vars_of(e):
                if e.degree(gens[v]) != 1: continue
                c = e.coeff_wrt(gens[v], 1)
                if not (c.is_ground and c != 0): continue
                cost = (len(e) - 1) * len(occ[v])
                if best is None or cost < best[0]: best = (cost, i, v, c)
        if best is None: break
        _, i, v, c = best; e = E[i]; x = gens[v]
        rest = e - c * x
        if rest.degree(x) > 0: occ[v].discard(i); continue          # (cannot happen for degree 1 in x)
        val = -rest * (1 / c.LC); subs.append((x, val))
        for j in list(occ[v]):
            old = E.pop(j, None)
            if old is None: continue
            for w in vars_of(old): occ[w].discard(j)
            new = old.compose(x, val)
            if new:
                E[j] = new
                for w in vars_of(new): occ[w].add(j)
        occ.pop(v, None)
    return sorted(set(E.values()), key=str), subs
