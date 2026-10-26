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

def _radical_factors(e, maxterms):
    """Distinct non-constant irreducible factors of e over ℚ (monomials directly; factor_list only when e has at most
    maxterms terms), or None when e was not factored."""
    if len(e) == 1:
        (mon, _), = e.items()
        R = e.ring
        return [R.gens[i] for i, k in enumerate(mon) if k]
    if len(e) > maxterms: return None
    hit = _FCACHE.get(e)
    if hit is not None: return hit
    R = e.ring
    vi = [i for i, g in enumerate(e.degrees()) if g]                  # factor in ℚ[vars of e]: sympy's dense
    key = tuple(vi)                                                   # factorisation is over every ring generator
    if key not in _SUBRINGS: _SUBRINGS[key] = ring(','.join(str(R.gens[i]) for i in vi), QQ)[0]
    Rs = _SUBRINGS[key]
    es = Rs.from_dict({tuple(m[i] for i in vi): c for m, c in e.items()})
    _, fl = es.factor_list()
    out = []
    for f, _ in fl:
        if f.is_ground: continue
        full = [0] * R.ngens
        out.append(R.from_dict({tuple(_scatter(full, vi, m)): c for m, c in f.items()}))
    _FCACHE[e] = out
    return out

_FCACHE, _SUBRINGS = {}, {}
def _scatter(full, vi, m):
    out = list(full)
    for i, k in zip(vi, m): out[i] = k
    return out

def branch_presolve(S, batch=None, maxterms=None, log=None):
    """v2.57: staged_presolve with case splits. The stratum is covered by finitely many leaves, each (remaining
    equations, substitutions) as returned by staged_presolve; their varieties' union is the stratum.
    After each presolve, every equation that factors (monomials always; others with ≤ maxterms terms) is replaced by
    its radical part (a product of distinct irreducible factors), and if it has ≥ 2 factors f₁⋯f_k the branch splits
    into branches with f_i added. The zero set is the union of the branches' zero sets, so the irreducible
    components of the stratum are among the maximal irreducible components of the leaves.
    Branches are disjoint: splitting f₁⋯f_k = 0 gives branches f_i = 0, f_1 ≠ 0, …, f_{i−1} ≠ 0, and a known-nonzero
    factor is divided out of later equations (a stronger equation, so every leaf variety still lies in the stratum).
    Returns a list of (remaining equations, substitutions, nonzero polynomials); the stratum is the union of the leaf
    varieties V(remaining) (the inequations only prune: the union with them is already the whole stratum)."""
    R, P0, ctx = S['R'], S['P'], S['ctx']
    d = ctx['d']; M = np.asarray(ctx['M'], dtype=np.int64); DL = np.asarray(ctx['DL'], dtype=np.int64)
    nzM = [(i, j, k, int(M[i, j, k])) for i in range(d) for j in range(d) for k in range(d) if M[i, j, k]]
    batch = batch or int(os.environ.get('FOLD_BATCH', str(d)))
    maxterms = maxterms or int(os.environ.get('FOLD_FACTOR_TERMS', '12'))
    PS = presolve_markowitz if os.environ.get('PRESOLVE_RULE') == 'markowitz' else presolve_fast
    pairs = [(i, j) for i in range(d) for j in range(d)]
    def apply(P_, new):
        for x, val in new: P_ = [[e.compose(x, val) for e in row] for row in P_]
        return P_
    def mul(x, y):
        out = [R(0)] * d
        for i, j, k, c in nzM:
            if x[i] and y[j]: out[k] += c * x[i] * y[j]
        return out
    def matvec(A, v): return [sum((A[r][c] * v[c] for c in range(d) if A[r][c] and v[c]), R(0)) for r in range(d)]
    mon = lambda f: f * (1 / f.LC)
    def settle(eqs, P_, subs, nz):
        """presolve + radical replacement + division by known-nonzero factors, until stable. Returns
        ('leaf', eqs, P_, subs, nz), ('split', eqs, P_, subs, nz, factors) or ('empty',). nz: set of monic
        polynomials known to be nonzero on this branch."""
        eqs = set(eqs)
        while True:
            eqs.discard(R(0))
            if any(e.is_ground for e in eqs): return ('empty',)
            pend, new = PS(eqs, R)
            subs = subs + new; P_ = apply(P_, new)
            if any(e.is_ground and e != 0 for e in pend): return ('empty',)
            if new:
                nz2 = set()
                for f in nz:
                    for x, val in new: f = f.compose(x, val)
                    if f == 0: return ('empty',)
                    if not f.is_ground: nz2.add(mon(f))
                nz = nz2
            changed = bool(new); best = None; out = set()
            for e in pend:
                fs = _radical_factors(e, maxterms)
                if fs is None: out.add(e); continue
                fs = [mon(f) for f in fs]
                kept = [f for f in fs if f not in nz]
                if not kept: return ('empty',)
                if len(kept) == 1:
                    f = kept[0]
                    if f != mon(e): changed = True
                    out.add(f)
                else:
                    if len(kept) < len(fs): changed = True
                    prod = kept[0]
                    for f in kept[1:]: prod = prod * f
                    out.add(prod)
                    key = (len(kept), sum(len(f) for f in kept))
                    if best is None or key < best[0]: best = (key, kept, prod)
            eqs = out
            if not changed:
                if best is None: return ('leaf', sorted(eqs, key=str), P_, subs, nz)
                return ('split', sorted(eqs - {best[2]}, key=str), P_, subs, nz, best[1])   # f_i = 0 implies e = 0
    leaves = []
    stack = [(set(_equations(P0, R, d, nzM, DL, ('idem', 'delta'))), P0, [], 0, frozenset())]   # (eqs, P, subs, next pair, nz)
    nodes = 0
    while stack:
        eqs, P_, subs, b0, nz = stack.pop(); nodes += 1
        if log and os.environ.get('FOLD_TRACE'):
            log(f"node {nodes}: pair {b0}, {len(eqs)} eqs, {sum(len(e) for e in eqs)} terms, {len(R.gens) - len(subs)} unknowns")
        r = settle(eqs, P_, subs, nz)
        if r[0] == 'empty': continue
        if r[0] == 'split':
            _, eqs, P_, subs, nz, fs = r
            fs = sorted(fs, key=lambda f: (len(f), str(f)))
            for k in range(len(fs) - 1, -1, -1):      # branch k: f_k = 0 and f_0, …, f_{k-1} ≠ 0 (disjoint cover)
                stack.append((set(eqs) | {fs[k]}, P_, subs, b0, frozenset(nz) | frozenset(fs[:k])))
            continue
        _, eqs, P_, subs, nz = r
        if b0 >= len(pairs):
            leaves.append((eqs, subs, sorted(nz, key=str)))
            if log: log(f"leaf {len(leaves)}: {len(eqs)} equations, {len(R.gens) - len(subs)} unknowns, {len(nz)} nonzero; nodes {nodes}, stack {len(stack)}")
            continue
        Qm = [[(R(1) if r_ == c else R(0)) - P_[r_][c] for c in range(d)] for r_ in range(d)]
        new = set(eqs)
        for i, j in pairs[b0:b0 + batch]:
            u = mul([P_[r_][i] for r_ in range(d)], [P_[r_][j] for r_ in range(d)])
            w = mul([Qm[r_][i] for r_ in range(d)], [Qm[r_][j] for r_ in range(d)])
            Pu = matvec(P_, u); Pw = matvec(P_, w)
            for r_ in range(d): new.add(Pu[r_] - u[r_]); new.add(Pw[r_])
        stack.append((new, P_, subs, b0 + batch, nz))
    return leaves
