"""v2.50 (exploratory, not pre-registered): polynomial invariants of 𝔥′ on V (7-dim, linear) and on L = D¹ (6-dim).
𝔥′ acts on the affine-linear space by the affine vector fields ξ_D(P) = [D, P] (D(A) ⊆ A). F is 𝔥′-invariant iff
ξ_D F = 0 for every basis element D. Solves for all invariant polynomials of degree ≤ DEG in the free coordinates.
Usage: python3 analysis/v2_50_invariants.py REPS.pkl REP_INDEX LINES_LOG T V_PRIME L_PRIME [DEG]"""
import os, sys, io, itertools, contextlib
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
rf, ri, logf, T, vprime, lprime = sys.argv[1:7]; DEG = int(sys.argv[7]) if len(sys.argv) > 7 else 3
src = open(os.path.join(HERE, 'v2_50_curve_slice.py')).read().split("# (33a)")[0]
sys.argv = ['v2_50_curve_slice.py', rf, ri, logf, T, lprime, '/dev/null']
G = {'__name__': 'cs', '__file__': os.path.join(HERE, 'v2_50_curve_slice.py')}
with contextlib.redirect_stdout(io.StringIO()):
    exec(src, G)
d, P0, h, Z, BAm, RK, Binv, BKm, nA, nK = (G[k] for k in ('d', 'P0', 'h', 'Z', 'BAm', 'RK', 'Binv', 'BKm', 'nA', 'nK'))
zs = list(Z)
def param(primef, pfx):
    gens = [sp.sympify(x.replace('^', '**'), locals={str(z): z for z in zs})
            for x in open(primef).read().split('\n') if x.strip() and not x.startswith(('dim', 'comp'))
            for x in x.split(',') if x.strip()]
    sol = sp.solve(gens, zs, dict=True)[0]
    free = [z for z in zs if z not in sol]
    ws = sp.symbols(f'{pfx}0:{len(free)}')
    par = {z: sp.expand((sol[z] if z in sol else z).subs(dict(zip(free, ws)))) for z in zs}
    return free, ws, par
from sympy.polys.matrices import DomainMatrix
def span_basis(mats):
    rows = [list(m) for m in mats]
    R = DomainMatrix.from_list_sympy(len(rows), len(rows[0]), rows).convert_to(sp.QQ).rref()[0].to_Matrix()
    return [sp.Matrix(d, d, list(R.row(i))) for i in range(R.rows) if any(R.row(i))]
hd = span_basis([A_ * B_ - B_ * A_ for A_ in h for B_ in h])
print(f"dim 𝔥′ = {len(h)}, dim [𝔥′, 𝔥′] = {len(hd)}")
def fields_of(alg, P, free):
    out = []
    for D in alg:
        Zdot = (Binv[:nA, :] * (D * P - P * D) * BKm).applyfunc(sp.expand)
        out.append([Zdot[zs.index(z) // nK, zs.index(z) % nK] for z in free])
    return out
def invariants(fields, ws, deg, const=False):
    mons = ([sp.Integer(1)] if const else []) + [sp.Mul(*m) for k in range(1, deg + 1) for m in itertools.combinations_with_replacement(ws, k)]
    cs = sp.symbols(f'c0:{len(mons)}')
    F = sum(c * m for c, m in zip(cs, mons))
    eqs = []
    for fld in fields:
        e = sp.expand(sum(fld[i] * sp.diff(F, ws[i]) for i in range(len(ws))))
        if e != 0: eqs += sp.Poly(e, *ws).coeffs()
    (vec,) = list(sp.linsolve(eqs, cs))
    params = sorted(set().union(*[x.free_symbols for x in vec]) & set(cs), key=str)
    return [sp.expand(sum(vc.subs({q: (1 if q == pq else 0) for q in params}) * m for vc, m in zip(vec, mons))) for pq in params]
for tag, primef, pfx in (("V", vprime, 'w'), ("L", lprime, 'u')):
    free, ws, par = param(primef, pfx)
    ZP = Z.subs(par); P = P0 + BAm * ZP * RK
    fields = []
    for D in h:
        Zdot = (Binv[:nA, :] * (D * P - P * D) * BKm).applyfunc(sp.expand)
        fields.append([Zdot[zs.index(z) // nK, zs.index(z) % nK] for z in free])   # components along the free coords
    mons = [sp.Mul(*m) for k in range(1, DEG + 1) for m in itertools.combinations_with_replacement(ws, k)]
    cs = sp.symbols(f'c0:{len(mons)}')
    F = sum(c * m for c, m in zip(cs, mons))
    eqs = []
    for fld in fields:
        e = sp.expand(sum(fld[i] * sp.diff(F, ws[i]) for i in range(len(ws))))
        if e != 0: eqs += sp.Poly(e, *ws).coeffs()
    # semi-invariants: invariants of the derived algebra [𝔥′, 𝔥′] (computed below), then the 𝔥′-action on them
    sol = sp.linsolve(eqs, cs)
    (vec,) = list(sol)
    params = sorted(set().union(*[x.free_symbols for x in vec]) & set(cs), key=str)
    basis = [sp.factor(sp.expand(sum(vc.subs({q: (1 if q == pq else 0) for q in params}) * m for vc, m in zip(vec, mons))))
             for pq in params]
    print(f"[{tag}] free coordinates {[str(z) for z in free]} = {list(ws)}; invariants of degree ≤ {DEG}: {len(basis)}")
    for b in basis: print(f"    {b}")
    # semi-invariants of degree ≤ DEG: [𝔥′, 𝔥′]-invariants, and the matrices of 𝔥′ on that space
    fd = fields_of(hd, P, free)
    sinv = invariants(fd, ws, DEG, const=True)
    print(f"[{tag}] [𝔥′, 𝔥′]-invariants of degree ≤ {DEG} (with constants): {len(sinv)}")
    for b in sinv: print(f"    {sp.factor(b)}")
    if sinv:
        cvars = sp.symbols(f'k0:{len(sinv)}')
        for k_, fld in enumerate(fields):
            imgs = [sp.expand(sum(fld[i] * sp.diff(b, ws[i]) for i in range(len(ws)))) for b in sinv]
            # express each image in the basis sinv
            Mrows = []
            for im in imgs:
                sol_ = sp.solve(sp.Poly(im - sum(c * b for c, b in zip(cvars, sinv)), *ws).coeffs(), cvars, dict=True)
                Mrows.append([sol_[0].get(c, c) for c in cvars] if sol_ else None)
            print(f"    action of basis element {k_} of 𝔥′ on them (rows = images): {Mrows}")
    # (36) characters of candidate semi-invariants (read off the matrices above) and the resulting rational invariant
    # σ: the sign with which u1 enters the semi-invariant b − 1 (A: +1, B: −1; mirror), chosen so its character is constant
    def is_semi(sv): return all(not sp.simplify(sum(fld[i] * sp.diff(sv, ws[i]) for i in range(len(ws))) / sv).free_symbols for fld in fields)
    SIG = next(sg for sg in (1, -1) if tag != "L" or is_semi(ws[0] * ws[5] + sg * ws[1] - 1))
    cand = {"V": {'x6-1': ws[2] - ws[1] - 1, 'x8+1': ws[4] - ws[5] + ws[6] + 1, 'w0': ws[0]} if tag == "V" else {},
            "L": {'b-1': ws[0] * ws[5] + SIG * ws[1] - 1, 'u0(c+1)': ws[0] * (ws[3] - ws[4] + ws[5] + 1), 'u0': ws[0],
                  'c+1': ws[3] - ws[4] + ws[5] + 1}}[tag]
    for nm_, sv in cand.items():
        if nm_ == 'b-1': nm_ += f' (σ = {SIG})'
        ch = [sp.simplify(sum(fld[i] * sp.diff(sv, ws[i]) for i in range(len(ws))) / sv) for fld in fields]
        print(f"    semi-invariant {nm_} = {sv}: character {ch}")
    if tag == "L":
        kappa = cand['u0(c+1)'] / cand['b-1']
        ok = all(sp.simplify(sum(fld[i] * sp.diff(kappa, ws[i]) for i in range(len(ws)))) == 0 for fld in fields)
        print(f"    κ = {kappa}: 𝔥′-invariant: {ok}")
