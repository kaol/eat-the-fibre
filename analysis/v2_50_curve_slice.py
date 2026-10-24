"""v2.50: the curve-slice argument on the divisor component D¹ (an affine-linear 6-space L), exactly over ℚ.

(33a) L ⊂ V: the fold equations vanish identically on the parametrisation of L (from the Singular prime).
(33b) 𝔤_L = {D ∈ Der(kQ) : D(A) ⊆ A, [D, P] ∈ L⃗ for all P ∈ L}; 𝔥′ ⊆ 𝔤_L; rank [𝔤_L, p].
(33c) rank [Der(kQ), p] versus dim Der(kQ).
(33d) dim(L⃗ ∩ [Der(kQ), p]).
(33e) coordinate directions v of L⃗ outside [Der(kQ), p]; writes p and an integer v for the Lean line C = p + s·v.
(33f) (33c), (33d) again at p + 2v.
Usage: python3 analysis/v2_50_curve_slice.py REPS.pkl REP_INDEX LINES_LOG T PRIME_FILE OUT_JSON"""
import os, sys, io, re, json, contextlib
import sympy as sp
from sympy.polys.matrices import DomainMatrix
HERE = os.path.dirname(os.path.abspath(__file__))
rf, ri, logf, TV, primef, outj = sys.argv[1], sys.argv[2], sys.argv[3], sp.Rational(sys.argv[4]), sys.argv[5], sys.argv[6]
# reuse v2_47_divisor.py up to the incidence variety: gives d, P0, h (𝔥′), M, eqs (the fold equations of V), Z, BA, ...
src = open(os.path.join(HERE, 'v2_47_divisor.py')).read().split("nh = len(h)")[0]
sys.argv = ['v2_47_divisor.py', rf, ri, logf, '/dev/null']
G = {'__name__': 'div', '__file__': os.path.join(HERE, 'v2_47_divisor.py')}
with contextlib.redirect_stdout(io.StringIO()):
    exec(src, G)
d, P0, h, M, eqs, Z = G['d'], G['P0'], G['h'], G['M'], list(G['eqs']), G['Z']
BAm, RK, Binv, BK, nA, nK = G['BAm'], G['RK'], G['Binv'], G['BK'], G['nA'], G['nK']
BKm = sp.Matrix.hstack(*BK)
QQ = sp.QQ
def rank_rows(rows): return DomainMatrix.from_list_sympy(len(rows), len(rows[0]), rows).convert_to(QQ).rank() if rows else 0
def nullsp(rows, n):
    Nm = DomainMatrix.from_list_sympy(len(rows), n, rows).convert_to(QQ).nullspace().to_Matrix()
    return [list(Nm.row(k)) for k in range(Nm.rows)]
zs = list(Z)
# --- L from the prime (affine-linear generators)
gens = [sp.sympify(x.replace('^', '**'), locals={str(z): z for z in zs})
        for x in open(primef).read().replace('\n', '').split(',') if x.strip()]
assert all(sp.Poly(g_, *zs).total_degree() <= 1 for g_ in gens), "prime is not affine-linear"
sol = sp.solve(gens, zs, dict=True); assert len(sol) == 1; sol = sol[0]
free = [z for z in zs if z not in sol]
us = sp.symbols(f'u0:{len(free)}')
par = {z: (sol[z] if z in sol else z) for z in zs}
par = {z: sp.expand(e.subs(dict(zip(free, us)))) for z, e in par.items()}
print(f"L: {len(gens)} affine-linear forms, dim {len(free)} (free {free})", flush=True)
# (33a)
bad = [e for e in eqs if sp.expand(e.subs(par)) != 0]
print(f"(33a) fold equations not vanishing on L: {len(bad)} of {len(eqs)}", flush=True)
ZL = Z.subs(par)                                   # nA × nK matrix, affine in u
def mat_of_Z(Zm): return BAm * Zm * RK
Ldir = [mat_of_Z(ZL.diff(u)) for u in us]          # basis of L⃗ (matrices)
# the point p on the v2.44 line, as a point of L
Nd = [__import__('ast').literal_eval(m.group(1)) for m in re.finditer(r"N(?: \(nonzero entries\):| =) (\{.*\})", open(logf).read())][0]
N = sp.zeros(d, d)
for (i, j), v in Nd.items(): N[i, j] = sp.Rational(v)
p = P0 + TV * N
Zp = Binv[:nA, :] * (p - P0) * BKm
up = sp.solve([ZL[i, j] - Zp[i, j] for i in range(nA) for j in range(nK)], us, dict=True)
assert len(up) == 1, "p is not on L"; up = up[0]
print(f"p = P0 + ({TV})·N lies on L at u = {[up.get(u, u) for u in us]}", flush=True)
# Der(kQ)
Mi = [[[int(M[i, j, m]) for m in range(d)] for j in range(d)] for i in range(d)]
rows = []
for i in range(d):
    for j in range(d):
        for m in range(d):
            row = [0] * (d * d)
            for k_ in range(d):
                if Mi[i][j][k_]: row[m * d + k_] += Mi[i][j][k_]
            for r_ in range(d):
                if Mi[r_][j][m]: row[r_ * d + i] -= Mi[r_][j][m]
                if Mi[i][r_][m]: row[r_ * d + j] -= Mi[i][r_][m]
            if any(row): rows.append(row)
Der = [sp.Matrix(d, d, v) for v in nullsp(rows, d * d)]
I = sp.eye(d); Q0 = I - P0
# (33b) 𝔤_L: x ∈ ℚ^{dim Der} with D = Σ x_k Der_k, (1 − P0) D B_A = 0 and the L⃗-normal forms of [D, P(u)] vanish for all u
xs = sp.symbols(f'x0:{len(Der)}')
Dx = sum((xs[k] * Der[k] for k in range(len(Der))), sp.zeros(d, d))
PL = (P0 + mat_of_Z(ZL)).applyfunc(sp.expand)
lin_forms = [sp.Poly(g_, *zs).as_expr() - sp.Poly(g_, *zs).as_expr().subs({z: 0 for z in zs}) for g_ in gens]   # linear parts
conds = list(Q0 * Dx * BAm)
Zdot = Binv[:nA, :] * (Dx * PL - PL * Dx) * BKm
for lf in lin_forms:
    e = sp.expand(lf.subs({zs[i * nK + j]: Zdot[i, j] for i in range(nA) for j in range(nK)}))
    conds += [c for c in sp.Poly(e, *us).coeffs()] if e != 0 else []
crow = [[sp.expand(c).coeff(x) for x in xs] for c in conds if sp.expand(c) != 0]
gL = [sum((v[k] * Der[k] for k in range(len(Der))), sp.zeros(d, d)) for v in nullsp(crow, len(Der))]
def in_span(vecs, w): return rank_rows([list(x) for x in vecs] + [list(w)]) == rank_rows([list(x) for x in vecs])
h_in = all(in_span(gL, Dh) for Dh in h)
def report(P, tag):
    bD = [list(D * P - P * D) for D in Der]
    rD = rank_rows(bD)
    rg = rank_rows([list(D * P - P * D) for D in gL])
    Lv = [list(X) for X in Ldir]
    inter = len(Lv) + rD - rank_rows(bD + Lv)
    out = [k for k, X in enumerate(Ldir) if not in_span([sp.Matrix(d, d, r) for r in bD], X)]
    print(f"[{tag}] rank [𝔤_L, P] = {rg}; rank [Der, P] = {rD} (dim Der {len(Der)}); dim(L⃗ ∩ [Der, P]) = {inter}; "
          f"directions of L⃗ outside [Der, P]: {[str(free[k]) for k in out]}", flush=True)
    return out
# (33c′) the common centraliser of L in Der(kQ): rank [Der, P] ≤ dim Der − dim c_L on all of L
PLb = PL.subs({u: 0 for u in us})
cL_rows = []
for X in [PLb] + Ldir:
    for e in (Dx * X - X * Dx).applyfunc(sp.expand):
        if e != 0: cL_rows.append([e.coeff(x) for x in xs])
cL = [sum((v[k] * Der[k] for k in range(len(Der))), sp.zeros(d, d)) for v in nullsp(cL_rows, len(Der))]
print(f"(33c′) common centraliser of L in Der: dim {len(cL)}, so rank [Der, P] ≤ {len(Der) - len(cL)} on L; "
      f"basis (nonzero entries): {[{(i, j): str(D[i, j]) for i in range(d) for j in range(d) if D[i, j] != 0} for D in cL]}", flush=True)
print(f"(33b) dim 𝔤_L = {len(gL)}; dim 𝔥′ = {len(h)}; 𝔥′ ⊆ 𝔤_L: {h_in}", flush=True)
out = report(p, "33b-d at p")
# (33e) integer direction for the Lean line
k0 = out[0]
V = Ldir[k0]; den = sp.ilcm(*[sp.fraction(x)[1] for x in V]) if any(V) else 1
V = (V * den).applyfunc(sp.nsimplify)
assert all(x.is_integer for x in p) and all(x.is_integer for x in V)
report(p + 2 * V, "33f at p + 2v")
json.dump({'d': d, 'p': [[int(x) for x in p.row(i)] for i in range(d)], 'v': [[int(x) for x in V.row(i)] for i in range(d)],
           'direction': str(free[k0])}, open(outj, 'w'))
print(f"(33e) C = p + s·v with v = {den}·∂/∂{free[k0]} written to {outj}", flush=True)
# (33g) random points of L: is rank [Der, P] = dim Der attained (then the open set U₁ of maximal orbit dimension is
# {rank = dim Der}), and is L⃗ ⊄ [Der, P] there?
import random
random.seed(int(os.environ.get('SEED', '50')))
for k in range(int(os.environ.get('NRAND', '4'))):
    uv = {u: sp.Integer(random.randint(-3, 3)) for u in us}
    report(PL.subs(uv), f"33g random u = {[int(uv[u]) for u in us]}")
# (33h) generic rank of [Der, P] on L, over the function field ℚ(u): a nonzero x(u) with Σ x_k(u) [Der_k, P(u)] = 0
# identically means rank ≤ dim Der − 1 everywhere on L. Computed as the nullspace over ℚ(u) of the (d², dim Der) matrix.
if os.environ.get('SYMB'):
    from sympy.polys.domains import QQ as _QQ
    K = _QQ.frac_field(*us)
    cols = [list((D * PL - PL * D).applyfunc(sp.expand)) for D in Der]
    Mt = [[cols[k][e] for k in range(len(Der))] for e in range(d * d)]
    Mt = [r for r in Mt if any(x != 0 for x in r)]
    DM = DomainMatrix.from_list_sympy(len(Mt), len(Der), Mt).convert_to(K)
    ns = DM.nullspace().to_Matrix()
    print(f"(33h) rank of [Der, P(u)] over ℚ(u) = {len(Der) - ns.rows}; generic centraliser dimension {ns.rows}", flush=True)
    for r_ in range(ns.rows):
        print(f"      x(u) = {[sp.factor(x) for x in ns.row(r_)]}", flush=True)
