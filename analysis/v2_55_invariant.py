"""v2.55: an explicit census invariant on a component where the pipeline fails (no open orbit).

Reuses v2_55_pipeline.py up to the component loop. For component COMP: solve the prime's generators for the dependent
coordinates in terms of an independent set (assumes the component is a graph over that set, as in the control's
failures). Then compute the vector fields of 𝔥_Z = Der_Δ ⊕ E_{W_Z} in the independent coordinates, and search for
- polynomial invariants of degree ≤ DEG, and
- ratios of semi-invariants of degree ≤ DEG with the same character (an invariant rational function).
A nonconstant invariant proves that the connected census group has no open orbit on Z: a one-parameter (or larger)
family of orbits.
Usage: python3 analysis/v2_55_invariant.py QUIVER BLOCK COMP OUT_PREFIX [DEG]"""
import os, sys, io, itertools, contextlib
HERE = os.path.dirname(os.path.abspath(__file__))
quiver, block, comp, outp = sys.argv[1:5]; DEG = int(sys.argv[5]) if len(sys.argv) > 5 else 2
src = open(os.path.join(HERE, 'v2_55_pipeline.py')).read().split("for i_, dim_, gens, indep, nf in comps:")[0]
sys.argv = ['v2_55_pipeline.py', quiver, block, outp]
G = {'__name__': 'pl', '__file__': os.path.join(HERE, 'v2_55_pipeline.py')}
with contextlib.redirect_stdout(io.StringIO()):
    exec(src, G)
sp, comps, Sv, vs, cl, nz, DerD, Ec, d, nullsp = (G[k] for k in ('sp', 'comps', 'Sv', 'vs', 'cl', 'nz', 'DerD', 'Ec', 'd', 'nullsp'))
PrE = sp.Matrix([[e.as_expr() for e in row] for row in G['Pr']])   # (the pipeline keeps P in ring arithmetic since v2.57)
i_, dim_, gens, indep, nf = [c for c in comps if c[0] == int(comp)][0]
n_ = len(nz); wv = sp.symbols(f'w0:{len(cl)}'); lin = []
for k in range(n_):
    e = sp.expand(sum(wv[c] * nf[c * n_ + k] for c in range(len(cl))))
    if e != 0:
        for co in sp.Poly(e, *Sv).coeffs(): lin.append([sp.expand(co).coeff(w) for w in wv])
W = nullsp(lin, len(cl))
h = DerD + [sum((w[c] * Ec[c] for c in range(len(cl))), sp.zeros(d, d)) for w in W]
free = [Sv[j] for j in range(len(vs)) if indep[j] == 1]
dep = [v for v in Sv if v not in free]
sol = sp.solve(gens, dep, dict=True)
assert len(sol) == 1, f"not a graph over the independent set: {len(sol)} solutions"
sol = sol[0]
# coordinates absent from the prime are free too
extra = [v for v in dep if v not in sol]; free += extra
par = {v: sol.get(v, v) for v in Sv}
P_Z = PrE.subs(par)
print(f"{quiver} V = {block}, component {i_}: dim {dim_}, coordinates {free}; W_Z = {[list(w) for w in W]}; dim 𝔥_Z = {len(h)}")
# the fields: [h, P] = Σ_f ξ_h^f ∂P/∂f along Z; solve for ξ (P_Z is an embedding of Z, the Jacobian has full column rank)
J = sp.Matrix([[sp.diff(P_Z[a, b], f) for f in free] for a in range(d) for b in range(d)])
fields = []
for D in h:
    rhs = sp.Matrix(list((D * P_Z - P_Z * D).applyfunc(sp.expand)))
    xi = sp.symbols(f'xi0:{len(free)}')
    eq = (J * sp.Matrix(xi) - rhs).applyfunc(sp.expand)
    s_ = sp.solve(list(eq), xi, dict=True)
    assert s_, "[h, P] is not tangent to Z"
    fields.append([sp.simplify(s_[0].get(x, 0)) for x in xi])
def act(fld, F): return sp.expand(sum(fld[i] * sp.diff(F, free[i]) for i in range(len(free))))
mons = [sp.Integer(1)] + [sp.Mul(*m) for k in range(1, DEG + 1) for m in itertools.combinations_with_replacement(free, k)]
cs = sp.symbols(f'c0:{len(mons)}'); F = sum(c * m for c, m in zip(cs, mons))
eqs = []
for fld in fields:
    e = sp.expand(sp.numer(sp.together(act(fld, F))))
    if e != 0: eqs += sp.Poly(e, *free).coeffs()
(vec,) = list(sp.linsolve(eqs, cs))
params = sorted(set().union(*[x.free_symbols for x in vec]) & set(cs), key=str)
inv = [sp.factor(sum(vc.subs({q: (1 if q == pq else 0) for q in params}) * m for vc, m in zip(vec, mons))) for pq in params]
print(f"polynomial invariants of degree ≤ {DEG}: {[str(f) for f in inv]}")
nonconst = [f for f in inv if f.free_symbols]
if nonconst:
    f = nonconst[0]
    ok = all(sp.simplify(act(fld, f)) == 0 for fld in fields)
    print(f"invariant {f}: annihilated by every field of 𝔥_Z: {ok}")
# semi-invariants: polynomials invariant under the derived algebra [𝔥_Z, 𝔥_Z], then a simultaneous eigenbasis of 𝔥_Z
# acting on them (via a generic combination); equal characters give invariant ratios.
from sympy.polys.matrices import DomainMatrix
br = [A_ * B_ - B_ * A_ for A_ in h for B_ in h]
rows_ = [list(m) for m in br if any(m)]
hd = []
if rows_:
    Rr = DomainMatrix.from_list_sympy(len(rows_), d * d, rows_).convert_to(sp.QQ).rref()[0].to_Matrix()
    hd = [sp.Matrix(d, d, list(Rr.row(i))) for i in range(Rr.rows) if any(Rr.row(i))]
def fields_of(alg):
    out = []
    for D in alg:
        rhs = sp.Matrix(list((D * P_Z - P_Z * D).applyfunc(sp.expand)))
        xi = sp.symbols(f'xi0:{len(free)}')
        s_ = sp.solve(list((J * sp.Matrix(xi) - rhs).applyfunc(sp.expand)), xi, dict=True)
        out.append([sp.simplify(s_[0].get(x, 0)) for x in xi])
    return out
fd = fields_of(hd)
eqs = []
for fld in fd:
    e = sp.expand(sp.numer(sp.together(act(fld, F))))
    if e != 0: eqs += sp.Poly(e, *free).coeffs()
(vec,) = list(sp.linsolve(eqs, cs))
params = sorted(set().union(*[x.free_symbols for x in vec]) & set(cs), key=str)
Vb = [sp.expand(sum(vc.subs({q: (1 if q == pq else 0) for q in params}) * m for vc, m in zip(vec, mons))) for pq in params]
print(f"dim [𝔥, 𝔥] = {len(hd)}; [𝔥, 𝔥]-invariants of degree ≤ {DEG}: {len(Vb)}")
kv = sp.symbols(f'kk0:{len(Vb)}')
def matrix_of(fld):
    Mx = []
    for b in Vb:
        im = sp.expand(sp.numer(sp.together(act(fld, b))))
        s_ = sp.solve(sp.Poly(im - sum(c * v for c, v in zip(kv, Vb)), *free).coeffs(), kv, dict=True)
        if not s_: return None
        Mx.append([s_[0].get(c, 0) for c in kv])
    return sp.Matrix(Mx).T
mats = [matrix_of(fld) for fld in fields]
if Vb and all(m is not None for m in mats):
    import random as _r
    gen = sum((_r.randint(1, 97) * m for m in mats), sp.zeros(len(Vb), len(Vb)))
    semis = []
    for ev, mult, vecs in gen.eigenvects():
        for v in vecs:
            f = sp.factor(sum(v[i] * Vb[i] for i in range(len(Vb))))
            chi = [sp.simplify(sp.together(act(fld, f) / f)) for fld in fields]
            if all(not c.free_symbols for c in chi): semis.append((f, tuple(chi)))
    print(f"semi-invariants (constant characters): {len(semis)}")
    for f, chi in semis: print(f"    {f}: character {list(chi)}")
    bych = {}
    for f, chi in semis:
        if f.free_symbols: bych.setdefault(chi, []).append(f)
    for chi, fs in bych.items():
        if len(fs) >= 2:
            kappa = sp.factor(fs[0] / fs[1])
            ok = all(sp.simplify(act(fld, kappa)) == 0 for fld in fields)
            print(f"INVARIANT κ = {kappa} (character {list(chi)}): annihilated by every field of 𝔥_Z: {ok}")
    # multiplicative relations: integer n with Σ n_i χ_i = 0 among the nonconstant semi-invariants give κ = Π f_i^{n_i}
    base = [(f, chi) for f, chi in semis if f.free_symbols and sp.factor_list(f)[1] and len(sp.factor_list(f)[1]) == 1
            and sp.factor_list(f)[1][0][1] == 1]                       # irreducible ones
    if len(base) >= 2:
        Cm = sp.Matrix([list(chi) for _, chi in base]).T                 # columns: characters
        for v in Cm.nullspace():
            den = sp.ilcm(*[sp.fraction(x)[1] for x in v]); v = [int(x * den) for x in v]
            kappa = sp.Mul(*[f ** n for (f, _), n in zip(base, v)])
            ok = all(sp.simplify(act(fld, kappa)) == 0 for fld in fields)
            print(f"INVARIANT κ = {kappa}: annihilated by every field of 𝔥_Z: {ok}; on Z: {sp.simplify(kappa)}")
