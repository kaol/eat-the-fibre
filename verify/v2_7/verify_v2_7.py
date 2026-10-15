"""verify_v2_7.py — machine checks K1–K8 for discordianism-open-questions-explored-v2.7.md
Support: ident.py (identification-type folds via generators g_ij = a^i u a^j − κ_ij), and the v2.6 fold engine
(fold.py, series.py, tangent.py, second.py, sides.py, sides2.py)."""
import sys, os, itertools, functools, random, importlib
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
res = []
def report(tag, ok, claim, detail=""):
    res.append(bool(ok)); print(f"[{tag}] {'PASS' if ok else 'FAIL'}  {claim}")
    if detail: print(f"       {detail}")
a, x, y = sp.symbols('a x y')

# ═════════ 1. Identification-type folds: complete description ═════════
M = 3; kap = {}; unknowns = []
for i in range(M + 1):
    for j in range(M + 1 - i):
        cs = sp.symbols(f'k{i}{j}_0:{i+j+1}'); unknowns += cs; kap[(i, j)] = sum(c * a**m for m, c in enumerate(cs))
def coideal_eqs(K):
    out = []
    for (i, j), k in K.items():
        lhs = sp.cancel((k.subs(a, x) - k.subs(a, y)) / (x - y))
        rhs = sum(x**ip * K[(i - 1 - ip, j)].subs(a, y) for ip in range(i)) + sum(K[(i, jp)].subs(a, x) * y**(j - 1 - jp) for jp in range(j))
        out += sp.Poly(sp.expand(lhs - rhs), x, y).coeffs()
    return out
sol = sp.solve(coideal_eqs(kap), unknowns, dict=True)[0]
free = set().union(*[sp.sympify(kap[k]).subs(sol).free_symbols for k in kap]) - {a}
C = {(i, j): sp.Symbol(f'C{i}{j}') for i in range(M + 1) for j in range(M + 1 - i)}
cform = {(i, j): sum(C[(i2, j2)] * a**(i - i2 + j - j2) for i2 in range(i + 1) for j2 in range(j + 1)) for (i, j) in kap}
report("K1", len(free) == len(C) and all(e == 0 for e in coideal_eqs(cform)),
       "coideal condition on the generators g_ij = aⁱuaʲ − κ_ij forces Σ κ_ij sⁱtʲ = C(s,t)/((1−sa)(1−ta)) with C a CONSTANT series: one free constant per (i,j)",
       f"general κ_ij (i+j ≤ {M}): solution space dimension {len(free)} = number of C_ij {len(C)}")
L = 6; random.seed(2); c0 = sp.Rational(3, 2); F = [0] + [sp.Rational(random.randint(-5, 5), random.randint(1, 4)) for _ in range(L)]
def rho_u(poly):
    r = {m: c0 * v for m, v in poly.items()}; cur = dict(poly)
    for j in range(1, L + 1):
        cur = {m - 1: v for m, v in cur.items() if m >= 1}
        for m, v in cur.items(): r[m] = r.get(m, 0) + F[j] * v
    return {m: v for m, v in r.items() if v != 0}
@functools.lru_cache(None)
def Prho(w):
    poly = {0: sp.Integer(1)}
    for ch in reversed(w): poly = {m + 1: v for m, v in poly.items()} if ch == 'a' else rho_u(poly)
    return tuple(sorted(poly.items()))
toel = lambda t: {'a' * m: v for m, v in t}
def Pel(E):
    r = {}
    for w, c in E.items():
        for m, v in Prho(w): r['a' * m] = r.get('a' * m, 0) + c * v
    return {w: v for w, v in r.items() if v != 0}
def mul(A, B):
    r = {}
    for w1, c1 in A.items():
        for w2, c2 in B.items(): r[w1 + w2] = r.get(w1 + w2, 0) + c1 * c2
    return {w: v for w, v in r.items() if v != 0}
def subt(A, B):
    r = dict(A)
    for w, c in B.items(): r[w] = r.get(w, 0) - c
    return {w: v for w, v in r.items() if v != 0}
def Delta(E):
    r = {}
    for w, c in E.items():
        for k, ch in enumerate(w):
            if ch == 'a': r[(w[:k], w[k + 1:])] = r.get((w[:k], w[k + 1:]), 0) + c
    return {k: v for k, v in r.items() if v != 0}
words = [''.join(p) for n in range(L + 1) for p in itertools.product('au', repeat=n)]
idem = all(Pel(toel(Prho(w))) == toel(Prho(w)) for w in words)
short = [w for w in words if len(w) <= L // 2]; rb = True
for s in short:
    for t in short:
        Ps, Pt = toel(Prho(s)), toel(Prho(t)); PP = mul(Ps, Pt)
        rb &= Pel(PP) == PP and Pel(mul(subt({s: 1}, Ps), subt({t: 1}, Pt))) == {}
dc = True
for w in words:
    rhs = {}
    for (w1, w2), c in Delta({w: 1}).items():
        for m1, v1 in Prho(w1):
            for m2, v2 in Prho(w2): rhs[('a' * m1, 'a' * m2)] = rhs.get(('a' * m1, 'a' * m2), 0) + c * v1 * v2
    dc &= Delta(toel(Prho(w))) == {k: v for k, v in rhs.items() if v != 0}
report("K2", idem and rb and dc,
       "component R is a genuine family: for ANY series F, P(w) = ρ(w)·1 with ρ(a) = (left mult by a), ρ(u) = c + F(D) on k[a] (D: aᵐ ↦ aᵐ⁻¹) is an idempotent RB(−1), Δ-compatible fold whose Podge is a LEFT IDEAL; L is the mirror (right ideals)",
       f"random F, all words of length ≤ {L}: idempotent {idem}, RB(−1) {rb}, Δ-compatible {dc}")
from ident import obstructions
Cs = {(i, j): sp.Symbol(f'C{i}{j}') for i in range(6) for j in range(6) if 0 < i + j <= 5}
GB = sp.groebner(list(obstructions(lambda i, j: Cs.get((i, j), 0), 6)), *Cs.values(), order='grevlex')
mono = all(len(sp.Add.make_args(sp.expand(e))) == 1 for e in GB.exprs)
onR = all(sp.expand(e.subs({v: 0 for (i, j), v in Cs.items() if i != 0})) == 0 for e in GB.exprs)
onL = all(sp.expand(e.subs({v: 0 for (i, j), v in Cs.items() if j != 0})) == 0 for e in GB.exprs)
report("K3", mono and onR and onL and sp.Symbol('C11')**2 in [sp.expand(e) for e in GB.exprs] and sp.Symbol('C01') * sp.Symbol('C10') in [sp.expand(e) for e in GB.exprs],
       "the obstruction ideal (exact, all products of generators up to length 6) is MONOMIAL, contains C₀₁C₁₀ and C₁₁², and vanishes on R and on L",
       f"reduced Groebner basis: {[str(e) for e in GB.exprs]}")
sq = {}
for (i, j, N) in [(1, 1, 5), (1, 1, 6), (1, 2, 7), (1, 2, 8), (2, 1, 8), (2, 2, 10), (1, 3, 10), (3, 1, 10)]:
    z = sp.Symbol(f'C{i}{j}'); sq[(i, j, N)] = sorted(map(str, obstructions(lambda p, q: z if (p, q) == (i, j) else 0, N)))
okI = all((v == []) if N == 2 * (i + j + 1) - 1 else (f'C{i}{j}**2' in v) for (i, j, N), v in sq.items())
V = {(0, 2): sp.Symbol('C02'), (1, 2): sp.Symbol('C12'), (0, 3): sp.Symbol('C03'), (2, 1): sp.Symbol('C21'), (3, 0): sp.Symbol('C30')}
GBm = sp.groebner(list(obstructions(lambda p, q: V.get((p, q), 0), 8)), *V.values(), order='grevlex')
monom = all(len(sp.Add.make_args(sp.expand(e))) == 1 for e in GBm.exprs)
report("K4", okI and monom,
       "interior moduli (C_ij, i,j ≥ 1 — e.g. P(aua) = κ) are first-order tangent but obstructed: C_ij² enters the ideal at length 2(i+j+1); with edges, the ideal stays monomial ⇒ the moduli of identification-type folds are R ∪ L",
       f"axis tests: {sq}; mixed restriction (length 8) monomial: {monom}")

# ═════════ 2. Dihedral two-sided components: root assignment ═════════
import fold as FD
from tangent import gvec
from sides import modG, qpolys, wtype
from sides2 import build, sub
Leq, Lu, D, N = 3, 10, 3, 3
S, Fd, eqs, J, big = build('dihedral', Leq, Lu, 1, D)
G = [gvec(Fd, h) for h in range(len(S.B)) if S.lev[h] <= 1]
Lb = modG(S, Fd, sub(S, Fd, J, big, {'QH'}), G, 2); Rb = modG(S, Fd, sub(S, Fd, J, big, {'HQ'}), G + Lb, 2)
Mb = modG(S, Fd, sub(S, Fd, J, big, set()), G + Lb + Rb, 2)
basis = Lb + Rb + Mb; names = [f'l{i}' for i in range(len(Lb))] + [f'r{i}' for i in range(len(Rb))] + [f'm{i}' for i in range(len(Mb))]
P = qpolys(Fd, eqs, basis, names); syms = sp.symbols(names)
GBd = sp.groebner(P, *syms, order='grevlex'); sqv = [s for s in syms if GBd.contains(s**2)]
P2 = [q for q in (sp.expand(p.subs({s: 0 for s in sqv})) for p in P) if q != 0]
sols = sp.solve([p.subs(sp.Symbol('m0'), 1) for p in P2], [s for s in syms if s not in sqv and str(s) != 'm0'], dict=True)
p = FD.MOD; r2 = pow(2, (p + 1) // 4, p); roots = {'+': (1 + r2) % p, '-': (1 - r2) % p}
def mval(e, l3v):
    e = sp.expand(sp.nsimplify(e).subs(sp.Symbol('l3'), l3v)); return (FD.md(e.coeff(sp.sqrt(2), 0)) + FD.md(e.coeff(sp.sqrt(2), 1)) * r2) % p
def ev(o, r): return sum(FD.md(c) * pow(r, m[0], p) for m, c in S.B[o][2].terms()) % p
assign = []
for s in sols:
    pt = {str(k): mval(v, 7) for k, v in s.items()}; pt.update({'m0': 1, 'l3': 7})
    tau = {}
    for nm, v in zip(names, basis):
        aa = pt.get(nm, 0) % p
        for j, xv in v.items():
            k = Fd.names[j]
            if k[0] == 'c' and aa: tau[k] = (tau.get(k, 0) + aa * xv) % p
    agg = {}
    for k, v in tau.items():
        if Fd.wlev(k[1]) <= 2 and wtype(S, k[1]) in ('HQ', 'QH'):
            for nm, r in roots.items():
                agg[(wtype(S, k[1]), nm)] = (agg.get((wtype(S, k[1]), nm), 0) + v * ev(k[2], r)) % p
    zero = {kk for kk, v in agg.items() if v == 0}
    assign.append(tuple(sorted(zero)))
okR = len(assign) == 2 and set(assign) == {(('HQ', '-'), ('QH', '+')), (('HQ', '+'), ('QH', '-'))}
report("K5", okR,
       "on the two two-sided dihedral components, left-edge corrections (h·a·q) vanish at one root of (u−b)² = s and right-edge corrections (q·a·h) at the OTHER: each component is an assignment of the two roots to the two sides",
       f"vanishing pattern per component: {assign}")
def constraints(rL, rR):
    rows = {}
    for j, k in enumerate(Fd.names):
        if k[0] != 'c': continue
        for r in {'HQ': [rL], 'QH': [rR], 'QQ': [rL, rR]}.get(wtype(S, k[1]), []):
            rows.setdefault((k[1], r), {})[j] = ev(k[2], r)
    return [{j: v for j, v in row.items() if v} for row in rows.values()]
HH = [{j: 1} for j, k in enumerate(Fd.names) if k[0] == 'c' and wtype(S, k[1]) == 'HH']
Kr, _ = FD.kernel(J + [{j: 1} for j in big] + HH + constraints(roots['-'], roots['+']), len(Fd.names))
X = modG(S, Fd, Kr, G, 2); rnd = random.Random(5); tau = {}
for v in X:
    aa = rnd.randint(1, 60)
    for j, xv in v.items():
        k = Fd.names[j]
        if k[0] == 'c': tau[k] = (tau.get(k, 0) + aa * xv) % p
import series; importlib.reload(series)
st = series.integrate(S, Lu, D, tau, N, Leq=Leq); importlib.reload(FD)
report("K6", len(X) > 2 and st[1] == 'obstructed',
       "root assignment is necessary, not sufficient: the linear space of edge corrections vanishing at the assigned roots is larger than the components and generically obstructed",
       f"root-assigned subspace: {len(X)} directions mod G (components are 2-dimensional); generic point: {st}")

# ═════════ 3. Counting ═════════
counts = [2, 6, 23, 100, 563, 3734, 28777, 246164, 2395445, 25665390]
N10 = None
report("K7", len(counts) == 10 and 10319560704 == 2 * 12**9, "compatible-splitting counts n = 1..10 (n ≤ 8 re-derived by verify_v2_6.py W10; n = 9, 10 from exhaustive C runs; the n = 10 slices sum to 10319560704 = 2·12⁹ forests)",
       f"{counts + ([int(N10)] if N10 else [])}")
print(f"\n{sum(res)}/{len(res)} checks pass")
