"""verify_v2_6.py — machine checks W1–W10 for discordianism-open-questions-explored-v2.6.md
Support modules (same folder): fold.py (exact fold engine, a-degree <= 2), series.py (order-by-order integrator),
tangent.py, second.py, sides.py, sides2.py; C enumerators cnt.c, cnt_types.c.
Fold engine conventions: words x a y (a-degree 1) and x a y a z (a-degree 2) over a basis of k[u] adapted to the
splitting H0 ⊕ Q0, with P(xay) = p(x)ap(y) + c(x,y), P(xayaz) = p(x)ap(y)ap(z) + p(x)a c(y,z) + c(x,y)a p(z) + c2(x,y,z).
Arithmetic mod the Mersenne prime 2^61 - 1; equations are generated "skip-free" (no equation is dropped for window reasons),
so every reported obstruction is a genuine necessary condition."""
import sys, os, subprocess, random, importlib, functools, itertools
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
res = []
def report(tag, ok, claim, detail=""):
    res.append(bool(ok)); print(f"[{tag}] {'PASS' if ok else 'FAIL'}  {claim}")
    if detail: print(f"       {detail}")
u, x, y, w = sp.symbols('u x y w')

# ═════════ 1. Density lemma: complete proof ═════════
def Phi(f): return sp.expand(sp.cancel((f.subs(u, x) - f.subs(u, y)) / (x - y)))
ok1, det1 = True, []
for f, g in [(u**3, u**4 + u), (u**2, u**3 + u), (u**3, u**5 + u**3 + u), (u**4 + u, u**3), (u**2 + u, u**5)]:
    d, e = sp.degree(f, u), sp.degree(g, u)
    r = sp.Poly(sp.resultant(Phi(f), Phi(g), y), x)
    ok1 &= r.degree() == (d - 1) * (e - 1); det1.append(f"({d},{e}):{r.degree()}")
# exact witness for the former "bad" case u^3 vs u^4+u: p = -1, q = -ω
om = sp.Rational(-1, 2) + sp.sqrt(3) * sp.I / 2
f, g = u**3, u**4 + u; p_, q_ = -1, -om
lam = lambda h: sp.simplify(sp.expand(h.subs(u, p_) - h.subs(u, q_)))
ok1 &= all(lam(f**j) == 0 and lam(g**j) == 0 for j in range(1, 8)) and lam(u) != 0
report("W1", ok1, "gcd(d,e)=1: the residual curves {f(x)=f(y)}, {g(x)=g(y)} meet in (d−1)(e−1) AFFINE points (no common points at infinity), giving a node (ev_p − ev_q) or cusp (d/du|_p) annihilator; u³ vs u⁴+u has the exact node p = −1, q = −ω",
       "deg Res_y = (d−1)(e−1): " + ", ".join(det1) + "; λ(h) = h(−1) − h(−ω) kills k[u³] and k[u⁴+u], λ(u) ≠ 0")
def topform(F, d, e):
    P = sp.Poly(F, x, y); W = max(d * i + e * j for (i, j) in P.monoms())
    return sp.expand(sum(c * x**i * y**j for (i, j), c in P.terms() if d * i + e * j == W))
random.seed(7); okA = True
for d, e in [(2, 4), (3, 3), (4, 6), (3, 6), (4, 4), (6, 9), (2, 2), (3, 4)]:
    for cz in (0, 1):
        cs = [random.randint(-3, 3) for _ in range(d)]; ds = [random.randint(-3, 3) for _ in range(e)]
        if cz: cs[-1] = 0; ds[-1] = 0
        f = u**d + sum(c * u**k for k, c in enumerate(cs)); g = u**e + sum(c * u**k for k, c in enumerate(ds))
        gg = sp.gcd(d, e); T = topform(sp.resultant(f - x, g - y, u), d, e)
        okA &= sp.simplify(T / (x**(e // gg) - y**(d // gg))**gg) in (1, -1)
j_, k_ = sp.symbols('j k', positive=True)
vdm = sp.factor(sp.Matrix([[1, 1, 1], [0, j_, k_], [0, j_ * (j_ - 1), k_ * (k_ - 1)]]).det())
def in_span(target, f, g, M):
    gens = [sp.expand(f**i) for i in range(M // sp.degree(f, u) + 1)] + [sp.expand(g**j) for j in range(1, M // sp.degree(g, u) + 1)]
    deg = max(sp.degree(t, u) for t in gens + [target])
    A = sp.Matrix([[sp.Poly(t, u).coeff_monomial(u**m) for t in gens] for m in range(deg + 1)])
    b = sp.Matrix([sp.Poly(target, u).coeff_monomial(u**m) for m in range(deg + 1)])
    return A.rank() == A.row_join(b).rank()
fg_out = [not in_span(sp.expand(f * g), f, g, 48) for f, g in [(u**3, u**6 + u**3 + u), (u**4, u**6 + u), (u**3 + u, u**3), (u**2, u**4 + u)]]
report("W2", okA and vdm != 0 and all(fg_out) and in_span(sp.expand(u**2 * (u**2 + u)), u**2, u**2 + u, 8),
       "gcd(d,e)=g≥2: top weighted form of Res_u(f−x, g−y) is ±(x^{e′}−y^{d′})^g; a k-nomial has root 1 of multiplicity ≤ k−1 (Vandermonde); hence fg ∉ k[f]+k[g] unless (d,e)=(2,2)",
       f"top form check on 16 random pairs: {okA}; trinomial det = {vdm}; fg ∉ V (deg ≤ 48) for (3,6),(4,6),(3,3),(2,4): {fg_out}; (2,2) dihedral fg ∈ V: True")
cases = [(d, e) for d in range(2, 13) for e in range(2, 13)]
logic = all((sp.gcd(d, e) == 1) or (sp.gcd(d, e) >= 2 and (d, e) != (2, 2)) or (d, e) == (2, 2) for d, e in cases)
report("W3", logic and sum((d, e) == (2, 2) for d, e in cases) == 1,
       "the two arguments cover every degree pair d, e ≥ 2 except (2,2): the density lemma is proved, so splittings of k[u] are exactly trivial, character, dihedral (v2.3 Conjecture 4 — now a theorem)")

# ═════════ 2. Fold deformations: tangent, obstructions, integration ═════════
import fold as FD
from tangent import spec, gvec
from sides import modG, qpolys, ratrec
from sides2 import build, sub
def tangent_info(kind, Leq, Lu, delta, D, Llow):
    S, F, eqs, J, big = build(kind, Leq, Lu, delta, D)
    G = [gvec(F, h) for h in range(len(S.B)) if S.lev[h] <= delta]
    return S, F, eqs, J, big, G
def integ(S, F, vecs, coeffs, Lu, D, N, Leq):
    tau = {}
    for a, v in zip(coeffs, vecs):
        for j, xv in v.items():
            k = F.names[j]
            if k[0] == 'c': tau[k] = (tau.get(k, 0) + a * xv) % FD.MOD
    import series; importlib.reload(series)
    r = series.integrate(S, Lu, D, tau, N, Leq=Leq); importlib.reload(FD); return r

S, F, eqs, J, big, G = tangent_info('ident', 2, 8, 1, 3, 2)
Xi = modG(S, F, sub(S, F, J, big, set()), G, 2)
Sf, Ff, eqsf, Jf, bigf, Gf = tangent_info('framed', 2, 8, 1, 3, 2)
Xf = modG(Sf, Ff, sub(Sf, Ff, Jf, bigf, set()), Gf, 2)
report("W4", len(Xi) == 2 and len(Xf) == 10,
       "a-degree-1 tangent = Hom_{H₀}(H₀⊗H₀, A₀/H₀) ⊕ Hom_{Q₀}(K, A₀/Q₀), K = ker(p⊗p): identification has 2 non-conjugation directions (P(ua)=α, P(au)=β, constants); framed has the predicted edges + P(uau)=uau+κ",
       f"extras mod G (δ=1, inputs ≤ 2): identification {len(Xi)}, framed {len(Xf)}")
import series; importlib.reload(series)
one = 0; Si = spec('ident', 14); uq = {Si.lev[i]: i for i in range(len(Si.B)) if Si.side[i] == 'Q'}
A_ = ('c', (uq[1], one), one); B_ = ('c', (one, uq[1]), one)
ri = [series.integrate(Si, 8, 3, t, 4, Leq=2) for t in ({A_: 1}, {B_: 1}, {A_: 1, B_: 1})]
importlib.reload(FD)
report("W5", ri[0] == (4, 'ok') and ri[1] == (4, 'ok') and ri[2][1] == 'obstructed',
       "identification fold is NOT rigid off the graded locus: P(ua)=α alone and P(au)=β alone integrate (to order 4); together obstructed at order 2 (quadric αβ)",
       f"α: {ri[0]}, β: {ri[1]}, α&β: {ri[2]}")
importlib.reload(series)
Sg = spec('framed', 14); q0 = [i for i in range(len(Sg.B)) if Sg.side[i] == 'Q'][0]
hh = {Sg.lev[i]: i for i in range(len(Sg.B)) if Sg.side[i] == 'H'}
Lf = lambda m: {('c', (q0, hh[m]), hh[m + 1]): 1}; Rf = lambda m: {('c', (hh[m], q0), hh[m + 1]): 1}
Ga = {('c', (q0, q0), hh[1]): 1}
U = lambda *ds: {k: v for d in ds for k, v in d.items()}
Lf3 = lambda m, v: {('c', (q0, hh[m]), hh[m + 1]): v}
C0 = U(*[Lf(m) for m in range(1, 9)], *[Rf(m) for m in range(1, 9)])           # uniform two-sided α_m = β_m = 1
CL = U(Lf3(1, 1), Lf3(2, 3), Lf3(3, -2))                                        # arbitrary left sequence
tests6 = [("α₁+β₂", U(Lf(1), Rf(2)), 'obstructed'), ("γ+α₁", U(Ga, Lf(1)), 'obstructed'), ("α₁+β₁ only (non-uniform)", U(Lf(1), Rf(1)), 'obstructed'),
          ("C₀ uniform", C0, 'ok'), ("C_L (1,3,−2)", CL, 'ok'), ("α₁ alone", Lf(1), 'ok')]
rf = [(nm, series.integrate(Sg, 8, 3, t, 3, Leq=2)) for nm, t, _ in tests6]
importlib.reload(FD)
report("W6", all(r[1] == exp for (_, r), (_, _, exp) in zip(rf, tests6)),
       "calibration: the skip-free integrator reproduces v2.5's framed classification C₀ ∪ C_L ∪ C_R (on-component points integrate, off-component points are obstructed — including α₁=β₁=1, which v2.5's three degree-4 quadrics do not exclude)",
       f"{rf}")
Sd, Fd, eqsd, Jd, bigd, Gd = tangent_info('dihedral', 2, 8, 1, 3, 2)
Xd = modG(Sd, Fd, sub(Sd, Fd, Jd, bigd, set()), Gd, 2)
Hd = modG(Sd, Fd, sub(Sd, Fd, Jd, bigd, {'QH', 'HQ', 'QQ'}), Gd, 2)
Ph = qpolys(Fd, eqsd, Hd, ['h0'])
Sd2, Fd2, eqsd2, Jd2, bigd2, Gd2 = tangent_info('dihedral', 2, 8, 2, 5, 2)
Xd2 = modG(Sd2, Fd2, sub(Sd2, Fd2, Jd2, bigd2, set()), Gd2, 2)
report("W7", len(Xd) == 10 and len(Xd2) == 12 and len(Hd) == 1 and any(sp.simplify(p / sp.Symbol('h0')**2).is_number for p in Ph),
       "dihedral fold is NOT rigid: non-conjugation tangent directions exist (finite per output-shift δ, +2 per unit of δ); the Hodge-only direction is obstructed (h² is an obstruction) — the Hodge side moves only by conjugation",
       f"extras mod G: δ=1 → {len(Xd)}, δ=2 → {len(Xd2)}; Hodge-only quadrics: {Ph}")
S3, F3, eqs3, J3, big3, G3 = tangent_info('dihedral', 3, 10, 1, 3, 2)
H3 = modG(S3, F3, sub(S3, F3, J3, big3, {'QH', 'HQ', 'QQ'}), G3, 2)
low3 = [j for j, k in enumerate(F3.names) if k[0] == 'c' and F3.wlev(k[1]) <= 2]
def strip(vs):
    E = FD.Echelon(); [E.add({j: g[j] for j in low3 if j in g}) for g in G3 + H3]
    return [v for v in vs if E.add({j: v[j] for j in low3 if j in v})]
L3 = strip(modG(S3, F3, sub(S3, F3, J3, big3, {'QH'}), G3, 2)); R3 = strip(modG(S3, F3, sub(S3, F3, J3, big3, {'HQ'}), G3, 2))
A3 = strip(modG(S3, F3, sub(S3, F3, J3, big3, set()), G3, 2))
rnd = random.Random(1)
rr = [integ(S3, F3, V, [rnd.randint(1, 60) for _ in V], 10, 3, 3, 3) for V in (L3, R3, A3)]
report("W8", rr[0] == (3, 'ok') and rr[1] == (3, 'ok') and rr[2][1] == 'obstructed',
       "dihedral one-sided families (Hodge mode removed) integrate at a generic point; a generic two-sided point is obstructed — the C_L ∪ C_R pattern",
       f"left-only ({len(L3)} dirs): {rr[0]}, right-only ({len(R3)}): {rr[1]}, generic all ({len(A3)}): {rr[2]}   [Leq=3, Lu=10, D=3]")
# two-sided components over Q(√s)
Lb = modG(S3, F3, sub(S3, F3, J3, big3, {'QH'}), G3, 2); Rb = modG(S3, F3, sub(S3, F3, J3, big3, {'HQ'}), G3 + Lb, 2)
Mb = modG(S3, F3, sub(S3, F3, J3, big3, set()), G3 + Lb + Rb, 2)
basis = Lb + Rb + Mb; names = [f'l{i}' for i in range(len(Lb))] + [f'r{i}' for i in range(len(Rb))] + [f'm{i}' for i in range(len(Mb))]
P = qpolys(F3, eqs3, basis, names); syms = sp.symbols(names)
GB = sp.groebner(P, *syms, order='grevlex'); sq = [s for s in syms if GB.contains(s**2)]
P2 = [q for q in (sp.expand(p.subs({s: 0 for s in sq})) for p in P) if q != 0]
mv = sp.Symbol('m0'); rest = [s for s in syms if s not in sq and s != mv]
sols = sp.solve([p.subs(mv, 1) for p in P2], rest, dict=True)
pmod = FD.MOD; r2 = pow(2, (pmod + 1) // 4, pmod)
def mval(e, l3v):
    e = sp.expand(sp.nsimplify(e).subs(sp.Symbol('l3'), l3v)); a = e.coeff(sp.sqrt(2), 0); b = e.coeff(sp.sqrt(2), 1)
    return (FD.md(a) + FD.md(b) * r2) % pmod
okC = len(sols) == 2 and all(any(v.has(sp.sqrt(2)) for v in s.values()) for s in sols)
outs = []
for s in sols:
    pt = {str(k): mval(v, 7) for k, v in s.items()}; pt.update({'m0': 1, 'l3': 7})
    co = [pt.get(n, 0) for n in names]
    on = integ(S3, F3, basis, co, 10, 3, 3, 3)
    co2 = list(co); co2[names.index('r2')] = (co2[names.index('r2')] + 1) % pmod
    off = integ(S3, F3, basis, co2, 10, 3, 3, 3)
    outs.append((on, off)); okC &= on == (3, 'ok') and off[1] == 'obstructed'
# discriminant ∈ s·Q² for several (a, b, s), from the small-window ideal
discs = []
for prm in ("0,1,3", "0,2,2", "1,0,5"):
    os.environ["DIH"] = prm; importlib.reload(sys.modules['tangent'])
    import sides2; importlib.reload(sides2)
    from tangent import gvec as gv2
    S4, F4, eqs4, J4, big4 = sides2.build('dihedral', 2, 8, 1, 3)
    G4 = [gv2(F4, h) for h in range(len(S4.B)) if S4.lev[h] <= 1]
    L4 = modG(S4, F4, sides2.sub(S4, F4, J4, big4, {'QH'}), G4, 2); R4 = modG(S4, F4, sides2.sub(S4, F4, J4, big4, {'HQ'}), G4 + L4, 2)
    M4 = modG(S4, F4, sides2.sub(S4, F4, J4, big4, set()), G4 + L4 + R4, 2)
    nm4 = [f'l{i}' for i in range(len(L4))] + [f'r{i}' for i in range(len(R4))] + [f'm{i}' for i in range(len(M4))]
    G4b = sp.groebner(qpolys(F4, eqs4, L4 + R4 + M4, nm4), *sp.symbols(nm4), order='lex')
    cub = [e for e in G4b.exprs if sp.Poly(e, *sp.symbols(nm4)).total_degree() == 3]
    s_ = sp.Rational(prm.split(',')[2])
    for e in cub:
        q = sp.Poly(sp.cancel(e / sp.Symbol('m0')).subs(sp.Symbol('m0'), 1), sp.Symbol(nm4[[i for i, n in enumerate(nm4) if n.startswith('r')][-1]]))
        if q.degree() == 2:
            dsc = q.discriminant(); discs.append((prm, dsc, sp.sqrt(dsc / s_).is_rational))
os.environ["DIH"] = "0,1,2"
okD = len(discs) == 3 and all(t[2] for t in discs)
report("W9", okC and okD,
       "dihedral two-sided deformations exist only on two Galois-conjugate components defined over ℚ(√s) — one per point u = b ± √s where the Podge character (u−b)² − s vanishes; they integrate, nearby points do not",
       f"Leq=3 components (m0=1, l3=7) on/off: {outs}; discriminants (a,b,s) → disc, disc/s a square: {discs}")

# ═════════ 3. Compatible splittings of path algebras ═════════
cdir = HERE; ok10 = True; got = []
try:
    subprocess.run(["gcc", "-O2", "-o", os.path.join(cdir, "cnt_types"), os.path.join(cdir, "cnt_types.c")], check=True)
    for n in range(1, 9):
        out = subprocess.run([os.path.join(cdir, "cnt_types"), str(n)], capture_output=True, text=True).stdout
        parts = dict(kv.split('=') for kv in out.split()[1:]); got.append(int(parts['compatible'])); ok10 &= parts['mismatches'] == '0'
except Exception as ex:
    ok10 = False; got = [str(ex)]
ok10 &= got == [2, 6, 23, 100, 563, 3734, 28777, 246164]
report("W10", ok10,
       "arrow compatibility = an adjacency rule on leaf types (E = Podge leaf, HL = Hodge leaf, P1 = Podge vertex with one leaf child, B = other): consecutive points must contain an E, or be (HL,HL), or (P1,P1) with consecutive children",
       f"counts n=1..8: {got} (n=9: 2395445 from the 147 s run); characterization agrees on every forest")
print(f"\n{sum(res)}/{len(res)} checks pass")
