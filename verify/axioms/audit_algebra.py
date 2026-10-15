"""Audit of the algebraic claims in Le Chat's v1.1–v1.3 (Q4 / NQ-1 / NQ-2 / OQ-1 / OQ-2 / OQ-8).
Free algebra k<a,b> as {word: coeff}; f = Aneristic fold (w ↦ w if w ∈ a*, else 0).
Tags [Rn] are this audit's checks."""
from itertools import product
from collections import defaultdict
import sympy as sp

res = []
def report(tag, ok, claim, detail=""):
    res.append(ok); print(f"[{tag}] {'PASS' if ok else 'FAIL'}  {claim}" + (f"\n       {detail}" if detail else ""))

# ---------- free algebra ----------
def clean(d): return {w: c for w, c in d.items() if c != 0}
def add(*ds):
    r = defaultdict(int)
    for d in ds:
        for w, c in d.items(): r[w] += c
    return clean(r)
def sc(k, d): return clean({w: k * c for w, c in d.items()})
def mul(x, y):
    r = defaultdict(int)
    for u, c in x.items():
        for v, e in y.items(): r[u + v] += c * e
    return clean(r)
def f(x): return clean({w: c for w, c in x.items() if set(w) <= {'a'}})
W = lambda w: {w: 1}
words = [''] + [''.join(p) for n in (1, 2, 3) for p in product('ab', repeat=n)] + ['aaaa']
basis = [W(w) for w in words]                       # Le Chat's 16-word corpus
eq = lambda x, y: add(x, sc(-1, y)) == {}

L = lambda x, y: mul(x, f(y))       # ⊣
R = lambda x, y: mul(f(x), y)       # ⊢
# Loday, Dialgebras (LNM 1763), Def. 2.1: ⊣,⊢ associative and
#   x⊣(y⊣z) = x⊣(y⊢z),  (x⊢y)⊣z = x⊢(y⊣z),  (x⊣y)⊢z = (x⊢y)⊢z
loday = {
 'A1 (x⊣y)⊣z=x⊣(y⊣z)': lambda x, y, z: eq(L(L(x, y), z), L(x, L(y, z))),
 'A2 x⊣(y⊣z)=x⊣(y⊢z)': lambda x, y, z: eq(L(x, L(y, z)), L(x, R(y, z))),
 'A3 (x⊢y)⊣z=x⊢(y⊣z)': lambda x, y, z: eq(L(R(x, y), z), R(x, L(y, z))),
 'A4 (x⊣y)⊢z=(x⊢y)⊢z': lambda x, y, z: eq(R(L(x, y), z), R(R(x, y), z)),
 'A5 (x⊢y)⊢z=x⊢(y⊢z)': lambda x, y, z: eq(R(R(x, y), z), R(x, R(y, z))),
}
T3 = list(product(basis, repeat=3))
bad = {k: sum(not p(*t) for t in T3) for k, p in loday.items()}
report("R1", all(v == 0 for v in bad.values()),
       "free k<a,b>: (x f(y), f(x) y) satisfies all five Loday axioms on all 4096 corpus triples (v1 [V11] stands)",
       str(bad))
bogus = lambda x, y, z: eq(R(L(x, y), z), L(x, R(y, z)))   # Le Chat's 'middle exchange'
cx = mul(f(mul(W(''), f(W('')))), W('b')), L(W(''), R(W(''), W('b')))
report("R2", sum(not bogus(*t) for t in T3) > 0 and cx == ({'b': 1}, {}),
       "Le Chat's withdrawn axiom '(x⊣y)⊢z = x⊣(y⊢z)' is not a Loday axiom; its counterexample (ε,ε,b): b vs 0 reproduces",
       f"(ε,ε,b): {cx}")

br = lambda x, y: add(L(x, y), sc(-1, R(y, x)))            # [x,y] = x⊣y − y⊢x
rightL = lambda x, y, z: eq(br(br(x, y), z), add(br(br(x, z), y), br(x, br(y, z))))
leftL = lambda x, y, z: eq(br(x, br(y, z)), add(br(br(x, y), z), br(y, br(x, z))))
nr = sum(not rightL(*t) for t in T3); nl = sum(not leftL(*t) for t in T3)
a, b = W('a'), W('b')
report("R3", nr == 0 and nl > 0 and not leftL(b, a, a),
       "[x,y]=xf(y)−f(y)x: RIGHT Leibniz 0/4096 bad; Le Chat's (b,a,a) witness fails only the LEFT identity (v1 [V12] stands)",
       f"right-Leibniz bad={nr}, left-Leibniz bad={nl}")

# ---------- dendriform family: exact variety ----------
# Generic test: noncommutative symbols p_i ∈ im f, n_i ∈ ker f; f kills any monomial containing an n.
p1, p2, p3, n1, n2, n3 = sp.symbols('p1 p2 p3 n1 n2 n3', commutative=False)
Pset = {p1, p2, p3}
def fz(e):
    e = sp.expand(e); out = 0
    for t in sp.Add.make_args(e):
        nc = [s for s in t.free_symbols if not s.is_commutative]   # ignore scalar b_i
        if all(s in Pset for s in nc): out += t
    return out
B = sp.symbols('b1:5')
def prec(x, y): return sp.expand(x * fz(y))
def succ(x, y, bb=B): return sp.expand(bb[0]*x*y + bb[1]*fz(x)*y + bb[2]*x*fz(y) + bb[3]*fz(x)*fz(y))
X, Y, Z = p1 + n1, p2 + n2, p3 + n3
star = lambda x, y: sp.expand(prec(x, y) + succ(x, y))
D = [sp.expand(prec(prec(X, Y), Z) - prec(X, star(Y, Z))),
     sp.expand(prec(succ(X, Y), Z) - succ(X, prec(Y, Z))),
     sp.expand(succ(star(X, Y), Z) - succ(X, succ(Y, Z)))]
def coeffs(e):
    d = defaultdict(int)
    for t in sp.Add.make_args(e):
        c, nc = t.args_cnc(); d[sp.Mul(*nc)] += sp.Mul(*c)
    return [sp.factor(v) for v in d.values() if sp.expand(v) != 0]
eqs = [c for e in D for c in coeffs(e)]
sols = sp.solve(eqs, B, dict=True)
print("       D1–D3 coefficient equations:", sorted(set(map(str, eqs))))
print("       solutions:", sols)
# Le Chat's parametrisation (b1..b4): line1 (−1−t,1+t,t,−t); line2 (−t,1+t,t,−1−t); (0,1,−1,0); 0
t = sp.symbols('t')
fam = [(-1-t, 1+t, t, -t), (-t, 1+t, t, -1-t), (0, 1, -1, 0), (0, 0, 0, 0)]
inV = all(all(sp.simplify(e.subs(dict(zip(B, m)))) == 0 for e in eqs) for m in fam)
dbl = (1, -1, -1, 1)
dblbad = not all(sp.expand(e.subs(dict(zip(B, dbl)))) == 0 for e in eqs)
report("R4", inV and dblbad,
       "dendriform completions of ≼=xf(y): Le Chat's two lines hold for ALL t (symbolic), plus two points; doublefold excluded",
       "so '34 members, t∈[−8,7]' is a search-range artifact: the family is two affine lines over k plus two points")

# Baxter star in blocks: x⋆y = f(x)f(y) − (x−f(x))(y−f(y))
can = (-1, 1, 0, 0)
st = sp.expand(prec(X, Y) + succ(X, Y, can))
report("R5", sp.expand(st - (p1*p2 - n1*n2)) == 0,
       "canonical Baxter star is x⋆y = pq − nm: (A,⋆) ≅ im f × ker f (n ↦ −n), so ann[−,−]⋆ = Z(im f) ⊕ Z(ker f)")

# ---------- UT_n ----------
def E(n, i, j):
    M = sp.zeros(n, n); M[i, j] = 1; return M
def UTbasis(n): return [(i, j) for i in range(n) for j in range(i, n)]
def diagp(M): return sp.diag(*[M[i, i] for i in range(M.shape[0])])
def ann_UT(n):
    xs = sp.symbols(f'x0:{n*(n+1)//2}')
    Xm = sp.zeros(n, n)
    for k, (i, j) in enumerate(UTbasis(n)): Xm[i, j] = xs[k]
    sstar = lambda P, Q: P*diagp(Q) + diagp(P)*Q - P*Q
    eqs = []
    for (i, j) in UTbasis(n):
        Yb = E(n, i, j); eqs += list(sstar(Xm, Yb) - sstar(Yb, Xm))
    sol = sp.solve([e for e in eqs if e != 0], xs, dict=True)[0]
    free = [(i, j) for k, (i, j) in enumerate(UTbasis(n)) if xs[k] not in sol]
    return free
annok = all(ann_UT(n) == sorted([(i, i) for i in range(n)] + [(0, n-1)]) for n in (3, 4, 5))
report("R6", annok, "UT_n, n=3,4,5: ⋆-commutator annihilator = diag ⊕ k·E_1n (= Z(n_n)); Le Chat's 'E14 fingerprint' is the centre of n_4, general n")

# naive dialgebra on UT_4 (correct axioms)
n = 4
UB = [E(n, i, j) for (i, j) in UTbasis(n)]
Lm = lambda x, y: x*diagp(y); Rm = lambda x, y: diagp(x)*y
lod = [lambda x, y, z: Lm(Lm(x, y), z) - Lm(x, Lm(y, z)), lambda x, y, z: Lm(x, Lm(y, z)) - Lm(x, Rm(y, z)),
       lambda x, y, z: Lm(Rm(x, y), z) - Rm(x, Lm(y, z)), lambda x, y, z: Rm(Lm(x, y), z) - Rm(Rm(x, y), z),
       lambda x, y, z: Rm(Rm(x, y), z) - Rm(x, Rm(y, z))]
ub = [sum(any(v != 0 for v in ax(x, y, z)) for x in UB for y in UB for z in UB) for ax in lod]
report("R7", ub == [0]*5, "UT_4: (xf(y), f(x)y) satisfies all five Loday axioms (Le Chat's [0,12,0,0,12] used wrong forms)", str(ub))

# ---------- ε-bialgebras ----------
# Aguiar's canonical ε-bialgebra on a path algebra: Δ = delete one arrow (letter).
def tens_add(*ds):
    r = defaultdict(int)
    for d in ds:
        for k, c in d.items(): r[k] += c
    return {k: c for k, c in r.items() if c}
def Dl(x):  # letter deletion, on words
    r = defaultdict(int)
    for w, c in x.items():
        for i in range(len(w)): r[(w[:i], w[i+1:])] += c
    return dict(r)
def eps_rule(Delta, x, y):
    lhs = Delta(mul(x, y))
    r1 = {(u, v + vv): c*e for (u, v), c in Delta(x).items() for vv, e in y.items()}
    r2 = {(uu + u, v): c*e for (u, v), c in Delta(y).items() for uu, e in x.items()}
    return tens_add(lhs, {k: -c for k, c in tens_add(r1, r2).items()}) == {}
def coassoc(Delta, x):
    l = defaultdict(int); r = defaultdict(int)
    for (u, v), c in Delta(x).items():
        for (u1, u2), e in Delta(W(u)).items(): l[(u1, u2, v)] += c*e
        for (v1, v2), e in Delta(W(v)).items(): r[(u, v1, v2)] += c*e
    return {k: c for k, c in l.items() if c} == {k: c for k, c in r.items() if c}
ok_eps = all(eps_rule(Dl, x, y) for x in basis for y in basis)
ok_co = all(coassoc(Dl, x) for x in basis)
def Ddec(x):   # Le Chat's 'left-rooted deconcatenation' (splits i < |w|)
    r = defaultdict(int)
    for w, c in x.items():
        for i in range(len(w)): r[(w[:i], w[i:])] += c
    return dict(r)
dec_eps = all(eps_rule(Ddec, x, y) for x in basis for y in basis)
dec_co = sum(not coassoc(Ddec, x) for x in basis if x != W(''))
report("R8", ok_eps and ok_co and dec_eps and dec_co == 15,
       "k<a,b> IS an ε-bialgebra via letter deletion (coassoc + ε-rule); deconcatenation satisfies the ε-rule but is not coassociative (15/15)",
       "so v1.1 [W11] ('free algebra is an ε-bialgebra with deconcatenation') was false, and v1.3's 'barren of coalgebra' question is answered: Aguiar's own example")
fmorph = all(Dl(f(x)) == {k: c for k, c in ((k, c) for (u, v), c in Dl(x).items()
              for k in [(u, v)] if set(u + v) <= {'a'})} for x in basis)
witness = (Dl(f(W('ab'))), {(u, v): c for (u, v), c in Dl(W('ab')).items() if set(u + v) <= {'a'}})
report("R9", not fmorph, "the fold f is NOT a coalgebra morphism for the genuine ε-coproduct",
       f"x=ab: Δf(x)={witness[0]}  (f⊗f)Δ(x)={witness[1]}")

# UT_4 arrow deletion: Δ(E_ij) = Σ_{i≤m<j} E_im ⊗ E_{m+1,j}
def DUT(i, j): return {((i, m), (m+1, j)): 1 for m in range(i, j)}
def mUT(p, q): return (p[0], q[1]) if p[1] == q[0] else None
UTb = UTbasis(4)
def eps_UT(x, y):
    xy = mUT(x, y); lhs = DUT(*xy) if xy else {}
    rhs = defaultdict(int)
    for (u, v), c in DUT(*x).items():
        vy = mUT(v, y)
        if vy: rhs[(u, vy)] += c
    for (u, v), c in DUT(*y).items():
        xu = mUT(x, u)
        if xu: rhs[(xu, v)] += c
    return lhs == {k: c for k, c in rhs.items() if c}
def co_UT(x):
    l = {(u1, u2, v) for (u, v) in DUT(*x) for (u1, u2) in DUT(*u)}
    r = {(u, v1, v2) for (u, v) in DUT(*x) for (v1, v2) in DUT(*v)}
    return l == r
report("R10", all(eps_UT(x, y) for x in UTb for y in UTb) and all(co_UT(x) for x in UTb) and DUT(0, 3) != {},
       "UT_4 = path algebra of A_4 carries a nonzero coassociative ε-coproduct (arrow deletion); v1.3's 'only Δ≡0' is an ansatz artifact (degree-0 only)")

# quasitriangularity obstruction by path-length degree: Δ has degree −1, Δ_r(x)=x·r−r·x has degree ≥ 0.
# Exhaustive small confirmation on UT_3: no r ∈ UT_3⊗UT_3 with Δ_r = Δ on basis.
n3 = UTbasis(3)
cs = {(p, q): sp.Symbol(f"r_{p[0]}{p[1]}_{q[0]}{q[1]}") for p in n3 for q in n3}
eqs = []
for x in n3:
    tgt = defaultdict(int)
    for k, c in DUT(*x).items(): tgt[k] += c
    got = defaultdict(int)
    for (u, v), c in cs.items():
        xu = mUT(x, u)
        if xu: got[(xu, v)] += c         # x·r
        vx = mUT(v, x)
        if vx: got[(u, vx)] -= c         # −r·x
    for k in set(tgt) | set(got): eqs.append(got[k] - tgt[k])
report("R11", sp.solve(eqs, list(cs.values()), dict=True) == [],
       "UT_3 with arrow deletion: Δ = Δ_r has no solution r (degree obstruction: Δ lowers path length, Δ_r cannot)")

print(f"\n{sum(res)}/{len(res)} checks pass")
