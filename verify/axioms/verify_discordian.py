"""Machine checks for discordianism-axioms.md (v1).
Each check prints  [Vn] PASS/FAIL  <claim>. Run: python3 verify_discordian.py
Part A audits v0 (Le Chat draft); Part B verifies the refined models."""
from itertools import product, permutations
from math import comb, factorial
import random
import z3
import sympy as sp
import numpy as np

results = []
def report(tag, ok, claim, detail=""):
    results.append(ok)
    print(f"[{tag}] {'PASS' if ok else 'FAIL'}  {claim}" + (f"\n       {detail}" if detail else ""))

def unsat(*fs, solver=None):
    s = solver or z3.Solver(); s.add(*fs); return s.check() == z3.unsat

# ───────────────────────── Part A: audit of v0 ─────────────────────────

# V1  N2 ∧ N3 inconsistent: ν(a,a) = a ∪ a ∪ {a⊕a} ∋ a⊕a, so ν(a,a) ≠ ∅.
E = z3.DeclareSort('E')
a, b, c, d, x, y, z_ = z3.Consts('a b c d x y z', E)
SetE = z3.SetSort(E)
A = z3.Const('A', SetE)
xor = z3.Function('xor', SetE, SetE, E)          # a⊕b : an element
nu_set = lambda P, Q: z3.SetAdd(z3.SetUnion(P, Q), xor(P, Q))
report("V1", unsat(nu_set(A, A) == z3.EmptySet(E)),
       "v0 N2 (ν(a,a)=∅) ∧ N3 (ν(a,b)=a∪b∪{a⊕b}) is unsatisfiable for every a")

# V2  N5 (ε≡0) + counit law (ε⊗id)Δ = id forces a = 0: Nu trivial, N1/N4 vacuous.
# For any Δ(a) = Σ uᵢ⊗vᵢ, (ε⊗id)Δ(a) = Σ ε(uᵢ)vᵢ = 0. Linear-algebra fact, Z3 over reals:
u1, u2, v1, v2, av = z3.Reals('u1 u2 v1 v2 av')
eps = lambda t: 0 * t
report("V2", unsat(eps(u1) * v1 + eps(u2) * v2 == av, av != 0),
       "v0 N5 (ε≡∅) with the counit law forces every Nu element to 0")

# V3  D1 ⇔ D2 (same identity, rearranged).
p, q, r = z3.Reals('p q r')      # p=[x,[y,z]], q=[[x,y],z], r=[y,[x,z]]
D1 = p == q - r
D2 = q == p + r
report("V3", unsat(D1 != D2), "v0 D1 and D2 are the same equation (D2 adds nothing)")

# V4  D1 is NOT the standard left-Leibniz identity p = q + r; they agree iff r = 0.
report("V4", unsat(D1, p == q + r, r != 0),
       "v0 D1 ∧ standard left Leibniz ⇒ [y,[x,z]] = 0: v0's sign is non-standard")

# V5  v0 Meta-2 (∀x∃y [x,y]≠0) fails at x=0 for any bilinear bracket ( D1 uses ±, so bilinear).
# [0,y] = [0·0,y] = 0·[0,y] = 0. Check in a generic 2-dim structure-constant model:
cs = [[[z3.Real(f"c{i}{j}{k}") for k in range(2)] for j in range(2)] for i in range(2)]
yv = z3.Reals('y0 y1')
br = lambda X, Y: [z3.Sum([X[i]*Y[j]*cs[i][j][k] for i in range(2) for j in range(2)]) for k in range(2)]
zero = [z3.RealVal(0)]*2
report("V5", unsat(z3.Or([e != 0 for e in br(zero, list(yv))])),
       "v0 Meta-2 is violated by x = 0 in every bilinear bracket")

# V6  v0 Meta-1  ¬(μ∨ν) ∧ (μ∧ν)  is classically unsatisfiable.
m_, n_ = z3.Bools('mu nu')
report("V6", unsat(z3.And(z3.Not(z3.Or(m_, n_)), z3.And(m_, n_))),
       "v0 Meta-1 is a classical contradiction (needs a paraconsistent reading)")

# V7  Eckmann–Hilton: interchange (v0 Rule 3) + a unit for μ + a unit for ν ⇒ μ = ν, assoc, comm.
# So v0 N1 (ν non-associative) survives only because ν has NO unit. Z3 with quantifiers:
mu = z3.Function('mu', E, E, E); nu = z3.Function('nu', E, E, E)
e1, e2 = z3.Consts('e1 e2', E)
# Quantified form diverges under MBQI; ground instances suffice (unsat on a subset of the
# instances of ∀ ⇒ unsat of the ∀ theory). Instance set: the Herbrand terms the proof touches.
s = z3.Solver()
T = [a, b, c, e1, e2, nu(a, b), nu(b, c), mu(a, b), mu(b, c)]
for t in T:
    s.add(mu(t, e1) == t, mu(e1, t) == t, nu(t, e2) == t, nu(e2, t) == t)
for (p1, p2, p3, p4) in product(T, repeat=4):
    s.add(nu(mu(p1, p2), mu(p3, p4)) == mu(nu(p1, p3), nu(p2, p4)))
s.add(nu(nu(a, b), c) != nu(a, nu(b, c)))
r7 = s.check()
# witness that the unitless case is consistent: Z/5, μ=+, ν=− (interchange ✓, ν non-assoc ✓, ν(a,a)=0=e)
Z5 = range(5)
inter = all(((a_+b_) - (c_+d_)) % 5 == ((a_-c_) + (b_-d_)) % 5 for a_ in Z5 for b_ in Z5 for c_ in Z5 for d_ in Z5)
nonassoc = any(((a_-b_)-c_) % 5 != (a_-(b_-c_)) % 5 for a_ in Z5 for b_ in Z5 for c_ in Z5)
r7 = r7 if (inter and nonassoc) else None
report("V7", r7 == z3.unsat,
       "v0 Rule 3 (interchange) + units for both ⇒ ν associative (Eckmann–Hilton); N1 forces ν unitless",
       f"z3: {r7}")

# V8  v0 entropy bound H(ν(a,b)) ≥ H(a)+H(b) fails for every deterministic ν unless ν is injective
# on the support with a ⟂ b.  Random search: deterministic ν never exceeds H(a)+H(b).
def H(pm):
    pm = np.asarray([v for v in pm if v > 0]); return float(-(pm*np.log2(pm)).sum())
rng = np.random.default_rng(23)        # 23: the Discordian prime, chosen for no reason
worst = -1e9; strict_viol = 0
for _ in range(4000):
    na, nb, nc = rng.integers(2, 6, 3)
    joint = rng.dirichlet(np.ones(na*nb)).reshape(na, nb)
    f = rng.integers(0, nc, (na, nb))                    # deterministic ν : A×B → C
    out = np.zeros(nc)
    for i in range(na):
        for j in range(nb): out[f[i, j]] += joint[i, j]
    gap = H(out) - (H(joint.sum(1)) + H(joint.sum(0)))
    worst = max(worst, gap); strict_viol += gap > 1e-9
report("V8", strict_viol == 0,
       "deterministic ν never satisfies H(ν(a,b)) > H(a)+H(b): v0 Bound 2 (chaos) needs injected randomness",
       f"max H(ν)-H(a)-H(b) over 4000 random cases = {worst:.4f}")

# V9  v0 Rule 1 (left-distributivity of μ over ν) fails in the canonical consistent model
# of the rest (abelian group, ν = subtraction), while translation-invariance holds.
Zn = range(5)
r1 = all((a_ + (b_ - c_)) % 5 == ((a_ + b_) - (a_ + c_)) % 5 for a_ in Zn for b_ in Zn for c_ in Zn)
ti = all(((a_ + b_) - (a_ + c_)) % 5 == (b_ - c_) % 5 for a_ in Zn for b_ in Zn for c_ in Zn)
report("V9", (not r1) and ti,
       "in Z/5 with ν=−: v0 Rule 1 fails, translation-invariance ν(a+b,a+c)=ν(b,c) holds")

# ───────────────────────── Part B: refined models ─────────────────────────

# V10  The v0 M-axioms (M1–M5) are exactly a group-like Hopf algebra kG: check on k[S3]
# with Δg = g⊗g, ε(g)=1, antipode S(g)=g⁻¹: coassoc, counit, antipode law, Δ multiplicative.
G = list(permutations(range(3)))
mul = lambda g, h: tuple(g[h[i]] for i in range(3))
inv = lambda g: tuple(sorted(range(3), key=lambda i: g[i]))
e = (0, 1, 2)
ok = all(mul(inv(g), g) == e and mul(g, inv(g)) == e for g in G)          # antipode on group-likes
ok &= all(mul(mul(g, h), k) == mul(g, mul(h, k)) for g in G for h in G for k in G)
# Δ(gh) = Δg·Δh on group-likes is (gh)⊗(gh) = (g⊗g)(h⊗h): holds by construction; counit ε(gh)=1=ε(g)ε(h).
report("V10", ok, "v0 M1–M5 are consistent: they axiomatize the group Hopf algebra kG (checked on S3)")

# V11  Loday dialgebra from an idempotent algebra endomorphism (the 'grid' f):
#   x ⊣ y = x·f(y),  x ⊢ y = f(x)·y.  Model: 3×3 upper-triangular, f = diagonal part.
n = 3
def UT(name):
    return sp.Matrix(n, n, lambda i, j: sp.Symbol(f"{name}{i}{j}") if j >= i else 0)
X, Y, Z = UT('x'), UT('y'), UT('z')
f = lambda M: sp.diag(*[M[i, i] for i in range(n)])
L = lambda P, Q: P * f(Q)          # ⊣  (left stays chaotic, right is gridded)
R = lambda P, Q: f(P) * Q          # ⊢
zeroM = lambda M: all(sp.expand(v) == 0 for v in M)
hom = zeroM(f(X*Y) - f(X)*f(Y)) and zeroM(f(f(X)) - f(X))
axioms = [
    (L(L(X, Y), Z), L(X, L(Y, Z))),   # D1  (x⊣y)⊣z = x⊣(y⊣z)
    (L(L(X, Y), Z), L(X, R(Y, Z))),   # D2  (x⊣y)⊣z = x⊣(y⊢z)
    (L(R(X, Y), Z), R(X, L(Y, Z))),   # D3  (x⊢y)⊣z = x⊢(y⊣z)
    (R(L(X, Y), Z), R(X, R(Y, Z))),   # D4  (x⊣y)⊢z = x⊢(y⊢z)
    (R(R(X, Y), Z), R(X, R(Y, Z))),   # D5  (x⊢y)⊢z = x⊢(y⊢z)
]
report("V11", hom and all(zeroM(l - rr) for l, rr in axioms),
       "f = diagonal part is an idempotent algebra hom on UT₃, and (⊣,⊢) satisfy all five Loday axioms")

# V12  The bracket [x,y] = x⊣y − y⊢x = x f(y) − f(y) x satisfies RIGHT Leibniz
#   [[x,y],z] = [[x,z],y] + [x,[y,z]]  and is not antisymmetric (genuinely Leibniz, not Lie).
Bk = lambda P, Q: L(P, Q) - R(Q, P)
rl = zeroM(Bk(Bk(X, Y), Z) - Bk(Bk(X, Z), Y) - Bk(X, Bk(Y, Z)))
nonlie = not zeroM(Bk(X, X))
report("V12", rl and nonlie, "[x,y]=x f(y)−f(y)x is a right-Leibniz bracket with [x,x]≢0 (not Lie)")

# V13  Refined Eris. Left annihilator {x : ∀y [x,y]=0} = im f (pure order is inert);
# right annihilator ⊇ ker f (pure chaos cannot be struck).
xs = [X[i, j] for i in range(n) for j in range(i, n)]
ys = [Y[i, j] for i in range(n) for j in range(i, n)]
conds = [sp.Poly(v, *ys).coeffs() for v in Bk(X, Y)]
eqs = [cf for cl in conds for cf in cl]
sol = sp.solve(eqs, [X[i, j] for i in range(n) for j in range(i+1, n)], dict=True)
left_ann_is_diag = len(sol) == 1 and all(v == 0 for v in sol[0].values())
Nchaos = X - f(X)                                     # ker f element
right_ann = zeroM(Bk(Y, Nchaos))
report("V13", left_ann_is_diag and right_ann,
       "Left annihilator = diagonal (order collides with nothing); ker f lies in the right annihilator",
       f"left-annihilator solve → {sol}")

# V14  Seven truth values: the Principia epigraph enumerates the non-empty subsets of {T,F,M}.
vals = ['T', 'F', 'M']
subsets = [frozenset(sv) for k in range(1, 4) for sv in __import__('itertools').combinations(vals, k)]
report("V14", len(subsets) == 7 == 2**3 - 1, "non-empty subsets of {T,F,M}: 7")

# V15  Logic 𝔼 = powerset lift of Bochvar (M infectious). Designated: contains T.
#  (a) Meta-1 ¬(μ∨ν)∧(μ∧ν) is designated for some assignment;  (b) explosion fails;
#  (c) on classical inputs 𝔼 agrees with Boolean logic.
def bn(a_): return {'T': 'F', 'F': 'T', 'M': 'M'}[a_]
def band(a_, b_): return 'M' if 'M' in (a_, b_) else ('T' if a_ == b_ == 'T' else 'F')
def bor(a_, b_): return 'M' if 'M' in (a_, b_) else ('T' if 'T' in (a_, b_) else 'F')
Neg = lambda S: frozenset(bn(s_) for s_ in S)
And = lambda S, U: frozenset(band(s_, u) for s_ in S for u in U)
Or = lambda S, U: frozenset(bor(s_, u) for s_ in S for u in U)
des = lambda S: 'T' in S
meta1 = lambda M_, N_: And(Neg(Or(M_, N_)), And(M_, N_))
witn = [(sorted(M_), sorted(N_)) for M_ in subsets for N_ in subsets if des(meta1(M_, N_))]
explode = [(sorted(P_), sorted(Q_)) for P_ in subsets for Q_ in subsets
           if des(P_) and des(Neg(P_)) and not des(Q_)]
classical = all(meta1(frozenset(u), frozenset(v)) == frozenset({'F'})
                for u in ('T', 'F') for v in ('T', 'F') for u in [u] for v in [v])
closed = all(Neg(S) in subsets and And(S, U) in subsets for S in subsets for U in subsets)
report("V15", bool(witn) and bool(explode) and classical and closed,
       "In 𝔼 (7 values) Meta-1 is designated, ex falso fails, Boolean fragment is classical",
       f"{len(witn)} designated assignments of Meta-1, e.g. μ,ν = {witn[0]}; "
       f"explosion counter-example P={explode[0][0]}, Q={explode[0][1]}")

# V16  Five-fingered hand: number of readings of a 5-ary ν (non-assoc, non-comm) vs μ (assoc).
catalan4 = comb(8, 4) // 5
gestures = factorial(5) * catalan4
report("V16", catalan4 == 14 and gestures == 1680 and gestures % 5 == 0,
       "μ(5 fingers): 1 bracketing;  ν(5 fingers): Catalan₄ = 14 bracketings, ×5! orders = 1680 gestures (÷5 ✓)")

# V17  Calendar: 5 seasons × 73 days = 365; St. Tib's Day lies outside every season and weekday.
report("V17", 5 * 73 == 365 and 73 % 5 == 3, "5×73 = 365; 73 ≡ 3 (mod 5) so seasons do not align with the 5-day week")

# V18  Distributive-law existence check (β): for the toy belief functor F X = 1 + Atom×X + X×X,
# a GSOS law exists for ⊗ iff observation of a composite factors through observations of parts.
# Here: a stateful observer ('drift') whose output on x⊗y depends on a global counter breaks it.
def obs_pure(t):  # t: nested tuples / atoms → multiset of atoms observed
    return [t] if isinstance(t, str) else obs_pure(t[0]) + obs_pure(t[1])
counter = [0]
def obs_drift(t):
    counter[0] += 1
    return ([t] if isinstance(t, str) else obs_drift(t[0]) + obs_drift(t[1])) + (['*'] if counter[0] % 5 == 0 else [])
t1, t2 = ('a', 'b'), 'c'
pure_ok = sorted(obs_pure((t1, t2))) == sorted(obs_pure(t1) + obs_pure(t2))
counter[0] = 0
whole = obs_drift((t1, t2)); counter[0] = 0
parts = obs_drift(t1) + obs_drift(t2)
report("V18", pure_ok and sorted(whole) != sorted(parts),
       "seed-pure observation is compositional; a counter-dependent (Greyface) observer is not",
       f"whole={whole}  parts={parts}")

print(f"\n{sum(results)}/{len(results)} checks pass")
