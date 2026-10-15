"""Checks for claims new in discordianism-axioms v2 (tags [Nn]).
Re-uses the 𝔼 evaluator from audit_logic_topology.py without running its checks."""
from itertools import product, combinations, chain
from math import log, ceil, exp
import numpy as np

src = open(__file__.replace('verify_v2_new.py', 'audit_logic_topology.py')).read()
exec(src.split('# S1 ')[0])            # definitions only: ev, var, formulas, assigns, cons, VALS, LPV, CLV, …
res = []
def report(tag, ok, claim, detail=""):
    res.append(ok); print(f"[{tag}] {'PASS' if ok else 'FAIL'}  {claim}" + (f"\n       {detail}" if detail else ""))

at2 = ['p', 'q']; F2 = formulas(at2, 2); F1 = formulas(at2, 1)
def subsets(G): return chain.from_iterable(combinations(G, k) for k in range(len(G) + 1))
def clsat(G, atoms): return any(all('T' in ev(g, v) for g in G) for v in assigns(atoms, CLV))

# Weak Kleene (Bochvar) 3-valued: values {T},{F},{M} as singletons of 𝔼 — the lift restricted to singletons IS weak Kleene.
WK = [frozenset('T'), frozenset('F'), frozenset('M')]
PWK = lambda G, phi: cons(G, phi, at2, WK, lambda S: S != frozenset('F'))     # designated {T, M}
B3  = lambda G, phi: cons(G, phi, at2, WK, lambda S: S == frozenset('T'))     # designated {T}
CL  = lambda G, phi: cons(G, phi, at2, CLV, lambda S: 'T' in S)

# N1 Ciuni–Carrara: Γ ⊨_PWK φ ⟺ ∃Γ'⊆Γ. Var(Γ') ⊆ Var(φ) ∧ Γ' ⊨_CL φ   (PWK = left variable-inclusion companion of CL)
bad = n = 0
for G in [(g,) for g in F2] + list(combinations(F1, 2)):
    for phi in F2:
        n += 1
        cc = any(set().union(set(), *map(var, Gp)) <= var(phi) and CL(list(Gp), phi) for Gp in subsets(G))
        bad += PWK(list(G), phi) != cc
report("N1", bad == 0, f"Ciuni–Carrara characterization of PWK matches the weak-Kleene matrix on {n} (Γ,φ) pairs")

# N2 Bochvar B3 = right companion of CL: Γ ⊨_B3 φ ⟺ Γ CL-unsat ∨ (Γ ⊨_CL φ ∧ Var(φ) ⊆ Var(Γ))
bad = n = 0
for G in [(g,) for g in F2] + list(combinations(F1, 2)):
    for phi in F2:
        n += 1
        rc = (not clsat(G, at2)) or (CL(list(G), phi) and var(phi) <= set().union(*map(var, G)))
        bad += B3(list(G), phi) != rc
report("N2", bad == 0, f"B3 consequence = CL-unsat ∨ (CL ∧ right containment) on {n} pairs; 𝔼 = the same with LP for CL and no antitheorems")

# N3 Channel crystallization (exponential weights on a FIXED grid set of 256 words)
def minper(w): return next(p for p in range(1, len(w)+1) if all(w[i] == w[i+p] for i in range(len(w)-p)))
Wd = [''.join(t) for t in product('01', repeat=8)]
An = np.array([minper(w) <= 3 for w in Wd]); Fv = np.array([w.count('1') == 5 for w in Wd])
def run(sel, beta=0.8, T=200):
    lw = np.zeros(256)
    for _ in range(T): lw += beta * sel
    p = np.exp(lw - lw.max()); p /= p.sum()
    return p[An].sum(), float(-(p[p > 0] * np.log(p[p > 0])).sum())
neutral = run(np.zeros(256)); aner = run(An.astype(float)); eris = run((~An).astype(float)); five = run(Fv.astype(float))
ok = (abs(neutral[0] - 10/256) < 1e-12 and abs(neutral[1] - log(256)) < 1e-9 and aner[0] > 1 - 1e-9 and abs(aner[1] - log(10)) < 1e-6
      and eris[0] < 1e-9 and abs(eris[1] - log(246)) < 1e-6 and abs(five[1] - log(56)) < 1e-6 and An.sum() == 10 and Fv.sum() == 56)
report("N3", ok, "fixed grid set, flowing measure: neutral stays at ln256 nats; biased channels crystallize onto exactly their selected sets (ln10, ln246, ln56)",
       f"P(A): neutral {neutral[0]:.4f}, aneristic {aner[0]:.4f}, eristic {eris[0]:.2e}")

# N4 Wrong-basin probability is monotone non-decreasing in step size s (exact discrete gambler's ruin)
x0, top, bot = log(10/246), log(99), -log(99)
def p_wrong(s, th):
    rho = (1 - th) / th; k1 = ceil((x0 - bot) / s - 1e-12); k2 = ceil((top - x0) / s - 1e-12)
    return (rho**k1 - rho**(k1 + k2)) / (1 - rho**(k1 + k2))
def p_diff(s, th):            # diffusion limit: drift (2θ−1)s, variance s² per step
    a, b, lam = x0 - bot, top - x0, 2 * (2*th - 1) / s
    return (exp(-lam*a) - exp(-lam*(a+b))) / (1 - exp(-lam*(a+b)))
mono_diff = True; worst_drop = 0.0
for th in (0.55, 0.6, 0.75, 0.9):
    S = np.linspace(0.05, 6, 600)
    pd = [p_diff(s, th) for s in S]; mono_diff &= all(b >= a - 1e-15 for a, b in zip(pd, pd[1:]))
    pl = [p_wrong(s, th) for s in S]
    worst_drop = max(worst_drop, max(a - b for a, b in zip(pl, pl[1:])))
# lattice formula P(k1,k2) = ρ^k1 (1−ρ^k2)/(1−ρ^(k1+k2)): ↓ in k1 (near barrier farther), ↑ in k2 (far barrier farther)
Pk = lambda k1, k2, r: r**k1 * (1 - r**k2) / (1 - r**(k1 + k2))
from fractions import Fraction as Q          # exact: floats saturate at ρ^k ≈ 1e-17
part = all(Pk(k1 + 1, k2, r) < Pk(k1, k2, r) < Pk(k1, k2 + 1, r)
           for r in (Q(1, 10), Q(1, 3), Q(2, 3), Q(9, 11)) for k1 in range(1, 30) for k2 in range(1, 30))
report("N4", mono_diff and part,
       "P(wrong basin) increases with step size in the diffusion limit; on the lattice it falls in k₁ (near-barrier steps) and rises in k₂",
       f"so the lattice curve is NOT monotone in general: largest decrease found {worst_drop:.3f} (both k's shrink as the step grows)")

# N5 kG (group-likes) violates Aguiar's ε-rule: Δ(gh)=gh⊗gh vs g⊗gh + gh⊗h  (checked on Z/5 as formal tensors)
viol = any({(g*h % 5, g*h % 5)} != {(g, g*h % 5), (g*h % 5, h)} for g in range(1, 5) for h in range(1, 5))
report("N5", viol, "group Hopf algebras are not ε-bialgebras: Drinfeld vs Aguiar quasitriangularity are different questions")

# N6 ⋆-annihilator strata over the completion family (free algebra, x in span of words ≤3, y over words ≤3)
import sympy as sp
Wk = [''] + [''.join(t) for n_ in (1, 2, 3) for t in product('ab', repeat=n_)]
def fw(w): return set(w) <= {'a'}
def starw(u, v, b):      # ⋆ = x f(y) + b1 xy + b2 f(x)y + b3 x f(y) + b4 f(x)f(y), on words → {word: coeff}
    out = {}
    def put(w, c): out[w] = out.get(w, 0) + c
    if fw(v): put(u + v, 1)
    put(u + v, b[0])
    if fw(u): put(u + v, b[1])
    if fw(v): put(u + v, b[2])
    if fw(u) and fw(v): put(u + v, b[3])
    return out
def ann(b):
    cs = sp.symbols(f'c0:{len(Wk)}'); eqs = {}
    for y in Wk:
        for i, x in enumerate(Wk):
            for w, c in starw(x, y, b).items(): eqs[(y, w)] = eqs.get((y, w), 0) + c * cs[i]
            for w, c in starw(y, x, b).items(): eqs[(y, w)] = eqs.get((y, w), 0) - c * cs[i]
    sol = sp.solve([e for e in eqs.values() if e != 0], cs, dict=True)
    sol = sol[0] if sol else {}
    return sorted(Wk[i] for i, c in enumerate(cs) if c not in sol)      # free coefficients = annihilator support
L1 = lambda t: (-1 - t, 1 + t, t, -t); L2 = lambda t: (-t, 1 + t, t, -1 - t)
imf = ['', 'a', 'aa', 'aaa']
got = {'L1(t), t≠−1': all(ann(L1(t)) == imf for t in (0, 1, 2, -2)), 'L1(−1)': ann(L1(-1)) == sorted(Wk),
       'L2(t)': all(ann(L2(t)) == [''] for t in (0, 1, -1, 2)), 'points': ann((0, 0, 0, 0)) == [] and ann((0, 1, -1, 0)) == []}
report("N6", all(got.values()), "⋆-annihilator strata: line 1 → im f (all of A at t=−1), line 2 → scalars, the two points → 0 (words ≤3)", str(got))

print(f"\n{sum(res)}/{len(res)} checks pass")
