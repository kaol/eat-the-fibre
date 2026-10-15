"""Audit of Le Chat's v1.1–v1.3 claims on 𝔼 (Q5/NQ-3/NQ-4/OQ-5), topology (Q3/NQ-5/OQ-3),
culture kernels (Q1/NQ-6/OQ-7) and the Arrow-fiber witness (Q2). Tags [Sn]."""
from itertools import product, combinations
from math import log, ceil
import random

res = []
def report(tag, ok, claim, detail=""):
    res.append(ok); print(f"[{tag}] {'PASS' if ok else 'FAIL'}  {claim}" + (f"\n       {detail}" if detail else ""))

# ───────── 𝔼 ─────────
VALS = [frozenset(s) for k in (1, 2, 3) for s in combinations('TFM', k)]
bn = {'T': 'F', 'F': 'T', 'M': 'M'}
def band(a, b): return 'M' if 'M' in (a, b) else ('T' if a == b == 'T' else 'F')
def bor(a, b):  return 'M' if 'M' in (a, b) else ('T' if 'T' in (a, b) else 'F')
def ev(phi, v):              # phi: atom str | ('~',φ) | ('&',φ,ψ) | ('|',φ,ψ)
    if isinstance(phi, str): return v[phi]
    if phi[0] == '~': return frozenset(bn[x] for x in ev(phi[1], v))
    A, B = ev(phi[1], v), ev(phi[2], v)
    op = band if phi[0] == '&' else bor
    return frozenset(op(x, y) for x in A for y in B)
def var(phi): return {phi} if isinstance(phi, str) else set().union(*map(var, phi[1:]))
def show(phi):
    if isinstance(phi, str): return phi
    if phi[0] == '~': return '¬' + show(phi[1])
    return f"({show(phi[1])}{'∧' if phi[0]=='&' else '∨'}{show(phi[2])})"
def formulas(atoms, k):
    F = {0: list(atoms)}
    for n in range(1, k + 1):
        F[n] = [('~', g) for g in F[n-1]]
        for i in range(n):
            for g in F[i]:
                for h in F[n-1-i]:
                    F[n] += [('&', g, h), ('|', g, h)]
    return [g for n in range(k + 1) for g in F[n]]
def assigns(atoms, vals):
    for t in product(vals, repeat=len(atoms)): yield dict(zip(atoms, t))
def cons(G, phi, atoms, vals, des):
    return all(des(ev(phi, v)) for v in assigns(atoms, vals) if all(des(ev(g, v)) for g in G))
desE = lambda S: 'T' in S
LPV = [frozenset('T'), frozenset('F'), frozenset('TF')]
CLV = [frozenset('T'), frozenset('F')]
E = lambda G, phi, at: cons(G, phi, at, VALS, desE)
LP = lambda G, phi, at: cons(G, phi, at, LPV, desE)
CL = lambda G, phi, at: cons(G, phi, at, CLV, desE)

# S1  [W18] counterexample: disjunctive syllogism as a single premise (3 connectives: outside Le Chat's ≤2 corpus)
at2 = ['p', 'q']
g = ('&', ('|', 'p', 'q'), ('~', 'p'))
report("S1", CL([g], 'q', at2) and var('q') <= var(g) and not E([g], 'q', at2),
       "v1.1 [W18] is false: ((p∨q)∧¬p) ⊨_CL q with Var ⊆, but ⊭_𝔼 q (countermodel p={T,F}, q={F})")
g2 = ('&', ('&', 'p', ('~', 'p')), 'q')
report("S2", CL([g2], ('&', 'q', ('~', 'q')), at2) and not E([g2], ('&', 'q', ('~', 'q')), at2),
       "second counterexample: (p∧¬p)∧q ⊨_CL q∧¬q, Var ⊆, ⊭_𝔼 (local explosion fails)")

# S3  Theorem: Γ ⊨_𝔼 φ ⟺ Γ ⊨_LP φ ∧ Var(φ) ⊆ Var(Γ). Exhaustive single premise ≤3 connectives / 2 atoms
F3 = formulas(at2, 3)
mism = 0; nCLred = 0; n = 0
for g in F3:
    for phi in F3:
        n += 1
        e = E([g], phi, at2)
        lp = LP([g], phi, at2) and var(phi) <= var(g)
        mism += e != lp
        if lp and not CL([g], phi, at2): nCLred += 1
report("S3", mism == 0 and nCLred == 0,
       f"single premise, ≤3 connectives, 2 atoms ({n} pairs): ⊨_𝔼 = LP ∩ right-containment; CL filter never binds",
       f"mismatches={mism}")
# multi-premise sample, 3 atoms
at3 = ['p', 'q', 'r']; F2 = formulas(at3, 2)
random.seed(23); bad = 0; m = 0
for _ in range(4000):
    G = random.sample(F2, random.choice([1, 2, 3])); phi = random.choice(F2); m += 1
    vg = set().union(*map(var, G))
    bad += E(G, phi, at3) != (LP(G, phi, at3) and var(phi) <= vg)
report("S4", bad == 0, f"multi-premise, 3 atoms, {m} random (Γ,φ): ⊨_𝔼 = LP ∩ containment; Le Chat's CL filter is redundant (LP ⊆ CL)")

# S5  7 values collapse to 4: designation depends only on the classical part of atom values.
four = {frozenset(): 'M', frozenset('T'): 'T', frozenset('F'): 'F', frozenset('TF'): 'B'}
def collapse(v): return {k: (s - {'M'}) or frozenset('M') for k, s in v.items()}
ok5 = all(desE(ev(phi, v)) == desE(ev(phi, collapse(v))) for phi in formulas(at2, 2) for v in assigns(at2, VALS))
report("S5", ok5, "𝔼's designation factors through {T,F,B,M}: 𝔼 is the 4-valued matrix (LP on T,F,B + infectious M, D={T,B})")

# S6  coNP-hardness reduction: γ CL-unsat ⟺ γ ⊨_𝔼 ⋁_{p∈Var γ}(p∧¬p)   (Le Chat [W47]) — corpus ≤3 connectives
def glutdisj(g):
    vs = sorted(var(g)); d = ('&', vs[0], ('~', vs[0]))
    for p in vs[1:]: d = ('|', d, ('&', p, ('~', p)))
    return d
ok6 = all((not any(desE(ev(g, v)) for v in assigns(at2, CLV))) == E([g], glutdisj(g), at2) for g in F3)
report("S6", ok6, "coNP-hardness reduction verified on all ≤3-connective 2-atom formulas (claim [W47] stands)")

# S7  Meta-1: 4/49 designated (reproduces v1 [V15], Le Chat [W24])
m1 = sum(desE(ev(('&', ('~', ('|', 'm', 'n')), ('&', 'm', 'n')), {'m': a, 'n': b})) for a in VALS for b in VALS)
report("S7", m1 == 4, "Meta-1 designated in exactly 4/49 assignments")

# ───────── topology ─────────
def minper(w): return next(p for p in range(1, len(w)+1) if all(w[i] == w[i+p] for i in range(len(w)-p)))
W8 = [''.join(t) for t in product('01', repeat=8)]
dist = {}
for w in W8: dist[minper(w)] = dist.get(minper(w), 0) + 1
report("S8", dist == {1: 2, 2: 2, 3: 6, 4: 12, 5: 30, 6: 48, 7: 82, 8: 74} and sum(v for k, v in dist.items() if k <= 5) == 52,
       "minimal-period table [W10] and the 52 count [W7] reproduce", str(dist))
def hcost(w):
    best = 10**9
    for lu in range(len(w)):
        for lv in range(1, len(w) - lu + 1):
            u, v = w[:lu], w[lu:lu+lv]
            s = u + v * (len(w) // lv + 2)
            if s[:len(w)] == w: best = min(best, lu + 2*lv - 1)
    return best
hd = {}
for w in W8: hd[hcost(w)] = hd.get(hcost(w), 0) + 1
cheaper = sum(hcost(w) < 2*minper(w) - 1 for w in W8)
report("S9", hd == {1: 2, 2: 2, 3: 6, 4: 10, 5: 26, 6: 46, 7: 90, 8: 74} and cheaper == 236,
       "hand-off cost table [W36] and '236 cheaper' reproduce", f"{hd}, cheaper={cheaper}")
# S10  'every infinite ω-regular language is uncountable' (v1.1 Q3) is false: a*b^ω is countably infinite.
#  Check: a Büchi automaton for a*b^ω accepts exactly a^n b^ω (lasso enumeration up to n=6).
def buchi_acc(u, v):          # states 0 (reading a's), 1 (reading b's, accepting); deterministic
    st = 0
    for ch in u:
        st = {(0,'a'):0,(0,'b'):1,(1,'b'):1}.get((st, ch))
        if st is None: return False
    seen = []
    while (st,) not in seen:
        seen.append((st,))
        for ch in v:
            st = {(0,'a'):0,(0,'b'):1,(1,'b'):1}.get((st, ch))
            if st is None: return False
    return st == 1
acc = {(u, v) for n in range(7) for u in [''.join(t) for k in range(n+1) for t in product('ab', repeat=k)] for v in ['a','b','ab'] if buchi_acc(u, v)}
report("S10", all(set(v) == {'b'} and u.rstrip('b') == 'a'*len(u.rstrip('b')) for u, v in acc) and len(acc) > 6,
       "a*b^ω is ω-regular and countably infinite: v1.1's reason for 'not ω-regular' is false (conclusion still true via Calbrix–Nivat–Podelski)")

# ───────── culture kernels ─────────
# [W35]: wrong-basin probability is gambler's ruin: start x0 = ln(10/246), barriers ±ln 99, step s = β ln(θ/(1−θ)).
th = 0.75; rho = (1 - th) / th; x0 = log(10/246); top = log(99); bot = -log(99)
def p_true(beta):
    s = beta * log(th / (1 - th)); k1 = ceil((x0 - bot) / s - 1e-12); k2 = ceil((top - x0) / s - 1e-12)
    return (1 - rho**k1) / (1 - rho**(k1 + k2))
lechat = {0.2: 1.000, 0.4: 0.995, 0.8: 0.888, 1.2: 0.888, 1.6: 0.672, 2.0: 0.672}
pred = {b: round(p_true(b), 3) for b in lechat}
report("S11", all(abs(pred[b] - lechat[b]) < 0.03 for b in lechat),
       "[W35] curve = exact discrete gambler's ruin (1−ρ^k₁)/(1−ρ^(k₁+k₂)), k = ⌈distance/step⌉ — Dubins–Savage timid play, not a new effect",
       f"exact: {pred}\n       Le Chat MC: {lechat}")
report("S12", abs(log(10/246) + 3.2028) < 1e-3 and abs(log(99) - log(10/246) - 7.798) < 1e-3,
       "[W33] intercept −3.203 = ln(10/246) and [W34] constant 7.80 = ln 99 − ln(10/246): logistic-in-β is exact by construction (R² = 1 tautological)")

# ───────── Arrow witness ─────────
F_ = lambda s: ((s[0]+1) % 4, s[1], s[2]); H_ = lambda s: (s[0], s[1] ^ (s[0] == 3), s[2]); G_ = lambda s: (s[0], s[1], s[2] ^ s[1] ^ 1)
def run(seq, s):
    for h in seq: s = h(s)
    return s
st = [(c, i, b) for c in range(4) for i in (0, 1) for b in (0, 1)]
FG = all(run([F_, G_], s) == run([G_, F_], s) for s in st)
r1, r2 = run([F_, H_, G_], (2, 0, 0)), run([F_, G_, H_], (2, 0, 0))
report("S13", FG and (r1, r2) == ((3, 1, 0), (3, 1, 1)),
       "[W5–W6] reproduce, but F and G touch disjoint coordinates: the 'fiber Arrow' is plain commutativity, not the Arrow condition",
       f"FHG={r1}, FGH={r2}")

# S14  𝔼 = Deutsch's S_fde (N infectious, FDE/LP on T,B,F, D={T,B}); tables as in Petrukhin, EPTCS 358 (2022)
NEG = {'T': 'F', 'F': 'T', 'B': 'B', 'N': 'N'}
ORD = {'F': 0, 'B': 1, 'T': 2}
def sand(a, b): return 'N' if 'N' in (a, b) else min(a, b, key=ORD.get)
def sor(a, b):  return 'N' if 'N' in (a, b) else max(a, b, key=ORD.get)
def sev(phi, v):
    if isinstance(phi, str): return v[phi]
    if phi[0] == '~': return NEG[sev(phi[1], v)]
    return (sand if phi[0] == '&' else sor)(sev(phi[1], v), sev(phi[2], v))
def Sfde(G, phi, atoms):
    D = {'T', 'B'}
    return all(sev(phi, v) in D for v in assigns(atoms, 'TFBN') if all(sev(g, v) in D for g in G))
F2b = formulas(at2, 2)
mm = sum(E([g], phi, at2) != Sfde([g], phi, at2) for g in F3 for phi in F2b)
random.seed(5); mm2 = 0
for _ in range(3000):
    G = random.sample(F2, random.choice([1, 2])); phi = random.choice(F2)
    mm2 += E(G, phi, at3) != Sfde(G, phi, at3)
report("S14", mm == 0 and mm2 == 0, "𝔼 coincides with Deutsch's S_fde (single ≤3/≤2-connective pairs + 3000 random multi-premise, 3 atoms)",
       f"mismatches {mm}, {mm2}")
print(f"\n{sum(res)}/{len(res)} checks pass")
