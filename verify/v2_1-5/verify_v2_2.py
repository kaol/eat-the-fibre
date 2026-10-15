"""Checks for 'Open Questions Explored (v2.2)'. Tags [Un]."""
from itertools import product, combinations, permutations
from collections import defaultdict
from fractions import Fraction as Q
from math import comb, factorial
import random
import numpy as np
import sympy as sp

res = []
def report(tag, ok, claim, detail=""):
    res.append(ok); print(f"[{tag}] {'PASS' if ok else 'FAIL'}  {claim}" + (f"\n       {detail}" if detail else ""))

# ═════════ Q2′  Folds compatible with letter deletion, beyond endomorphisms ═════════
# Work in the letters a and u = a − b of k<a,b> = k<a,u>. Letter deletion: Δa = 1⊗1, Δu = 0, ε-rule ⇒ Δ deletes one 'a'.
def words(n, al='au'): return [''.join(t) for k in range(n + 1) for t in product(al, repeat=k)]
def Dl(w): return [(w[:i], w[i+1:]) for i, ch in enumerate(w) if ch == 'a']
def clean(d): return {k: v for k, v in d.items() if v != 0}
def mul(x, y):
    r = defaultdict(int)
    for u, c in x.items():
        for v, e in y.items(): r[u + v] += c * e
    return clean(r)
def ap(P, x):
    r = defaultdict(int)
    for w, c in x.items():
        for v, e in P[w].items(): r[v] += c * e
    return clean(r)
def forced(P1, Pu, N):
    """Bigraded (a-degree, u-degree) P with Δ∘P = (P⊗P)∘Δ, given P(1) and P on k[u].
    Δ is injective on every bidegree with a-degree ≥ 1 (ker Δ = k[u]), so P is forced recursively."""
    P = {'': {'': 1} if P1 else {}}
    free_syms = 0
    for w in sorted(words(N), key=len)[1:]:
        if 'a' not in w: P[w] = {w: 1} if Pu else {}; continue
        T = defaultdict(int)
        for l, r in Dl(w):
            for x, c in P[l].items():
                for y, e in P[r].items(): T[(x, y)] += c * e
        T = clean(T)
        blk = [v for v in words(len(w)) if len(v) == len(w) and v.count('a') == w.count('a')]
        cs = sp.symbols(f'c0:{len(blk)}'); D = defaultdict(int)
        for v, c in zip(blk, cs):
            for k in Dl(v): D[k] += c
        sol = sp.solve([D[k] - T.get(k, 0) for k in set(D) | set(T)], cs, dict=True)
        if not sol: return None, 0
        s = sol[0]; free_syms += sum(1 for c in cs if c not in s)
        P[w] = clean({v: s.get(c, c) for v, c in zip(blk, cs)})
    return P, free_syms
def RB(P, x, y):          # P(x)P(y) = P(P(x)y + xP(y) − xy)
    px, py = ap(P, {x: 1}), ap(P, {y: 1})
    s = defaultdict(int)
    for k, v in mul(px, {y: 1}).items(): s[k] += v
    for k, v in mul({x: 1}, py).items(): s[k] += v
    s[x + y] -= 1
    return mul(px, py) == ap(P, clean(s))
N = 5; W2 = words(2)
sols = {}
for P1, Pu in product((1, 0), repeat=2):
    P, fr = forced(P1, Pu, N)
    ok = P is not None and fr == 0 and all(ap(P, ap(P, {w: 1})) == ap(P, {w: 1}) for w in words(N)) \
         and all(RB(P, x, y) for x in W2 for y in W2)
    sols[(P1, Pu)] = (ok, P)
report("U1", all(v[0] for v in sols.values()),
       "bigraded Δ-compatible idempotent RB(−1) operators: exactly four, each forced uniquely by (P(1), P|k[u]) ∈ {0,1}² — 0, id, identification fold, framed fold",
       "uniqueness: Δ injective off k[u]; P|k[u] all-or-nothing because {j : u^j ∈ A₁} and its complement are both additive semigroups")
import re
framed = re.compile(r'^u+(au+)*$')                 # u⁺(a u⁺)*: every a isolated and framed by u's
Pfr = sols[(0, 1)][1]
desc = all((Pfr[w] == {w: 1}) == bool(framed.match(w)) and Pfr[w] in ({w: 1}, {}) for w in words(N))
# direct check of the word-partition description up to length 7
W7 = words(7)
A1 = {w for w in W7 if framed.match(w)}; A2 = set(W7) - A1
sub1 = all(x + y in A1 for x in A1 for y in A1 if len(x + y) <= 7)
sub2 = all(x + y in A2 for x in A2 for y in A2 if len(x + y) <= 7)
co1 = all(l in A1 and r in A1 for w in A1 for l, r in Dl(w))
co2 = all(not (l in A1 and r in A1) for w in A2 for l, r in Dl(w))
report("U2", desc and sub1 and sub2 and co1 and co2,
       "the framed fold f□: Hodge = span u⁺(au⁺)* (every a isolated between u's), Podge = the rest incl. Void; both subalgebras; Δ(Hodge) ⊆ Hodge⊗Hodge, Podge a coideal (words ≤ 7)")
report("U3", ap(Pfr, {'uau': 1}) == {'uau': 1} and ap(Pfr, {'a': 1}) == {} and ap(Pfr, {'': 1}) == {},
       "f□ is not multiplicative (f□(uau) = uau ≠ f□(u)f□(a)f□(u) = 0) and f□(1) = 0: Void belongs to Podge")

# UT_n with arrow deletion: monomial (basis-subset) splittings
def ut_check(n):
    B = [(i, j) for i in range(n) for j in range(i, n)]
    m = lambda p, q: (p[0], q[1]) if p[1] == q[0] else None
    D = lambda i, j: [((i, k), (k + 1, j)) for k in range(i, j)]
    out = []
    for r in range(len(B) + 1):
        for S in combinations(B, r):
            S = set(S); C = set(B) - S
            if any((pq := m(p, q)) is not None and pq not in S for p in S for q in S): continue
            if any((pq := m(p, q)) is not None and pq not in C for p in C for q in C): continue
            if any(not (l in S and r_ in S) for x in S for l, r_ in D(*x)): continue
            if any(l in S and r_ in S for x in C for l, r_ in D(*x)): continue
            out.append(S)
    return out
def vertex_support(n):     # A₁ = paths whose every vertex lies in V′, for each V′ ⊆ vertices
    out = []
    for r in range(n + 1):
        for V in combinations(range(n), r):
            out.append({(i, j) for i in range(n) for j in range(i, n) if all(k in V for k in range(i, j + 1))})
    return out
key = lambda L: sorted(sorted(s) for s in L)
ut = {n: ut_check(n) for n in (2, 3, 4, 5)}
report("U4", all(key(ut[n]) == key(vertex_support(n)) for n in ut),
       "UT_n (arrow deletion): the monomial Δ-compatible splittings are exactly the 2ⁿ vertex supports — Hodge = paths that never leave V′ (n = 2..5)",
       f"counts: {[len(ut[n]) for n in ut]}")

# ═════════ Q1′  Expressivity moves to the index sets ═════════
# Weighted-automaton bound: two sources with n1, n2 states agree on all prefixes iff they agree on all prefixes of length ≤ n1+n2−1.
def rows(n):
    out = [tuple(Q(int(i == j)) for j in range(n)) for i in range(n)]
    out += [tuple(Q(1, 2) if j in pr else Q(0) for j in range(n)) for pr in combinations(range(n), 2)]
    return out
def rand_src(rng):
    n = rng.choice([1, 2, 3]); R = rows(n)
    return n, [rng.choice(R) for _ in range(n)], rng.choice(R), [rng.choice([Q(0), Q(1, 2), Q(1)]) for _ in range(n)]
def law(src, L):
    n, T, init, emit = src
    vecs = {'': list(init)}; out = {}
    for _ in range(L):
        nv = {}
        for w, vec in vecs.items():
            for s in (0, 1):
                e = [vec[i] * (emit[i] if s else 1 - emit[i]) for i in range(n)]
                out[w + str(s)] = sum(e)
                nv[w + str(s)] = [sum(e[i] * T[i][j] for i in range(n)) for j in range(n)]
        vecs = nv
    return out
rng = random.Random(7); worst = 0; viol = 0; tested = 0
for _ in range(6000):
    s1, s2 = rand_src(rng), rand_src(rng)
    L1, L2 = law(s1, 8), law(s2, 8)
    diff = [len(w) for w in L1 if L1[w] != L2[w]]
    if not diff: continue
    tested += 1; first = min(diff); bound = s1[0] + s2[0] - 1
    worst = max(worst, first - bound); viol += first > bound
report("U5", viol == 0, "inequivalent finite-state sources always disagree on some prefix of length ≤ n₁+n₂−1 (equivalence decidable from finitely many prefixes)",
       f"{tested} inequivalent random pairs, max(first disagreement − bound) = {worst}")

# ═════════ Q3′  The 𝔼 calculus: locality no-go, interpolation ═════════
_src = open(__file__.replace('verify_v2_2.py', 'audit_logic_topology.py')).read()
_ns = {}; exec(_src.split('# S1 ')[0], _ns)
ev, var, formulas, assigns, VALS = (_ns[k] for k in ('ev', 'var', 'formulas', 'assigns', 'VALS'))
def E(G, phi, atoms): return all('T' in ev(phi, v) for v in assigns(atoms, VALS) if all('T' in ev(g, v) for g in G))
def leaves(L, R):
    for i, g in enumerate(L):
        if not isinstance(g, str) and g[0] == '&': return leaves(L[:i] + [g[1], g[2]] + L[i+1:], R)
        if not isinstance(g, str) and g[0] == '|': return leaves(L[:i] + [g[1]] + L[i+1:], R) + leaves(L[:i] + [g[2]] + L[i+1:], R)
    for i, g in enumerate(R):
        if not isinstance(g, str) and g[0] == '|': return leaves(L, R[:i] + [g[1], g[2]] + R[i+1:])
        if not isinstance(g, str) and g[0] == '&': return leaves(L, R[:i] + [g[1]] + R[i+1:]) + leaves(L, R[:i] + [g[2]] + R[i+1:])
    return [(frozenset(L), frozenset(R))]
pq = ('|', 'p', 'q')
lv_valid, lv_invalid = set(leaves([pq], [pq])), set(leaves(['p'], [pq]))
report("U6", E([pq], pq, ['p', 'q']) and not E(['p'], pq, ['p', 'q']) and lv_invalid <= lv_valid,
       "no leaf-local axiom set works with the invertible LP rules: the leaves of the INVALID p ⇒ p∨q are a subset of the leaves of the VALID p∨q ⇒ p∨q",
       f"valid leaves {sorted(map(lambda t: (sorted(t[0]), sorted(t[1])), lv_valid))}; invalid leaves {sorted(map(lambda t: (sorted(t[0]), sorted(t[1])), lv_invalid))}")
at3 = ['p', 'q', 'r']; F2 = formulas(at3, 2); random.seed(3); ci = n_ = 0
for _ in range(3000):
    G = random.sample(F2, random.choice([1, 2])); phi = random.choice(F2)
    if E(G, phi, at3):
        n_ += 1; VG = set().union(*map(var, G))
        ci += var(phi) <= VG & var(phi) and E(G, phi, at3) and E([phi], phi, at3)
report("U7", n_ > 0 and ci == n_, f"Craig interpolation is trivial in 𝔼: the conclusion is its own interpolant ({n_} valid instances)")

# ═════════ Q5′  The noisy bureaucracy ═════════
forms = ['X', 'Y', 'Z']; lanes = ('fast', 'mid', 'slow'); structs = list(permutations(lanes))
def traffic(s, q=0.6):
    fast = forms[s.index('fast')]; return [q if f == fast else (1 - q) / 2 for f in forms]
def policy_dist(counts, efficient):
    """Rank forms by observed count (desc for efficient, asc for fair), ties broken uniformly; return law over structures."""
    key = (lambda i: -counts[i]) if efficient else (lambda i: counts[i])
    orders = [o for o in permutations(range(3)) if all(key(o[k]) <= key(o[k+1]) for k in range(2))]
    d = defaultdict(float)
    for o in orders:
        lane_of = {o[0]: 'fast', o[1]: 'mid', o[2]: 'slow'}
        d[tuple(lane_of[i] for i in range(3))] += 1 / len(orders)
    return d
def chain(M, efficient):
    T = np.zeros((6, 6))
    for a, s in enumerate(structs):
        p = traffic(s)
        for x in range(M + 1):
            for y in range(M + 1 - x):
                z = M - x - y
                pr = factorial(M) / (factorial(x) * factorial(y) * factorial(z)) * p[0]**x * p[1]**y * p[2]**z
                for s2, w in policy_dist((x, y, z), efficient).items(): T[a, structs.index(s2)] += pr * w
    return T
def spectrum(T): return sorted(np.linalg.eigvals(T).real, reverse=True)
eff = {M: spectrum(chain(M, True)) for M in (2, 5, 10, 20, 40)}
fai = {M: spectrum(chain(M, False)) for M in (2, 5, 10, 20, 40)}
lam2 = [eff[M][1] for M in eff]; lamneg = [fai[M][-1] for M in fai]
stat = np.linalg.matrix_power(chain(20, True), 4000)[0]
report("U8", all(np.diff(lam2) > 0) and lam2[-1] > 0.95 and abs(stat.reshape(-1).sum() - 1) < 1e-9,
       "efficient policy under sampling noise: the second eigenvalue climbs toward +1 as the traffic window M grows — lock-in becomes metastable switching",
       f"λ₂ for M = 2,5,10,20,40: {[round(x, 4) for x in lam2]}")
Tf = chain(40, False)
fast_of = lambda s: s.index('fast')
stay = np.mean([sum(Tf[a, b] for b, s2 in enumerate(structs) if fast_of(s2) == fast_of(s)) for a, s in enumerate(structs)])
report("U9", all(np.diff(lamneg) < 0) and abs(lamneg[-1] + 0.5) < 0.01 and stay < 0.01,
       "fair policy under noise: λ_min → −½, not −1; the fast lane never stays and jumps uniformly to another form (a walk on K₃), so the 2-cycle was a tie-breaking artifact",
       f"λ_min for M = 2,5,10,20,40: {[round(x, 4) for x in lamneg]}; P(fast form unchanged), M=40: {stay:.2e}")

# ═════════ Q4′  Empirical anchor: Glauber culture = random-field threshold model ═════════
# P(adopt) = ½(1+tanh(β(Jm+h)))  ⟺  adopt iff φᵢ + Jm + h > 0 with φᵢ ~ Logistic(0, 1/(2β)) (zero-temperature RFIM, Michard–Bouchaud)
xs = np.linspace(-5, 5, 2001); beta = 1.7
glaub = 0.5 * (1 + np.tanh(beta * xs)); logistic_cdf_tail = 1 / (1 + np.exp(-2 * beta * xs))
rngn = np.random.default_rng(0); phis = rngn.logistic(0, 1 / (2 * beta), 400000)
mc = [np.mean(phis + x > 0) for x in (-1.0, -0.3, 0.0, 0.4, 1.2)]
exact = [0.5 * (1 + np.tanh(beta * x)) for x in (-1.0, -0.3, 0.0, 0.4, 1.2)]
report("U10", np.allclose(glaub, logistic_cdf_tail) and np.allclose(mc, exact, atol=3e-3),
       "the v2.1 Glauber culture is exactly a zero-temperature random-field threshold model with logistic propensities of scale 1/(2β): K_c = 1 ⟺ J_c ∝ propensity spread")

print(f"\n{sum(res)}/{len(res)} checks pass")
