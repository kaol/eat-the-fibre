"""Checks for 'Open Questions Explored (v2.3)'. Tags [Wn] are reserved by Le Chat; this round uses [Xn]."""
from itertools import product, combinations, permutations
from collections import defaultdict
from math import factorial, cos, sin, pi
import cmath
import numpy as np
import sympy as sp

res = []
def report(tag, ok, claim, detail=""):
    res.append(ok); print(f"[{tag}] {'PASS' if ok else 'FAIL'}  {claim}" + (f"\n       {detail}" if detail else ""))
u = sp.symbols('u')

# ═════════ 1. Splittings of k[u] ═════════
# Lemma: in any splitting k[u] = B₁ ⊕ B₂ (subalgebras), 1 ∈ B₁ or 1 ∈ B₂.
#   1 = e₁ + e₂; if e₁ = t nonconstant then t − t² ∈ tk[t] ⊆ B₁ and t − t² = (1−t) − (1−t)² ∈ (1−t)k[1−t] ⊆ B₂.
t = sp.symbols('t')
report("X1", sp.expand((t - t**2) - ((1 - t) - (1 - t)**2)) == 0,
       "unit lemma: 1 lies in one summand (t − t² is simultaneously a positive-power polynomial in e₁ and in e₂ = 1 − e₁)")
# Dihedral splitting: k[u] = k[u²] ⊕ (u²+u)·k[u²+u]
def dihedral_ok(a, b, s, m):
    f, g = (u - a)**2, (u - b)**2
    gens = [sp.expand(f**i) for i in range(m + 1)] + [sp.expand((g - s) * g**j) for j in range(m)]
    M = sp.Matrix([[sp.Poly(p, u).coeff_monomial(u**k) for k in range(2 * m + 1)] for p in gens])
    return M.shape == (2 * m + 1, 2 * m + 1) and M.rank() == 2 * m + 1
cases = [(0, sp.Rational(-1, 2), sp.Rational(1, 4), 8), (0, 3, 0, 8), (1, -2, 5, 7), (0, sp.Rational(-1, 2), 0, 10)]
report("X2", all(dihedral_ok(*c) for c in cases),
       "DIHEDRAL splittings exist: k[u] = k[(u−a)²] ⊕ ((u−b)²−s)·k[(u−b)²] for a ≠ b, any s — basis check in every degree ≤ 2m",
       "e.g. k[u²] ⊕ (u²+u)k[u²+u]; the two reflections u ↦ 2a−u, u ↦ 2b−u generate D∞, whose invariants are only k")
def covers(f, g, N=10, Mdeg=30):
    gens = [sp.expand(f**i) for i in range(Mdeg // sp.degree(f, u) + 1)] + [sp.expand(g**j) for j in range(1, Mdeg // sp.degree(g, u) + 1)]
    D = max(sp.degree(p, u) for p in gens)
    M = sp.Matrix([[sp.Poly(p, u).coeff_monomial(u**k) for k in range(D + 1)] for p in gens]).T
    r = M.rank()
    return all(M.row_join(sp.Matrix([1 if k == j else 0 for k in range(D + 1)])).rank() == r for j in range(N + 1))
pairs = [(u**2, u**2 + u, True), (u**2, (u - 3)**2, True), (u**2, u**3 + u, False), (u**3, (u + 1)**3, False),
         (u**3, u**3 + u**2 + u, False), (u**2, u**4 + u, False), (u**2 + u, u**3, False), (u**2, u**2, False)]
report("X3", all(covers(f, g) == exp for f, g, exp in pairs),
       "density evidence: k[f] + k[g] reaches all of k[u] (cancellations allowed to degree 30) only for two quadratics with distinct centres")

# ═════════ 2. Deconstructing folds on k<a,b>: run-wise extension + conjugation ═════════
def P_on_ku(kind, n):
    if kind == 'dihedral':          # B₁ = k[u²], B₂ = (u²+u)k[u²+u]
        g = u**2 + u; es = sp.symbols(f'e0:{n+1}'); rs = sp.symbols(f'r1:{n+2}')
        expr = sum(es[i] * u**(2*i) for i in range(n + 1)) + sum(rs[j-1] * g**j for j in range(1, n + 2)) - u**n
        s = sp.solve(sp.Poly(sp.expand(expr), u).all_coeffs(), list(es) + list(rs), dict=True)[0]
        return sp.expand(sum(s.get(es[i], 0) * u**(2*i) for i in range(n + 1)).subs({x: 0 for x in es + rs}))
    if kind == 'char':              # B₁ = k, B₂ = (u−c): evaluation at c (c = 2)
        return sp.Integer(2)**n
    if kind == 'framed': return sp.Integer(0) if n == 0 else u**n
    if kind == 'swapdih':           # B₁ = (u²+u)k[u²+u], B₂ = k[u²]  (roles exchanged: 1 ∈ B₂)
        return sp.expand(u**n - P_on_ku('dihedral', n))
cache = {}
def Pk(kind, n):
    if (kind, n) not in cache: cache[(kind, n)] = P_on_ku(kind, n)
    return cache[(kind, n)]
def poly_to_words(p):
    p = sp.Poly(sp.expand(p), u); return {'u' * m[0]: c for m, c in zip(p.monoms(), p.coeffs())} if not p.is_zero else {}
def clean(d): return {k: v for k, v in d.items() if v != 0}
def mulw(x, y):
    r = defaultdict(int)
    for a_, c in x.items():
        for b_, e in y.items(): r[a_ + b_] += c * e
    return clean(r)
def runs(w): return [len(r) for r in w.split('a')]
def Phat(kind, x):                  # run-wise extension: P on each maximal u-run, letters a kept
    out = defaultdict(int)
    for w, c in x.items():
        acc = {'': 1}
        for k, r in enumerate(runs(w)):
            if k: acc = mulw(acc, {'a': 1})
            acc = mulw(acc, poly_to_words(Pk(kind, r)))
        for k_, v in acc.items(): out[k_] += c * v
    return clean(out)
def Dl(x):
    r = defaultdict(int)
    for w, c in x.items():
        for i, ch in enumerate(w):
            if ch == 'a': r[(w[:i], w[i+1:])] += c
    return clean(r)
def tensP(kind, T):
    r = defaultdict(int)
    for (l, rr), c in T.items():
        for l2, c2 in Phat(kind, {l: 1}).items():
            for r2, c3 in Phat(kind, {rr: 1}).items(): r[(l2, r2)] += c * c2 * c3
    return clean(r)
def addd(*ds):
    r = defaultdict(int)
    for d in ds:
        for k, v in d.items(): r[k] += v
    return clean(r)
def neg(d): return {k: -v for k, v in d.items()}
W3 = [''.join(t_) for n_ in range(4) for t_ in product('au', repeat=n_)]
W2 = [w for w in W3 if len(w) <= 2]
def fold_ok(kind):
    idem = all(Phat(kind, Phat(kind, {w: 1})) == Phat(kind, {w: 1}) for w in W3)
    rb = all(mulw(Phat(kind, {x: 1}), Phat(kind, {y: 1})) ==
             Phat(kind, addd(mulw(Phat(kind, {x: 1}), {y: 1}), mulw({x: 1}, Phat(kind, {y: 1})), {x + y: -1}))
             for x in W2 for y in W2)
    dc = all(Dl(Phat(kind, {w: 1})) == tensP(kind, Dl({w: 1})) for w in W3)
    return idem, rb, dc
outs = {k: fold_ok(k) for k in ('dihedral', 'swapdih', 'char', 'framed')}
report("X4", all(all(v) for v in outs.values()),
       "run-wise extension: for ANY splitting of k[u] (character, framed, dihedral, swapped dihedral), P̂ = P on each u-run is an idempotent RB(−1) fold commuting with letter deletion",
       f"{outs}")
# Conjugation: φ_{c,g} (a ↦ a + g(u), u ↦ c) equals θ_g ∘ P̂_char ∘ θ_g⁻¹ with θ_g: a ↦ a + g(u)
def subst(img, x):
    out = defaultdict(int)
    for w, c in x.items():
        acc = {'': 1}
        for ch in w: acc = mulw(acc, img[ch])
        for k, v in acc.items(): out[k] += c * v
    return clean(out)
gpoly = {'u': -2, 'uu': 1}                         # g(u) = u² − 2u, g(2) = 0
theta = {'a': addd({'a': 1}, gpoly), 'u': {'u': 1}}
theta_inv = {'a': addd({'a': 1}, neg(gpoly)), 'u': {'u': 1}}
phi = {'a': addd({'a': 1}, gpoly), 'u': {'': 2}}
conj_ok = all(subst(theta, Phat('char', subst(theta_inv, {w: 1}))) == subst(phi, {w: 1}) for w in W3)
report("X5", conj_ok,
       "v2.1's endomorphism folds are conjugates θ_g P̂ θ_g⁻¹ of the character run-wise extension by Δ-automorphisms θ_g: a ↦ a + g(u)")

# ═════════ 3. Path algebras UT_n: vertex splittings + arrow compatibility ═════════
def ut_folds(n):
    B = [(i, j) for i in range(n) for j in range(i, n)]; idx = {b: k for k, b in enumerate(B)}; N = len(B)
    def mulv(x, y):
        r = [0] * N
        for i, p in enumerate(B):
            if x[i] == 0: continue
            for j, q in enumerate(B):
                if y[j] != 0 and p[1] == q[0]: r[idx[(p[0], q[1])]] += x[i] * y[j]
        return r
    Dl_ = lambda i, j: [((i, k), (k + 1, j)) for k in range(i, j)]
    subsets = [frozenset(s) for r in range(1, n + 1) for s in combinations(range(n), r)]
    fams = [[]]
    for s in subsets:
        fams += [f + [s] for f in fams if all(not (s & t_) for t_ in f)]
    ind = lambda s: sp.Matrix([1 if v in s else 0 for v in range(n)])
    splits = []
    for F1 in fams:
        for F2 in fams:
            if len(F1) + len(F2) == n and sp.Matrix.hstack(*[ind(s) for s in F1 + F2]).rank() == n: splits.append((F1, F2))
    seen = set(); uniq = []
    for F1, F2 in splits:
        kk = (frozenset(F1), frozenset(F2))
        if kk not in seen: seen.add(kk); uniq.append((F1, F2))
    splits = uniq
    found = set(); compat = set()
    E = lambda b: [1 if k == idx[b] else 0 for k in range(N)]
    for F1, F2 in splits:
        Mx = sp.Matrix.hstack(*[ind(s) for s in F1 + F2]); proj = Mx * sp.diag(*([1] * len(F1) + [0] * len(F2))) * Mx.inv()
        key = (frozenset(F1), frozenset(F2))
        Pe = [proj[:, v] for v in range(n)]
        if all(Pe[i].is_zero_matrix or Pe[i+1].is_zero_matrix or any(
               all(Pe[i][w] == 0 for w in range(n) if w != k) and all(Pe[i+1][w] == 0 for w in range(n) if w != k + 1)
               for k in range(n - 1)) for i in range(n - 1)):
            compat.add(key)
        unknowns = []; cols = []
        for b in B:
            if b[0] == b[1]:
                v = [0] * N
                for w in range(n): v[idx[(w, w)]] = proj[w, b[0]]
                cols.append(v)
            else:
                xs = sp.symbols(f'x_{b[0]}{b[1]}_0:{N}'); unknowns += xs; cols.append(list(xs))
        P = lambda v: [sum(cols[k][r] * v[k] for k in range(N)) for r in range(N)]
        eqs = []
        for b in B:
            eqs += [a - c for a, c in zip(P(P(E(b))), P(E(b)))]
            for c in B:
                Px, Py = P(E(b)), P(E(c))
                inner = [a + bb - cc for a, bb, cc in zip(mulv(Px, E(c)), mulv(E(b), Py), mulv(E(b), E(c)))]
                eqs += [a - r for a, r in zip(mulv(Px, Py), P(inner))]
            Px = P(E(b)); T = defaultdict(int)
            for k, coef in enumerate(Px):
                p = B[k]
                if p[0] < p[1]:
                    for l, r in Dl_(*p): T[(l, r)] += coef
            S = defaultdict(int)
            for l, r in (Dl_(*b) if b[0] < b[1] else []):
                Pl, Pr = P(E(l)), P(E(r))
                for i1 in range(N):
                    for j1 in range(N): S[(B[i1], B[j1])] += Pl[i1] * Pr[j1]
            for k_ in set(T) | set(S): eqs.append(T[k_] - S[k_])
        eqs = [sp.expand(e) for e in eqs if sp.expand(e) != 0]
        if sp.solve(eqs, unknowns, dict=True): found.add(key)
    return len(splits), found, compat
s2, f2, c2 = ut_folds(2); s3, f3, c3 = ut_folds(3)
report("X6", f2 == c2 and f3 == c3 and len(f2) == 6 and len(f3) == 23,
       "UT_n deconstructing folds = splittings of the vertex algebra kⁿ satisfying ARROW COMPATIBILITY, each extending to a fold (exhaustive, n = 2, 3)",
       f"n=2: {len(f2)} of {s2} splittings; n=3: {len(f3)} of {s3}. Compatibility: for each arrow i→i+1, P(eᵢ)=0, P(eᵢ₊₁)=0, or (P(eᵢ),P(eᵢ₊₁)) ∝ (e_k, e_(k+1))")

# ═════════ 4. Borel question: semialgebraic termination ═════════
def lfp(p, iters=20000):            # least fixed point of x = (1−p) + p x²  (return of a p-biased walk)
    x = 0.0
    for _ in range(iters): x = (1 - p) + p * x * x
    return x
ps = [0.1, 0.3, 0.45, 0.5, 0.55, 0.7, 0.9]
report("X7", all(abs(lfp(p) - min(1, (1 - p) / p)) < 2e-3 for p in ps),
       "termination probabilities of pushdown/one-counter grids are least fixed points of polynomial systems, hence SEMIALGEBRAIC in the coin parameter (here min(1,(1−p)/p), kinked at ½)")

# ═════════ 5. Bureaucracies: symmetry decides whether there can be a clock ═════════
def macro_chain(k, M, rule, q=0.6):
    """Markov chain on 'which form holds the fast lane'. Traffic: fast form q, others (1−q)/(k−1).
    rule(counts, current) -> distribution over next fast form."""
    T = np.zeros((k, k))
    def comps(M, k):
        if k == 1: yield (M,); return
        for x in range(M + 1):
            for rest in comps(M - x, k - 1): yield (x,) + rest
    for cur in range(k):
        p = [q if f == cur else (1 - q) / (k - 1) for f in range(k)]
        for c in comps(M, k):
            pr = factorial(M) * np.prod([p[f]**c[f] / factorial(c[f]) for f in range(k)])
            for nxt, w in rule(c, cur).items(): T[cur, nxt] += pr * w
    return T
def argset(c, best):
    m = best(c); S = [f for f in range(len(c)) if c[f] == m]; return {f: 1 / len(S) for f in S}
most = lambda c, cur: argset(c, max)
least = lambda c, cur: argset(c, min)
def second(c, cur):
    vals = sorted(set(c), reverse=True); v = vals[1] if len(vals) > 1 else vals[0]
    S = [f for f in range(len(c)) if c[f] == v]; return {f: 1 / len(S) for f in S}
def successor(c, cur):                 # fast lane to the cyclic successor of the most-filed form
    d = argset(c, max); return {(f + 1) % len(c): w for f, w in d.items()}
sym_ok = True; worst = {}
for k in (3, 4, 5):
    for name, rule in (('most', most), ('least', least), ('second', second)):
        ev = np.linalg.eigvals(macro_chain(k, 12 if k < 5 else 8, rule))
        sym_ok &= np.max(np.abs(ev.imag)) < 1e-9 and ev.real.min() >= -1 / (k - 1) - 1e-9
        worst[(k, name)] = round(ev.real.min(), 4)
report("X8", sym_ok,
       "S_k-symmetric policies: the fast-lane chain commutes with all permutations, so its spectrum is {1, 1−kβ} ⊂ [−1/(k−1), 1], all real: no clock is possible",
       f"min eigenvalues: {worst}")
clock = {}
clk_ok = True
for k in (3, 4, 5):
    mods = []
    for M in (4, 8, 12):
        if k == 5 and M > 8: continue
        ev = np.linalg.eigvals(macro_chain(k, M, successor))
        z = max((e for e in ev if abs(e.imag) > 1e-9), key=abs)
        mods.append(round(abs(z), 4)); ang = abs(cmath.phase(z))
        clk_ok &= abs(ang - 2 * pi / k) < 1e-6
    clk_ok &= all(b > a for a, b in zip(mods, mods[1:]))
    clock[k] = mods
report("X9", clk_ok,
       "a C_k-equivariant 'successor of the most-filed' policy is circulant: complex eigenvalues at exactly e^{±2πi/k}, modulus climbing toward 1 with the window — a data-driven clock",
       f"|λ| by window M: {clock}")

print(f"\n{sum(res)}/{len(res)} checks pass")
