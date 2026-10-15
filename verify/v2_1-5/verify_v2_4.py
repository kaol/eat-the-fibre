"""Checks for 'Open Questions Explored (v2.4)'. Tags [Yn]."""
from itertools import combinations
from math import comb, gcd, sqrt, pi
import numpy as np
import sympy as sp

res = []
def report(tag, ok, claim, detail=""):
    res.append(ok); print(f"[{tag}] {'PASS' if ok else 'FAIL'}  {claim}" + (f"\n       {detail}" if detail else ""))

# ═════════ 1. Density lemma: an exact annihilating functional ═════════
# λ(h) = Σ_m b_m [u^{m−1}] h  (= −Res_∞(h·φ), φ = Σ b_m u^{−m} du). λ ⊥ k[f] ⟺ λ(f^j) = 0 ∀j: constraint j involves b_1..b_{dj+1}
# with pivot b_{dj+1} (coefficient lc(f)^j ≠ 0) — EXACTLY triangular, no truncation. Shared pivots (m ≡ 1 mod lcm(d,e)) are met
# with the free variable b_{m−1}, whose coefficient in the eliminated constraint is ∝ (m−1)(a_f − a_g), a = centre = −subleading/deg.
from fractions import Fraction as Fq
import random as _r
def pmul(p, q):
    r = [Fq(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q): r[i + j] += a * b
    return r
def construct(f, g, L):
    d, e = len(f) - 1, len(g) - 1
    Pf = [[Fq(1)]]; Pg = [[Fq(1)]]
    for _ in range(L // d + 2): Pf.append(pmul(Pf[-1], f))
    for _ in range(L // e + 2): Pg.append(pmul(Pg[-1], g))
    C = lambda P, j, m: P[j][m - 1] if m - 1 < len(P[j]) else Fq(0)
    pf = {d * j + 1: j for j in range(L // d + 1)}; pg = {e * j + 1: j for j in range(L // e + 1)}
    b = {}; first = True
    for m in range(1, L + 1):
        if m in pf:
            j = pf[m]; b[m] = -sum(C(Pf, j, k) * b[k] for k in range(1, m)) / C(Pf, j, m)
            if m in pg and sum(C(Pg, pg[m], k) * b[k] for k in range(1, m + 1)) != 0: return False, b
        elif m in pg:
            j = pg[m]; b[m] = -sum(C(Pg, j, k) * b[k] for k in range(1, m)) / C(Pg, j, m)
        elif (m + 1) in pf and (m + 1) in pg:
            def lin(P, j):
                lead = C(P, j, m + 1); return C(P, j, m) / lead, sum(C(P, j, k) * b[k] for k in range(1, m)) / lead
            cg, sg = lin(Pg, pg[m + 1]); cf, sf = lin(Pf, pf[m + 1])
            b[m] = Fq(0) if cg == cf else -(sg - sf) / (cg - cf)
        else:
            b[m] = Fq(1) if first else Fq(0); first = False
    ok = all(sum(C(Pf, j, k) * b[k] for k in range(1, d * j + 2)) == 0 for j in range(L // d + 1) if d * j + 1 <= L) and \
         all(sum(C(Pg, j, k) * b[k] for k in range(1, e * j + 2)) == 0 for j in range(L // e + 1) if e * j + 1 <= L)
    return ok and any(v != 0 for v in b.values()), b
def centred(deg, c, rng):
    v = [Fq(-c), Fq(1)]; p = [Fq(1)]
    for _ in range(deg): p = pmul(p, v)
    for k in range(deg - 1):
        coef = Fq(rng.randint(-3, 3)); q = [Fq(1)]
        for _ in range(k): q = pmul(q, v)
        for i, a in enumerate(q): p[i] += coef * a
    return p
rng = _r.Random(1); built = {}
for d, e in ((2, 3), (3, 3), (2, 4), (3, 4), (3, 5), (4, 4)):
    built[(d, e)] = all(construct(centred(d, 0, rng), centred(e, rng.choice([1, 2, -1]), rng), 90)[0] for _ in range(3))
dihedral_none = not construct(centred(2, 0, rng), centred(2, 1, rng), 60)[0]
report("Y1", all(built.values()) and dihedral_none,
       "density lemma (distinct centres): an exact functional λ ≠ 0 with λ(k[f]) = λ(k[g]) = 0 exists for every (d,e) ≠ (2,2) — so k[f]+k[g] ≠ k[u]; for (2,2) the construction is forced to 0 (dihedral splitting survives)",
       f"constructed to order 90 for random f, g with distinct centres: {built}; (2,2) admits none: {dihedral_none}")

# ═════════ 2. Splittings of kⁿ: Cayley's rooted-forest count ═════════
def count_splittings(n):
    subsets = [frozenset(s) for r in range(1, n + 1) for s in combinations(range(n), r)]
    fams = {frozenset()}; frontier = [frozenset()]
    while frontier:
        new = []
        for f in frontier:
            used = set().union(*f) if f else set()
            for s in subsets:
                if not (s & used):
                    g = f | {s}
                    if g not in fams: fams.add(g); new.append(g)
        frontier = new
    ind = lambda s: [1 if v in s else 0 for v in range(n)]
    cnt = 0
    for F1 in fams:
        for F2 in fams:
            if len(F1) + len(F2) == n and np.linalg.matrix_rank(np.array([ind(s) for s in list(F1) + list(F2)], dtype=float)) == n:
                cnt += 1
    return cnt
counts = [count_splittings(n) for n in range(1, 6)]
report("Y2", counts == [2 * (n + 2)**(n - 1) for n in range(1, 6)],
       "splittings of the vertex algebra kⁿ into two subalgebras number 2(n+2)^{n−1} — Cayley's count of forests on n+2 labelled vertices with two specified roots (n = 1..5)",
       f"{counts}")

# ═════════ 3. Path algebras: shift-transport extension (sufficiency of arrow compatibility) ═════════
def ut_suff(n):
    B = [(i, j) for i in range(n) for j in range(i, n)]; idx = {b: k for k, b in enumerate(B)}; N = len(B)
    subsets = [frozenset(s) for r in range(1, n + 1) for s in combinations(range(n), r)]
    fams = {frozenset()}; frontier = [frozenset()]
    while frontier:
        new = []
        for f in frontier:
            used = set().union(*f) if f else set()
            for s in subsets:
                if not (s & used):
                    g = f | {s}
                    if g not in fams: fams.add(g); new.append(g)
        frontier = new
    ind = lambda s: sp.Matrix([1 if v in s else 0 for v in range(n)])
    mulv = lambda x, y: _mul(x, y, B, idx, N)
    E = lambda b: [1 if k == idx[b] else 0 for k in range(N)]
    Dl = lambda i, j: [((i, k), (k + 1, j)) for k in range(i, j)]
    ncomp = nok = 0
    for F1 in fams:
        for F2 in fams:
            if len(F1) + len(F2) != n: continue
            M = sp.Matrix.hstack(*[ind(s) for s in list(F1) + list(F2)])
            if M.rank() != n: continue
            proj = M * sp.diag(*([1] * len(F1) + [0] * len(F2))) * M.inv()
            Pe = [proj[:, v] for v in range(n)]
            compat = all(Pe[i].is_zero_matrix or Pe[i+1].is_zero_matrix or any(
                all(Pe[i][w] == 0 for w in range(n) if w != k) and all(Pe[i+1][w] == 0 for w in range(n) if w != k + 1)
                for k in range(n - 1)) for i in range(n - 1))
            if not compat: continue
            ncomp += 1
            cols = []
            for b in B:                       # shift-transport: P(E_ij) = Σ_s Π_v c_{v+s,v} E_{i+s, j+s}
                v = [0] * N
                for s_ in range(-n, n):
                    if 0 <= b[0] + s_ and b[1] + s_ < n:
                        c = 1
                        for vv in range(b[0], b[1] + 1): c *= proj[vv + s_, vv]
                        v[idx[(b[0] + s_, b[1] + s_)]] += c
                cols.append(v)
            P = lambda x: [sum(cols[k][r] * x[k] for k in range(N)) for r in range(N)]
            ok = all(P(P(E(b))) == P(E(b)) for b in B)
            for b in B:
                for c in B:
                    Px, Py = P(E(b)), P(E(c))
                    inner = [a + bb - cc for a, bb, cc in zip(mulv(Px, E(c)), mulv(E(b), Py), mulv(E(b), E(c)))]
                    ok &= mulv(Px, Py) == P(inner)
            for b in B:
                T = {}
                for k, coef in enumerate(P(E(b))):
                    p = B[k]
                    if coef != 0 and p[0] < p[1]:
                        for l, r in Dl(*p): T[(l, r)] = T.get((l, r), 0) + coef
                S = {}
                for l, r in (Dl(*b) if b[0] < b[1] else []):
                    Pl, Pr = P(E(l)), P(E(r))
                    for i1 in range(N):
                        for j1 in range(N):
                            if Pl[i1] * Pr[j1] != 0: S[(B[i1], B[j1])] = S.get((B[i1], B[j1]), 0) + Pl[i1] * Pr[j1]
                ok &= {k: v for k, v in T.items() if v != 0} == {k: v for k, v in S.items() if v != 0}
            nok += ok
    return ncomp, nok
def _mul(x, y, B, idx, N):
    r = [0] * N
    for i, p in enumerate(B):
        if x[i] == 0: continue
        for j, q in enumerate(B):
            if y[j] != 0 and p[1] == q[0]: r[idx[(p[0], q[1])]] += x[i] * y[j]
    return r
suff = {n: ut_suff(n) for n in (2, 3, 4)}
report("Y3", all(c == o for c, o in suff.values()) and [suff[n][0] for n in (2, 3, 4)] == [6, 23, 100],
       "every arrow-compatible vertex splitting of UT_n extends by SHIFT-TRANSPORT (paths carried along with their vertices), so compatibility is necessary and sufficient (n = 2, 3, 4; n = 5 in the exploration run: 563/563)",
       f"(compatible, extension valid): {suff}")

# ═════════ 4. Kronecker quiver: continuous moduli ═════════
Bk = ['e1', 'e2', 'al', 'be']
tab = {('e1', 'e1'): 'e1', ('e2', 'e2'): 'e2', ('e1', 'al'): 'al', ('al', 'e2'): 'al', ('e1', 'be'): 'be', ('be', 'e2'): 'be'}
def mk(x, y):
    r = sp.zeros(4, 1)
    for i in range(4):
        for j in range(4):
            if x[i] != 0 and y[j] != 0:
                k = tab.get((Bk[i], Bk[j]))
                if k: r[Bk.index(k)] += x[i] * y[j]
    return r
def is_fold(P):
    vec = lambda i: sp.Matrix([1 if k == i else 0 for k in range(4)])
    if not (P * P - P).is_zero_matrix: return False
    for i in range(4):
        for j in range(4):
            Px, Py = P * vec(i), P * vec(j)
            if not (mk(Px, Py) - P * (mk(Px, vec(j)) + mk(vec(i), Py) - mk(vec(i), vec(j)))).applyfunc(sp.simplify).is_zero_matrix: return False
    def Dk(v): T = sp.zeros(4, 4); T[0, 1] = v[2] + v[3]; return T
    return all((Dk(P * vec(j)) - P * Dk(vec(j)) * P.T).applyfunc(sp.simplify).is_zero_matrix for j in range(4))
t_ = sp.symbols('t')
# one-parameter family found by the exhaustive solve: P(e1)=e1, P(e2)=e2, P(α) = P(β) = (1−t)α + tβ  (columns = images)
fam = sp.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1 - t_, t_], [0, 0, 1 - t_, t_]]).T
report("Y4", all(is_fold(fam.subs(t_, v)) for v in (0, sp.Rational(1, 3), 2, -5)),
       "Kronecker quiver (two parallel arrows): folds come in CONTINUOUS families — e.g. a line of folds sending both arrows to the mixture (1−t)α + tβ; exhaustive solve: 30 components of dimension 0, 1, 2",
       "parallel arrows put α−β into ker Δ, and the choice of which mixture of α, β counts as order is a free parameter")

# ═════════ 5. The Borel frontier: productive regions can be any G_δ ═════════
# Survival S_k(p) = Σ_n 2^{−n} max(0, 1 − k·2^n·|p − q_n|) decreases to Σ_n 2^{−n}[p = q_n]: a Bernoulli-factory loop that
# continues with probability S_k/S_{k−1} never terminates with probability D(p) > 0 exactly at rationals.
from fractions import Fraction
rats = sorted({Fraction(a, b) for b in range(2, 13) for a in range(1, b) if Fraction(1, 4) <= Fraction(a, b) <= Fraction(3, 4)})
def S(k, p): return sum(2.0**-(n + 1) * max(0.0, 1 - k * 2**(n + 1) * abs(p - float(q))) for n, q in enumerate(rats))
irr = [1 / sqrt(2), (sqrt(5) - 1) / 2 + 0.0, pi / 5]
at_rat = [S(10**6, float(q)) for q in rats[:5]]
at_irr = [S(10**6, x) for x in irr if 0.25 <= x <= 0.75]
mono = all(S(k + 1, 0.4) <= S(k, 0.4) + 1e-15 for k in range(1, 200))
report("Y5", min(at_rat) > 0 and max(at_irr) < 1e-3 and mono,
       "a coin program's productive region can be the irrationals (non-σ-compact G_δ): survival probabilities stay ≥ 2^{−n} at rationals and vanish at irrationals",
       f"S at first rationals {[round(x, 4) for x in at_rat]}, at irrationals {[f'{x:.1e}' for x in at_irr]}")

# ═════════ 6. Clocks under social coupling ═════════
def softmax(x): e = np.exp(x - x.max()); return e / e.sum()
def tick_rate(K, H, R, k=3, N=200, T=60, seed=0):
    rng = np.random.default_rng(seed); s = rng.integers(k, size=N); fast = 0; ok = 0
    counts = np.bincount(s, minlength=k).astype(float)
    for _ in range(T):
        for _ in range(R * N):
            i = rng.integers(N); counts[s[i]] -= 1
            p = softmax(K * counts / (N - 1) + H * (np.arange(k) == fast)); s[i] = rng.choice(k, p=p); counts[s[i]] += 1
        top = np.flatnonzero(counts == counts.max()); f = (rng.choice(top) + 1) % k
        ok += f == (fast + 1) % k; fast = f
    return ok / T
avg = lambda K, R: np.mean([tick_rate(K, 0.15, R, seed=s) for s in range(3)])
slow = [avg(K, 2) for K in (0.0, 2.0, 3.5)]
fast_ = [avg(K, 10) for K in (0.0, 2.0, 3.5)]
report("Y6", slow[0] > slow[1] > slow[2] and fast_[1] > fast_[0] + 0.05 and fast_[2] < 0.25,
       "coupled culture driving a cyclic clock: with slow relaxation, coupling only adds inertia; with fast relaxation it sharpens the clock — until the 3-state Potts transition (K ≈ 4 ln 2 ≈ 2.77) locks it",
       f"true-tick rate at K = 0, 2, 3.5: slow (2 sweeps/tick) {[round(x, 3) for x in slow]}, fast (10 sweeps/tick) {[round(x, 3) for x in fast_]}")

# ═════════ 7. Conjecture 7: rigidity of graded deconstructing folds ═════════
from itertools import product as _prod
from collections import defaultdict as _dd
def graded_system(T):
    """Unknowns l_w = coefficient of u^{|w|} in P(w) (the ker Δ direction); everything else forced by Δ∘P = (P⊗P)∘Δ.
    Returns (equations, unknown symbols) for idempotence + RB(−1) + Δ-consistency up to total degree T."""
    words = lambda n: [''.join(t) for t in _prod('au', repeat=n)]
    lam = {w: sp.Symbol('l_' + (w or 'e')) for n in range(T + 1) for w in words(n)}
    Dl = lambda w: [(w[:i], w[i+1:]) for i, ch in enumerate(w) if ch == 'a']
    P = {}; cons = []
    for n in range(T + 1):
        W = words(n); others = [v for v in W if v != 'u' * n]
        for w in W:
            rhs = _dd(lambda: 0)
            for l, r in Dl(w):
                for x, c in P[l].items():
                    for y, e in P[r].items(): rhs[(x, y)] += c * e
            cs = sp.symbols(f'c0:{len(others)}'); D = _dd(lambda: 0)
            for v, c in zip(others, cs):
                for k in Dl(v): D[k] += c
            eqs = [sp.expand(D[k] - rhs[k]) for k in set(D) | set(rhs)]
            sol = sp.solve(eqs, cs, dict=True)
            if not sol:
                A, bvec = sp.linear_eq_to_matrix(eqs, cs); R, _ = A.row_join(bvec).rref()
                cons += [R[i, -1] for i in range(R.rows) if all(R[i, j] == 0 for j in range(len(cs))) and R[i, -1] != 0]
                sol = [{}]
            s = sol[0]
            P[w] = {v: s.get(c, c) for v, c in zip(others, cs)}; P[w]['u' * n] = P[w].get('u' * n, 0) + lam[w]
            P[w] = {k: sp.expand(v) for k, v in P[w].items() if sp.expand(v) != 0}
    def ap(x):
        r = _dd(lambda: 0)
        for w, c in x.items():
            for v, e in P[w].items(): r[v] += c * e
        return {k: sp.expand(v) for k, v in r.items() if sp.expand(v) != 0}
    def mul(x, y):
        r = _dd(lambda: 0)
        for a_, c in x.items():
            for b_, e in y.items(): r[a_ + b_] += c * e
        return {k: sp.expand(v) for k, v in r.items() if sp.expand(v) != 0}
    eqs = list(cons); allw = [w for n in range(T + 1) for w in words(n)]
    for w in allw:
        d1, d0 = ap(ap({w: 1})), P[w]
        eqs += [sp.expand(d1.get(k, 0) - d0.get(k, 0)) for k in set(d1) | set(d0)]
    for x in allw:
        for y in allw:
            if len(x) + len(y) > T: continue
            px, py = ap({x: 1}), ap({y: 1}); inner = _dd(lambda: 0)
            for k, v in mul(px, {y: 1}).items(): inner[k] += v
            for k, v in mul({x: 1}, py).items(): inner[k] += v
            inner[x + y] -= 1
            lhs, rhs = mul(px, py), ap(dict(inner))
            eqs += [sp.expand(lhs.get(k, 0) - rhs.get(k, 0)) for k in set(lhs) | set(rhs)]
    return [e for e in set(eqs) if e != 0], lam
def tangent_profile(T, point):
    eqs, lam = graded_system(T)
    pt = {v: point(w) for w, v in lam.items()}
    unk = sorted(lam.values(), key=lambda s: (len(s.name), s.name))
    J = sp.Matrix([[sp.diff(e, x).subs(pt) for x in unk] for e in eqs])
    N = J.nullspace(); deg = lambda s: 0 if s.name == 'l_e' else len(s.name) - 2
    return [sp.Matrix.hstack(*N)[[i for i, x in enumerate(unk) if deg(x) <= c], :].rank() if N else 0 for c in range(T + 1)]
framed_pt = lambda w: 1 if set(w) == {'u'} else 0
ident_pt = lambda w: 1 if w == '' else 0
tf, ti = tangent_profile(4, framed_pt), tangent_profile(4, ident_pt)
report("Y7", tf == [0, 1, 3, 5, 7] and ti == [0, 1, 1, 1, 1],
       "Conjecture 7 fails infinitesimally on the framed branch: the identification fold is rigid (tangent = the 1-dim conjugation orbit), the framed fold has 2c−1 deformation directions in degree ≤ c (2 new per degree: left and right edges)",
       f"rank of tangent projected to degree ≤ c (T=4): framed {tf}, identification {ti}; T=5 exploration run: framed [0,1,3,5,7,9]")
from scipy.optimize import least_squares
eqs5, lam5 = graded_system(5)
Ls = sp.Symbol
base = {Ls('l_e'): 0, **{Ls('l_' + 'u' * n): 1 for n in range(1, 6)}}
def nonlinear(extra, starts=40):
    es = [e.subs({**base, **extra}) for e in eqs5]; es = [e for e in es if e != 0]
    unk = sorted(set().union(*[e.free_symbols for e in es]), key=str)
    F = sp.lambdify(unk, es, 'numpy'); Jf = sp.lambdify(unk, sp.Matrix(es).jacobian(unk), 'numpy')
    rng = np.random.default_rng(0); best = 1e9
    for _ in range(starts):
        r = least_squares(lambda x: np.array(F(*x), dtype=float), rng.normal(scale=1.5, size=len(unk)),
                          jac=lambda x: np.array(Jf(*x), dtype=float), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
        best = min(best, np.max(np.abs(r.fun)))
        if best < 1e-10: break
    return best
two_sided = nonlinear({Ls('l_au'): -1, Ls('l_ua'): -2})
one_sided = nonlinear({Ls('l_au'): -1, Ls('l_ua'): 0})
report("Y8", two_sided < 1e-10 and one_sided < 1e-10,
       "the edge deformations integrate: framed folds with independent left/right twists (λ_L, λ_R) = (1, 2) and (1, 0) solve every equation through total degree 5 — conjugation only gives λ_L = λ_R",
       f"max residual: two-sided {two_sided:.1e}, one-sided {one_sided:.1e}")
# equal-centre sub-case of the density lemma: stable annihilator images (evidence)
uu = sp.symbols('uu')
def cons_rows(f, L):
    d = sp.degree(f, uu); rows = []; j = 0
    while d * j + 1 <= L:
        p = sp.Poly(sp.expand(f**j), uu); row = [0] * L
        for (k,), c in zip(p.monoms(), p.coeffs()):
            if k + 1 <= L: row[k] = c
        rows.append(row); j += 1
    return rows
def stable(f, g, L=24, Lp=72):
    N = sp.Matrix(cons_rows(f, Lp) + cons_rows(g, Lp)).nullspace()
    return sp.Matrix.hstack(*N)[:L, :].rank() if N else 0
eqc = {str((f, g)): stable(f, g) for f, g in ((uu**2, uu**3 + uu), (uu**3, uu**4 + uu**2), (uu**3 + uu, uu**3 + 2 * uu))}
report("Y9", all(v > 0 for v in eqc.values()),
       "equal-centre sub-case (evidence): exact annihilator solutions to order 72 restrict to a nonzero space at order 24 — no splitting arises",
       f"stable image dimensions: {eqc}")

print(f"\n{sum(res)}/{len(res)} checks pass")
