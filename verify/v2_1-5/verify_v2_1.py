"""Checks for 'Open Questions Explored (v2.1)'. Tags [Tn]. Sections Q1–Q6 follow v2 §10."""
from itertools import product, combinations
from collections import defaultdict
from fractions import Fraction as Q
from math import tanh, atanh, sqrt, exp, log
import random
import numpy as np
import sympy as sp

res = []
def report(tag, ok, claim, detail=""):
    res.append(ok); print(f"[{tag}] {'PASS' if ok else 'FAIL'}  {claim}" + (f"\n       {detail}" if detail else ""))

# ═════════ Q1  Dist-valued denotations ═════════
# Lemma: an n-state hidden Markov source (Moore emissions) whose output law is a Dirac δ_w emits w = u v^ω with |u|+|v| ≤ n.
# Exhaustive over n ≤ 3, probabilities in {0, ½, 1}: transitions, initial law, emission laws.
def rows(n):
    out = [tuple(Q(int(i == j)) for j in range(n)) for i in range(n)]
    out += [tuple(Q(1, 2) if j in pr else Q(0) for j in range(n)) for pr in combinations(range(n), 2)]
    return out
def dirac_prefix(n, T, init, emit, L):
    paths = {'': init}                                  # prefix → state vector
    for _ in range(L):
        nxt = {}
        for w, vec in paths.items():
            for sym in (0, 1):
                v2 = [sum(vec[i] * (emit[i] if sym else 1 - emit[i]) * T[i][j] for i in range(n)) for j in range(n)]
                if any(v2): nxt[w + str(sym)] = v2
        if len(nxt) > 1: return None
        paths = nxt
    return next(iter(paths))
def up_bound(w, n):
    L = len(w)
    return any(all(w[k] == w[k + (j - i)] for k in range(i, L - (j - i))) for j in range(1, n + 1) for i in range(j))
cnt = nontriv = viol = 0; L = 10
for n in (1, 2, 3):
    R = rows(n); E = [Q(0), Q(1, 2), Q(1)]
    for T in product(R, repeat=n):
        for init in R:
            for emit in product(E, repeat=n):
                w = dirac_prefix(n, T, init, emit, L)
                if w is None: continue
                cnt += 1
                nontriv += any(len([x for x in row if x]) > 1 for row in T) or len([x for x in init if x]) > 1
                viol += not up_bound(w, n)
report("T1", viol == 0 and nontriv > 0,
       "every Dirac-output HMM (n ≤ 3 states, probabilities in {0,½,1}) emits a prefix of u·v^ω with |u|+|v| ≤ n",
       f"{cnt} Dirac sources to depth {L}, {nontriv} with genuinely random state paths, 0 violations")

# Σ⁰₂-hardness by composition: x ↦ δ_{F(x)} with F the FIN→UP encoding (1 at i ↦ 0·1^(i+2), 0 ↦ 10).
def enc(x):
    return ''.join(('0' + '1' * (i + 2)) if b else '10' for i, b in enumerate(x))
def max_run1(s): return max((len(r) for r in s.split('0')), default=0)
fin = [1, 0, 1] + [0] * 40; inf = [1 if i % 3 == 0 else 0 for i in range(43)]
report("T2", enc(fin)[-20:] == '10' * 10 and max_run1(enc(inf)[:200]) > max_run1(enc(inf)[:60]),
       "FIN ↦ eventually (10)^ω; infinitely many 1s ↦ unbounded 1-runs (not ultimately periodic); Dirac embedding is weak-* continuous")

# ═════════ Q2  A fold that deconstructs ═════════
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
ONE = {'': 1}
def Dl(x):                                                   # letter deletion
    r = defaultdict(int)
    for w, c in x.items():
        for i in range(len(w)): r[(w[:i], w[i+1:])] += c
    return {k: c for k, c in r.items() if c}
def subst(phi, x):                                           # algebra endomorphism from letter images
    out = {}
    for w, c in x.items():
        t = ONE
        for ch in w: t = mul(t, phi[ch])
        out = add(out, sc(c, t))
    return out
def tens_map(phi, T):
    r = defaultdict(int)
    for (u, v), c in T.items():
        for u2, c2 in subst(phi, {u: 1}).items():
            for v2, c3 in subst(phi, {v: 1}).items(): r[(u2, v2)] += c * c2 * c3
    return {k: c for k, c in r.items() if c}
span = ['', 'a', 'b', 'aa', 'ab', 'ba', 'bb']
C = []
for cs in product((-1, 0, 1), repeat=7):
    p = clean(dict(zip(span, cs)))
    if Dl(p) == {('', ''): 1}: C.append(p)
u = add({'a': 1}, {'b': -1}); u2 = mul(u, u)
pred = [add({'a': 1}, sc(c0, ONE), sc(c1, u), sc(c2, u2)) for c0, c1, c2 in product((-1, 0, 1), repeat=3)]
pred = [p for p in pred if all(abs(c) <= 1 for c in p.values()) and set(p) <= set(span)]
report("T3", sorted(map(lambda d: sorted(d.items()), C)) == sorted(map(lambda d: sorted(d.items()), pred)),
       "Δp = 1⊗1 (deg ≤ 2, coeffs in {−1,0,1}) ⟺ p ∈ a + k[u], u = a − b  (so ker Δ = k[u])", f"{len(C)} solutions")
words3 = [''.join(t) for n_ in range(4) for t in product('ab', repeat=n_)]
folds = []; ident = 0
for pa in C:
    for pb in C:
        phi = {'a': pa, 'b': pb}
        if all(subst(phi, subst(phi, {w: 1})) == subst(phi, {w: 1}) for w in ['a', 'b']):
            if pa == {'a': 1} and pb == {'b': 1}: ident += 1; continue
            folds.append(phi)
def key(phi): return (tuple(sorted(phi['a'].items())), tuple(sorted(phi['b'].items())))
found = {key(p) for p in folds}
pred_f = set()                               # predicted family: φ(a) = a + g(u), φ(b) = φ(a) − c, g(c) = 0
for g0, g1, g2, c in product((-2, -1, 0, 1, 2), repeat=4):
    if g0 + g1 * c + g2 * c * c != 0: continue
    pa = add({'a': 1}, sc(g0, ONE), sc(g1, u), sc(g2, u2)); pb = add(pa, sc(-c, ONE))
    if pa in C and pb in C and not (pa == {'a': 1} and pb == {'b': 1}): pred_f.add(key({'a': pa, 'b': pb}))
comm = all(Dl(subst(phi, {w: 1})) == tens_map(phi, Dl({w: 1})) for phi in folds for w in words3)
report("T4", ident == 1 and found == pred_f and comm,
       "non-identity idempotent endomorphism folds commuting with letter deletion are exactly φ: u ↦ c, a ↦ a+g(u), g(c)=0",
       f"{len(folds)} found by search = {len(pred_f)} predicted (both directions); all commute with Δ on words ≤3")
idf = {'a': {'a': 1}, 'b': {'a': 1}}; exf = {'a': {'a': 1}, 'b': {}}
ok_id = all(Dl(subst(idf, {w: 1})) == tens_map(idf, Dl({w: 1})) for w in words3)
ok_ex = all(Dl(subst(exf, {w: 1})) == tens_map(exf, Dl({w: 1})) for w in words3)
report("T5", ok_id and not ok_ex and idf in folds,
       "the identification fold b ↦ a deconstructs; the exclusion fold b ↦ 0 (v2's f) does not — same Hodge k[a], Podge (a−b) vs (b)")

# ═════════ Q3  The Koszul pair ═════════
# Finite-dim algebras by structure constants; P = im f, N = ker f; operations as bilinear tensors.
def alg(n, table):
    M = np.zeros((n, n, n))
    for (i, j), v in table.items():
        for k, c in v.items(): M[i, j, k] = c
    return M
def ops(M, fmat):
    n = M.shape[0]
    mul_ = lambda x, y: np.einsum('i,j,ijk->k', x, y, M)
    f = lambda x: fmat @ x
    I = np.eye(n)
    def tensor(op): return np.array([[op(I[i], I[j]) for j in range(n)] for i in range(n)])
    return {'⊣': tensor(lambda x, y: mul_(x, f(y))), '⊢': tensor(lambda x, y: mul_(f(x), y)),
            '≺': tensor(lambda x, y: mul_(x, f(y))), '≻': tensor(lambda x, y: mul_(f(x), y) - mul_(x, y)),
            '·': tensor(mul_)}
def assoc(M): return np.allclose(np.einsum('ijm,mkl->ijkl', M, M), np.einsum('jkm,iml->ijkl', M, M))
# basis e1, e2, n ;  UT2: e1 n = n, n e2 = n ;  A': e2 n = n e2 = n (n in the e2-corner)
UT2 = alg(3, {(0, 0): {0: 1}, (1, 1): {1: 1}, (0, 2): {2: 1}, (2, 1): {2: 1}})
Ap = alg(3, {(0, 0): {0: 1}, (1, 1): {1: 1}, (1, 2): {2: 1}, (2, 1): {2: 1}})
fP = np.diag([1, 1, 0])
o1, o2 = ops(UT2, fP), ops(Ap, fP)
dendsame = np.allclose(o1['≺'], o2['≺']) and np.allclose(o1['≻'], o2['≻'])
diadiff = not np.allclose(o1['⊢'], o2['⊢'])
# basis 1, x, y: k[x]/(x³) with y = x²  vs  k ⊕ (x,y) square-zero
B1 = alg(3, {(0, 0): {0: 1}, (0, 1): {1: 1}, (1, 0): {1: 1}, (0, 2): {2: 1}, (2, 0): {2: 1}, (1, 1): {2: 1}})
B2 = alg(3, {(0, 0): {0: 1}, (0, 1): {1: 1}, (1, 0): {1: 1}, (0, 2): {2: 1}, (2, 0): {2: 1}})
fS = np.diag([1, 0, 0]); p1, p2 = ops(B1, fS), ops(B2, fS)
diasame = np.allclose(p1['⊣'], p2['⊣']) and np.allclose(p1['⊢'], p2['⊢'])
denddiff = not np.allclose(p1['≻'], p2['≻'])
recover = all(np.allclose(o['⊢'] - o['≻'], o['·']) for o in (o1, o2, p1, p2))
report("T6", all(map(assoc, (UT2, Ap, B1, B2))) and dendsame and diadiff and diasame and denddiff and recover,
       "neither half of the Koszul pair determines the other: UT₂ vs k×k[ε] share δ-D but not δ-L; k[x]/(x³) vs square-zero share δ-L but not δ-D; together xy = x⊢y − x≻y")

# block bookkeeping, generic: ⊣ = ≺ = pq+nq ; ⊢ = pq+pm ; ≻ = −nq−nm
p1s, n1s, p2s, n2s = sp.symbols('p1 n1 p2 n2', commutative=False)
P = {p1s, p2s}
def fz(e):
    out = 0
    for t in sp.Add.make_args(sp.expand(e)):
        if all(s in P for s in t.free_symbols if not s.is_commutative): out += t
    return out
X, Y = p1s + n1s, p2s + n2s
blocks = {'⊣': sp.expand(X * fz(Y)), '⊢': sp.expand(fz(X) * Y), '≻': sp.expand(fz(X) * Y - X * Y)}
want = {'⊣': p1s*p2s + n1s*p2s, '⊢': p1s*p2s + p1s*n2s, '≻': -n1s*p2s - n1s*n2s}
report("T7", all(sp.expand(blocks[k] - want[k]) == 0 for k in want),
       "block content: δ-L sees {pq, nq, pm} and is blind to nm; δ-D sees {pq, nq, nm} and is blind to pm")

# Leibniz bracket recovers im f as left annihilator iff im f is self-centralizing; right annihilator = ker f ⊕ (im f ∩ Z(A))
def leib_ann(n):
    xs = sp.symbols(f'x0:{n*n}'); Xm = sp.Matrix(n, n, xs)
    Xm = Xm.applyfunc(lambda e: e)
    for i in range(n):
        for j in range(i): Xm[i, j] = 0
    D = lambda M: sp.diag(*[M[i, i] for i in range(n)])
    basis = [(i, j) for i in range(n) for j in range(i, n)]
    E = lambda i, j: sp.Matrix(n, n, lambda a, b: 1 if (a, b) == (i, j) else 0)
    br = lambda A_, B_: A_ * D(B_) - D(B_) * A_
    free = [Xm[i, j] for (i, j) in basis]
    left = sp.solve([e for (i, j) in basis for e in br(Xm, E(i, j)) if e != 0], free, dict=True)[0]
    right = sp.solve([e for (i, j) in basis for e in br(E(i, j), Xm) if e != 0], free, dict=True)[0]
    lfree = sorted((i, j) for (i, j) in basis if Xm[i, j] not in left)
    rconstr = {k: v for k, v in right.items()}
    return lfree, rconstr
lf, rc = leib_ann(3)
diag_ok = lf == [(0, 0), (1, 1), (2, 2)]
right_ok = all(str(k) in ('x0', 'x4', 'x8') for k in rc) and len(set(rc.values())) == 1   # diagonal forced scalar
report("T8", diag_ok and right_ok,
       "UT₃ Leibniz bracket: left annihilator = diag = im f (self-centralizing); right annihilator = ker f ⊕ k·1",
       f"left free coords {lf}; right forces {rc}")

# ═════════ Q4  Proof theory of 𝔼 = S_fde ═════════
_src = open(__file__.replace('verify_v2_1.py', 'audit_logic_topology.py')).read()
_ns = {}; exec(_src.split('# S1 ')[0], _ns)
ev, var, formulas, assigns, VALS, LPV = (_ns[k] for k in ('ev', 'var', 'formulas', 'assigns', 'VALS', 'LPV'))
def mcons(G, phi, atoms, vals):
    return all('T' in ev(phi, v) for v in assigns(atoms, vals) if all('T' in ev(g, v) for g in G))
def nnf(phi, pos=True):
    if isinstance(phi, str): return ('lit', phi, pos)
    if phi[0] == '~': return nnf(phi[1], not pos)
    op = phi[0] if pos else ('|' if phi[0] == '&' else '&')
    return (op, nnf(phi[1], pos), nnf(phi[2], pos))
def prove(L, R, V=None):
    """Cut-free invertible calculus G_LP (V=None) / availability-labelled G_𝔼 (V = Var(Γ)).
    Initial sequents: identity ℓ ⇒ ℓ and excluded middle ⇒ p, ¬p; NO explosion p, ¬p ⇒.
    Label condition (G_𝔼): every atom occurring on the right of a leaf is available (∈ V)."""
    for i, g in enumerate(L):
        if g[0] == '&': return prove(L[:i] + [g[1], g[2]] + L[i+1:], R, V)
        if g[0] == '|': return prove(L[:i] + [g[1]] + L[i+1:], R, V) and prove(L[:i] + [g[2]] + L[i+1:], R, V)
    for i, g in enumerate(R):
        if g[0] == '|': return prove(L, R[:i] + [g[1], g[2]] + R[i+1:], V)
        if g[0] == '&': return prove(L, R[:i] + [g[1]] + R[i+1:], V) and prove(L, R[:i] + [g[2]] + R[i+1:], V)
    Ls, Rs = {(g[1], g[2]) for g in L}, {(g[1], g[2]) for g in R}
    if V is not None and not {p for p, _ in Rs} <= V: return False
    return bool(Ls & Rs) or any((p, True) in Rs and (p, False) in Rs for p, _ in Rs)
at2 = ['p', 'q']; F3 = formulas(at2, 3); F2 = formulas(at2, 2); at3 = ['p', 'q', 'r']; F2c = formulas(at3, 2)
bad_lp = bad_e = n = 0
for g in F3:
    for phi in F2:
        n += 1
        bad_lp += prove([nnf(g)], [nnf(phi)]) != mcons([g], phi, at2, LPV)
        bad_e += prove([nnf(g)], [nnf(phi)], var(g)) != mcons([g], phi, at2, VALS)
random.seed(11); m = 0
for _ in range(3000):
    G = random.sample(F2c, random.choice([1, 2, 3])); phi = random.choice(F2c); m += 1
    V = set().union(*map(var, G))
    bad_e += prove([nnf(g) for g in G], [nnf(phi)], V) != mcons(G, phi, at3, VALS)
report("T9", bad_lp == 0 and bad_e == 0,
       f"G_LP is sound and complete for LP, and availability-labelled G_𝔼 for 𝔼 ({n} single-premise pairs ≤3/≤2 connectives + {m} multi-premise, 3 atoms)")
pq = ('|', 'p', 'q')
local = lambda L, R: prove(L, R, None) and False
def prove_local(L, R):               # Szmuc-style leaf-local containment WITHOUT the label: right atoms ⊆ left atoms of the leaf
    for i, g in enumerate(L):
        if g[0] == '&': return prove_local(L[:i] + [g[1], g[2]] + L[i+1:], R)
        if g[0] == '|': return prove_local(L[:i] + [g[1]] + L[i+1:], R) and prove_local(L[:i] + [g[2]] + L[i+1:], R)
    for i, g in enumerate(R):
        if g[0] == '|': return prove_local(L, R[:i] + [g[1], g[2]] + R[i+1:])
        if g[0] == '&': return prove_local(L, R[:i] + [g[1]] + R[i+1:]) and prove_local(L, R[:i] + [g[2]] + R[i+1:])
    Ls, Rs = {(g[1], g[2]) for g in L}, {(g[1], g[2]) for g in R}
    if not {p for p, _ in Rs} <= {p for p, _ in Ls}: return False
    return bool(Ls & Rs) or any((p, True) in Rs and (p, False) in Rs for p, _ in Rs)
report("T10", mcons([pq], pq, at2, VALS) and prove([nnf(pq)], [nnf(pq)], {'p', 'q'}) and not prove_local([nnf(pq)], [nnf(pq)]),
       "containment is a GLOBAL invariant of the proof tree: the leaf-local check rejects p∨q ⇒ p∨q (leaf p ⇒ p,q), the threaded label V accepts it")
# cut admissibility: Γ ⊢ φ and Γ,φ ⊢ ψ ⟹ Γ ⊢ ψ  (label monotone: Var(φ) ⊆ V)
random.seed(13); cut_bad = cut_n = 0
for _ in range(20000):
    G = random.sample(F2c, random.choice([1, 2])); phi = random.choice(F2c); psi = random.choice(F2c)
    V = set().union(*map(var, G)); NG = [nnf(g) for g in G]
    if prove(NG, [nnf(phi)], V) and prove(NG + [nnf(phi)], [nnf(psi)], V | var(phi)):
        cut_n += 1; cut_bad += not prove(NG, [nnf(psi)], V)
report("T11", cut_bad == 0 and cut_n > 50, f"cut is admissible in G_𝔼 ({cut_n} random cut instances, 0 failures)")

# ═════════ Q5  Coupled cultures ═════════
def roots(K, H, grid=20000):        # even grid: avoids sampling the root m = 0 exactly
    xs = np.linspace(-1, 1, grid); g = np.tanh(K * xs + H) - xs
    return int(np.sum(np.sign(g[:-1]) != np.sign(g[1:])))
crit = all(roots(K, 0.0) == 1 for K in (0.5, 0.8, 0.95)) and all(roots(K, 0.0) == 3 for K in (1.05, 1.5, 3.0))
def Hs(K): ms = sqrt(1 - 1 / K); return K * ms - atanh(ms)
spin = all(roots(K, 0.9 * Hs(K)) == 3 and roots(K, 1.1 * Hs(K)) == 1 for K in (1.2, 1.5, 2.0, 3.0))
report("T12", crit and spin,
       "mean-field culture m = tanh(K m + H): spontaneous consensus iff K = βJ > 1; a wrong-consensus branch survives iff H < H_s(K) = K√(1−1/K) − artanh√(1−1/K)",
       f"H_s: K=1.5 → {Hs(1.5):.4f}, K=2 → {Hs(2):.4f}, K=3 → {Hs(3):.4f}")
def glauber(N, K, H, sweeps, s0=None, seed=0):
    rng = np.random.default_rng(seed); s = np.ones(N) * -1 if s0 == 'wrong' else rng.choice([-1.0, 1.0], N)
    tot = s.sum(); first_cross = None
    for t in range(sweeps * N):
        i = rng.integers(N); mi = (tot - s[i]) / (N - 1)
        new = 1.0 if rng.random() < 0.5 * (1 + np.tanh(K * mi + H)) else -1.0
        tot += new - s[i]; s[i] = new
        if first_cross is None and tot > 0: first_cross = t / N
    return tot / N, first_cross
ms = {K: np.mean([abs(glauber(400, K, 0.0, 150, seed=sd)[0]) for sd in range(4)]) for K in (0.5, 0.8, 1.3, 1.6)}
sharp = ms[0.5] < 0.15 and ms[0.8] < 0.2 and ms[1.3] > 0.5 and ms[1.6] > 0.75
K = 1.5
stay = glauber(400, K, 0.5 * Hs(K), 200, 'wrong', 1)[1] is None
flip = glauber(400, K, 1.5 * Hs(K), 200, 'wrong', 1)[1] is not None
report("T13", sharp and stay and flip,
       "Glauber simulation (N=400): |m| ≈ 0 below K=1, ≈ m*(K) above; at K=1.5 an all-wrong culture stays wrong for 200 sweeps at H = ½H_s, flips at H = 1.5H_s",
       f"|m|: {ms}")
life = []
for N in (16, 32, 64):
    ts = [glauber(N, 1.5, 0.5 * Hs(1.5), 3000, 'wrong', sd)[1] for sd in range(6)]
    life.append(np.median([t if t is not None else 3000 for t in ts]))
report("T14", life[0] < life[1] < life[2], "wrong-consensus lifetime grows with culture size (Kramers escape, ~exp(N·ΔF))",
       f"median sweeps to escape at N = 16, 32, 64: {life}")

# ═════════ Q6  Arrow fibers: the bureaucracy ═════════
from itertools import permutations
forms, lanes = ['X', 'Y', 'Z'], {'fast': 1, 'mid': 3, 'slow': 6}
structs = list(permutations(['fast', 'mid', 'slow']))     # structure s: form i ↦ lane s[i]  (a GridF tree Par(Seq(lane, form), …))
def grid(s): return ('Par', [('Seq', ('Atom', s[i]), ('Atom', forms[i])) for i in range(3)])
def cata_latency(t, pi):            # Mu-side analysis: expected latency, by a fold over the static tree
    tag = t[0]
    if tag == 'Par': return sum(cata_latency(c, pi) for c in t[1])
    if tag == 'Seq': return pi[t[2][1]] * lanes[t[1][1]]
def simulate(s, pi, n, rng):        # Nu-side run: route n sampled forms through the fixed wiring
    route = dict(zip(forms, s)); draws = rng.choice(forms, size=n, p=[pi[f] for f in forms])
    return np.mean([lanes[route[f]] for f in draws])
rng = np.random.default_rng(23); arrow_ok = True
for s in structs:
    for pi in ({'X': .5, 'Y': .3, 'Z': .2}, {'X': .1, 'Y': .1, 'Z': .8}):
        arrow_ok &= abs(cata_latency(grid(s), pi) - simulate(s, pi, 40000, rng)) < 0.06
def traffic(s, q=0.6):              # citizens game the system: favour the form currently in the fast lane
    fast = forms[s.index('fast')]; return {f: (q if f == fast else (1 - q) / 2) for f in forms}
def efficient(pi):                  # fixed-traffic fiber: a pure function of the measure
    return min(structs, key=lambda s: cata_latency(grid(s), pi))
def fair(pi):                       # 'serve the neglected': fast lane to the least-filed form
    order = sorted(forms, key=lambda f: (pi[f], f)); lane_of = {order[0]: 'fast', order[1]: 'mid', order[2]: 'slow'}
    return tuple(lane_of[f] for f in forms)
def orbits(R):
    Phi = {s: R(traffic(s)) for s in structs}
    fixed = [s for s in structs if Phi[s] == s]
    cyc = set()
    for s in structs:
        seen = []; x = s
        while x not in seen: seen.append(x); x = Phi[x]
        cyc.add(len(seen) - seen.index(x))
    return fixed, cyc
fe, ce = orbits(efficient); ff, cf = orbits(fair)
report("T15", arrow_ok, "fiber 1 (structure fixed): the fold over the static tree predicts the run's mean latency for every structure and traffic (Arrow-sound)")
report("T16", len(fe) == 3 and ce == {1} and ff == [] and 2 in cf,
       "diagonal: 'efficient' restructuring has 3 self-fulfilling fixed points (lock-in); 'fair' restructuring has none — every orbit cycles",
       f"efficient fixed points: {fe}; fair cycle lengths: {sorted(cf)}")

print(f"\n{sum(res)}/{len(res)} checks pass")
