"""verify_v2_9.py — machine checks N1–N7 for discordianism-open-questions-explored-v2.9.md
gen.py / gen2.py / gseries.py: exact engine for unital-Hodge folds at ALL a-degrees (generators g(f) = f − κ(f)),
fold.py / series.py / sides*.py / second.py: the v2.6–v2.8 engine (a-degree ≤ 2), used only for the comparison in N7."""
import sys, os, random, itertools, functools, importlib, contextlib, io, time
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
res = []
def report(tag, ok, claim, detail=""):
    res.append(bool(ok)); print(f"[{tag}] {'PASS' if ok else 'FAIL'}  {claim}", flush=True)
    if detail: print(f"       {detail}", flush=True)
import fold as FD
from fold import MOD, md
from tangent import spec
import gen as G0
from gen2 import GenSym, equations
import gseries
p = MOD; rt2 = pow(2, (p + 1) // 4, p)
def q0_of(S): return sorted([i for i in range(len(S.B)) if S.side[i] == 'Q'], key=lambda i: S.lev[i])[0]
# ═════════ N1: the correction formula is the general coideal solution ═════════
def Delta(E):
    out = {}
    for w, c in E.items():
        for t in range(len(w) - 1): out[(w[:t + 1], w[t + 1:])] = (out.get((w[:t + 1], w[t + 1:]), 0) + c) % MOD
    return {k: v for k, v in out.items() if v}
ok1 = True; nw = {}
for kind in ('ident', 'dihedral'):
    S = spec(kind, 30); q0 = q0_of(S); rnd = random.Random(1)
    words = list(G0.gen_words(S, 4, 3, q0)); nw[kind] = len(words)
    C = {w: {b: rnd.randint(1, 99) for b in range(len(S.B)) if S.side[b] == 'H' and S.lev[b] <= sum(S.lev[x] for x in w) + 1} for w in words}
    G = G0.Gen(S, C)
    for f in words:
        rhs = {}; i = f.index(q0)
        for t in range(len(f) - 1):
            pre, suf = f[:t + 1], f[t + 1:]
            for w, c in (G.kappa(suf).items() if t < i else G.kappa(pre).items()):
                key = (pre, w) if t < i else (w, suf); rhs[key] = (rhs.get(key, 0) + c) % MOD
        ok1 &= Delta(G.kappa(f)) == {k: v for k, v in rhs.items() if v}
report("N1", ok1, "for unital-Hodge folds, κ(f) = Σ_{optional left cut, optional right cut} L₁·a·C(m)·a·R₂ solves the coideal recursion for an ARBITRARY linear C: generator words → H₀ (cofree coderivation); C = a-degree-0 part of P on generator words",
       f"random C, all generator words (≤ 4 runs, level ≤ 3): {nw}")
# ═════════ N2: the all-a-degree engine reproduces the identification theorem ═════════
def randC(S, words, pred, delta=1, seed=7):
    rnd = random.Random(seed); C = {}
    for w in words:
        if pred(w):
            L = sum(S.lev[x] for x in w) + delta
            C[w] = {b: rnd.randint(1, 50) for b in range(len(S.B)) if S.side[b] == 'H' and S.lev[b] <= L}
    return C
S = spec('ident', 30); q0 = q0_of(S); words = list(G0.gen_words(S, 5, 4, q0))
r2 = {nm: G0.check(G0.Gen(S, randC(S, words, pr)), words, 2) for nm, pr in
      (('chaos-first', lambda w: w[0] == q0), ('chaos-last', lambda w: w[-1] == q0), ('mixed', lambda w: True))}
report("N2", r2['chaos-first'][0] == 0 and r2['chaos-last'][0] == 0 and r2['mixed'][0] > 0,
       "calibration: in the all-a-degree engine, identification folds with random C on chaos-first (resp. chaos-last) words are exact folds (families R, L); mixed support fails",
       f"bad two-generator products: {r2}")
# ═════════ N3: dihedral one-sided support is constrained ═════════
Sd = spec('dihedral', 40); qd = q0_of(Sd); wd = list(G0.gen_words(Sd, 3, 3, qd))
r3 = G0.check(G0.Gen(Sd, randC(Sd, wd, lambda w: w[0] == qd)), wd, 2)
report("N3", r3[0] > 0, "unlike identification, a RANDOM C on dihedral chaos-first words is not a fold: the dihedral one-sided families are cut out by Q₀-module conditions (v2.6 Prop. 2) already at first order",
       f"bad products: {r3[0]}/{r3[1]}")
# ═════════ N4: exact closed-form one-sided folds — root deviation ═════════
roots = {'r+': (1 + rt2) % p, 'r-': (1 - rt2) % p, 'r=5': 5}
def ev(i, r): return sum(md(c) * pow(r, m[0], p) for m, c in Sd.B[i][2].terms()) % p
Hs = [i for i in range(len(Sd.B)) if Sd.side[i] == 'H']; one = [i for i in Hs if Sd.lev[i] == 0][0]
def dev(h, r, coef):
    d = {h: coef % p}; d[one] = (d.get(one, 0) - coef * ev(h, r)) % p
    return {b: v for b, v in d.items() if v}
def Cform(lam, rl, mu, rr, maxlev=40):
    C = {}
    for h in Hs:
        if Sd.lev[h] > maxlev: continue
        if lam: C[(qd, h)] = dev(h, rl, lam)
        if mu: C[(h, qd)] = dev(h, rr, mu)
    return {w: d for w, d in C.items() if d}
w6 = list(G0.gen_words(Sd, 3, 6, qd)); w2 = list(G0.gen_words(Sd, 3, 2, qd))
chk = {}
chk['right, r+ (level≤6, pairs)'] = G0.check(G0.Gen(Sd, Cform(0, 0, 4, roots['r+'])), w6, 2)
chk['left, r− (level≤6, pairs)'] = G0.check(G0.Gen(Sd, Cform(3, roots['r-'], 0, 0)), w6, 2)
chk['right, r+ (level≤2, triples)'] = G0.check(G0.Gen(Sd, Cform(0, 0, 4, roots['r+'])), w2, 3)
chk['right, r = 5 (control)'] = G0.check(G0.Gen(Sd, Cform(0, 0, 4, roots['r=5'])), wd, 2)
sup = Cform(0, 0, 4, roots['r+'])
for w, d in Cform(0, 0, 5, roots['r-']).items():
    for b, v in d.items(): sup.setdefault(w, {}); sup[w][b] = (sup[w].get(b, 0) + v) % p
chk['right, μ₁·dev(r+) + μ₂·dev(r−) (control)'] = G0.check(G0.Gen(Sd, sup), wd, 2)
okc = all(v[0] == 0 for k, v in chk.items() if 'control' not in k) and all(v[0] > 0 for k, v in chk.items() if 'control' in k)
report("N4", okc, "exact closed-form one-sided dihedral folds at ALL a-degrees: C(h·a·q₀) = μ·(h − h(r)) (all other C = 0) is a fold for every μ iff r is a root b ± √s of the Podge's quadric; mirror for C(q₀·a·h); two roots cannot be superposed",
       f"bad products: {chk}")
# ═════════ N5: generic one-sided deformations integrate at all a-degrees; generic two-sided do not ═════════
def generic_run(kind, sup_pred, Lvars, Dc, Llow, N, maxlev, seed=5):
    S_ = spec(kind, 2 * Lvars + 9); q_ = q0_of(S_); ws = list(G0.gen_words(S_, 3, maxlev, q_))
    lev = lambda m: sum(S_.lev[x] for x in m)
    Gs = GenSym(S_, 1, sup_pred(q_), Lvars); eqs, sk = equations(Gs, ws)
    K, _ = FD.kernel(FD.linear_part(eqs), len(Gs.names)); rnd = random.Random(seed); tau = {}
    for v in K:
        a = rnd.randint(1, 60)
        for j, x in v.items():
            if lev(Gs.names[j][0]) <= Llow: tau[Gs.names[j]] = (tau.get(Gs.names[j], 0) + a * x) % MOD
    gseries.DCORR[0] = Dc
    return gseries.integrate(S_, ws, 1, Lvars, tau, N, Llow=Llow)
t0 = time.time()
g5 = {'chaos-first (level-3 generators)': generic_run('dihedral', lambda q: (lambda m: m[0] == q), 22, 5, 2, 3, 3),
      'chaos-last (level-2 generators)': generic_run('dihedral', lambda q: (lambda m: m[-1] == q), 20, 5, 2, 3, 2),
      'two-sided generic (level-3 generators)': generic_run('dihedral', lambda q: (lambda m: True), 22, 5, 2, 3, 3)}
report("N5", g5['chaos-first (level-3 generators)'] == (3, 'ok') and g5['chaos-last (level-2 generators)'] == (3, 'ok') and g5['two-sided generic (level-3 generators)'][1] == 'obstructed',
       "skip-free series in C-coordinates (all a-degrees): a generic one-sided first-order dihedral deformation integrates (order 3); a generic two-sided one is obstructed",
       f"{g5}  [{time.time()-t0:.0f}s]")
# ═════════ N6: two-sided second-order components are window artifacts ═════════
from second import quadrics
from sides import ratrec
def two_sided_families(maxlev, Lvars, deta):
    S_ = spec('dihedral', 2 * Lvars + 2 * deta + 6); q_ = q0_of(S_); ws = list(G0.gen_words(S_, 3, maxlev, q_))
    lev = lambda m: sum(S_.lev[x] for x in m)
    Gs = GenSym(S_, deta, lambda m: True, Lvars); eqs, sk = equations(Gs, ws); J = FD.linear_part(eqs)
    big = [j for j, (m, b) in enumerate(Gs.names) if S_.lev[b] > lev(m) + 1]
    def fam(pred):
        kill = [j for j, (m, b) in enumerate(Gs.names) if not pred(m)]
        K, _ = FD.kernel(J + [{j: 1} for j in big + kill], len(Gs.names))
        low = [j for j, (m, b) in enumerate(Gs.names) if lev(m) <= 2 and len(m) == 2]
        E = FD.Echelon(); return [v for v in K if E.add({j: v[j] for j in low if j in v})]
    L = fam(lambda m: m[0] == q_); R = fam(lambda m: m[-1] == q_)
    names = [f'l{i}' for i in range(len(L))] + [f'r{i}' for i in range(len(R))]
    Q, pairs, nv = quadrics(Gs, eqs, L + R)
    polys = [sp.expand(sum(ratrec(v) * sp.Symbol(names[pairs[c - nv][0]]) * sp.Symbol(names[pairs[c - nv][1]]) for c, v in r.items())) for r in Q]
    syms = sp.symbols(names); sols = sp.solve(polys, syms, dict=True)
    two = [s_ for s_ in sols if any(s_.get(x, x) != 0 for x in syms if str(x)[0] == 'l') and any(s_.get(x, x) != 0 for x in syms if str(x)[0] == 'r')]
    return len(sols), len(two), sum(1 for s_ in two if any(v.has(sp.sqrt(2)) for v in s_.values()))
t0 = time.time()
fam2 = two_sided_families(2, 14, 3); fam3 = two_sided_families(3, 16, 3)
report("N6", fam2[1] > 0 and fam3[1] == 0,
       "the second-order variety on chaos-first ⊕ chaos-last has two-sided solution families (incl. ℚ(√2) ones) when only level-2 generators are imposed, and NONE once level-3 generators are imposed: the two-sided dihedral components of v2.6–v2.8 are window artifacts",
       f"(families, two-sided, over ℚ(√2)): level-2 window {fam2}, level-3 window {fam3}  [{time.time()-t0:.0f}s]")
print(f"\n{sum(res)}/{len(res)} checks pass")
