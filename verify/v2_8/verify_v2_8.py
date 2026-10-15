"""verify_v2_8.py — machine checks L1–L7 for discordianism-open-questions-explored-v2.8.md"""
import sys, os, importlib, random, contextlib, io
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
res = []
def report(tag, ok, claim, detail=""):
    res.append(bool(ok)); print(f"[{tag}] {'PASS' if ok else 'FAIL'}  {claim}")
    if detail: print(f"       {detail}")
# ═════════ 1. Identification folds = R ∪ L: proof steps ═════════
with contextlib.redirect_stdout(io.StringIO()):
    from feq import Phi, E, mul, add, Cf
s, sig, tau, t = sp.symbols('s sigma tau t')
def cores(lhs, rhs):
    d = add((1, lhs), (-1, rhs)); out = {}
    for k, v in d.items():
        n, _ = sp.fraction(sp.together(v)); n = sp.expand(n)
        out[k] = sp.expand(sp.Mul(*[f**m for f, m in sp.factor_list(n)[1] if f.has(Cf)])) if n != 0 else 0
    return out
D1 = cores(Phi(E(s), E(sig), E(tau), E(t)), Phi(E(s), {0: 1}, mul(E(sig), E(tau)), E(t)))
D2 = cores(Phi(E(s), E(sig), E(tau), E(t)), Phi(E(s), mul(E(sig), E(tau)), {0: 1}, E(t)))
nz1 = [v for v in D1.values() if v != 0]; nz2 = [v for v in D2.values() if v != 0]
single = len(nz1) == 2 and sp.expand(nz1[0] - nz1[1]) == 0 or sp.expand(nz1[0] + nz1[1]) == 0
report("L1", single and D1.get(sig, 0) == 0 and D1.get(tau, 0) == 0,
       "two-u consistency (P(x u z₁z₂ u y) independent of the split) reduces, in coordinates e_r = (1−ra)⁻¹ where B(e_s,e_t) = C(s,t)e_se_t, to ONE functional equation N(s,σ,τ,t) = 0, bilinear in values of C")
st1 = sp.factor(nz1[0].subs(t, 0).subs(tau, 0)); st2 = sp.factor(nz2[0].subs(s, 0).subs(sig, 0))
tgt1 = (Cf(s, 0) - Cf(s, sig)) * (Cf(s, 0) - Cf(sig, 0)); tgt2 = (Cf(0, t) - Cf(0, tau)) * (Cf(0, t) - Cf(tau, t))
r1 = sp.simplify(st1 / tgt1); r2 = sp.simplify(st2 / tgt2)
report("L2", r1.free_symbols <= {s, sig} and not r1.has(Cf) and r2.free_symbols <= {t, tau} and not r2.has(Cf),
       "edge dichotomy: t=τ=0 gives (C(s,0)−C(s,σ))(C(s,0)−C(σ,0)) = 0 and the mirror gives (C(0,t)−C(0,τ))(C(0,t)−C(τ,t)) = 0; in the domain k[[·]] each forces C ∈ L or C(s,0) const, resp. C ∈ R or C(0,t) const",
       f"Eq.1 core = {r1}·(…)(…);  Eq.2 core = {r2}·(…)(…)")
k = sp.Function('k'); q = sp.Function('q'); qp = sp.Function('qp')
def En(n, qq): return sig**n * (t - s) * qq(sig, t) * qq(s, sig) + t**n * (s - sig) * qq(sig, t) * qq(s, t) + s**n * (sig - t) * qq(s, sig) * qq(s, t)
st3 = sp.factor(nz1[0].replace(Cf, lambda x_, y_: x_ * y_ * k(x_, y_)).subs(tau, 0).doit())
r3 = sp.simplify(st3 / En(2, k))
desc = all(sp.expand(En(n, q).subs(sig, 0) - s * t * q(s, t) * (t**(n - 1) * q(0, t) - s**(n - 1) * q(s, 0))) == 0 and
           sp.expand(En(n, lambda x_, y_: x_ * y_ * qp(x_, y_)) - s * sig * t * En(n + 1, qp)) == 0 for n in range(2, 10))
report("L3", not r3.has(k) and desc,
       "interior case (C(s,0)=C(0,t)=c=0, C = st·k): the equation is E₂(k) = 0; descent: σ=0 forces tⁿ⁻¹q(0,t) = sⁿ⁻¹q(s,0) ⇒ st | q, and E_n(st·q′) = sσt·E_{n+1}(q′) ⇒ lowest homogeneous part of k vanishes ⇒ k = 0.  THEOREM: identification-type folds are exactly R ∪ L",
       f"Eq.1 core / E₂(k) = {r3}; descent identities n=2..9: {desc}")
# ═════════ 2. Dihedral: module hierarchy of the components ═════════
import fold as FD
from tangent import gvec
from sides import modG, qpolys
from sides2 import build, sub
Leq, Lu, D, N = 3, 10, 3, 3
S, F, eqs, J, big = build('dihedral', Leq, Lu, 1, D)
G = [gvec(F, h) for h in range(len(S.B)) if S.lev[h] <= 1]; p = FD.MOD
tok = lambda v: {F.names[j]: x for j, x in v.items() if F.names[j][0] == 'c'}
import series
def solve(tau):
    importlib.reload(series); st = series.integrate(S, Lu, D, tau, N, Leq=Leq); return st, series.LAST[0]
def absorbers(known, side, orders=(0, 1, 2), corrected=True):
    series.K_ORDER[0] = 3; Fk = series.SFold(S, Lu, D, known)
    lab = lambda w: '·a·'.join(f"{S.side[i]}{S.lev[i]}" for i in w)
    cands = [(lab(w), {w: FD.ONE}) for w in FD.monomials(S, 1, 1)] + [(f"{S.side[b]}{S.lev[b]}", {(b,): FD.ONE}) for b in range(len(S.B)) if S.lev[b] <= 2]
    nc = len(cands); M = max(orders); rows = {}
    for w in [w for d in (0, 1) for w in FD.monomials(S, 1, d)]:
        try:
            Pw = Fk.P({w: FD.ONE}); qel = Fk.sub({w: FD.ONE}, Pw)
            outs = [Fk.P(Fk.mul(x, qel) if side == 'left' else Fk.mul(qel, x)) for _, x in cands]
        except FD.OutOfWindow: continue
        for ci, r in enumerate(outs):
            for wd, cf in r.items():
                for m in orders:
                    for kk in (range(m + 1) if corrected else [0]):
                        v = cf.get((m - kk, None), 0) % p
                        if v: rows.setdefault((w, wd, m), {})[kk * nc + ci] = v
    K, _ = FD.kernel(list(rows.values()), nc * (M + 1 if corrected else 1))
    E0 = FD.Echelon(); out = []
    for kv in K:
        v0 = {i: x for i, x in kv.items() if i < nc}
        if E0.add(v0): out.append({cands[i][0]: x for i, x in v0.items()})
    importlib.reload(FD); return out
def letter(X):  # does the absorber space contain the bare letter a = H0·a·H0 (possibly with Q-terms)?
    E = FD.Echelon(); idx = {}
    for v in X:
        E.add({idx.setdefault(k_, len(idx)): x for k_, x in v.items()})
    if 'H0·a·H0' not in idx: return False
    r = E.reduce({idx['H0·a·H0']: 1}); return not r
def core(X):
    E = FD.Echelon(); idx = {}
    for v in X: E.add({idx.setdefault(k_, len(idx)): x for k_, x in v.items()})
    need = ['H0·a·Q1', 'Q1·a·H0', 'Q1']
    return all(n in idx and not E.reduce({idx[n]: 1}) for n in need)
def profile(tau):
    st, kn = solve(tau); Lx = absorbers(kn, 'left'); st, kn = solve(tau); Rx = absorbers(kn, 'right')
    return st, (letter(Lx), letter(Rx), core(Lx), core(Rx))
prof = {'run-wise': profile({})}
rnd = random.Random(1)
for nm, kt in (("left-only", {'QH'}), ("right-only", {'HQ'})):
    X = modG(S, F, sub(S, F, J, big, kt), G, 2); H = modG(S, F, sub(S, F, J, big, {'QH', 'HQ', 'QQ'}), G, 2)
    low = [j for j, k_ in enumerate(F.names) if k_[0] == 'c' and F.wlev(k_[1]) <= 2]
    E_ = FD.Echelon(); [E_.add({j: g[j] for j in low if j in g}) for g in G + H]
    X = [v for v in X if E_.add({j: v[j] for j in low if j in v})]; tau = {}
    for v in X:
        aa = rnd.randint(1, 60)
        for k_, xv in tok(v).items(): tau[k_] = (tau.get(k_, 0) + aa * xv) % p
    prof[nm] = profile(tau)
Lb = modG(S, F, sub(S, F, J, big, {'QH'}), G, 2); Rb = modG(S, F, sub(S, F, J, big, {'HQ'}), G + Lb, 2)
Mb = modG(S, F, sub(S, F, J, big, set()), G + Lb + Rb, 2)
basis = Lb + Rb + Mb; names = [f'l{i}' for i in range(len(Lb))] + [f'r{i}' for i in range(len(Rb))] + [f'm{i}' for i in range(len(Mb))]
P_ = qpolys(F, eqs, basis, names); syms = sp.symbols(names)
GB = sp.groebner(P_, *syms, order='grevlex'); sq = [x_ for x_ in syms if GB.contains(x_**2)]
P2 = [qq for qq in (sp.expand(pp.subs({x_: 0 for x_ in sq})) for pp in P_) if qq != 0]
sols = sp.solve([pp.subs(sp.Symbol('m0'), 1) for pp in P2], [x_ for x_ in syms if x_ not in sq and str(x_) != 'm0'], dict=True)
rt2 = pow(2, (p + 1) // 4, p)
def mval(e, l3v):
    e = sp.expand(sp.nsimplify(e).subs(sp.Symbol('l3'), l3v)); return (FD.md(e.coeff(sp.sqrt(2), 0)) + FD.md(e.coeff(sp.sqrt(2), 1)) * rt2) % p
taus_pm = []
for i, s_ in enumerate(sols):
    pt = {str(k_): mval(v, 7) for k_, v in s_.items()}; pt.update({'m0': 1, 'l3': 7}); tau = {}
    for nm, v in zip(names, basis):
        aa = pt.get(nm, 0) % p
        for k_, xv in tok(v).items():
            if aa: tau[k_] = (tau.get(k_, 0) + aa * xv) % p
    taus_pm.append(tau); prof[f"two-sided {i}"] = profile(tau)
expect = {'run-wise': (True, True, True, True), 'left-only': (False, True, True, True), 'right-only': (True, False, True, True),
          'two-sided 0': (False, False, True, True), 'two-sided 1': (False, False, True, True)}
report("L4", all(prof[k_][1] == v for k_, v in expect.items()) and all(prof[k_][0][1] == 'ok' for k_ in prof),
       "module hierarchy (to order 2, absorbers allowed ε-corrections): every component is a bimodule over B_core = ⟨Q₀, aQ₀, Q₀a⟩; the letter a is absorbed on both sides (run-wise), on one side (one-sided families), on neither (C_±)",
       "(letter-left, letter-right, core-left, core-right): " + ", ".join(f"{k_}: {v[1]}" for k_, v in prof.items()))
# first-order root letters on C_±
r_p, r_m = (1 + rt2) % p, (1 - rt2) % p
def first_order_letters(tau, side):
    st, kn = solve(tau); X = absorbers(kn, side, orders=(0, 1), corrected=False)
    found = []
    for v in X:
        if 'H1·a·H0' in v or 'H0·a·H1' in v:
            key = 'H1·a·H0' if side == 'left' else 'H0·a·H1'
            c = v.get('H0·a·H0', 0) * pow(v[key], p - 2, p) % p
            found.append('r+' if (c + r_p * r_p) % p == 0 else ('r-' if (c + r_m * r_m) % p == 0 else str(c)))
    return found
fl = [(first_order_letters(tt, 'left'), first_order_letters(tt, 'right')) for tt in taus_pm]
def has_root_letter_order2(tau):
    out = []
    for side, key in (('left', 'H1·a·H0'), ('right', 'H0·a·H1')):
        st, kn = solve(tau); X = absorbers(kn, side)
        out.append(any(key in v for v in X))
    return out
o2 = [has_root_letter_order2(tt) for tt in taus_pm]
okfl = all(not any(x_) for x_ in o2) and len(fl) == 2 and all(len(a_) == 1 and len(b_) == 1 and a_[0] != b_[0] for a_, b_ in fl) and fl[0][0] == fl[1][1]
report("L5", okfl,
       "at FIRST order, each C_± component absorbs a root-tagged letter on each side — (u² − r²)·a (+Q-terms) on the left, a·(u² − r′²) on the right, with {r, r′} = {b ± √s} swapped between the two components; at second order (corrections allowed) no root-tagged letter survives — only B_core",
       f"roots of (left letter, right letter) per component: {fl}; root letter present at order 2 (left, right): {o2}")
print(f"\n{sum(res)}/{len(res)} checks pass")
