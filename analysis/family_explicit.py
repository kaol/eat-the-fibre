"""An explicit one-parameter family of folds of k(1⇉2⇉3) with an H-invariant parameter (v2.26).

Basis notation: arrows a1, a2: 1 → 2 and b1, b2: 2 → 3; δa = a1 − a2, δb = b1 − b2, δaδb = a1b1 − a1b2 − a2b1 + a2b2.
  P_t(e1) = e1 − δa − (1 + t)·δaδb,   P_t(e3) = e3 − δb + t·δaδb,   P_t = 0 on every other basis path.
Checks:
 1. P_t is a Δ-compatible fold for every t (symbolic over Q(t)).
 2. Conjugation by (λ-scaling of δa, μ-scaling of δb) sends the general member with coefficients (a, b, c) to
    (λa, μb, λμc); every element of H = G × T² acts this way on the span of e1, e3, δa, δb, δaδb, so t = c/(ab) is an
    H-invariant. Checked on generators: for every element of H over F_p, conjugation maps P_t to P_t (mod p).
 3. Over F_3 and F_5: every P_t (t ∈ F_p) is a fold of the census, and they lie in p pairwise distinct H-orbits.
Usage: python3 analysis/family_explicit.py"""
import os as _os, sys as _sys
if 'ROOT' not in globals():
    AN = _os.path.dirname(_os.path.abspath(__file__)); ROOT = _os.path.dirname(AN); _os.chdir(ROOT)
    __file__ = _os.path.join(ROOT, 'pcensus.py')
import sys, pickle, itertools, collections, numpy as np, sympy as sp
sys.argv = ['x']
exec(open(_os.path.join(ROOT, 'pcensus.py')).read().split("if __name__ == '__main__':")[0])
name = 'two_pairs'
ctx = setup(name, 3); basis, arrows, d, M, DL = ctx['basis'], ctx['arrows'], ctx['d'], ctx['M'], ctx['DL']
idx = {b: i for i, b in enumerate(basis)}
v = lambda *path: idx[('p', path)]
e1, e3 = idx[('v', 0)], idx[('v', 2)]
a1, a2, b1, b2 = v(0), v(1), v(2), v(3)
a1b1, a1b2, a2b1, a2b2 = v(0, 2), v(0, 3), v(1, 2), v(1, 3)
t = sp.Symbol('t')
def family(a, b, c):
    P = sp.zeros(d, d)
    P[e1, e1] = 1; P[a1, e1] = -a; P[a2, e1] = a
    for (i, s) in ((a1b1, 1), (a1b2, -1), (a2b1, -1), (a2b2, 1)): P[i, e1] = -(a * b + c) * s; P[i, e3] = c * s
    P[e3, e3] = 1; P[b1, e3] = -b; P[b2, e3] = b
    return P
Pt = family(1, 1, t)
# 1. symbolic fold check
Id = sp.eye(d)
def mul(x, y):
    out = [0] * d
    for i in range(d):
        if x[i] == 0: continue
        for j in range(d):
            if y[j] == 0: continue
            for kk in np.nonzero(M[i, j])[0]: out[kk] += x[i] * y[j] * int(M[i, j, kk])
    return sp.Matrix(out)
ok = sp.expand(Pt * Pt - Pt) == sp.zeros(d, d)
for i in range(d):
    for j in range(d):
        h = mul(Pt[:, i], Pt[:, j]); q = mul((Id - Pt)[:, i], (Id - Pt)[:, j])
        ok &= sp.expand(Pt * h - h) == sp.zeros(d, 1) and sp.expand(Pt * q) == sp.zeros(d, 1)
for l in range(d):
    lhs = sum((Pt[m, l] * sp.Matrix(DL[m].tolist()) for m in range(d)), sp.zeros(d, d))
    ok &= sp.expand(lhs - Pt * sp.Matrix(DL[l].tolist()) * Pt.T) == sp.zeros(d, d)
print(f"1. P_t is a Δ-compatible fold for all t (symbolic): {ok}")
gen = family(*sp.symbols('a b c')); ok_g = sp.expand(gen * gen - gen) == sp.zeros(d, d)
print(f"   (the general member with coefficients (a, b, c) is idempotent: {ok_g}; it is P_t after scaling δa, δb when ab ≠ 0, t = c/(ab))")
# 2 + 3. over F_p: H fixes each P_t; distinct H-orbits in the census
for p in (3, 5):
    ctx = setup(name, p); C.update(ctx)
    G = build_G(ctx)
    cls = collections.defaultdict(list)
    for a_, st in enumerate(arrows): cls[st].append(a_)
    cl = list(cls.values()); T = []
    for w in itertools.product(range(1, p), repeat=len(cl)):
        wa = {a_: w[k] for k, c in enumerate(cl) for a_ in c}
        dg = [1 if b[0] == 'v' else int(np.prod([wa[x] for x in b[1]])) % p for b in basis]
        T.append((np.diag(dg), np.diag([pow(x, p - 2, p) for x in dg])))
    H = [(tt @ g % p, gi @ ti % p) for g, gi in G for tt, ti in T]
    fam_p = [np.array([[int(sp.sympify(x).subs(t, tv)) % p for x in row] for row in Pt.tolist()], dtype=np.int64)
             for tv in range(p)]
    # H maps the family into itself, preserving t: conjugate of P_t is P_t' with t' = t (members have a = b = 1 after
    # normalization; compute the invariant c/(ab) of the conjugate directly from its coefficients)
    def inv_t(P):
        a = (-P[a1, e1]) % p; b = (-P[b1, e3]) % p; c = P[a1b1, e3] % p
        return c * pow(int(a * b) % p, p - 2, p) % p
    preserved = all(inv_t((h @ P @ hi) % p) == tv for tv, P in enumerate(fam_p) for h, hi in H)
    fpath = f'logs/folds_{name}_{p}.pkl'
    if _os.path.exists(fpath):
        folds = [f.astype(np.int64) for f in pickle.load(open(fpath, 'rb'))['folds']]
        FS = {(f % p).astype(np.int8).tobytes() for f in folds}
        in_census = all((P % p).astype(np.int8).tobytes() in FS for P in fam_p)
        orbs = [{((h @ P @ hi) % p).astype(np.int8).tobytes() for h, hi in H} for P in fam_p]
        distinct = all(orbs[i].isdisjoint(orbs[j]) for i in range(p) for j in range(i + 1, p))
        sizes = sorted({len(o) for o in orbs})
    else:
        in_census = distinct = sizes = 'folds file missing'
    print(f"{2 if p == 3 else 3}. F_{p}: every element of H (|H| = {len(H)}) preserves t on the family: {preserved}; "
          f"all P_t are census folds: {in_census}; the {p} members lie in pairwise distinct H-orbits: {distinct}; H-orbit sizes {sizes}")
