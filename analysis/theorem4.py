"""Theorem 4 (v2.28): consecutive parallel classes give a one-parameter family of folds that vanish on every vertex, on
every acyclic quiver (no hypothesis on the outer ends).

Let A: u → v and B: v → w be parallel classes with |A|, |B| ≥ 2 (Q acyclic). δ_a = a₁ − a₂, δ_b = b₁ − b₂,
δ_aδ_b = a₁b₁ − a₁b₂ − a₂b₁ + a₂b₂. For x, y, z, w′ ∈ k with x·w′ = y·z − 1 put
    P(a₂) = x·e_u − δ_a + y·δ_aδ_b,   P(b₂) = z·(e_u + e_v) − δ_b + w′·δ_aδ_b,   P = 0 on every other basis path.
Invariant: t = y·z, read off the image (see t_of). Normal form (z = w′ = 1, x = t − 1):
    P_t(a₂) = (t − 1)e_u − δ_a + t·δ_aδ_b,   P_t(b₂) = e_u + e_v − δ_b + δ_aδ_b.
Checks:
 1. symbolic: the general member over Q(x, y, z) (w′ = (yz − 1)/x) and the normal form over Q(t) are folds, on
    several quivers, including ones where Theorem 2's hypothesis fails;
 2. over F_p (1⇉2⇉3 + 1→3 and 1⇉2⇉3, p = 3, 5): t_of(P_t) = t; every conjugate h P_t h⁻¹ (h ∈ G × all class
    weightings) that is again a fold has the same t; with a --save-reps dump, the p members lie in pairwise distinct
    local-torus orbits.
Usage: python3 analysis/theorem4.py"""
import os as _os, sys as _sys
if 'ROOT' not in globals():
    AN = _os.path.dirname(_os.path.abspath(__file__)); ROOT = _os.path.dirname(AN); _os.chdir(ROOT)
    __file__ = _os.path.join(ROOT, 'pcensus.py')
import sys, collections, itertools, pickle, numpy as np, sympy as sp
ARGS = list(sys.argv); sys.argv = ['x']
_MAIN = "if __name__ == '__main__'" + ":"                                 # (not spelled out: this file is split on it too)
exec(open(_os.path.join(ROOT, 'pcensus.py')).read().split(_MAIN)[0])
x, y, z, t = sp.symbols('x y z t')

def members(name):
    """for every consecutive pair of classes (A: u→v, B: v→w), the general member and the normal form"""
    nv, arrows = QUIVERS[name]; basis, M, DL = build(nv, arrows); d = len(basis); idx = {b: i for i, b in enumerate(basis)}
    cls = collections.defaultdict(list)
    for a, st in enumerate(arrows): cls[st].append(a)
    out = []
    for (u, v), A in cls.items():
        for (v2, w), B in cls.items():
            if v2 != v or len(A) < 2 or len(B) < 2: continue
            a1, a2, b1, b2 = A[0], A[1], B[0], B[1]
            def vec(eu=0, ev=0, da=0, db=0, dadb=0):
                V = sp.zeros(d, 1); V[idx[('v', u)]] += eu; V[idx[('v', v)]] += ev
                V[idx[('p', (a1,))]] += da; V[idx[('p', (a2,))]] -= da; V[idx[('p', (b1,))]] += db; V[idx[('p', (b2,))]] -= db
                for (p_, q_), sg in (((a1, b1), 1), ((a1, b2), -1), ((a2, b1), -1), ((a2, b2), 1)):
                    V[idx[('p', (p_, q_))]] += sg * dadb
                return V
            def P_of(xx, yy, zz, ww):
                P = sp.zeros(d, d)
                P[:, idx[('p', (a2,))]] = vec(eu=xx, da=-1, dadb=yy)
                P[:, idx[('p', (b2,))]] = vec(eu=zz, ev=zz, db=-1, dadb=ww)
                return P
            cols = (idx[('v', u)], idx[('p', (a1,))], idx[('p', (b1,))], idx[('p', (a1, b1))])
            out.append(((u + 1, v + 1, w + 1), P_of(x, y, z, (y * z - 1) / x), P_of(t - 1, t, 1, 1), cols))
    return basis, M, DL, d, arrows, out

def is_fold(P, M, DL, d):
    Id = sp.eye(d)
    def mul(xv, yv):
        o = [0] * d
        for i in range(d):
            if xv[i] == 0: continue
            for j in range(d):
                if yv[j] == 0: continue
                for kk in np.nonzero(M[i, j])[0]: o[kk] += xv[i] * yv[j] * int(M[i, j, kk])
        return sp.Matrix(o)
    zero = lambda E: all(sp.simplify(e) == 0 for e in E)
    ok = zero(P * P - P)
    for i in range(d):
        for j in range(d):
            h = mul(P[:, i], P[:, j]); q = mul((Id - P)[:, i], (Id - P)[:, j])
            ok = ok and zero(P * h - h) and zero(P * q)
    for l in range(d):
        lhs = sum((P[m, l] * sp.Matrix(DL[m].tolist()) for m in range(d)), sp.zeros(d, d))
        ok = ok and zero(lhs - P * sp.Matrix(DL[l].tolist()) * P.T)
    return bool(ok)

def t_of(P, p, cols):
    """the invariant of the image: X, Y ∈ Im P with (δ_a, δ_b)-coordinates (−1, 0), (0, −1) (the coefficients of a₁ and
    b₁); t = (a₁b₁-coefficient of X)·(e_u-coefficient of Y). None if Im P is not of this shape."""
    cu, ca1, cb1, ca1b1 = cols
    F = Field(p); Bc = F.rref((P % p).T)[0]                               # rows: a basis of Im P
    if len(Bc) != 2: return None
    Mab = Bc[:, [ca1, cb1]] % p
    if F.rank(Mab) != 2: return None
    XY = (F.inv(Mab) @ Bc) % p                                            # rows with (a₁, b₁)-coordinates e₁, e₂
    X, Y = (-XY[0]) % p, (-XY[1]) % p
    return int(X[ca1b1] * Y[cu] % p)

def fp_checks(name, p):
    ctx = setup(name, p); C.update(ctx); basis, arrows, d = ctx['basis'], ctx['arrows'], ctx['d']
    _, M, DL, _, _, out = members(name); _, _, Pt, cols = out[0]
    key = lambda A_: (A_ % p).astype(np.int8).tobytes()
    fam = [np.array([[int(sp.sympify(e).subs(t, tv)) % p for e in row] for row in Pt.tolist()], dtype=np.int64) for tv in range(p)]
    G = build_G(ctx); cls = collections.defaultdict(list)
    for a_, st in enumerate(arrows): cls[st].append(a_)
    cl = list(cls.values()); T = []
    for wv in itertools.product(range(1, p), repeat=len(cl)):
        wa = {a_: wv[k] for k, c in enumerate(cl) for a_ in c}
        dg = [1 if b[0] == 'v' else int(np.prod([wa[x_] for x_ in b[1]])) % p for b in basis]
        T.append((np.diag(dg), np.diag([pow(x_, p - 2, p) for x_ in dg])))
    DLp = DL.astype(np.int64)
    dcompat = lambda P: ((np.einsum('ml,mab->lab', P, DLp) - np.einsum('ai,lij,bj->lab', P, DLp, P)) % p == 0).all()
    ts = [t_of(P, p, cols) for P in fam]
    pres, nconj = True, 0
    for tv, P in enumerate(fam):
        for g, gi in G:
            for t_, ti in T:
                Q = (t_ @ g % p) @ P @ (gi @ ti % p) % p
                if not dcompat(Q): continue
                nconj += 1; pres &= t_of(Q, p, cols) == tv
    msg = f"t read off P_t is t: {ts == list(range(p))}; all {nconj} fold conjugates under G × T keep t: {pres}"
    rp = f'logs/reps_{name}_{p}.pkl'
    if _os.path.exists(rp):
        D = pickle.load(open(rp, 'rb')); gk = {key(R_.astype(np.int64)): i_ for i_, R_ in enumerate(D['greps'])}
        comp = [D['loc'][next(gk[key(g @ P @ gi)] for g, gi in G if key(g @ P @ gi) in gk)] for P in fam]
        msg += f"; the {p} members lie in pairwise distinct local-torus orbits: {len(set(comp)) == p}"
    return msg

if __name__ == '__main__':
    print("1. symbolic fold check (general member over Q(x, y, z); normal form over Q(t)):")
    for q in ('two_pairs', 'two_pairs_short', 'pair_triple', 'two_pairs_pre', 'two_pairs_post', 'two_pairs_dense'):
        basis, M, DL, d, arrows, out = members(q)
        for (u, v, w), Pg, Pt, _ in out:
            print(f"   {q:16s} classes {u}⇉{v}⇉{w}: general member a fold: {is_fold(Pg, M, DL, d)}; "
                  f"normal form: {is_fold(Pt, M, DL, d)}", flush=True)
    print("2. over F_p:")
    for name in ('two_pairs_short', 'two_pairs'):
        for p in (3, 5):
            print(f"   {name} F_{p}: " + fp_checks(name, p), flush=True)
