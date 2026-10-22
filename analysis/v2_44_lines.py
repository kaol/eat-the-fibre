"""v2.44: at a lifted fold P₀ of a (b, e) = (1, 0) orbit, find the tangent directions of Folds that are not
derivation directions, and test the idempotent lines P₀ + t·N with N = (1 − P₀)XP₀ (image moves, kernel fixed) or
N = P₀X(1 − P₀) (kernel moves, image fixed), X ranging over a basis of T_P Folds plus small combinations.

Linear algebra is done modulo a 7-digit prime PR (no int64 overflow; ranks agree with ℚ for all but finitely many
primes). A line is a fold curve iff the fold equations, polynomials of degree ≤ 3 in t, vanish; that is tested at 5
values of t mod PR. Candidate N are rationally reconstructed and then verified exactly over ℚ[t] with sympy; their
𝔽_q-points are located in the merged censuses.
Usage: python3 analysis/v2_44_lines.py REPS_p.pkl REP_INDEX [REPS_q.pkl ...]"""
import os, sys, io, contextlib, itertools
import numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
rf, ri = sys.argv[1], int(sys.argv[2]); others = sys.argv[3:]
PR = 9999991
_FC = {}
def load_fc(f):
    if f in _FC: return _FC[f]
    sys.argv = ['family_curves.py', f]
    g = {'__name__': 'fc', '__file__': os.path.join(HERE, 'family_curves.py')}
    with contextlib.redirect_stdout(io.StringIO()):
        exec(open(os.path.join(HERE, 'family_curves.py')).read().split("fams = []; nlift")[0], g)
    _FC[f] = g
    return g
g = load_fc(rf); d = g['d']
M = np.array(g['M'], np.int64); DL = np.array(g['DL'], np.int64)
sys.argv = ['fibre_rigidity.py']
fr = {'__name__': 'fr', '__file__': os.path.join(HERE, 'fibre_rigidity.py')}
exec(open(os.path.join(HERE, 'fibre_rigidity.py')).read().split("\ndef main():")[0], fr)
Pz = g['lifted'](ri); assert Pz is not None
P0 = Pz.astype(np.int64) % PR
C = fr['setup'](g['r']['name'], PR)
Cj = dict(d=d, M=M, DL=DL)
J = fr['jacobian'](P0, Cj, PR)
Tb = fr['nullspace'](J, PR)                                    # rows: basis of T_P (vectorised d×d)
br = np.array([((D @ P0 - P0 @ D) % PR).ravel() for D in C['Der']])
rbr = fr['rank'](br, PR)
print(f"mod {PR}: dim T = {len(Tb)}, dim [Der, P] = {rbr}, dim [Der,P] ∩ T = "
      f"{len(Tb) + rbr - fr['rank'](np.concatenate([Tb, br]), PR)}", flush=True)
I = np.eye(d, dtype=np.int64)
def fold_resid(P):
    """all fold equations at P mod PR (zero iff fold)"""
    Q = (I - P) % PR
    e1 = (P @ P - P) % PR
    h = np.einsum('ia,jb,ijm->mab', P, P, M) % PR
    e2 = (np.einsum('km,mab->kab', P, h) - h) % PR
    qq = np.einsum('ia,jb,ijm->mab', Q, Q, M) % PR
    e3 = np.einsum('km,mab->kab', P, qq) % PR
    e4 = (np.einsum('ml,mab->lab', P, DL) - np.einsum('ai,lij,bj->lab', P, DL, P)) % PR
    return any(x.any() for x in (e1, e2, e3, e4))
def line_ok(N):
    return all(not fold_resid((P0 + s * N) % PR) for s in (1, 2, 3, 5, 7))
cands = []
nT = len(Tb)
coeffs = [c for c in itertools.product(range(-1, 2), repeat=nT) if any(c)] if nT <= 9 else [tuple(int(i == j) for i in range(nT)) for j in range(nT)]
for c in coeffs:
    X = (np.array(c, dtype=np.int64) @ Tb % PR).reshape(d, d)
    for name, N in (("image moves", (I - P0) @ X @ P0 % PR), ("kernel moves", P0 @ X @ (I - P0) % PR)):
        if not N.any(): continue
        if fr['rank'](np.concatenate([br, N.reshape(1, -1)]), PR) == rbr: continue   # a derivation direction
        if line_ok(N): cands.append((name, c, N))
print(f"non-derivation idempotent lines that are fold curves mod {PR}: {len(cands)}", flush=True)
def ratrec(a):
    a %= PR; r0, r1, s0, s1 = PR, a, 0, 1
    while r1 * r1 * 2 > PR: qq = r0 // r1; r0, r1 = r1, r0 - qq * r1; s0, s1 = s1, s0 - qq * s1
    return sp.Rational(r1, s1)
t = sp.Symbol('t')
seen = set(); FAMS = []
for name, c, N in cands:
    # normalise: scale so the first nonzero entry is 1, then reconstruct over Q
    k = np.flatnonzero(N)[0]; inv = pow(int(N.flat[k]), PR - 2, PR)
    Nn = N * inv % PR
    key = Nn.tobytes()
    if key in seen: continue
    seen.add(key)
    Nq = sp.Matrix(d, d, [ratrec(int(x)) for x in Nn.ravel()])
    Pt = (sp.Matrix(Pz.tolist()) + t * Nq).applyfunc(sp.expand)
    exact = g['is_fold_family'](Pt) and (Pt * Pt - Pt).applyfunc(sp.expand) == sp.zeros(d, d)
    print(f"{name}, X-coefficients {c}: exact Δ-compatible and idempotent over Q[t]: {exact}; "
          f"denominators {sorted({x.q for x in Nq})}", flush=True)
    if not exact: continue
    hits = {}
    for f2 in [rf] + others:
        g2 = load_fc(f2); q = g2['p']
        if any(x.q % q == 0 for x in Nq): print(f"  F_{q}: denominator divisible by q"); continue
        hit = [g2['key2loc'].get(g2['key'](np.array(Pt.subs(t, s).applyfunc(lambda z: int(z.p * pow(z.q, -1, q)) % q).tolist(), dtype=np.int64) % q)) for s in range(q)]
        hits[q] = frozenset(hit)
        print(f"  F_{q}: points lie in {len(set(hit))} orbits; all census folds: {None not in hit}", flush=True)
    same = [k for k, h in enumerate(FAMS) if h == hits]
    print(f"  same orbits as line #{same[0]}" if same else f"  new line #{len(FAMS)}", flush=True)
    if not same: FAMS.append(hits)
    print("  N (nonzero entries):", {(i, j): str(Nq[i, j]) for i in range(d) for j in range(d) if Nq[i, j] != 0}, flush=True)
