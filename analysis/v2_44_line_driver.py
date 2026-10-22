"""v2.44: idempotent-line family search over many representatives (contexts loaded once).
For each merged orbit of the given signature group and (b, e) class at prime p (one liftable representative each):
non-derivation tangent directions X ∈ T_P Folds give lines P₀ + t·N with N = (1−P₀)XP₀ or P₀X(1−P₀). They are tested
as fold curves mod a 7-digit prime, reconstructed over ℚ, verified exactly (Δ-compatible, idempotent, image and kernel
closed over ℚ[t]) and located in the merged censuses at every prime. A line meeting q orbits at every prime is a
family; families are deduplicated by their orbit sets (two lines sharing all but ≤ 1 orbit per prime count as one).
Usage: python3 analysis/v2_44_line_driver.py SIG_SUBSTRING "(b, e)" MAX_REPS REPS_p.pkl REPS_q.pkl..."""
import os, sys, io, contextlib, itertools, collections
import numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
tg, cls, maxr = sys.argv[1], eval(sys.argv[2]), int(sys.argv[3]); files = sys.argv[4:]
PR = 9999991
def load_fc(f):
    sys.argv = ['family_curves.py', f]
    g = {'__name__': 'fc', '__file__': os.path.join(HERE, 'family_curves.py')}
    with contextlib.redirect_stdout(io.StringIO()):
        exec(open(os.path.join(HERE, 'family_curves.py')).read().split("fams = []; nlift")[0], g)
    return g
G = [load_fc(f) for f in files]; g = G[0]; p = g['p']; d = g['d']
M = np.array(g['M'], np.int64); DL = np.array(g['DL'], np.int64)
sys.argv = ['fibre_rigidity.py']
fr = {'__name__': 'fr', '__file__': os.path.join(HERE, 'fibre_rigidity.py')}
exec(open(os.path.join(HERE, 'fibre_rigidity.py')).read().split("\ndef main():")[0], fr)
Cp = fr['setup'](g['r']['name'], p); CP = fr['setup'](g['r']['name'], PR)
I = np.eye(d, dtype=np.int64)
def fold_resid(P):
    Q = (I - P) % PR
    if ((P @ P - P) % PR).any(): return True
    h = np.einsum('ia,jb,ijm->mab', P, P, M) % PR
    if ((np.einsum('km,mab->kab', P, h) - h) % PR).any(): return True
    qq = np.einsum('ia,jb,ijm->mab', Q, Q, M) % PR
    if (np.einsum('km,mab->kab', P, qq) % PR).any(): return True
    return bool(((np.einsum('ml,mab->lab', P, DL) - np.einsum('ai,lij,bj->lab', P, DL, P)) % PR).any())
def ratrec(a):
    a %= PR; r0, r1, s0, s1 = PR, a, 0, 1
    while r1 * r1 * 2 > PR: qq = r0 // r1; r0, r1 = r1, r0 - qq * r1; s0, s1 = s1, s0 - qq * s1
    return sp.Rational(r1, s1)
t = sp.Symbol('t')
def exact_fold(Pt):
    if not g['is_fold_family'](Pt) or (Pt * Pt - Pt).applyfunc(sp.expand) != sp.zeros(d, d): return False
    Qt = (sp.eye(d) - Pt).applyfunc(sp.expand)
    for a in range(d):
        for b in range(d):
            u = g['smul'](list(Pt[:, a]), list(Pt[:, b]))
            if any(sp.expand(x) for x in (Pt * sp.Matrix(u) - sp.Matrix(u))): return False
            w = g['smul'](list(Qt[:, a]), list(Qt[:, b]))
            if any(sp.expand(x) for x in (Pt * sp.Matrix(w))): return False
    return True
r, sig = g['r'], g['sig_of_loc']
orbits = []
seen_l = set()
for i, l in enumerate(r['loc']):
    if l in seen_l or tg not in str(sig[l]): continue
    if g['lifted'](i) is None: continue
    a = fr['at'](r['greps'][i], Cp, p)
    seen_l.add(l)
    if (a['b'], a['e']) == cls: orbits.append((l, i))
print(f"q={p}: {len(orbits)} orbits of class {cls} with a liftable representative", flush=True)
FAMS = []
for l, i in orbits[:maxr]:
    Pz = g['lifted'](i).astype(np.int64); P0 = Pz % PR
    Tb = fr['nullspace'](fr['jacobian'](P0, dict(d=d, M=M, DL=DL), PR), PR)
    br = np.array([((D @ P0 - P0 @ D) % PR).ravel() for D in CP['Der']]); rbr = fr['rank'](br, PR)
    nT = len(Tb)
    coeffs = [c for c in itertools.product(range(-1, 2), repeat=nT) if any(c) and next(x for x in c if x) > 0] if nT <= 9 else \
             [tuple(int(a_ == b_) for a_ in range(nT)) for b_ in range(nT)]
    seenN = set()
    for c in coeffs:
        X = (np.array(c, dtype=np.int64) @ Tb % PR).reshape(d, d)
        for name, N in (("image moves", (I - P0) @ X @ P0 % PR), ("kernel moves", P0 @ X @ (I - P0) % PR)):
            if not N.any() or fr['rank'](np.concatenate([br, N.reshape(1, -1)]), PR) == rbr: continue
            if any(fold_resid((P0 + s * N) % PR) for s in (1, 2, 3, 5, 7)): continue
            k = np.flatnonzero(N)[0]; Nn = N * pow(int(N.flat[k]), PR - 2, PR) % PR
            if Nn.tobytes() in seenN: continue
            seenN.add(Nn.tobytes())
            Nq = sp.Matrix(d, d, [ratrec(int(x)) for x in Nn.ravel()])
            Pt = (sp.Matrix(Pz.tolist()) + t * Nq).applyfunc(sp.expand)
            if not exact_fold(Pt): continue
            hits = {}
            for g2 in G:
                q = g2['p']
                if any(x.q % q == 0 for x in Nq): hits = None; break
                pts = [np.array(Pt.subs(t, s).applyfunc(lambda z: int(z.p * pow(z.q, -1, q)) % q).tolist(), dtype=np.int64) % q
                       for s in range(q)]
                hits[q] = frozenset(g2['key2loc'].get(g2['key'](P)) for P in pts)
            if hits is None or any(len(h) != q for q, h in hits.items()): continue
            dup = [k_ for k_, (h, _) in enumerate(FAMS) if all(len(h[q] & hits[q]) >= q - 1 for q in hits)]
            if dup: continue
            FAMS.append((hits, (l, i, name, {(a_, b_): str(Nq[a_, b_]) for a_ in range(d) for b_ in range(d) if Nq[a_, b_] != 0})))
            print(f"family #{len(FAMS) - 1}: from orbit {l} rep {i}, {name}; meets q orbits at q = {sorted(hits)}; "
                  f"N = {FAMS[-1][1][3]}", flush=True)
print(f"done: {len(FAMS)} families from {min(maxr, len(orbits))} representatives")
