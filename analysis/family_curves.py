"""v2.36: inner-curve family search on any census with two parallel classes of size 2.

For every local-torus orbit (one liftable G-orbit representative P₀ each; symmetric residues giving a fold over ℤ) and
every direction x = α δ_a + β δ_b (+ γ δ_aδ_b when the classes are consecutive), with small integer coefficients up to
sign and scale:
  P_t = Inn(1 + t·x) P₀ Inn(1 + t·x)⁻¹ is checked to be a fold for all t over ℚ[t] (Δ-compatibility as a polynomial
  identity; the rest is automatic by transport, v2.34), and its 𝔽_q-points are located in the census.
Local-torus orbits are located through G-orbits: the census G-orbit representatives carry their local-torus label, and
a fold is matched to its G-orbit by a minimum hash over G (G = column-sum-1 substitutions, exactly as in the census).
A curve meeting ≥ 3 local orbits is reported as a family.
Usage: python3 analysis/family_curves.py REPS.pkl
"""
import os, sys, pickle, itertools, collections, json
import numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
_FC_FILE = sys.argv[1]; sys.argv = ['x']
__file__ = os.path.join(HERE, 'family_signatures.py')
exec(open(__file__).read().split("\ndef main():")[0])

r = pickle.load(open(_FC_FILE, 'rb')); name, p = r['name'], r['p']
nv, arrows = QUIVERS[name]; basis, M, DL = build(nv, arrows)
M = np.array(M, np.int64); DL = np.array(DL, np.int64); d = len(basis)
idx = {b: i for i, b in enumerate(basis)}
classes = collections.defaultdict(list)
for ai, st in enumerate(arrows): classes[st].append(ai)
big = [c for c in classes.values() if len(c) == 2]
assert len(big) == 2, "needs exactly two parallel classes of size 2"
(a1, a2), (b1, b2) = big
v = lambda *path: idx[('p', path)]

# ---- G (column-sum-1 substitutions on the classes of size ≥ 2) as matrices, and a minimum-hash key
def gl1(m):
    for vals in itertools.product(range(p), repeat=m * m):
        A = np.array(vals, np.int64).reshape(m, m)
        if round(np.linalg.det(A)) % p and ((A.sum(0) - 1) % p == 0).all(): yield A
def subst(mats):
    img = {}
    for ai in range(len(arrows)):
        vec = np.zeros(d, np.int64); vec[v(ai)] = 1; img[ai] = vec
    for Mc, c in zip(mats, big):
        for jx, x in enumerate(c):
            vec = np.zeros(d, np.int64)
            for iy, y in enumerate(c): vec[v(y)] = Mc[iy, jx] % p
            img[x] = vec
    Phi = np.zeros((d, d), np.int64)
    for i, b in enumerate(basis):
        if b[0] == 'v': Phi[i, i] = 1; continue
        vec = img[b[1][0]]
        for x in b[1][1:]: vec = np.einsum('i,j,ijm->m', vec, img[x], M) % p
        Phi[:, i] = vec
    return Phi
G1 = list(gl1(2))
Gm = [(subst((A, B)), subst((np.array(sp.Matrix(A.tolist()).inv_mod(p).tolist(), np.int64),
                              np.array(sp.Matrix(B.tolist()).inv_mod(p).tolist(), np.int64)))) for A in G1 for B in G1]
H = np.array([g for g, _ in Gm], np.float64); Hi = np.array([gi for _, gi in Gm], np.float64)
W = np.random.default_rng(11).integers(1, 2 ** 52, size=d * d).astype(np.int64)
def key(P):
    Q = np.mod(H @ (P.astype(np.float64) % p) @ Hi, p).reshape(len(H), -1).astype(np.int64)
    return int((Q * W).sum(1).min())

loc = r['loc']
key2loc = {}
for i, P in enumerate(r['greps']): key2loc[key(P)] = loc[i]
assert len(key2loc) == len(r['greps']), "G-orbit keys collide"
locrep = {}
for i in range(len(r['greps'])): locrep.setdefault(loc[i], i)
sig_of_loc = {l: signature(r['greps'][i], M, p) for l, i in locrep.items()}

# ---- inner curves
t = sp.Symbol('t')
da = {v(a1): 1, v(a2): -1}; db = {v(b1): 1, v(b2): -1}
prod = np.einsum('i,j,ijm->m', np.array([da.get(q, 0) for q in range(d)]), np.array([db.get(q, 0) for q in range(d)]), M)
consecutive = bool(prod.any())
dirs = []
for a, b, c in itertools.product(range(-2, 3), repeat=3):
    if not consecutive and c: continue
    vec = (a, b, c)
    if vec == (0, 0, 0): continue
    first = next(z for z in vec if z)
    g = np.gcd.reduce([abs(z) for z in vec if z])
    if first < 0 or g != 1: continue
    dirs.append(vec)
one = [1 if b[0] == 'v' else 0 for b in basis]
def smul(x, y): return [sp.expand(sum(int(M[i, j, m]) * x[i] * y[j] for i in range(d) for j in range(d) if M[i, j, m])) for m in range(d)]
E = lambda j: [1 if i == j else 0 for i in range(d)]
DLs = [sp.Matrix(DL[m].tolist()) for m in range(d)]
def inner(vec):
    a, b, c = vec
    x = [t * (a * da.get(q, 0) + b * db.get(q, 0) + c * int(prod[q])) for q in range(d)]
    u = [one[q] + x[q] for q in range(d)]; uinv = one[:]; pw = one[:]
    for _ in range(d): pw = smul(pw, [-y for y in x]); uinv = [sp.expand(z + w) for z, w in zip(uinv, pw)]
    return (sp.Matrix([smul(smul(u, E(j)), uinv) for j in range(d)]).T,
            sp.Matrix([smul(smul(uinv, E(j)), u) for j in range(d)]).T)
INN = {vec: inner(vec) for vec in dirs}
def is_fold_family(Ps):
    for l in range(d):
        lhs = sum((Ps[m, l] * DLs[m] for m in range(d)), sp.zeros(d, d))
        if (lhs - Ps * DLs[l] * Ps.T).applyfunc(sp.expand) != sp.zeros(d, d): return False
    return True
def lifted(i):
    Pz = ((r['greps'][i].astype(np.int64) % p + p // 2) % p) - p // 2
    Q = np.eye(d, dtype=np.int64) - Pz
    ok = (Pz @ Pz == Pz).all()
    h = np.einsum('ia,jb,ijm->mab', Pz, Pz, M); ok &= (np.einsum('km,mab->kab', Pz, h) == h).all()
    ok &= not np.einsum('km,mab->kab', Pz, np.einsum('ia,jb,ijm->mab', Q, Q, M)).any()
    ok &= (np.einsum('ml,mab->lab', Pz, DL) == np.einsum('ai,lij,bj->lab', Pz, DL, Pz)).all()
    return Pz if ok else None

fams = []; nlift = 0; ncurves = 0
for l, i in sorted(locrep.items()):
    Pz = lifted(i)
    if Pz is None: continue
    nlift += 1
    P0 = sp.Matrix(Pz.tolist())
    for vec, (Inn, Innv) in INN.items():
        Ps = (Inn * P0 * Innv).applyfunc(sp.expand)
        if Ps == P0 or not is_fold_family(Ps): continue
        ncurves += 1
        hit = set()
        for sv in range(p):
            Pn = np.array(Ps.subs(t, sv).tolist(), dtype=np.int64) % p
            hit.add(key2loc.get(key(Pn)))
        assert None not in hit, "a curve point is not a census fold"
        if len(hit) >= 3 and not any(hit & f for f, _ in fams):
            fams.append((hit, (i, vec)))
            print(f"family: {len(hit)} local orbits, from rep {i} along {vec}, signature classes "
                  f"{sorted({sig_of_loc[x][0] for x in hit})}", flush=True)
grow_sigs = collections.Counter(sig_of_loc[x] for f, _ in fams for x in f)
print(json.dumps(dict(name=name, p=p, local_orbits=len(locrep), liftable=nlift, directions=len(dirs),
                      fold_curves=ncurves, families=len(fams),
                      family_signatures=sorted({str((s[0], s[2], s[3])) for s in grow_sigs}))))
