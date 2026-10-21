"""v2.36: closed forms for the families of 1⇉2⇉3.

For every growing signature group (v2.35) and every liftable member P0 (a census representative whose symmetric
residues give a fold over ℤ), and for x ∈ {δ_a, δ_b, δ_aδ_b}:
  1. P_s = Inn(1 + s·x) P0 Inn(1 + s·x)⁻¹ is checked to be a fold for all s, symbolically over ℚ[s]
     (image and kernel closure are automatic by transport, v2.34; Δ-compatibility is checked as a polynomial identity);
  2. its points over 𝔽_p are located in the census: local-torus orbit = (G × T²)-orbit on 1⇉2⇉3 (every fold admits every
     class weighting, v2.24), identified by a minimum hash over G × T²;
  3. a curve meeting ≥ 3 local orbits is a family (curves meeting 2 are degenerations: a generic orbit plus a special point); families are identified by their sets of local orbits.
Output: per signature group, the families found (number of local orbits met over 𝔽_p) and one closed form each.
Usage: python3 analysis/family_normal_forms.py REPS_two_pairs_p.pkl
"""
import os, sys, pickle, itertools, collections, json
import numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
_NF_FILE = sys.argv[1]; sys.argv = ['x']
__file__ = os.path.join(HERE, 'family_signatures.py')
exec(open(__file__).read().split("\ndef main():")[0])

r = pickle.load(open(_NF_FILE, 'rb')); name, p = r['name'], r['p']
assert name == 'two_pairs'
nv, arrows = QUIVERS[name]; basis, M, DL = build(nv, arrows)
M = np.array(M, np.int64); DL = np.array(DL, np.int64); d = len(basis)
idx = {b: i for i, b in enumerate(basis)}
pname = lambda b: 'e%d' % (b[1] + 1) if b[0] == 'v' else ''.join('ab'[x // 2] + str(x % 2 + 1) for x in b[1])
names = [pname(b) for b in basis]
v = lambda *path: idx[('p', path)]

# ---- H = G × T² as matrices
def gl(m):
    for vals in itertools.product(range(p), repeat=m * m):
        A = np.array(vals, np.int64).reshape(m, m)
        if round(np.linalg.det(A)) % p and ((A.sum(0) - 1) % p == 0).all(): yield A
cls = [[0, 1], [2, 3]]
def subst(Ms, lam):
    img = {}
    for (Mc, c), l in zip(zip(Ms, cls), lam):
        for jx, x in enumerate(c):
            vec = np.zeros(d, np.int64)
            for iy, y in enumerate(c): vec[v(y)] = Mc[iy, jx] * l % p
            img[x] = vec
    Phi = np.zeros((d, d), np.int64)
    for i, b in enumerate(basis):
        if b[0] == 'v': Phi[i, i] = 1; continue
        vec = img[b[1][0]]
        for x in b[1][1:]: vec = np.einsum('i,j,ijm->m', vec, img[x], M) % p
        Phi[:, i] = vec
    return Phi
Gc = list(gl(2))
H = np.array([subst((A, B), (la, lb)) for A in Gc for B in Gc for la in range(1, p) for lb in range(1, p)], np.float64)
Hi = np.array([np.array(sp.Matrix(h.astype(np.int64).tolist()).inv_mod(p).tolist(), np.int64) for h in H], np.float64)
W = np.random.default_rng(7).integers(1, 2 ** 52, size=d * d).astype(np.float64)
def key(P):
    Q = np.mod(H @ (P.astype(np.float64) % p) @ Hi, p).reshape(len(H), -1)
    hv = (Q.astype(np.int64) * W.astype(np.int64)).sum(1)               # wraps mod 2^64
    return int(hv.min())

# ---- census local orbits: key → loc label, and signature per loc
loc = r['loc']; locrep = {}
for i in range(len(r['greps'])): locrep.setdefault(loc[i], i)
key2loc = {}; sig_of_loc = {}
for l, i in locrep.items():
    key2loc[key(r['greps'][i])] = l
    sig_of_loc[l] = signature(r['greps'][i], M, p)
assert len(key2loc) == len(locrep), "local-orbit keys collide"

# ---- inner curves
s_ = sp.Symbol('s')
def ad_elems():
    da = [0] * d; da[v(0)] = 1; da[v(1)] = -1
    db = [0] * d; db[v(2)] = 1; db[v(3)] = -1
    dadb = [0] * d
    for (i, j), sg in (((0, 2), 1), ((0, 3), -1), ((1, 2), -1), ((1, 3), 1)): dadb[v(i, j)] = sg
    return {'δa': da, 'δb': db, 'δaδb': dadb}
X = ad_elems()
one = [1 if b[0] == 'v' else 0 for b in basis]
def smul(a, b): return [sp.expand(sum(int(M[i, j, m]) * a[i] * b[j] for i in range(d) for j in range(d) if M[i, j, m])) for m in range(d)]
E = lambda j: [1 if i == j else 0 for i in range(d)]
DLs = [sp.Matrix(DL[m].tolist()) for m in range(d)]
def inner(x):
    u = [one[i] + s_ * x[i] for i in range(d)]; uinv = one[:]; pw = one[:]
    for _ in range(3): pw = smul(pw, [-s_ * y for y in x]); uinv = [sp.expand(a + b) for a, b in zip(uinv, pw)]
    return (sp.Matrix([smul(smul(u, E(j)), uinv) for j in range(d)]).T,
            sp.Matrix([smul(smul(uinv, E(j)), u) for j in range(d)]).T)
INN = {k: inner(x) for k, x in X.items()}
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

growing = {}
for l, s in sig_of_loc.items(): growing.setdefault(s, []).append(l)
TARGET = json.loads(os.environ.get('NF_GROUPS', 'null'))
families = collections.defaultdict(dict)                 # signature → {frozenset(locs): (rep, x, Ps)}
for sig, locs in sorted(growing.items(), key=lambda kv: str(kv[0])):
    tag = f"{sig[0]}|{sig[2]}|{sig[3]}"
    if TARGET and tag not in TARGET: continue
    covered = set()
    for l in locs:
        if l in covered: continue
        i = locrep[l]; Pz = lifted(i)
        if Pz is None: continue
        P0 = sp.Matrix(Pz.tolist())
        for xn, (Inn, Innv) in INN.items():
            Ps = (Inn * P0 * Innv).applyfunc(sp.expand)
            if Ps == P0 or not is_fold_family(Ps): continue
            hit = set()
            for sv in range(p):
                Pn = np.array(Ps.subs(s_, sv).tolist(), dtype=np.int64) % p
                hit.add(key2loc.get(key(Pn), None))
            hit.discard(None)
            if len(hit) >= 3:
                fam = frozenset(hit)
                if not any(fam & f for f in families[tag]):
                    families[tag][fam] = (i, xn, Ps)
                covered |= hit
    print(f"{tag}: {len(locs)} local orbits, families found {len(families[tag])}, "
          f"sizes {sorted(len(f) for f in families[tag])}", flush=True)
pickle.dump({t: {tuple(sorted(f)): (i, xn, str(Ps)) for f, (i, xn, Ps) in fs.items()} for t, fs in families.items()},
            open(os.environ.get('NF_OUT', '/dev/null'), 'wb'))
