"""v2.60 (Program 3): the four kinds of folds of Q⁺ = Q + a leaf, for any pair of census quivers.

Q⁺ must extend Q by one vertex v (the last one) and one arrow ℓ (the last one), with Q's arrows first and in the
same order. J = paths through v (e_v and the paths containing ℓ), an ideal; kQ⁺ = kQ ⊕ J.
Per census class of Q⁺ (from --save-reps dumps): split (P(kQ) ⊆ kQ and P(J) ⊆ J), only P(kQ) ⊆ kQ, only P(J) ⊆ J,
neither; the compression π P ι checked as a fold of Q; for split classes λ = coefficient of e_v in P(e_v) and whether
P_J = 0 (the trivial lift P̄ ⊕ 0); fibre sizes over Q's classes for the classes whose compression is a fold.
Usage: python3 analysis/v2_60_leaf_pairs.py Q QPLUS q REPS_DIR"""
import os, sys, pickle, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
import numpy as np
import pcensus as pc
Qn, Qp, q, rd = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
class SafeU(pickle.Unpickler):
    OK = {('numpy._core.numeric', '_frombuffer'), ('numpy.core.numeric', '_frombuffer'), ('numpy', 'dtype'),
          ('numpy._core.multiarray', '_reconstruct'), ('numpy.core.multiarray', '_reconstruct'), ('numpy', 'ndarray')}
    def find_class(self, mod, name):
        if (mod, name) not in self.OK: raise pickle.UnpicklingError(f"refused global {mod}.{name}")
        return super().find_class(mod, name)
def reps(name):
    for fn in (os.path.join(rd, f'{name}_{q}.pkl'), os.path.join(os.path.dirname(HERE), 'logs', f'reps_{name}_{q}.pkl')):
        if os.path.exists(fn): return SafeU(open(fn, 'rb')).load()
    raise FileNotFoundError(name)
nv0, ar0 = pc.QUIVERS[Qn]; nv1, ar1 = pc.QUIVERS[Qp]
assert nv1 == nv0 + 1 and len(ar1) == len(ar0) + 1 and ar1[:len(ar0)] == ar0 and nv0 in ar1[-1]
leaf_kind = 'sink' if ar1[-1][1] == nv0 else 'source'
ell = len(ar0)
R1, R0 = reps(Qp), reps(Qn)
c1, c0 = pc.setup(Qp, q), pc.setup(Qn, q)
b1, b0 = c1['basis'], c0['basis']; d1, d0 = len(b1), len(b0)
Jx = [i for i, b in enumerate(b1) if (b[0] == 'v' and b[1] == nv0) or (b[0] == 'p' and ell in b[1])]
Qx = [b1.index(b) for b in b0]
assert sorted(Jx + Qx) == list(range(d1))
ev = b1.index(('v', nv0))
M0 = np.asarray(c0['M']); DL0 = np.asarray(c0['DL'])
def is_fold(P, M, DL, p):
    d = P.shape[0]; P = P % p
    if ((P @ P - P) % p).any(): return False
    Qm = (np.eye(d, dtype=np.int64) - P) % p
    for A, check in ((P, lambda u: (P @ u - u) % p), (Qm, lambda u: (P @ u) % p)):
        prod = np.einsum('ia,jb,ijk->abk', A, A, M) % p
        if check(prod.reshape(-1, d).T).any(): return False
    lhs = np.einsum('ml,mab->lab', P, DL) % p
    rhs = np.einsum('ax,lxy,by->lab', P, DL, P) % p
    return not ((lhs - rhs) % p).any()
key0 = {}
for P, lab in zip(R0['greps'], R0['loc']):
    P = P.astype(np.int64)
    for g, gi in pc.build_G(c0): key0[((g @ P @ gi) % q).astype(np.int8).tobytes()] = lab
cls = collections.defaultdict(collections.Counter); img = collections.defaultdict(set); lam = collections.defaultdict(set)
triv = set(); nofold = 0; h1 = collections.Counter()
for P, lab in zip(R1['greps'], R1['loc']):
    P = P.astype(np.int64) % q
    J = not P[np.ix_(Qx, Jx)].any(); Qs = not P[np.ix_(Jx, Qx)].any()
    PQ = P[np.ix_(Qx, Qx)]
    C = is_fold(PQ, M0, DL0, q)
    if C:
        k0 = key0.get(PQ.astype(np.int8).tobytes()); assert k0 is not None
        img[lab].add(k0)
    cls[lab][(J, Qs, C)] += 1
    if J and not Qs:                                   # H1: the J-component of P on kQ lies in k e_v
        h1[not np.delete(P[np.ix_(Jx, Qx)], Jx.index(ev), axis=0).any()] += 1
    if J and Qs:
        lam[lab].add(int(P[ev, ev]))
        if not P[np.ix_(Jx, Jx)].any(): triv.add(lab)
kind = collections.Counter(); bad = 0
for lab, c in cls.items():
    ks = set(c)
    if len(ks) > 1: bad += 1
    J, Qs, C = next(iter(ks))
    kind[('split' if J and Qs else 'kQ only' if Qs else 'J only' if J else 'neither') + (', compression a fold' if C else ', compression not a fold')] += 1
n0 = len(set(R0['loc'])); n1 = len(cls)
over = collections.defaultdict(set)
for lab in cls:
    if len(img[lab]) == 1: over[next(iter(img[lab]))].add(lab)
print(f"{Qn} → {Qp} ({leaf_kind} leaf), q = {q}: census classes {n0} → {n1}")
print("  kinds: " + "; ".join(f"{k}: {v}" for k, v in sorted(kind.items(), key=lambda kv: -kv[1])) + (f"; classes of mixed kind: {bad}" if bad else ""))
print(f"  split classes by λ: {dict(collections.Counter(frozenset(v) for v in lam.values()))}; trivial lifts P̄ ⊕ 0: {len(triv)}"
      f" (over distinct Q classes: {len(set(next(iter(img[l])) for l in triv))})")
print(f"  J-only G-orbit representatives with P(kQ) ⊆ kQ ⊕ k e_v: {h1[True]}, without: {h1[False]}")
print(f"  classes over several Q classes: {sum(1 for l in cls if len(img[l]) > 1)}; Q classes hit {len(over)} of {n0}; "
      f"fibre sizes {dict(sorted(collections.Counter(len(v) for v in over.values()).items()))}")
