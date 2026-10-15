"""Torus-reduced fast census of folds of k(1⇉2→3) over F_p. G ≅ Aff(1) acts on K = span(e1,e2,e3,δ) through
λ = a − b only (gδ = λδ; translations fix K pointwise), so fibres over kernel folds in one torus orbit are G-translates.
Census = Σ over torus-orbit representatives PK of |torus orbit| · |fibre(PK)|; G-orbits of folds = Σ over reps of
Stab_G(PK)-orbits on fibre(PK). Fibres are enumerated exactly as in q123f.py."""
import sys, itertools, time, collections
import numpy as np
exec(open("q123f.py").read().split("\nfolds = []; t0")[0])
blocks = [np.array(v, dtype=np.int64).reshape(2, 2) for v in itertools.product(range(p), repeat=4)]
Gb = [B for B in blocks if (B.sum(axis=0) % p == 1).all() and (int(B[0, 0] * B[1, 1] - B[0, 1] * B[1, 0]) % p)]
def L8(B):
    A = np.eye(d, dtype=np.int64); A[3:5, 3:5] = B; A[6:8, 6:8] = B; return A
inv = {B.tobytes(): next(C for C in blocks if ((B @ C) % p == np.eye(2)).all()) for B in Gb}
lam = lambda B: int(B[0, 0] - B[0, 1]) % p
def fibre(PK):
    out = []
    P = np.zeros((d, d), dtype=np.int64); PKfull = Kmat @ PK % p
    for i in range(3): P[:, i] = PKfull[:, i]
    Pdelta = PKfull[:, 3]
    La, Lb = lift(P, 3), lift(P, 5)
    if La is None or Lb is None: return out
    def set_a2(Pb): Pb[:, :, 4] = (Pb[:, :, 3] - Pdelta) % p
    for Pa in extend(P, 3, La, 4, set_a2):
        for Pab in extend(Pa, 5, Lb, 5):
            L6, L7 = lift(Pab, 6), lift(Pab, 7)
            if L6 is None or L7 is None: continue
            for P6 in extend(Pab, 6, L6, 6):
                for P7 in extend(P6, 7, L7, 7): out.append(P7.copy())
    return out
fib = collections.Counter()
t0 = time.time(); done = set(); N = 0; orbits = 0; sizes = collections.Counter(); nonfilt = 0; reps_done = 0
key = lambda A: (A % p).tobytes()
for PK in Kfolds:
    if key(PK) in done: continue
    torb = {}
    for l in range(1, p):
        gK = np.diag([1, 1, 1, l]); gKi = np.diag([1, 1, 1, pow(l, p - 2, p)])
        torb.setdefault(key(gK @ PK @ gKi), l)
    done |= set(torb)
    F = fibre(PK); reps_done += 1
    N += len(torb) * len(F); fib[len(F)] += len(torb)
    nonfilt += len(torb) * sum(1 for f in F if (f[6:, :6] % p).any())
    stab = [B for B in Gb if key(np.diag([1, 1, 1, lam(B)]) @ PK @ np.diag([1, 1, 1, pow(lam(B), p - 2, p)])) == key(PK)]
    pairs = [(L8(B), L8(inv[B.tobytes()])) for B in stab]
    seen = set()
    for P in F:
        if key(P) in seen: continue
        so = {key(g @ P @ gi) for g, gi in pairs}; seen |= so
        orbits += 1; sizes[len(so) * len(torb)] += 1        # G-orbit size = |Stab-orbit| · |torus orbit of PK|
    if reps_done % 20 == 0: print(f"  reps {reps_done}  kernel folds covered {len(done)}/{len(Kfolds)}  folds {N}  t={time.time()-t0:.0f}s", flush=True)
assert len(done) == len(Kfolds)
print(f"F_{p}: {N} folds (prediction 21p²+6p+12 = {21*p*p+6*p+12}), |G| = {len(Gb)}, orbits = {orbits}, "
      f"orbit sizes {dict(sorted(sizes.items()))}, non-length-filtered folds = {nonfilt}; {reps_done} torus reps, t={time.time()-t0:.0f}s")
print(f"F_{p}: fibre sizes over the {len(Kfolds)} kernel folds: {dict(sorted(fib.items()))}")
