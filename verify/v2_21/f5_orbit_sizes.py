"""Orbit-size distribution of G(F_p) = {column-sum-1 substitutions of α1, α2} ≅ Aff(1) on the folds stored in a census
checkpoint (from `q123v.py p --checkpoint FILE`). Prediction from the F_3 census plus stabilizer types over Q (stabcheck.py):
  23 orbits of size 1, 11 of size p−1 (unipotent stabilizer), 16 of size p (torus stabilizer), 21 of size p(p−1) (free),
hence N(p) = 21p² + 6p + 12 folds in 71 orbits.  Usage: python3 f5_orbit_sizes.py f5_checkpoint.pkl"""
import sys, pickle, itertools, collections
import numpy as np
ck = pickle.load(open(sys.argv[1] if len(sys.argv) > 1 else 'f5_checkpoint.pkl', 'rb'))
p, folds = ck['prime'], ck['folds']
assert ck['next_kernel'] == ck['kernel_count'], f"checkpoint incomplete: {ck['next_kernel']}/{ck['kernel_count']} kernel folds"
d = folds[0].shape[0]
def lift_block(B):
    A = np.eye(d, dtype=np.int64); A[3:5, 3:5] = B; A[6:8, 6:8] = B; return A
blocks = [np.array(v, dtype=np.int64).reshape(2, 2) for v in itertools.product(range(p), repeat=4)]
G = [B for B in blocks if (B.sum(axis=0) % p == 1).all() and int(round(np.linalg.det(B))) % p]
inv = {B.tobytes(): next(C for C in blocks if ((B @ C) % p == np.eye(2)).all()) for B in G}
pairs = [(lift_block(B), lift_block(inv[B.tobytes()])) for B in G]
seen, sizes = set(), collections.Counter()
for P in folds:
    if P.tobytes() in seen: continue
    orb = {((g @ P @ gi) % p).tobytes() for g, gi in pairs}
    seen |= orb; sizes[len(orb)] += 1
pred = collections.Counter({1: 23}); pred[p - 1] += 11; pred[p] += 16; pred[p * (p - 1)] += 21
print(f"F_{p}: {len(folds)} folds, |G| = {len(G)}, {sum(sizes.values())} orbits; sizes {dict(sorted(sizes.items()))}")
print(f"prediction: {21*p*p + 6*p + 12} folds; sizes {dict(sorted(pred.items()))}")
print("MATCH" if sizes == pred and len(folds) == 21*p*p + 6*p + 12 else "MISMATCH")
