"""Fast census of Δ-compatible folds of k(1⇉2→3) over F_p (same kernel fibration as q123v.py).
Each K-component (of α1, β, α1β, α2β in turn) is enumerated alone (p^4 candidates) and the masked necessary conditions
on the columns known so far are applied before the next one: p^8 + p^8 per kernel fold becomes ≈ 4·p^4·(survivors).
The last stage checks idempotency and both subalgebra conditions on all columns; Δ-compatibility holds by construction.
Usage: python3 q123f.py p  →  fold count, orbit count and orbit-size distribution under G = Aff(1)."""
import sys, itertools, time, collections
import numpy as np
exec(open("q123v.py").read().split("\nallk = ")[0])        # algebra, kernel folds, stage_ok
allk = np.array(list(itertools.product(range(p), repeat=4)), dtype=np.int64) @ np.array(Kb) % p   # all of K, (p^4, d)
def lift(Pcur, x_idx):
    T = delta(I[x_idx]); acc = np.zeros((d, d), dtype=np.int64)
    for a in range(d):
        for b in range(d):
            if T[a, b]: acc = (acc + T[a, b] * np.outer(Pcur[:, a], Pcur[:, b])) % p
    sol = solve(DC, acc.ravel())
    if sol is None: return None
    return sum((sol[k] * I[Cb[k]] for k in range(4)), np.zeros(d, dtype=np.int64)) % p
def extend(P, col, L, upto, extra=None):
    """All K-completions of column col (= L + k), keeping those passing the masked checks on columns 0..upto."""
    Pb = np.repeat(P[None], len(allk), axis=0)
    Pb[:, :, col] = (L + allk) % p
    if extra is not None: extra(Pb)
    return Pb[stage_ok(Pb, list(range(upto + 1)))]
folds = []; t0 = time.time(); fibre = collections.Counter()
for ki, PK in enumerate(Kfolds):
    P = np.zeros((d, d), dtype=np.int64); PKfull = Kmat @ PK % p
    for i in range(3): P[:, i] = PKfull[:, i]
    Pdelta = PKfull[:, 3]
    La, Lb = lift(P, 3), lift(P, 5)
    if La is None or Lb is None: fibre[0] += 1; continue
    n0 = len(folds)
    def set_a2(Pb): Pb[:, :, 4] = (Pb[:, :, 3] - Pdelta) % p          # α2 = α1 − δ
    for Pa in extend(P, 3, La, 4, set_a2):
        for Pab in extend(Pa, 5, Lb, 5):
            L6, L7 = lift(Pab, 6), lift(Pab, 7)
            if L6 is None or L7 is None: continue
            for P6 in extend(Pab, 6, L6, 6):
                for P7 in extend(P6, 7, L7, 7): folds.append(P7.copy())
    fibre[len(folds) - n0] += 1
print(f"F_{p}: algebra folds of K: {len(Kfolds)}; fibre-size histogram over them: {dict(sorted(fibre.items()))}; t={time.time()-t0:.0f}s", flush=True)
assert len({f.tobytes() for f in folds}) == len(folds)
nonfilt = sum(1 for f in folds if (f[6:, :6] % p).any())
blocks = [np.array(v, dtype=np.int64).reshape(2, 2) for v in itertools.product(range(p), repeat=4)]
Gb = [B for B in blocks if (B.sum(axis=0) % p == 1).all() and (int(B[0, 0] * B[1, 1] - B[0, 1] * B[1, 0]) % p)]
def L8(B):
    A = np.eye(d, dtype=np.int64); A[3:5, 3:5] = B; A[6:8, 6:8] = B; return A
inv = {B.tobytes(): next(C for C in blocks if ((B @ C) % p == np.eye(2)).all()) for B in Gb}
pairs = [(L8(B), L8(inv[B.tobytes()])) for B in Gb]
seen, sizes = set(), collections.Counter()
for P in folds:
    if P.tobytes() in seen: continue
    orb = {((g @ P @ gi) % p).tobytes() for g, gi in pairs}
    seen |= orb; sizes[len(orb)] += 1
print(f"F_{p}: {len(folds)} folds (prediction 21p²+6p+12 = {21*p*p+6*p+12}), |G| = {len(Gb)}, orbits = {sum(sizes.values())}, "
      f"orbit sizes {dict(sorted(sizes.items()))}, non-length-filtered folds = {nonfilt}")
