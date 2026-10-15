"""Vectorized version of q123.py (same algebra, same fibration), validated against it over F_3."""
import sys, itertools, time
import numpy as np
exec(open("q123.py").read().split("\nfolds = []\n")[0])
def bmul(X, Y):            # batched products: X, Y (N, d) -> (N, d)
    return np.einsum('ni,nj,ijk->nk', X, Y, M) % p
def bapply(Pb, V):          # Pb (N, d, d), V (N, d) -> (N, d)
    return np.einsum('nij,nj->ni', Pb, V) % p
def stage_ok(Pb, cols):
    ok = np.ones(Pb.shape[0], dtype=bool)
    for i in cols:
        v = Pb[:, :, i]
        known = ~v[:, max(cols) + 1:].any(axis=1) if max(cols) + 1 < d else np.ones(len(v), bool)
        ok &= ~((((bapply(Pb, v) - v) % p).any(axis=1)) & known)
    for i in cols:
        for j in cols:
            x, y = Pb[:, :, i], Pb[:, :, j]; h = bmul(x, y)
            within = ~h[:, max(cols) + 1:].any(axis=1) if max(cols) + 1 < d else np.ones(len(h), bool)
            bad = ((bapply(Pb, h) - h) % p).any(axis=1) & within
            ok &= ~bad
            q = bmul((I[:, i] - x) % p, (I[:, j] - y) % p)
            within = ~q[:, max(cols) + 1:].any(axis=1) if max(cols) + 1 < d else np.ones(len(q), bool)
            ok &= ~(((bapply(Pb, q)) % p).any(axis=1) & within)
    return ok
allk = np.array(list(itertools.product(range(p), repeat=4)), dtype=np.int64) @ np.array(Kb) % p   # (p^4, d)
folds = []
t0 = time.time()
for ki, PK in enumerate(Kfolds):
    P = np.zeros((d, d), dtype=np.int64); PKfull = Kmat @ PK % p
    for i in range(3): P[:, i] = PKfull[:, i]
    Pdelta = PKfull[:, 3]
    def lift(Pcur, x_idx):
        T = delta(I[x_idx]); acc = np.zeros((d, d), dtype=np.int64)
        for a in range(d):
            for b in range(d):
                if T[a, b]: acc = (acc + T[a, b] * np.outer(Pcur[:, a], Pcur[:, b])) % p
        sol = solve(DC, acc.ravel())
        if sol is None: return None
        return sum((sol[k] * I[Cb[k]] for k in range(4)), np.zeros(d, dtype=np.int64)) % p
    La = lift(P, 3); Lb = lift(P, 5)
    if La is None or Lb is None: continue
    nk = len(allk); idx = np.array(list(itertools.product(range(nk), repeat=2)))
    for chunk in np.array_split(idx, max(1, len(idx) // 40000)):
        Pb = np.repeat(P[None], len(chunk), axis=0)
        Pb[:, :, 3] = (La + allk[chunk[:, 0]]) % p; Pb[:, :, 4] = (Pb[:, :, 3] - Pdelta) % p
        Pb[:, :, 5] = (Lb + allk[chunk[:, 1]]) % p
        good = Pb[stage_ok(Pb, list(range(6)))]
        for Q1 in good:
            L6 = lift(Q1, 6); L7 = lift(Q1, 7)
            if L6 is None or L7 is None: continue
            Qb = np.repeat(Q1[None], len(idx), axis=0)
            Qb[:, :, 6] = (L6 + allk[idx[:, 0]]) % p; Qb[:, :, 7] = (L7 + allk[idx[:, 1]]) % p
            okb = stage_ok(Qb, list(range(8)))
            # Δ-compatibility is built in by the lifts; idempotency and both subalgebra conditions are checked in full above
            for Q2 in Qb[okb]: folds.append(Q2.copy())
    if ki % 8 == 0: print(f"  K-fold {ki}/{len(Kfolds)} folds {len(folds)} t={time.time()-t0:.0f}s", flush=True)
nonfilt = sum(1 for f in folds if (f[6:, :6] % p).any())
print(f"folds whose P raises path length on A_(≤1) (e.g. P(β) ∋ δβ): {nonfilt}")
exec("G = []" + open('q123.py').read().split("G = []")[1])
