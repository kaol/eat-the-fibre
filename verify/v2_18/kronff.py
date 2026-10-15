"""Exhaustive census of Δ-compatible folds of the generalized Kronecker algebra K_m over F_p, and of their orbits under
G = {g ∈ GL_m(F_p) : σ∘g = σ} (σ = sum of arrow coefficients), acting on the arrow span and fixing e1, e2.
Basis: e1, e2, α_1..α_m.  e_i e_i = e_i, e1 α = α = α e2;  Δ(α_i) = e1⊗e2, Δ(e_i) = 0.
Prediction (hand classification): 14m − 8 orbits."""
import sys, itertools
import numpy as np
m, p = int(sys.argv[1]), int(sys.argv[2]); d = m + 2
def mul(x, y):
    z = np.zeros(d, dtype=np.int64)
    z[0] = x[0] * y[0]; z[1] = x[1] * y[1]
    z[2:] = x[0] * y[2:] + x[2:] * y[1]
    return z % p
def rref(M):
    M = M.copy() % p; r = 0; piv = []
    for c in range(M.shape[1]):
        pr = next((i for i in range(r, M.shape[0]) if M[i, c]), None)
        if pr is None: continue
        M[[r, pr]] = M[[pr, r]]; M[r] = M[r] * pow(int(M[r, c]), p - 2, p) % p
        for i in range(M.shape[0]):
            if i != r and M[i, c]: M[i] = (M[i] - M[i, c] * M[r]) % p
        piv.append(c); r += 1
        if r == M.shape[0]: break
    return M[:r], piv
def subspaces(k):
    for piv in itertools.combinations(range(d), k):
        free = [(i, j) for i in range(k) for j in range(d) if j not in piv and j > piv[i]]
        for vals in itertools.product(range(p), repeat=len(free)):
            B = np.zeros((k, d), dtype=np.int64)
            for i, c in enumerate(piv): B[i, c] = 1
            for (i, j), v in zip(free, vals): B[i, j] = v
            yield B
def in_span(B, v):
    if B.shape[0] == 0: return not v.any()
    R, _ = rref(np.vstack([B, v])); return R.shape[0] == B.shape[0]
def subalg(B):
    return all(in_span(B, mul(x, y)) for x in B for y in B)
e1 = np.eye(d, dtype=np.int64)[0]; e2 = np.eye(d, dtype=np.int64)[1]
sig = lambda x: int(x[2:].sum() % p)
Hs, Qs = [], []
for k in range(d + 1):
    for B in subspaces(k):
        if not subalg(B): continue
        sigzero = all(sig(x) == 0 for x in B)
        if sigzero or (in_span(B, e1) and in_span(B, e2)): Hs.append(B)
        if sigzero or in_span(B, e1) or in_span(B, e2): Qs.append(B)
folds = []
for H in Hs:
    for Q in Qs:
        if H.shape[0] + Q.shape[0] != d: continue
        M = np.vstack([H, Q]) if H.shape[0] and Q.shape[0] else (H if Q.shape[0] == 0 else Q)
        if rref(M)[0].shape[0] != d: continue
        folds.append((H, Q))
# orbits under G
G = []
for vals in itertools.product(range(p), repeat=m * m):
    g = np.array(vals, dtype=np.int64).reshape(m, m)
    if any(int(g[:, j].sum() % p) != 1 for j in range(m)): continue     # σ(g v) = σ(v): column sums 1
    if rref(g)[0].shape[0] != m: continue
    G.append(g)
def act(g, B):
    if B.shape[0] == 0: return B
    C = B.copy(); C[:, 2:] = (B[:, 2:] @ g.T) % p; return C
key = lambda H, Q: (rref(H)[0].tobytes() if H.shape[0] else b'', rref(Q)[0].tobytes() if Q.shape[0] else b'')
seen = set(); orbits = 0
for H, Q in folds:
    k0 = key(H, Q)
    if k0 in seen: continue
    orbits += 1
    for g in G: seen.add(key(act(g, H), act(g, Q)))
print(f"K_{m} over F_{p}: {len(folds)} folds, |G| = {len(G)}, orbits = {orbits}, predicted 14m−8 = {14 * m - 8}")
