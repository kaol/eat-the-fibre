"""Census of Δ-compatible folds of k(1⇉2→3) over F_p by kernel fibration, and their orbits under
G = {g ∈ GL_2(F_p): column sums 1} acting on (α1, α2) (and on α1β, α2β), fixing e1, e2, e3, β.
Basis: 0 e1, 1 e2, 2 e3, 3 α1, 4 α2, 5 β, 6 α1β, 7 α2β."""
import sys, itertools
import numpy as np
p = int(sys.argv[1]); d = 8
# multiplication table: start/end vertices and concatenation
paths = {0: ('v', 0), 1: ('v', 1), 2: ('v', 2), 3: ('p', (3,)), 4: ('p', (4,)), 5: ('p', (5,)), 6: ('p', (3, 5)), 7: ('p', (4, 5))}
src = {0: 0, 1: 1, 2: 2, 3: 0, 4: 0, 5: 1, 6: 0, 7: 0}; tgt = {0: 0, 1: 1, 2: 2, 3: 1, 4: 1, 5: 2, 6: 2, 7: 2}
seq_idx = {(3,): 3, (4,): 4, (5,): 5, (3, 5): 6, (4, 5): 7}
M = np.zeros((d, d, d), dtype=np.int64)
for i in range(d):
    for j in range(d):
        if tgt[i] != src[j]: continue
        if paths[i][0] == 'v': M[i, j, j] = 1
        elif paths[j][0] == 'v': M[i, j, i] = 1
        else:
            s = paths[i][1] + paths[j][1]
            if s in seq_idx: M[i, j, seq_idx[s]] = 1
DL = np.zeros((d, d, d), dtype=np.int64)      # Δ(basis l) = Σ DL[l,a,b] a⊗b  (arrow deletion)
vert = lambda v: v
for l, (kind, s) in paths.items():
    if kind == 'v': continue
    for i, a in enumerate(s):
        pre = seq_idx[s[:i]] if i > 0 else src[l] if i == 0 and len(s[:i]) == 0 else None
        pre = seq_idx[s[:i]] if i > 0 else [0, 1, 2][[0, 0, 0, 0, 0, 1][a] if a < 5 else 1]
        a_src = {3: 0, 4: 0, 5: 1}[a]; a_tgt = {3: 1, 4: 1, 5: 2}[a]
        pre = seq_idx[s[:i]] if i > 0 else a_src
        suf = seq_idx[s[i + 1:]] if i + 1 < len(s) else a_tgt
        DL[l, pre, suf] += 1
def mul(x, y): return np.einsum('i,j,ijk->k', x, y, M) % p
def delta(x): return np.einsum('l,lab->ab', x, DL) % p
def rank(Mx):
    A = Mx.copy() % p; r = 0
    for c in range(A.shape[1]):
        pr = next((i for i in range(r, A.shape[0]) if A[i, c]), None)
        if pr is None: continue
        A[[r, pr]] = A[[pr, r]]; A[r] = A[r] * pow(int(A[r, c]), p - 2, p) % p
        for i in range(A.shape[0]):
            if i != r and A[i, c]: A[i] = (A[i] - A[i, c] * A[r]) % p
        r += 1
    return r
def solve(A, b):
    """one solution x of A x = b over F_p, or None"""
    A = A % p; b = b % p; m_, n_ = A.shape
    Ab = np.concatenate([A, b.reshape(-1, 1)], axis=1); r = 0; piv = []
    for c in range(n_):
        pr = next((i for i in range(r, m_) if Ab[i, c]), None)
        if pr is None: continue
        Ab[[r, pr]] = Ab[[pr, r]]; Ab[r] = Ab[r] * pow(int(Ab[r, c]), p - 2, p) % p
        for i in range(m_):
            if i != r and Ab[i, c]: Ab[i] = (Ab[i] - Ab[i, c] * Ab[r]) % p
        piv.append(c); r += 1
    if any(Ab[i, n_] for i in range(r, m_)): return None
    x = np.zeros(n_, dtype=np.int64)
    for i, c in enumerate(piv): x[c] = Ab[i, n_]
    return x
I = np.eye(d, dtype=np.int64)
Kb = [I[0], I[1], I[2], (I[3] - I[4]) % p]               # kernel basis
Cb = [3, 5, 6, 7]                                         # complement basis indices
DC = np.array([delta(I[c]).ravel() for c in Cb]).T       # Δ on the complement (injective)
# 1. algebra folds of K (4-dim, coordinates in Kb)
Kmat = np.array(Kb).T                                     # d×4
def kvec(c): return Kmat @ np.array(c) % p
def in_span(B, v):
    if len(B) == 0: return not (v % p).any()
    return rank(np.array(B)) == rank(np.array(list(B) + [v]))
def subspaces(n, k):
    for piv in itertools.combinations(range(n), k):
        free = [(i, j) for i in range(k) for j in range(n) if j not in piv and j > piv[i]]
        for vals in itertools.product(range(p), repeat=len(free)):
            B = np.zeros((k, n), dtype=np.int64)
            for i, c in enumerate(piv): B[i, c] = 1
            for (i, j), v in zip(free, vals): B[i, j] = v
            yield B
subalgs = []
for k in range(5):
    for B in subspaces(4, k):
        vecs = [kvec(r) for r in B]
        if all(in_span(vecs, mul(x, y)) for x in vecs for y in vecs): subalgs.append(B)
Kfolds = []
for H in subalgs:
    for Q in subalgs:
        if H.shape[0] + Q.shape[0] != 4: continue
        Bm = np.vstack([H, Q]) if H.shape[0] and Q.shape[0] else (H if Q.shape[0] == 0 else Q)
        if rank(Bm) != 4: continue
        # projection onto H along Q in K-coordinates
        Bt = Bm.T % p; inv = np.round(np.linalg.inv(Bt.astype(float)) * round(np.linalg.det(Bt.astype(float)))).astype(np.int64)
        det = int(round(np.linalg.det(Bt.astype(float)))) % p
        inv = inv * pow(det, p - 2, p) % p
        E = np.diag([1] * H.shape[0] + [0] * Q.shape[0])
        PK = Bt @ E @ inv % p
        Kfolds.append(PK)
print(f"algebra folds of K over F_{p}: {len(Kfolds)}", flush=True)
def stage1_ok(P):
    low = list(range(6))
    for i in low:
        v = P[:, i]
        if v[6:].any(): return False
        if ((P @ v - v) % p).any(): return False
    for i in low:
        for j in low:
            x, y = P[:, i], P[:, j]; h = mul(x, y)
            if not h[6:].any() and ((P @ h - h) % p).any(): return False
            q = mul((I[:, i] - x) % p, (I[:, j] - y) % p)
            if not q[6:].any() and ((P @ q) % p).any(): return False
    return True
def full_check(P):
    if ((P @ P - P) % p).any(): return False
    for i in range(d):
        for j in range(d):
            x, y = P[:, i], P[:, j]; h = mul(x, y)
            if ((P @ h - h) % p).any(): return False
            q = mul((I[:, i] - x) % p, (I[:, j] - y) % p)
            if ((P @ q) % p).any(): return False
    for j in range(d):
        lhs = delta(P[:, j]); rhs = np.einsum('ab,xa,yb->xy', np.reshape(DL[j], (d, d)), P, P) % p
        if ((lhs - rhs) % p).any(): return False
    return True
folds = []
Kcols = np.array(Kb).T
import time
t0 = time.time()
LIMIT = int(sys.argv[2]) if len(sys.argv) > 2 else None
for ki, PK in enumerate(Kfolds[:LIMIT]):
    if ki % 16 == 0: print(f'  K-fold {ki}/{len(Kfolds)}  folds so far {len(folds)}  t={time.time()-t0:.0f}s', flush=True)
    P = np.zeros((d, d), dtype=np.int64)
    # P on K: columns for e1,e2,e3 and δ; express α2 = α1 − δ later
    PKfull = Kmat @ PK % p                                 # d×4: P(Kb_i)
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
    # stage 1: α1, β (α2 = α1 − δ)
    base1 = {}
    P1 = P.copy()
    La = lift(P1, 3); Lb = lift(P1, 5)
    if La is None or Lb is None: continue
    for ka in itertools.product(range(p), repeat=4):
        for kb in itertools.product(range(p), repeat=4):
            Q1 = P1.copy()
            Q1[:, 3] = (La + kvec(ka)) % p; Q1[:, 4] = (Q1[:, 3] - Pdelta) % p; Q1[:, 5] = (Lb + kvec(kb)) % p
            if not stage1_ok(Q1): continue
            L6 = lift(Q1, 6); L7 = lift(Q1, 7)
            if L6 is None or L7 is None: continue
            # stage 2: K-components of α1β, α2β
            for k6 in itertools.product(range(p), repeat=4):
                for k7 in itertools.product(range(p), repeat=4):
                    Q2 = Q1.copy(); Q2[:, 6] = (L6 + kvec(k6)) % p; Q2[:, 7] = (L7 + kvec(k7)) % p
                    if full_check(Q2): folds.append(Q2)
print(f"folds of k(1⇉2→3) over F_{p}: {len(folds)}", flush=True)
G = []
for vals in itertools.product(range(p), repeat=4):
    g = np.array(vals, dtype=np.int64).reshape(2, 2)
    if (g.sum(axis=0) % p != 1).any() or rank(g) != 2: continue
    Ah = np.eye(d, dtype=np.int64); Ah[3:5, 3:5] = g; Ah[6:8, 6:8] = g
    G.append(Ah)
def inv_mod(A):
    n = A.shape[0]; Ab = np.concatenate([A % p, np.eye(n, dtype=np.int64)], axis=1)
    for c in range(n):
        pr = next(i for i in range(c, n) if Ab[i, c]); Ab[[c, pr]] = Ab[[pr, c]]
        Ab[c] = Ab[c] * pow(int(Ab[c, c]), p - 2, p) % p
        for i in range(n):
            if i != c and Ab[i, c]: Ab[i] = (Ab[i] - Ab[i, c] * Ab[c]) % p
    return Ab[:, n:]
Ginv = [inv_mod(g) for g in G]
seen = set(); orbits = 0
for P in folds:
    if P.tobytes() in seen: continue
    orbits += 1
    for g, gi in zip(G, Ginv): seen.add(((g @ P @ gi) % p).tobytes())
print(f"1⇉2→3 over F_{p}: {len(folds)} folds, |G| = {len(G)}, orbits = {orbits}")
