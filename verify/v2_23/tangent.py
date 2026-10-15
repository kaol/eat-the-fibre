"""Exact (mod a large prime) Zariski tangent dimension of the fold variety at the transversal Q-representatives, compared with
r_G and r_H. The fold equations: P² = P; P(x·y) = x·y for x, y columns of P; P((1−P)_i·(1−P)_j) = 0; ΔP = (P⊗P)Δ.
Their derivative in direction V is computed exactly by the product rule (all equations are polynomial in P)."""
import sys, pickle, collections
import numpy as np, sympy as sp
name, p = sys.argv[1], int(sys.argv[2])
sys.argv = ['x']
exec(open('tfam.py').read().split("shown = 0; seen_kind")[0].replace("name, p, nshow = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 3\nsys.argv = ['x']\n", ""))
R = 1000003
def mulb(X, Y): return np.einsum('ai,aj,ijk->ak', X, Y, M)          # columnwise products: X[:,i]*Y[:,j] for paired columns
def jac(P):
    I = np.eye(d, dtype=np.int64); rows = []
    for idx in range(d * d):
        V = np.zeros((d, d), dtype=np.int64); V.flat[idx] = 1
        e1 = (P @ V + V @ P - V).ravel()
        e2 = []; e3 = []
        for i in range(d):
            for j in range(d):
                h = np.einsum('i,j,ijk->k', P[:, i], P[:, j], M)
                dh = np.einsum('i,j,ijk->k', V[:, i], P[:, j], M) + np.einsum('i,j,ijk->k', P[:, i], V[:, j], M)
                e2.append(V @ h + P @ dh - dh)
                q = np.einsum('i,j,ijk->k', (I - P)[:, i], (I - P)[:, j], M)
                dq = -np.einsum('i,j,ijk->k', V[:, i], (I - P)[:, j], M) - np.einsum('i,j,ijk->k', (I - P)[:, i], V[:, j], M)
                e3.append(V @ q + P @ dq)
        e4 = []
        for l in range(d):
            e4.append((np.einsum('m,mab->ab', V[:, l], DL) - V @ DL[l] @ P.T - P @ DL[l] @ V.T).ravel())
        rows.append(np.concatenate([e1, np.concatenate(e2), np.concatenate(e3), np.concatenate(e4)]))
    J = np.array(rows).T % R
    # rank mod R
    A = J.copy(); r = 0
    for c in range(A.shape[1]):
        pr = next((i for i in range(r, A.shape[0]) if A[i, c]), None)
        if pr is None: continue
        A[[r, pr]] = A[[pr, r]]; A[r] = A[r] * pow(int(A[r, c]), R - 2, R) % R
        nz = np.nonzero(A[:, c])[0]
        for i in nz:
            if i != r: A[i] = (A[i] - A[i, c] * A[r]) % R
        r += 1
    return d * d - r
tab = collections.Counter(); seen = set()
for P in folds:
    Pq = qfold(P)
    if Pq is None: continue
    rows = [(D @ Pq - Pq @ D).ravel().tolist() for D in gens]
    rG = sp.Matrix(rows).rank(); rH = sp.Matrix(rows + [(E @ Pq - Pq @ E).ravel().tolist()]).rank()
    if rH <= rG: continue
    k = Pq.tobytes()
    if k in seen: continue
    seen.add(k); tab[(rG, rH, jac(Pq))] += 1
print(f"{name}: transversal lifted folds by (r_G, r_H, tangent dim): {sorted(tab.items())}", flush=True)
