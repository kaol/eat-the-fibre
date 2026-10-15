"""Invariants of the 71 G-orbits of folds of k(1⇉2→3) over F_3: rank of P (dim H), the vertex fold P|_{A0} pattern
(images of e1, e2, e3 under P0 = projection of P to A0 along arrows), F_3 orbit size, and whether P is length-filtered."""
import pickle, collections, numpy as np
D = pickle.load(open('folds_p3.pkl', 'rb')); p = 3
def rk(A):
    A = A.copy() % p; r = 0
    for c in range(A.shape[1]):
        pr = next((i for i in range(r, A.shape[0]) if A[i, c]), None)
        if pr is None: continue
        A[[r, pr]] = A[[pr, r]]; A[r] = A[r] * pow(int(A[r, c]), p - 2, p) % p
        for i in range(A.shape[0]):
            if i != r and A[i, c]: A[i] = (A[i] - A[i, c] * A[r]) % p
        r += 1
    return r
exec(open('stabcheck.py').read().split("D0, D1 = Ds")[0].split("exec(open('tangent123.py')")[0])  # nothing heavy
import itertools
blocks = [np.array(v).reshape(2, 2) for v in itertools.product(range(3), repeat=4)]
Gb = [B for B in blocks if (B.sum(0) % 3 == 1).all() and round(np.linalg.det(B)) % 3]
def L8(B): A = np.eye(8, dtype=np.int64); A[3:5, 3:5] = B; A[6:8, 6:8] = B; return A
inv = lambda B: next(C for C in blocks if ((B @ C) % 3 == np.eye(2)).all())
tab = collections.Counter()
for P in D['reps']:
    size = len({((L8(B) @ P @ L8(inv(B))) % 3).tobytes() for B in Gb})
    P0 = P[:3, :3] % 3                       # vertex part: coefficients of e1,e2,e3 in P(e_j)
    P0s = tuple(''.join('0+-'[x] for x in P0[:, j]) for j in range(3))
    tab[(rk(P), P0s, size, bool((P[6:, :6] % 3).any()))] += 1
print("(dim H, P0 columns e1|e2|e3 as signs over (e1,e2,e3), F_3 orbit size, non-filtered): orbits")
for k, v in sorted(tab.items()): print(" ", k, v)
print("vertex folds P0 occurring:", len({k[1] for k in tab}))
