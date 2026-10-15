"""Stabilizer types over Q versus F_3 orbit sizes, for 1⇉2→3.
G ≅ Aff(1) acting on δ (gδ = (a−b)δ). Lie(G) = {D_(a,b): α1 ↦ aδ, α2 ↦ bδ}; D_(a,b)(δ) = (a−b)δ.
Stabilizer Lie algebra at P = {(a,b): [D_(a,b), P] = 0}; a = b is the unipotent (translation) line, a ≠ b a torus.
Prediction of N(q) = 21q² + 6q + 12: F_3 orbit size 1 ↔ stab dim 2; 2 = q−1 ↔ unipotent line; 3 = q ↔ torus line; 6 ↔ stab 0."""
import pickle, collections
exec(open('tangent123.py').read().split("lifted = 0")[0])
from sage.all__sagemath_singular import QQ, matrix
data = pickle.load(open('folds_p3.pkl', 'rb'))
G3 = []
import itertools
for vals in itertools.product(range(3), repeat=4):
    g = np.array(vals, dtype=np.int64).reshape(2, 2)
    if (g.sum(axis=0) % 3 != 1).any() or round(np.linalg.det(g)) % 3 == 0: continue
    Ah = np.eye(d, dtype=np.int64); Ah[3:5, 3:5] = g; Ah[6:8, 6:8] = g; G3.append(Ah)
def inv3(A):
    B = next(B for B in (np.array(v, dtype=np.int64).reshape(2, 2) for v in itertools.product(range(3), repeat=4)) if ((A[3:5, 3:5] @ B) % 3 == np.eye(2)).all())
    Ah = np.eye(d, dtype=np.int64); Ah[3:5, 3:5] = B; Ah[6:8, 6:8] = B; return Ah
Gi = [inv3(g) for g in G3]; assert all(((g @ gi) % 3 == np.eye(d)).all() for g, gi in zip(G3, Gi))
D0, D1 = Ds   # D_(1,0), D_(0,1)
table = collections.Counter(); unlifted = collections.Counter()
for P3 in data['reps']:
    orb = {((g @ P3 @ gi) % 3).tobytes(): (g @ P3 @ gi) % 3 for g, gi in zip(G3, Gi)}
    kind = None
    for f in orb.values():
        P = matrix(QQ, [[(-1 if x == 2 else x) for x in row] for row in f.tolist()])
        if not is_fold(P): continue
        A = matrix(QQ, [(D0 * P - P * D0).list(), (D1 * P - P * D1).list()]).transpose()
        K = A.right_kernel()
        if K.dimension() == 1:
            a, b = K.basis()[0]; kind = 'unipotent' if a == b else 'torus'
        else: kind = f'stab dim {K.dimension()}'
        break
    if kind is None: unlifted[len(orb)] += 1
    else: table[(len(orb), kind)] += 1
print("lifted orbits (F_3 orbit size, Q stabilizer type):", sorted(table.items()))
print("orbits with no member lifting to Q, by F_3 size:", dict(unlifted))
