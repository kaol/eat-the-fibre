"""v2.60 (Program 3): how the folds of Q⁺ = Q + a leaf relate to the folds of Q, on census data.

Q⁺ = pairs_in_short_tail (Q₁) = Q₀ + the sink leaf ℓ: 2→4. kQ⁺ = kQ ⊕ J with J = the paths through the new vertex
(e₄ and every path ending in ℓ): J is a two-sided ideal, kQ a subalgebra, and kQ⁺/J ≅ kQ.
For every G-orbit representative P of Folds(𝔽_q Q⁺) (census --save-reps dump) this records
- J: P(J) ⊆ J (then P induces P̄ on kQ⁺/J = kQ, and P̄ is a fold: image and kernel map to subalgebras, Δ descends);
- Q: P(kQ) ⊆ kQ (then P restricts to kQ);
- C: the compression P_Q = π P ι (π: kQ⁺ → kQ killing J) is a fold of kQ (checked directly over 𝔽_q);
- the census class (local-torus orbit) of P, and for C-folds the census class of P_Q among Q's folds.
Aggregated per local-torus orbit of Q⁺ (a class is counted as J/Q/C if all its G-orbits are).
Usage: python3 analysis/v2_60_leaf_census.py q"""
import os, sys, pickle, collections, io
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
import numpy as np
import pcensus as pc
q = int(sys.argv[1])
class SafeU(pickle.Unpickler):
    OK = {('numpy._core.numeric', '_frombuffer'), ('numpy.core.numeric', '_frombuffer'), ('numpy', 'dtype'),
          ('numpy._core.multiarray', '_reconstruct'), ('numpy.core.multiarray', '_reconstruct'), ('numpy', 'ndarray')}
    def find_class(self, mod, name):
        if (mod, name) not in self.OK: raise pickle.UnpicklingError(f"refused global {mod}.{name}")
        return super().find_class(mod, name)
load = lambda fn: SafeU(open(fn, 'rb')).load()
R1 = load(os.path.join(os.path.dirname(HERE), 'logs', f'reps_pairs_in_short_tail_{q}.pkl'))
R0 = load(os.path.join(os.path.dirname(HERE), 'logs', f'reps_pairs_in_short_{q}.pkl'))
c1 = pc.setup('pairs_in_short_tail', q); c0 = pc.setup('pairs_in_short', q)
b1, b0 = c1['basis'], c0['basis']; d1, d0 = len(b1), len(b0)
Jx = [i for i, b in enumerate(b1) if (b[0] == 'v' and b[1] == 3) or (b[0] == 'p' and 5 in b[1])]
Qx = [b1.index(b) for b in b0]                       # Q's basis inside Q⁺'s, in Q's order
assert sorted(Jx + Qx) == list(range(d1))
M0 = np.asarray(c0['M']); DL0 = np.asarray(c0['DL'])
def is_fold(P, M, DL, p):
    d = P.shape[0]; P = P % p
    if ((P @ P - P) % p).any(): return False
    Q = (np.eye(d, dtype=np.int64) - P) % p
    for A, check in ((P, lambda u: (P @ u - u) % p), (Q, lambda u: (P @ u) % p)):
        prod = np.einsum('ia,jb,ijk->abk', A, A, M) % p          # products of columns of A
        if check(prod.reshape(-1, d).T).any(): return False
    lhs = np.einsum('ml,mab->lab', P, DL) % p                   # Σ_m P[m,l] DL[m]
    rhs = np.einsum('ax,lxy,by->lab', P, DL, P) % p             # P DL[l] Pᵀ
    return not ((lhs - rhs) % p).any()
key0 = {}                                            # Q's folds → census class, via its G-orbit representatives
# Q's census classes: the reps are G-orbit representatives; to label an arbitrary fold of Q we need its G-orbit, so
# enumerate each representative's G-orbit (G small) and map every element to the class label.
G0 = pc.build_G(c0)
for P, lab in zip(R0['greps'], R0['loc']):
    P = P.astype(np.int64)
    for g, gi in G0: key0[((g @ P @ gi) % q).astype(np.int8).tobytes()] = lab
stats = collections.Counter(); cls = collections.defaultdict(lambda: collections.Counter()); img = collections.defaultdict(set)
for P, lab in zip(R1['greps'], R1['loc']):
    P = P.astype(np.int64) % q
    J = not P[np.ix_(Qx, Jx)].any()                  # P(J) ⊆ J: no Q-component in P of J-vectors
    Qs = not P[np.ix_(Jx, Qx)].any()                 # P(kQ) ⊆ kQ
    PQ = P[np.ix_(Qx, Qx)]
    C = is_fold(PQ, M0, DL0, q)
    k0 = key0.get(PQ.astype(np.int8).tobytes()) if C else None
    assert (not C) or k0 is not None, "compression is a fold but not found among Q's census folds"
    cls[lab][(J, Qs, C)] += 1
    if C: img[lab].add(k0)
summary = collections.Counter()
for lab, c in cls.items():
    allJ = all(k[0] for k in c); allQ = all(k[1] for k in c); allC = all(k[2] for k in c); anyC = any(k[2] for k in c)
    summary[(allJ, allQ, allC, anyC)] += 1
n0 = len(set(R0['loc'])); n1 = len(cls)
print(f"q = {q}: Q⁺ census classes {n1}, Q census classes {n0}")
print("classes by (all P(J)⊆J, all P(kQ)⊆kQ, all compressions folds, some compression fold):")
for k, v in sorted(summary.items(), key=lambda kv: -kv[1]): print(f"  {k}: {v}")
# J-preserving classes: how many lie over each Q class (via the induced fold P̄ = compression when P(J) ⊆ J)
over = collections.defaultdict(set)
for lab, c in cls.items():
    if all(k[0] for k in c) and len(img[lab]) == 1: over[next(iter(img[lab]))].add(lab)
multi = sum(1 for lab in cls if len(img[lab]) > 1)
print(f"J-preserving classes with a single base class: {sum(len(v) for v in over.values())}; classes whose G-orbits map to several base classes: {multi}")
fib = collections.Counter(len(v) for v in over.values())
print(f"Q classes hit: {len(over)} of {n0}; fibre sizes (number of Q⁺ classes over one Q class): {dict(sorted(fib.items()))}")
