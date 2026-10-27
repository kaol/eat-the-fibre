"""v2.66: path reversal is a duality of fold varieties, checked exactly over ℤ.

For quivers Q, Q' with Q' ≅ Q^op: find a vertex bijection σ and arrow bijection α with α(a): σ(t(a)) → σ(s(a)), and the
basis bijection r(e_i) = e_σ(i), r(a₁…a_k) = α(a_k)…α(a₁). Check over ℤ, entry by entry:
  M'[r j, r i, r k] = M[i, j, k]          (r is an anti-isomorphism kQ → kQ'),
  DL'[r m, r b, r a] = DL[m, a, b]        (Δ' ∘ r = τ ∘ (r⊗r) ∘ Δ).
Then P ↦ rPr⁻¹ maps folds to folds over every ring (image/kernel subalgebras go to subalgebras, idempotency and the
Δ-identity are conjugated, τ commutes with P⊗P), and conjugates the census group and weightings of Q onto those of Q'.
With REPS (a census --save-reps pickle of Q, our own output), also maps every G-orbit representative and checks that it
is a fold of Q' over 𝔽_p.
Usage: python3 analysis/v2_66_duality.py Q Q' [REPS]"""
import sys, os, itertools, pickle
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pcensus as pc
Q, Qp = sys.argv[1], sys.argv[2]
nv, A = pc.QUIVERS[Q]; nv2, A2 = pc.QUIVERS[Qp]
assert nv == nv2 and len(A) == len(A2)
B1, M1, D1 = pc.build(nv, A); B2, M2, D2 = pc.build(nv2, A2)
d = len(B1); assert d == len(B2)
idx2 = {b: i for i, b in enumerate(B2)}
found = None
for sig in itertools.permutations(range(nv)):
    rev = [(sig[t], sig[s]) for s, t in A]                    # reversed, relabelled arrows of Q
    if sorted(rev) != sorted(A2): continue
    # arrow bijection: within each (s,t) slot any matching works up to the census group; take the order-preserving one
    used = set(); alpha = []
    for st in rev:
        j = next(j for j, st2 in enumerate(A2) if st2 == st and j not in used); used.add(j); alpha.append(j)
    r = []
    for b in B1:
        rb = ('v', sig[b[1]]) if b[0] == 'v' else ('p', tuple(alpha[a] for a in reversed(b[1])))
        if rb not in idx2: break
        r.append(idx2[rb])
    else:
        found = (sig, alpha, r); break
assert found, "no reversal isomorphism"
sig, alpha, r = found; r = np.array(r)
assert sorted(r) == list(range(d))
okM = all(M2[r[j], r[i], r[k]] == M1[i, j, k] for i in range(d) for j in range(d) for k in range(d))
okD = all(D2[r[m], r[b], r[a]] == D1[m, a, b] for m in range(d) for a in range(d) for b in range(d))
print(f"{Q} → {Qp}: σ = {sig}, α = {alpha}; anti-isomorphism {okM}; Δ' r = τ (r⊗r) Δ {okD}")
if len(sys.argv) > 3:
    R = pickle.load(open(sys.argv[3], 'rb')); p = R['p']
    assert R['basis'] == B1
    Rm = np.zeros((d, d), dtype=np.int64); Rm[r, np.arange(d)] = 1     # r as a permutation matrix
    ctx2 = pc.setup(Qp, p); bad = 0
    for P in R['greps']:
        P2 = (Rm @ np.asarray(P, dtype=np.int64) @ Rm.T) % p
        Pq = P2 % p
        ok = ((Pq @ Pq - Pq) % p == 0).all()
        h = np.einsum('ia,jb,ijk->abk', Pq, Pq, M2) % p; I = np.eye(d, dtype=np.int64)
        qq = np.einsum('ia,jb,ijk->abk', I - Pq, I - Pq, M2) % p
        ok &= ((np.einsum('xk,abk->abx', Pq, h) - h) % p == 0).all() and (np.einsum('xk,abk->abx', Pq, qq) % p == 0).all()
        ok &= ((np.einsum('ml,mab->lab', Pq, D2) - np.einsum('ai,lij,bj->lab', Pq, D2, Pq)) % p == 0).all()
        bad += not ok
    print(f"q = {p}: {len(R['greps'])} G-orbit representatives of {Q} mapped; non-folds of {Qp}: {bad}; "
          f"local-torus classes of {Q}: {len(set(R['loc']))}")
