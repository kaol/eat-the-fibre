"""v2.66: the collapse retraction kQ ⇄ kQ^s, checked exactly over ℤ.

Q^s is Q with every parallel class replaced by one arrow (QUIVERS[QS] must be that, arrows listed in class order).
π: kQ → kQ^s sends every arrow of a class to the class's arrow; ι: kQ^s → kQ sends it to the first arrow of the class.
Checked entry by entry over ℤ (structure constants M, coproduct DL of pcensus.build):
  π and ι are algebra maps, Δ^s π = (π⊗π) Δ, Δ ι = (ι⊗ι) Δ^s, π ι = id.
Then P^s ↦ ι P^s π is an injective map Folds(kQ^s) → Folds(kQ) over every commutative ring, preserving the vertex block
(proof in docs v2.66 §3). With REPS (census --save-reps of QS, our own output) it also maps every G-orbit representative
and checks that the image is a fold of Q over 𝔽_p.
Usage: python3 analysis/v2_66_retraction.py Q QS [REPS]"""
import sys, os, pickle
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pcensus as pc
Q, QS = sys.argv[1], sys.argv[2]
nv, A = pc.QUIVERS[Q]; nvs, As = pc.QUIVERS[QS]; assert nv == nvs
B, M, D = pc.build(nv, A); Bs, Ms, Ds = pc.build(nvs, As)
d, ds = len(B), len(Bs)
cls = {st: i for i, st in enumerate(As)}; assert len(cls) == len(As) and set(cls) == set(A)
first = {}
for i, st in enumerate(A): first.setdefault(st, i)
col = {b: i for i, b in enumerate(B)}; cols = {b: i for i, b in enumerate(Bs)}
def pmap(b): return b if b[0] == 'v' else ('p', tuple(cls[A[a]] for a in b[1]))
def imap(b): return b if b[0] == 'v' else ('p', tuple(first[As[a]] for a in b[1]))
Pi = np.zeros((ds, d), dtype=np.int64); Io = np.zeros((d, ds), dtype=np.int64)
for i, b in enumerate(B): Pi[cols[pmap(b)], i] = 1
for i, b in enumerate(Bs): Io[col[imap(b)], i] = 1
hom = lambda F, Ma, Mb: np.array_equal(np.einsum('ijk,lk->ijl', Ma, F), np.einsum('ai,bj,abl->ijl', F, F, Mb))
# Δ_b(F x) = (F⊗F) Δ_a(x): for a basis vector e_m of the source, Σ_l F[l,m] Db[l] = F Da[m] Fᵀ
def coalg(F, Da, Db): return all(np.array_equal(np.einsum('l,lab->ab', F[:, m], Db), F @ Da[m] @ F.T) for m in range(F.shape[1]))
res = dict(pi_hom=hom(Pi, M, Ms), iota_hom=hom(Io, Ms, M), pi_coalg=coalg(Pi, D, Ds), iota_coalg=coalg(Io, Ds, D),
           pi_iota_id=np.array_equal(Pi @ Io, np.eye(ds, dtype=np.int64)))
print(f"{QS} ⇄ {Q}: {res}")
if len(sys.argv) > 3:
    R = pickle.load(open(sys.argv[3], 'rb')); p = R['p']; assert R['basis'] == Bs
    I = np.eye(d, dtype=np.int64); bad = 0; vid = [i for i, b in enumerate(B) if b[0] == 'v']; vids = [i for i, b in enumerate(Bs) if b[0] == 'v']
    vb_ok = True
    for Ps in R['greps']:
        P = Io @ np.asarray(Ps, dtype=np.int64) @ Pi % p
        ok = ((P @ P - P) % p == 0).all()
        h = np.einsum('ia,jb,ijk->abk', P, P, M) % p; qq = np.einsum('ia,jb,ijk->abk', I - P, I - P, M) % p
        ok &= ((np.einsum('xk,abk->abx', P, h) - h) % p == 0).all() and (np.einsum('xk,abk->abx', P, qq) % p == 0).all()
        ok &= ((np.einsum('ml,mab->lab', P, D) - np.einsum('ai,lij,bj->lab', P, D, P)) % p == 0).all()
        bad += not ok
        vb_ok &= np.array_equal(P[np.ix_(vid, vid)] % p, np.asarray(Ps, dtype=np.int64)[np.ix_(vids, vids)] % p)
    print(f"q = {p}: {len(R['greps'])} G-orbit representatives of {QS} lifted; non-folds of {Q}: {bad}; vertex blocks kept: {vb_ok}")
