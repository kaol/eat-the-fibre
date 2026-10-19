"""Families of non-generic orbits, located through fixed-point loci (v2.26).

H = G × T (Lie G + class generators E_c; all weightings act, as for 1⇉2⇉3). For every H-orbit over F_p, take a member P
that lifts to an exact Q-fold. Over Q, exactly:
  Lie S   = {x ∈ Lie H : [x, P] = 0}             (stabilizer; dim S = dim H − r)
  C_S     = {X : [s, X] = 0 for s ∈ Lie S}       (commutant: the linear span of the S-fixed locus)
  Lie N   = {x ∈ Lie H : [x, Lie S] ⊆ Lie S}     (normalizer of S in H)
  t_S     = dim of the tangent space of Folds ∩ C_S at P (fold Jacobian restricted to C_S, rank mod 1000003)
  n       = dim N·P = rank{[x, P] : x ∈ Lie N}
Near P, the orbits whose stabilizer contains S are the points of Folds^S, and H·P ∩ Folds^S = N·P (connected
stabilizers). If dim_P Folds^S > n, infinitely many H-orbits with stabilizer ⊇ S pass near P: a family.
t_S = n proves there is none at P. For t_S > n the local dimension of Folds^S is measured numerically (ε-scaling,
solutions constrained to C_S).
Usage: python3 analysis/fixedlocus.py QUIVER p [tests per type=3] [exact|numeric|all]   (needs logs/folds_QUIVER_p.pkl;
  the exact stage saves its candidates to logs/fixedlocus_cands_QUIVER_p.pkl for the numeric stage)"""
import os as _os, sys as _sys
if 'ROOT' not in globals():
    AN = _os.path.dirname(_os.path.abspath(__file__)); ROOT = _os.path.dirname(AN); _os.chdir(ROOT)
    __file__ = _os.path.join(ROOT, 'pcensus.py')
import sys, pickle, itertools, collections, numpy as np, sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy import QQ as SQQ
from scipy.optimize import least_squares
name, p = sys.argv[1], int(sys.argv[2]); per_type = int(sys.argv[3]) if len(sys.argv) > 3 else 3
stage = sys.argv[4] if len(sys.argv) > 4 else 'all'          # 'exact' | 'numeric' | 'all'
CANDS = f'logs/fixedlocus_cands_{name}_{p}' + ('_local' if _os.environ.get('FL_LOCAL') == '1' else '') + ('_aut' if _os.environ.get('FL_AUT') == '1' else '') + '.pkl'
sys.argv = ['x']
exec(open(_os.path.join(ROOT, 'pcensus.py')).read().split("if __name__ == '__main__':")[0])
ctx = setup(name, p); C.update(ctx); F = ctx['F']; basis, arrows, d, M, DL = ctx['basis'], ctx['arrows'], ctx['d'], ctx['M'], ctx['DL']
folds = [f.astype(np.int64) for f in pickle.load(open(f'logs/folds_{name}_{p}.pkl', 'rb'))['folds']]
# ---- H = G × T and its orbits over F_p
G = build_G(ctx)
cls = collections.defaultdict(list)
for a, st in enumerate(arrows): cls[st].append(a)
cl = list(cls.values())
T = []
for w in itertools.product(range(1, p), repeat=len(cl)):
    wa = {a: w[k] for k, c in enumerate(cl) for a in c}
    dg = [1 if b[0] == 'v' else int(np.prod([wa[x] for x in b[1]])) % p for b in basis]
    T.append((np.diag(dg), np.diag([pow(x, p - 2, p) for x in dg])))
H = [(t @ g % p, gi @ ti % p) for g, gi in G for t, ti in T]
lab, reps, sizes = orbit_label(folds, (np.array([h for h, _ in H]), np.array([hi for _, hi in H])), p)
members = collections.defaultdict(list)
for P in folds: members[lab[(P % p).astype(np.int8).tobytes()]].append(P)
# ---- Lie H over Q
gens, E = lie_G(ctx)
if _os.environ.get('FL_AUT') == '1':      # v2.28: Lie(Aut_Δ) (all Δ-derivations) instead of Lie G; larger where an
    import runpy as _rp                   # arrow runs parallel to a product of parallel differences
    _T3 = _rp.run_path(_os.path.join(AN if 'AN' in globals() else _os.path.join(ROOT, 'analysis'), 'theorem3.py'), run_name='t3')
    gens = [np.array(m, dtype=np.int64) for m in _T3['delta_derivations'](name)]
Ec = [np.diag([0 if b[0] == 'v' else sum(1 for x in b[1] if x in c) for b in basis]).astype(np.int64) for c in cl]
LH = [sp.Matrix(x.tolist()) for x in gens + Ec]
def qfold(P):
    Pq = np.where(P == p - 1, -1, P).astype(np.int64); I = np.eye(d, dtype=np.int64)
    if (Pq @ Pq != Pq).any(): return None
    for i in range(d):
        for j in range(d):
            h = np.einsum('i,j,ijk->k', Pq[:, i], Pq[:, j], M); q = np.einsum('i,j,ijk->k', (I - Pq)[:, i], (I - Pq)[:, j], M)
            if (Pq @ h != h).any() or (Pq @ q).any(): return None
    for l in range(d):
        if (np.einsum('m,mab->ab', Pq[:, l], DL) != Pq @ DL[l] @ Pq.T).any(): return None
    return Pq
def jac_matrix(P):
    """exact integer Jacobian of the fold equations at P (rows: equations, columns: the d² entries)"""
    I = np.eye(d, dtype=np.int64); cols = []
    H2 = np.einsum('ai,bj,abk->ijk', P, P, M); Q = I - P; Q2 = np.einsum('ai,bj,abk->ijk', Q, Q, M)
    for idx in range(d * d):
        V = np.zeros((d, d), dtype=np.int64); V.flat[idx] = 1
        e1 = P @ V + V @ P - V
        dH = np.einsum('ai,bj,abk->ijk', V, P, M) + np.einsum('ai,bj,abk->ijk', P, V, M)
        e2 = np.einsum('kl,ijl->ijk', V, H2) + np.einsum('kl,ijl->ijk', P, dH) - dH
        dQ = -np.einsum('ai,bj,abk->ijk', V, Q, M) - np.einsum('ai,bj,abk->ijk', Q, V, M)
        e3 = np.einsum('kl,ijl->ijk', V, Q2) + np.einsum('kl,ijl->ijk', P, dQ)
        e4 = np.einsum('ml,mab->lab', V, DL) - np.einsum('xa,lab,yb->lxy', V, DL, P) - np.einsum('xa,lab,yb->lxy', P, DL, V)
        cols.append(np.concatenate([e1.ravel(), e2.ravel(), e3.ravel(), e4.ravel()]))
    return np.array(cols).T
R = 1000003
def rank_mod(A):
    A = A.copy() % R; r = 0
    for c in range(A.shape[1]):
        pr = next((i for i in range(r, A.shape[0]) if A[i, c]), None)
        if pr is None: continue
        A[[r, pr]] = A[[pr, r]]; A[r] = A[r] * pow(int(A[r, c]), R - 2, R) % R
        nz = np.nonzero(A[:, c])[0]; nz = nz[nz != r]
        if len(nz): A[nz] = (A[nz] - np.outer(A[nz, c], A[r])) % R
        r += 1
    return r
def int_basis(vecs):
    out = []
    for v in vecs:
        den = sp.ilcm(*[sp.fraction(x)[1] for x in v]) if v else 1
        out.append([int(x * den) for x in v])
    return np.array(out, dtype=np.int64)
Mf, DLf, If = M.astype(float), DL.astype(float), np.eye(d)
def res(x):
    P = x.reshape(d, d); H2 = np.einsum('ai,bj,abk->ijk', P, P, Mf); Q = If - P; Q2 = np.einsum('ai,bj,abk->ijk', Q, Q, Mf)
    return np.concatenate([(P @ P - P).ravel(), (np.einsum('kl,ijl->ijk', P, H2) - H2).ravel(),
                           np.einsum('kl,ijl->ijk', P, Q2).ravel(),
                           (np.einsum('ml,mab->lab', P, DLf) - np.einsum('xa,lab,yb->lxy', P, DLf, P)).ravel()])
rng = np.random.default_rng(11)
import torch
Mt, DLt, It = torch.tensor(M, dtype=torch.float64), torch.tensor(DL, dtype=torch.float64), torch.eye(d, dtype=torch.float64)
def res_t(xt):
    P = xt.reshape(d, d); H2 = torch.einsum('ai,bj,abk->ijk', P, P, Mt); Q = It - P; Q2 = torch.einsum('ai,bj,abk->ijk', Q, Q, Mt)
    return torch.cat([(P @ P - P).ravel(), (torch.einsum('kl,ijl->ijk', P, H2) - H2).ravel(),
                      torch.einsum('kl,ijl->ijk', P, Q2).ravel(),
                      (torch.einsum('ml,mab->lab', P, DLt) - torch.einsum('xa,lab,yb->lxy', P, DLt, P)).ravel()])
def local_dim_in(P0, Bc, samples=24):
    """numerical local dimension at P0 of {X ∈ P0 + span(Bc) : X is a fold} (ε-scaling, as in v2.12); exact Jacobians"""
    x0 = torch.tensor(P0.astype(float).ravel()); Bt = torch.tensor(Bc.astype(float).T)       # (d², k)
    f = lambda y: res_t(x0 + Bt @ torch.tensor(y)).numpy()
    jf = lambda y: torch.autograd.functional.jacobian(lambda yy: res_t(x0 + Bt @ yy), torch.tensor(y), vectorize=True, strategy='forward-mode').numpy()
    J = jf(np.zeros(Bt.shape[1])); U, s_, Vt = np.linalg.svd(J); rk = int((s_ > 1e-7 * s_[0]).sum()); N = Vt[rk:].T
    if N.shape[1] == 0: return 0
    svs = []
    for eps in (1e-4, 1e-5):
        diffs = []
        for _ in range(samples):
            sol = least_squares(f, eps * (N @ rng.normal(size=N.shape[1])), jac=jf, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=200)
            if np.linalg.norm(f(sol.x)) < 1e-11: diffs.append(Bc.astype(float).T @ sol.x)
        sv = np.zeros(N.shape[1])
        if diffs:
            v = np.linalg.svd(np.array(diffs), compute_uv=False)[:N.shape[1]] / eps; sv[:len(v)] = v
        svs.append(sv)
    return int(sum((a > 0.2) and (b > 0.2) and (0.25 < b / a < 4) for a, b in zip(svs[0], svs[1])))
def dm_null(rows):
    """exact nullspace over Q of an integer matrix (list of rows), as a list of integer vectors"""
    A = DomainMatrix([[SQQ(int(x)) for x in r] for r in rows], (len(rows), len(rows[0])), SQQ)
    N = A.nullspace().to_Matrix()
    out = []
    for i in range(N.rows):
        v = list(N.row(i)); den = sp.ilcm(*[sp.fraction(x)[1] for x in v]); out.append([int(x * den) for x in v])
    return out
def dm_rank(rows):
    if not rows: return 0
    return DomainMatrix([[SQQ(int(x)) for x in r] for r in rows], (len(rows), len(rows[0])), SQQ).rank()
LHn = [x.astype(np.int64) for x in gens + Ec]; I_d = np.eye(d, dtype=np.int64)
# v2.28, FL_LOCAL=1: the torus part of Lie H at P is the Lie algebra of the local torus, {E_w : P compatible with Δ_w}
# (a linear condition on w), instead of all class generators. Needed where class weightings are not symmetries.
LOCAL = _os.environ.get('FL_LOCAL') == '1'
idxb_ = {b: i for i, b in enumerate(basis)}
DLc = []
for c in cl:
    D_ = np.zeros((d, d, d), dtype=np.int64)
    for l, b in enumerate(basis):
        if b[0] == 'v': continue
        for k, x in enumerate(b[1]):
            if x not in c: continue
            pre = idxb_[('p', b[1][:k])] if k > 0 else idxb_[('v', arrows[x][0])]
            suf = idxb_[('p', b[1][k + 1:])] if k + 1 < len(b[1]) else idxb_[('v', arrows[x][1])]
            D_[l, pre, suf] += 1
    DLc.append(D_)
def lie_at(Pq):
    if not LOCAL: return LHn
    cols_ = [(np.einsum('ml,mab->lab', Pq, D_) - np.einsum('ai,lij,bj->lab', Pq, D_, Pq)).ravel() for D_ in DLc]
    W = dm_null([list(r) for r in np.array(cols_).T])                    # w with Σ_c w_c (...)_c = 0
    return [x.astype(np.int64) for x in gens] + [sum(w[k] * Ec[k] for k in range(len(cl))).astype(np.int64) for w in W]
tab = collections.Counter(); cands = []; nolift = 0
if stage in ('exact', 'all'):
  for oid in range(len(reps)):
    Pq = next((q for q in (qfold(P) for P in members[oid]) if q is not None), None)
    if Pq is None: nolift += 1; continue
    LHP = lie_at(Pq)
    br = [(x @ Pq - Pq @ x).ravel() for x in LHP]                              # [x_i, P]
    r = dm_rank([list(v) for v in br])
    kerS = dm_null([list(row) for row in np.array(br).T])                      # coefficient vectors c with Σ c_i [x_i,P] = 0
    S = [sum(c[i] * LHP[i] for i in range(len(LHP))) for c in kerS]
    if not S: tab[(sizes[oid], r, 'free')] += 1; continue
    # commutant: (I⊗s − sᵀ⊗I) vec(X) = 0, row-major vec: [s,X] = sX − Xs  ->  kron(s, I) − kron(I, sᵀ)
    Acs = np.vstack([np.kron(s_, I_d) - np.kron(I_d, s_.T) for s_ in S])
    Cb = np.array(dm_null([list(row) for row in Acs]), dtype=np.int64)
    # normalizer: unknowns (c, m): Σ c_i [x_i, s_a] − Σ_b m_ab s_b = 0
    cols = []
    for x in LHP: cols.append(np.concatenate([(x @ s_ - s_ @ x).ravel() for s_ in S]))
    for a_ in range(len(S)):
        for b_ in range(len(S)):
            v = np.zeros(len(S) * d * d, dtype=np.int64); v[a_ * d * d:(a_ + 1) * d * d] = -S[b_].ravel(); cols.append(v)
    Nk = dm_null([list(row) for row in np.array(cols).T])
    Nv = [sum(c[i] * LHP[i] for i in range(len(LHP))) for c in Nk]
    n = dm_rank([list((x @ Pq - Pq @ x).ravel()) for x in Nv]) if Nv else 0
    J = jac_matrix(Pq); tS = Cb.shape[0] - rank_mod(J @ Cb.T)
    key_ = (sizes[oid], r, len(S), Cb.shape[0], n, tS)
    tab[key_] += 1
    if tS > n: cands.append((key_, Pq, Cb))
  pickle.dump(dict(tab=tab, cands=cands, nolift=nolift, norbits=len(reps)), open(CANDS, 'wb'))
else:
  dd = pickle.load(open(CANDS, 'rb')); tab, cands, nolift = dd['tab'], dd['cands'], dd['nolift']
print(f"{name} F_{p}: {len(reps)} H-orbits, {nolift} with no member lifting to Q")
print("  (orbit size, r, dim S, dim C_S, dim N·P, tangent of Folds^S) : orbits")
for k_, v in sorted(tab.items(), key=lambda kv: str(kv[0])): print(f"   {k_}: {v}")
if stage in ('numeric', 'all'):
    ver = collections.Counter(); seen_type = collections.Counter(); nfam = ntest = 0
    print("  numerical local dimension of Folds^S where tangent > dim N·P (at most %d per type):" % per_type, flush=True)
    sizes_only = {int(x) for x in _os.environ.get('FL_SIZES', '').split(',') if x}     # optional filter on orbit size
    for key_, Pq, Cb in cands:
        if sizes_only and key_[0] not in sizes_only: continue
        if seen_type[key_] >= per_type: continue
        seen_type[key_] += 1; ld = local_dim_in(Pq, Cb); ntest += 1
        fam = ld > key_[4]; nfam += fam
        print(f"   {key_}: local dim of Folds^S = {ld} -> {'FAMILY' if fam else 'no family'}", flush=True)
    print(f"  tested {ntest} of {len(cands)} candidates; with a family through them (local dim Folds^S > dim N·P): {nfam}", flush=True)
