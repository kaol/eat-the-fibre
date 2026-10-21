"""v2.34: infinitesimal rigidity of behaviours and of behaviour fibres, at every census orbit representative.

At a fold P (columns = images), over 𝔽_p:
  T = T_P Folds          kernel of the linearised fold equations (idempotent, image closed, kernel closed, Δ)
  F = [Der, P] ∩ T       Der = derivations of kQ (Lie Aut(kQ)); [D, P] = DP − PD
  S = [Der_Δ, P] + [E_w, P] (w ∈ W_P)   Der_Δ = Δ-compatible derivations; E_w = class-weight derivations;
                                         W_P = {w : P compatible with Δ_w = Σ_c w_c Δ_c} (the local torus, v2.24)
  b = dim T − dim F      fold directions not produced by automorphisms (behaviour excess)
  e = dim F − dim S      automorphism directions that stay among folds without being census symmetries (fibre excess)
S ⊆ F ⊆ T is asserted at every point.
Usage: python3 analysis/fibre_rigidity.py REPS.pkl ...   (REPS from pcensus.py --save-reps; logs/v2_32/reps.tar.gz)
"""
import os, sys, pickle, collections, json, itertools
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_files = sys.argv[1:]; sys.argv = ['x']; _cwd = os.getcwd(); os.chdir(ROOT); __file__ = os.path.join(ROOT, 'pcensus.py')
exec(open(__file__).read().split("if __name__ == '__main" + "__':")[0])
os.chdir(_cwd)


def rref(A, p):
    """row echelon form mod p; returns (R, pivots)"""
    A = np.array(A, dtype=np.int64) % p
    piv = []; r = 0
    for c in range(A.shape[1]):
        if r == A.shape[0]: break
        nz = np.nonzero(A[r:, c])[0]
        if not len(nz): continue
        i = r + nz[0]
        if i != r: A[[r, i]] = A[[i, r]]
        A[r] = A[r] * pow(int(A[r, c]), p - 2, p) % p
        col = A[:, c].copy(); col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr): A[nzr] = (A[nzr] - np.outer(col[nzr], A[r])) % p
        piv.append(c); r += 1
    return A[:r], piv


def rank(A, p): return len(rref(A, p)[1]) if np.size(A) else 0


def nullspace(A, p):
    """basis (rows) of {x : A x = 0} mod p"""
    n = A.shape[1]
    R, piv = rref(A, p)
    free = [c for c in range(n) if c not in piv]
    N = np.zeros((len(free), n), dtype=np.int64)
    for k, f in enumerate(free):
        N[k, f] = 1
        for r, c in enumerate(piv): N[k, c] = (-R[r, f]) % p
    return N


def setup(name, p):
    nv, arrows = QUIVERS[name]
    basis, M, DL = build(nv, arrows)
    M = np.array(M, dtype=np.int64); DL = np.array(DL, dtype=np.int64); d = len(basis)
    idxb = {b: i for i, b in enumerate(basis)}
    classes = collections.defaultdict(list)
    for ai, st in enumerate(arrows): classes[st].append(ai)
    cl = list(classes.values()); cls_of = {a: k for k, c in enumerate(cl) for a in c}
    # per-class deletion Δ_c (as in pcensus.finish)
    DLc = np.zeros((len(cl), d, d, d), dtype=np.int64)
    for l, b in enumerate(basis):
        if b[0] == 'v': continue
        for k, x in enumerate(b[1]):
            pre = idxb[('p', b[1][:k])] if k > 0 else idxb[('v', arrows[x][0])]
            suf = idxb[('p', b[1][k + 1:])] if k + 1 < len(b[1]) else idxb[('v', arrows[x][1])]
            DLc[cls_of[x], l, pre, suf] += 1
    assert ((DLc.sum(0) - DL) % p == 0).all(), "per-class deletions do not sum to Δ"
    # derivations: D (d×d, column j = D e_j) with D(e_i e_j) = D(e_i) e_j + e_i D(e_j)
    # unknown vector x = D.ravel() (row-major: D[a, j] at a*d + j)
    rows = []
    for i in range(d):
        for j in range(d):
            for m in range(d):
                row = np.zeros(d * d, dtype=np.int64)
                for a in range(d): row[m * d + a] += M[i, j, a]          # (D (e_i e_j))_m = Σ_a D[m,a] M[i,j,a]
                for a in range(d): row[a * d + i] -= M[a, j, m]          # (D e_i · e_j)_m = Σ_a D[a,i] M[a,j,m]
                for a in range(d): row[a * d + j] -= M[i, a, m]          # (e_i · D e_j)_m
                if row.any(): rows.append(row % p)
    Der = nullspace(np.array(rows), p).reshape(-1, d, d)
    # Δ-compatible: Δ(D e_l) = (D⊗1 + 1⊗D) Δ e_l, i.e. Σ_m D[m,l] DL[m] = D DL[l] + DL[l] Dᵀ
    rows = []
    for l in range(d):
        for a_ in range(d):
            for b_ in range(d):
                row = np.zeros(d * d, dtype=np.int64)
                for m in range(d): row[m * d + l] += DL[m, a_, b_]
                for k in range(d): row[a_ * d + k] -= DL[l, k, b_]
                for k in range(d): row[b_ * d + k] -= DL[l, a_, k]
                if row.any(): rows.append(row % p)
    coef = nullspace(np.concatenate([np.array(rows) @ Der.reshape(len(Der), -1).T % p]), p) if len(Der) else np.zeros((0, 0))
    DerD = (coef @ Der.reshape(len(Der), -1) % p).reshape(-1, d, d) if len(coef) else np.zeros((0, d, d), np.int64)
    # class-weight derivations E_c (diagonal: path ↦ (#arrows of class c in it) · path)
    E = np.zeros((len(cl), d, d), dtype=np.int64)
    for i, b in enumerate(basis):
        if b[0] == 'p':
            for x in b[1]: E[cls_of[x], i, i] += 1
    return dict(d=d, M=M, DL=DL, DLc=DLc, Der=Der, DerD=DerD, E=E, ncl=len(cl))


def jacobian(P, C, p):
    """rows: linearised fold equations in the direction V (unknown V.ravel(), row-major)"""
    d, M, DL = C['d'], C['M'], C['DL']
    I = np.eye(d, dtype=np.int64); Q = (I - P) % p
    cols = []
    for idx in range(d * d):
        V = np.zeros((d, d), dtype=np.int64); V.flat[idx] = 1
        e1 = (P @ V + V @ P - V).ravel()
        h = np.einsum('ia,jb,ijm->mab', P, P, M)                              # P e_a · P e_b
        dh = np.einsum('ia,jb,ijm->mab', V, P, M) + np.einsum('ia,jb,ijm->mab', P, V, M)
        e2 = (np.einsum('km,mab->kab', V, h) + np.einsum('km,mab->kab', P, dh) - dh).ravel()
        q = np.einsum('ia,jb,ijm->mab', Q, Q, M)
        dq = -np.einsum('ia,jb,ijm->mab', V, Q, M) - np.einsum('ia,jb,ijm->mab', Q, V, M)
        e3 = (np.einsum('km,mab->kab', V, q) + np.einsum('km,mab->kab', P, dq)).ravel()
        e4 = (np.einsum('ml,mab->lab', V, DL) - np.einsum('ai,lij,bj->lab', V, DL, P)
              - np.einsum('ai,lij,bj->lab', P, DL, V)).ravel()
        cols.append(np.concatenate([e1, e2, e3, e4]) % p)
    J = np.array(cols).T
    return J                                        # all rows, aligned with residual_series


def at(P, C, p):
    d = C['d']; P = P.astype(np.int64) % p
    J = jacobian(P, C, p)
    T = d * d - rank(J, p)
    br = lambda D: ((D @ P - P @ D) % p).ravel()
    orbit = np.array([br(D) for D in C['Der']]) if len(C['Der']) else np.zeros((0, d * d), np.int64)
    # F = img ∩ ker J: dim = rank(img) − rank(J · img) restricted to a basis of img
    if len(orbit):
        Rb, piv = rref(orbit, p); B = Rb                                    # rows: basis of [Der, P]
        dimF = len(B) - rank((J @ B.T) % p, p)
    else: B = orbit; dimF = 0
    # W_P: w with Σ_m P[m,l] Δ_w[m] = P Δ_w[l] Pᵀ, linear in w
    DLc = C['DLc']
    cond = (np.einsum('ml,cmab->clab', P, DLc) - np.einsum('ai,clij,bj->clab', P, DLc, P)) % p   # (c, l, a, b)
    W = nullspace(cond.reshape(C['ncl'], -1).T % p, p)                                          # rows: w
    Ew = [np.einsum('c,cij->ij', w, C['E']) % p for w in W]
    S = [br(D) for D in C['DerD']] + [br(Ew_) for Ew_ in Ew]
    Smat = np.array(S) if S else np.zeros((0, d * d), np.int64)
    dimS = rank(Smat, p)
    if len(Smat):                                                           # sanity: S ⊆ F (in [Der,P] and in T)
        assert not ((J @ Smat.T) % p).any(), "census symmetry direction leaves T_P Folds"
        assert rank(np.concatenate([B, Smat]), p) == len(B), "census symmetry direction not in [Der, P]"
    return dict(T=T, F=dimF, S=dimS, b=T - dimF, e=dimF - dimS)


def main():
    for f in _files:
        r = pickle.load(open(f, 'rb')); name, p = r['name'], r['p']
        C = setup(name, p)
        hist = collections.Counter(); ex = []
        for i, P in enumerate(r['greps']):
            a = at(P, C, p)
            hist[(a['b'], a['e'])] += 1
            if a['e'] or a['b']: ex.append((i, a))
        out = dict(name=name, p=p, reps=len(r['greps']), dimDer=len(C['Der']), dimDerD=len(C['DerD']),
                   be_hist={f"b={b},e={e}": n for (b, e), n in sorted(hist.items())})
        print(json.dumps(out), flush=True)



# ---- second order (v2.34): does an excess direction [D, P] extend to a curve s ↦ exp(sD + s²D₂ + …)·P in Folds?
def residual_series(Ps, C, p):
    """fold equations evaluated on a truncated series P(s) = Ps[0] + s Ps[1] + s² Ps[2]; returns the s² coefficient"""
    d, M, DL = C['d'], C['M'], C['DL']; n = len(Ps)
    I = np.eye(d, dtype=np.int64)
    def mm(A, B):                                   # series product of matrix series
        out = np.zeros((n, d, d), np.int64)
        for i in range(n):
            for j in range(n - i): out[i + j] += A[i] @ B[j]
        return out % p
    def amul(A, B):                                 # columnwise algebra products: (A e_a)·(B e_b), series
        out = np.zeros((n, d, d, d), np.int64)
        for i in range(n):
            for j in range(n - i): out[i + j] += np.einsum('ia,jb,ijm->mab', A[i], B[j], M)
        return out % p
    def act(A, X):                                  # A applied to vector-valued series X (n, d, ...)
        out = np.zeros_like(X)
        for i in range(n):
            for j in range(n - i): out[i + j] += np.einsum('km,m...->k...', A[i], X[j])
        return out % p
    Qs = (-Ps) % p; Qs[0] = (I - Ps[0]) % p
    e1 = (mm(Ps, Ps) - Ps) % p
    h = amul(Ps, Ps); e2 = (act(Ps, h) - h) % p
    e3 = act(Ps, amul(Qs, Qs))
    e4 = np.zeros((n, d, d, d), np.int64)
    for i in range(n): e4[i] += np.einsum('ml,mab->lab', Ps[i], DL)
    for i in range(n):                              # P Δ Pᵀ is quadratic in P
        for j in range(n - i): e4[i + j] -= np.einsum('ai,lij,bj->lab', Ps[i], DL, Ps[j])
    e4 %= p
    return np.concatenate([e1[2].ravel(), e2[2].ravel(), e3[2].ravel(), e4[2].ravel()]) % p


def second_order(P, C, p):
    """the excess directions of F over S, tested at second order. Returns (e, number of excess classes [D] mod S that
    extend to second order), counted over all nonzero classes in F/S with entries in 𝔽_p (up to scalars)."""
    d = C['d']; P = P.astype(np.int64) % p
    J = jacobian(P, C, p)
    br = lambda D: ((D @ P - P @ D) % p)
    Der = C['Der']
    # D_F = {D ∈ Der : J [D,P] = 0}, as coefficient vectors on the Der basis
    A = np.array([(J @ br(D).ravel()) % p for D in Der]).T
    cF = nullspace(A, p) if len(Der) else np.zeros((0, 0), np.int64)
    DF = [np.einsum('k,kij->ij', c, Der) % p for c in cF]
    # S directions
    cond = (np.einsum('ml,cmab->clab', P, C['DLc']) - np.einsum('ai,clij,bj->clab', P, C['DLc'], P)) % p
    W = nullspace(cond.reshape(C['ncl'], -1).T % p, p)
    Sdir = [br(D).ravel() for D in C['DerD']] + [br(np.einsum('c,cij->ij', w, C['E']) % p).ravel() for w in W]
    Smat = np.array(Sdir) if Sdir else np.zeros((0, d * d), np.int64)
    rS = rank(Smat, p) if len(Smat) else 0
    # complement of S inside [D_F, P]: greedy basis
    comp = []; cur = Smat.copy() if len(Smat) else np.zeros((0, d * d), np.int64); r0 = rS
    for D in DF:
        v = br(D).ravel()
        test = np.concatenate([cur, v[None]]) if len(cur) else v[None]
        if rank(test, p) > r0: comp.append(D); cur = test; r0 += 1
    e = len(comp)
    if e == 0: return 0, 0
    inv2 = pow(2, p - 2, p)
    imgJ = np.array([(J @ br(D2).ravel()) % p for D2 in Der]).T                # J[D₂, P] over the Der basis
    rimg = rank(imgJ, p)
    ext = 0
    seen = set()
    for coeffs in itertools.product(range(p), repeat=e):
        if not any(coeffs): continue
        # projective representative: first nonzero coefficient 1
        k0 = next(i for i, x in enumerate(coeffs) if x); inv = pow(coeffs[k0], p - 2, p)
        cn = tuple((x * inv) % p for x in coeffs)
        if cn in seen: continue
        seen.add(cn)
        D = sum(c * Dk for c, Dk in zip(cn, comp)) % p
        V = br(D); W0 = (D @ V - V @ D) % p * inv2 % p                               # ½ [D, [D, P]]
        Ps = np.array([P, V, W0]) % p
        r2 = residual_series(Ps, C, p)
        # need D₂ with J [D₂, P] = −r2
        aug = np.concatenate([imgJ, ((-r2) % p)[:, None]], 1)
        if rank(aug, p) == rimg: ext += 1
    return e, ext


if __name__ == '__main__':
    main()
