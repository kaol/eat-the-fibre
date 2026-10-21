"""v2.35: locate the families of a census by Aut(kQ)-invariant signatures.

Every invariant used is preserved by conjugation with any algebra automorphism of kQ (it depends only on the
subalgebras im P, ker P and their position up to isomorphism), so it is constant on behaviours, hence on families
(v2.33). Local-torus orbits are grouped by signature; a signature whose group grows with q carries families.
Signature: module class (L, R, K, E), the vanishing pattern of the 18 typed words (v2.32), and isomorphism invariants of
the algebras im P and ker P (dimension, powers, annihilators, centre; for im P also idempotent and square-zero counts).
Usage: python3 analysis/family_signatures.py REPS_q1.pkl REPS_q2.pkl
"""
import os, sys, pickle, collections, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
_FAM_FILES = sys.argv[1:]; sys.argv = ["x"]
__file__ = os.path.join(HERE, 'completion_words.py')
exec(open(__file__).read().split("\ndef main():")[0])     # brings module_classes (analyse, …), vanishing, QUIVERS, build

def basis_cols(Mat, p):
    """column space basis mod p"""
    A = [list(map(int, r)) for r in (Mat.T % p).tolist()]; rows = []
    for r in A:
        for (pc, br) in rows:
            if r[pc]: f = r[pc]; r = [(x - f * y) % p for x, y in zip(r, br)]
        nz = [i for i, x in enumerate(r) if x]
        if nz:
            pc = nz[0]; inv = pow(r[pc], p - 2, p); r = [x * inv % p for x in r]
            rows = [(q, [(x - br[pc] * y) % p for x, y in zip(br, r)]) for q, br in rows] + [(pc, r)]
    return np.array([br for _, br in rows], dtype=np.int64).T if rows else np.zeros((Mat.shape[0], 0), np.int64)

def coords(Bs, v, p):
    d, r = Bs.shape
    # solve Bs x = v mod p via augmented elimination
    M = np.concatenate([Bs, v[:, None]], 1) % p; M = [list(map(int, x)) for x in M]; row = 0; piv = []
    for col in range(r):
        sel = [i for i in range(row, d) if M[i][col] % p]
        M[row], M[sel[0]] = M[sel[0]], M[row]; inv = pow(M[row][col], p - 2, p)
        M[row] = [x * inv % p for x in M[row]]
        for i in range(d):
            if i != row and M[i][col]: f = M[i][col]; M[i] = [(x - f * y) % p for x, y in zip(M[i], M[row])]
        row += 1
    return np.array([M[i][r] for i in range(r)], np.int64)

def sub_struct(M, Bs, p):
    r = Bs.shape[1]; T = np.zeros((r, r, r), np.int64)
    for a in range(r):
        for b in range(r):
            T[a, b] = coords(Bs, np.einsum('i,j,ijm->m', Bs[:, a], Bs[:, b], M) % p, p)
    return T

def rk(rows, p):
    return 0 if not len(rows) else basis_cols(np.array(rows).T, p).shape[1]

def alg_inv(T, p, enumerate_small=True):
    r = T.shape[0]
    if r == 0: return (0,)
    E = np.eye(r, dtype=np.int64)
    pw = [r]; cur = list(E)
    for _ in range(r + 1):
        nxt = [np.einsum('i,j,ijm->m', x, E[b], T) % p for x in cur for b in range(r)]
        dm = rk(nxt, p); pw.append(dm)
        if dm == 0 or dm == pw[-2]: break
        Bn = basis_cols(np.array(nxt).T, p); cur = [Bn[:, i] for i in range(Bn.shape[1])]
    null = lambda rows: r - rk(rows, p)
    annL = null([[T[a, b, m] for a in range(r)] for b in range(r) for m in range(r)])
    annR = null([[T[b, a, m] for a in range(r)] for b in range(r) for m in range(r)])
    cen = null([[(T[a, b, m] - T[b, a, m]) % p for a in range(r)] for b in range(r) for m in range(r)])
    out = (r, tuple(pw), annL, annR, cen)
    if False and enumerate_small and p ** r <= 200_000:   # counts are polynomials in q: not comparable across q
        X = np.array(list(itertools.product(range(p), repeat=r)), np.int64)
        sq = np.einsum('ni,nj,ijm->nm', X, X, T) % p
        out += (int((sq == X).all(1).sum()), int((sq == 0).all(1).sum()))
    return out

def signature(P, M, p):
    d = P.shape[0]; P = P.astype(np.int64) % p
    cl, _ = (None, None)
    Q = (np.eye(d, dtype=np.int64) - P) % p
    mulc = lambda U, W: np.einsum('ia,jb,ijm->mab', U, W, M) % p
    Pm = lambda T_: np.einsum('km,m...->k...', P, T_) % p
    z = lambda T_: not (T_ % p).any()
    L = z(Pm(mulc(Q, P))); R = z(Pm(mulc(P, Q)))
    K = z(np.einsum('ia,njl,inm->majl', Q, mulc(P, Q), M))
    Em = z(np.einsum('km,mab->kab', P, mulc(np.eye(d, dtype=np.int64), np.eye(d, dtype=np.int64))) - mulc(P, P))
    cls = ''.join(k for k, v in zip('LRKE', (L, R, K, Em)) if v) or '-'
    van = tuple(sorted(map(str, vanishing(P, M, p))))
    Im, Kr = basis_cols(P, p), basis_cols(Q, p)
    return (cls, van, alg_inv(sub_struct(M, Im, p), p), alg_inv(sub_struct(M, Kr, p), p, enumerate_small=False))

def census_sigs(f):
    r = pickle.load(open(f, 'rb')); name, p = r['name'], r['p']
    nv, arrows = QUIVERS[name]; basis, M, DL = build(nv, arrows); M = np.array(M, np.int64)
    loc = r['loc']; seen = {}; sig_of_loc = {}
    for i, P in enumerate(r['greps']):
        s = signature(P, M, p)
        if loc[i] in sig_of_loc: assert sig_of_loc[loc[i]] == s, "signature not constant on a local-torus orbit"
        sig_of_loc[loc[i]] = s
    return name, p, collections.Counter(sig_of_loc.values())

def main():
    res = [census_sigs(f) for f in _FAM_FILES]
    (n1, p1, c1), (n2, p2, c2) = res
    sigs = sorted(set(c1) | set(c2), key=lambda s: (s[0], str(s)))
    grow = [(s, c1[s], c2[s]) for s in sigs if c1[s] != c2[s]]
    print(json.dumps(dict(name=n1, q=[p1, p2], local_orbits=[sum(c1.values()), sum(c2.values())],
                          signatures=[len(c1), len(c2)], growing=len(grow))))
    tot = 0
    for s, a, b in grow:
        slope = (b - a) / (p2 - p1); tot += slope
        print(f"  class {s[0]:5} im {s[2]}  ker {s[3]}  local orbits {a} → {b}  (slope {slope:g})")
    print("  total slope", tot)

if __name__ == '__main__':
    main()
