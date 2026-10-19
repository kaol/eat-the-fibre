"""v2.32: the module class and dendriform completion variety of every fold orbit in a census.

For a fold P (column convention: P[:, i] = P(e_i)), with Q = 1 − P:
  L  ker·im ⊆ ker      P (Q x · P y) = 0          (Completion.LeftMod)
  R  im·ker ⊆ ker      P (P x · Q y) = 0          (Completion.RightMod)
  K  ker·im·ker = 0    Q x · (P y · Q z) = 0      (Completion.KIK)
  E  P multiplicative  P(xy) = Px · Py            (M3)
Each is invariant under conjugation by algebra automorphisms (G, the class torus, local torus moves), so it is an
invariant of every orbit notion of the census.

The completion variety V(P) ⊆ F_p⁴ is the set of b with ≺ = xP(y), ≻ = b1 xy + b2 P(x)y + b3 xP(y) + b4 P(x)P(y)
dendriform. It is computed exactly: the three axioms are polynomials of degree ≤ 2 in b with tensor coefficients,
evaluated at all p⁴ points. It is compared with the strata of `lean/EatTheFibre/Discordia/Completion.lean`.

Usage: python3 analysis/module_classes.py REPS.pkl [...]   (REPS from `pcensus.py … --save-reps`)
"""
import os, sys, pickle, itertools, collections, json
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_argv = sys.argv; sys.argv = ['x']; _cwd = os.getcwd(); os.chdir(ROOT); __file__ = os.path.join(ROOT, 'pcensus.py')
exec(open(os.path.join(ROOT, 'pcensus.py')).read().split("if __name__ == '__main" + "__':")[0])
sys.argv = _argv; os.chdir(_cwd)


def strata(p):
    """the predicted strata as sets of tuples mod p (Completion.lean)"""
    m = lambda *v: tuple(x % p for x in v)
    S = {'baxter': {m(-1, 1, 0, 0)},
         'line1': {m(-1 - s, 1 + s, s, -s) for s in range(p)},
         'point0': {m(0, 0, 0, 0)},
         'line2': {m(-t, 1 + t, t, -1 - t) for t in range(p)},
         'pointB': {m(0, 1, -1, 0)},
         'line0': {m(c, -c, -c, c) for c in range(p)},
         'lineB': {m(c, 1 - c, -1 - c, c) for c in range(p)}}
    return S


def analyse(P, M, p):
    d = P.shape[0]
    P = P.astype(np.int64) % p
    Q = (np.eye(d, dtype=np.int64) - P) % p
    I = np.eye(d, dtype=np.int64)
    mulc = lambda U, W: np.einsum('ia,jb,ijm->mab', U, W, M) % p          # products of column sets
    L = not (np.einsum('km,mab->kab', P, mulc(Q, P)) % p).any()
    R = not (np.einsum('km,mab->kab', P, mulc(P, Q)) % p).any()
    PQ = mulc(P, Q)                                                       # (m, j, l): P e_j · Q e_l
    K = not (np.einsum('ia,njl,inm->majl', Q, PQ, M) % p).any()
    E = not ((np.einsum('km,mab->kab', P, mulc(I, I)) - mulc(P, P)) % p).any()
    # bilinear tensors T[i, j, m]: e_i ∘ e_j = Σ_m T[i,j,m] e_m
    Tm = M % p
    prec = np.einsum('kj,ikm->ijm', P, M) % p
    S = [Tm, np.einsum('ki,kjm->ijm', P, M) % p, prec, np.einsum('ki,lj,klm->ijm', P, P, M) % p]
    o1 = lambda A, B: np.einsum('ijn,nlm->ijlm', B, A) % p                # A(B(x, y), z)
    o2 = lambda A, B: np.einsum('jln,inm->ijlm', B, A) % p                # A(x, B(y, z))
    mons = [()] + [(k,) for k in range(4)] + [(k, l) for k in range(4) for l in range(k, 4)]
    idx = {mo: n for n, mo in enumerate(mons)}
    C = np.zeros((len(mons), 3, d, d, d, d), dtype=np.int64)
    C[idx[()], 0] += o1(prec, prec) - o2(prec, prec)
    for k in range(4):
        C[idx[(k,)], 0] -= o2(prec, S[k])
        C[idx[(k,)], 1] += o1(prec, S[k]) - o2(S[k], prec)
        C[idx[(k,)], 2] -= o1(S[k], prec)
        for l in range(4):
            C[idx[tuple(sorted((k, l)))], 2] += o2(S[k], S[l]) - o1(S[k], S[l])
    C = C.reshape(len(mons), -1) % p
    C = C[:, C.any(0)]
    pts = np.array(list(itertools.product(range(p), repeat=4)), dtype=np.int64)
    V = np.stack([np.prod(pts[:, list(mo)], axis=1) if mo else np.ones(len(pts), np.int64) for mo in mons])
    ok = ~((C.T @ V) % p).any(0)
    var = {tuple(int(x) for x in pts[i]) for i in np.nonzero(ok)[0]}
    return dict(L=L, R=R, K=K, E=E), var


def signature(var, p):
    S = strata(p)
    contained = tuple(n for n, s in S.items() if s <= var)
    union = set().union(*(S[n] for n in contained)) if contained else set()
    return contained, len(var - union)


def main():
    for f in sys.argv[1:]:
        r = pickle.load(open(f, 'rb'))
        name, p = r['name'], r['p']
        nv, arrows = QUIVERS[name]
        basis, M, DL = build(nv, arrows)
        M = np.array(M, dtype=np.int64)
        assert [tuple(b) for b in basis] == [tuple(b) for b in r['basis']]
        per_orbit, comp_cls = [], collections.defaultdict(set)
        for i, P in enumerate(r['greps']):
            cl, var = analyse(P, M, p)
            sig, extra = signature(var, p)
            key = ''.join(k for k in 'LRKE' if cl[k]) or '-'
            per_orbit.append((key, sig, extra, len(var)))
            comp_cls[r['loc'][i]].add((key, sig, extra))
        assert all(len(v) == 1 for v in comp_cls.values()), "class not constant on local-torus components"
        # consistency with the Lean theorems: predicted strata are present
        for key, sig, extra, n in per_orbit:
            assert 'baxter' in sig
            if 'L' in key: assert {'line1', 'point0'} <= set(sig), (key, sig)
            if 'L' in key and 'R' in key: assert {'line2', 'pointB'} <= set(sig), (key, sig)
            if 'L' in key and 'K' in key: assert 'line0' in sig, (key, sig)
            if 'L' in key and 'R' in key and 'K' in key: assert 'lineB' in sig, (key, sig)
        cG = collections.Counter((k, s, e) for k, s, e, _ in per_orbit)
        cL = collections.Counter(next(iter(v)) for v in comp_cls.values())
        out = dict(name=name, p=p, orbits_G=len(per_orbit), orbits_loc=len(comp_cls),
                   classes=[dict(cls=k, strata=list(s), extra=e, G=cG[(k, s, e)], loc=cL[(k, s, e)])
                            for (k, s, e) in sorted(cG, key=lambda x: -cG[x])])
        print(json.dumps(out, ensure_ascii=False))


if __name__ == '__main__':
    main()
