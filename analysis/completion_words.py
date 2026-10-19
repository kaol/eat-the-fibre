"""v2.32: is the completion variety of a fold determined by which typed words vanish?

Write x = p + n (p ∈ im P, n ∈ ker P). The three dendriform axioms, expanded, are combinations of degree-3 words in
the pieces of x, y, z (in that order), where two adjacent pieces may be grouped under P as an opaque image element
U = P(ab) (only for mixed ab; P(ii) = ii, P(kk) = 0). There are 16 such typed words:
    abc (a, b, c ∈ {i, k}),  U(ab)·c and a·U(bc) (ab, bc mixed).
`vanishing(P, M, p)` computes which of them vanish identically for a census fold (plus L, R: U_ki = 0, U_ik = 0).
`predicted(vanish, p)` treats every other word as free and returns the set of b ∈ F_p⁴ satisfying the axioms.
Usage: python3 analysis/completion_words.py REPS.pkl ...   (compares predicted with actual for every orbit)
"""
import os, sys, pickle, itertools, collections, functools
import numpy as np
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
_files = sys.argv[1:]; sys.argv = ['x']
__file__ = os.path.join(HERE, 'module_classes.py'); exec(open(__file__).read().split("def main")[0])

B = sp.symbols('b1:5')

# ---- symbolic side: words over letters ('p',v) ('n',v) ('U',(l1,l2)); elements are dict word -> coeff
def typ(l):
    if l[0] == 'p': return 'i'
    if l[0] == 'n': return 'k'
    return 'U' + typ(l[1][0]) + typ(l[1][1])
def add(*es):
    r = collections.defaultdict(lambda: 0)
    for e in es:
        for w, c in e.items(): r[w] += c
    return {w: sp.expand(c) for w, c in r.items() if sp.expand(c) != 0}
def sc(a, e): return {w: a * c for w, c in e.items()}
def mul(e, f):
    r = collections.defaultdict(lambda: 0)
    for (w1, c1), (w2, c2) in itertools.product(e.items(), f.items()): r[w1 + w2] += c1 * c2
    return dict(r)
def make(vanish):
    def P(e):
        out = {}
        for w, c in e.items():
            if len(w) == 1: o = {w: 1} if typ(w[0]) == 'i' else {}
            elif len(w) == 2:
                t = typ(w[0]) + typ(w[1])
                if t == 'ii': o = {w: 1}
                elif t == 'kk': o = {}
                elif ('U' + t) in vanish: o = {}
                else: o = {(('U', w),): 1}
            else: raise ValueError(w)
            out = add(out, sc(c, o))
        return out
    return P
def predicted_eqs(vanish):
    P = make(vanish)
    el = lambda v: {(('p', v),): 1, (('n', v),): 1}
    X, Y, Z = el(0), el(1), el(2)
    prec = lambda x, y: mul(x, P(y))
    succ = lambda x, y: add(sc(B[0], mul(x, y)), sc(B[1], mul(P(x), y)), sc(B[2], mul(x, P(y))), sc(B[3], mul(P(x), P(y))))
    star = lambda x, y: add(prec(x, y), succ(x, y))
    eqs = set()
    for lhs, rhs in [(prec(prec(X, Y), Z), prec(X, star(Y, Z))), (prec(succ(X, Y), Z), succ(X, prec(Y, Z))),
                     (succ(X, succ(Y, Z)), succ(star(X, Y), Z))]:
        for w, c in add(lhs, sc(-1, rhs)).items():
            if tuple(typ(l) for l in w) in vanish: continue
            eqs.add(sp.expand(c))
    return [e for e in eqs if e != 0]
@functools.lru_cache(None)
def predicted(vanish, p):
    eqs = predicted_eqs(vanish)
    pts = list(itertools.product(range(p), repeat=4))
    fs = [sp.lambdify(B, e, 'math') for e in eqs]
    return frozenset(pt for pt in pts if all(int(f(*pt)) % p == 0 for f in fs))

# ---- numeric side
def vanishing(P, M, p):
    d = P.shape[0]; P = P.astype(np.int64) % p; Q = (np.eye(d, dtype=np.int64) - P) % p
    mats = {'i': P, 'k': Q}
    def mulc(U, W):                                                       # (d, A...) × (d, B...) → (d, A..., B...)
        R = np.einsum('ia,jb,ijm->mab', U.reshape(d, -1), W.reshape(d, -1), M) % p
        return R.reshape((d,) + U.shape[1:] + W.shape[1:])
    Pm = lambda T: np.einsum('km,m...->k...', P, T) % p
    z = lambda T: not (T % p).any()
    v = set()
    for a, b in [('k', 'i'), ('i', 'k')]:
        if z(Pm(mulc(mats[a], mats[b]))): v.add('U' + a + b)
    for a, b, c in itertools.product('ik', repeat=3):
        if z(mulc(mulc(mats[a], mats[b]), mats[c])): v.add((a, b, c))
    for a, b in [('k', 'i'), ('i', 'k')]:
        for c in 'ik':
            if z(mulc(Pm(mulc(mats[a], mats[b])), mats[c])): v.add(('U' + a + b, c))
            if z(mulc(mats[c], Pm(mulc(mats[a], mats[b])))): v.add((c, 'U' + a + b))
    return frozenset(v)

def main():
    for f in _files:
        r = pickle.load(open(f, 'rb')); name, p = r['name'], r['p']; nv, arrows = QUIVERS[name]
        basis, M, DL = build(nv, arrows); M = np.array(M, dtype=np.int64)
        agree, bad, pats = 0, [], collections.Counter()
        for i, P in enumerate(r['greps']):
            _, var = analyse(P, M, p)
            van = vanishing(P, M, p)
            pats[van] += 1
            if predicted(van, p) == frozenset(var): agree += 1
            else: bad.append((i, sorted(map(str, van)), len(var), len(predicted(van, p))))
        print(f"{name} p={p}: {agree}/{len(r['greps'])} orbits agree; {len(pats)} vanishing patterns", flush=True)
        for b in bad[:5]: print('   mismatch', b)

if __name__ == '__main__':
    main()
