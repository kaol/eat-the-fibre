"""Explicit family from a fixed-point locus (v2.26).

Takes a candidate saved by fixedlocus.py (an H-orbit representative P over Q with stabilizer S, commutant basis B of
C_S, and a numerical family through P). Writes X = P + Σ y_i B_i, restricts all fold equations to C_S, and decomposes
the variety with Singular (minimal associated primes). Reports every component through P (y = 0): its dimension, and
whether the orbit of N = N_H(S) through P fills it. A component of dimension > dim N·P carries a family of H-orbits with
stabilizer S.
Usage: python3 analysis/family.py QUIVER p KEY_INDEX [which]   (KEY_INDEX: position in the sorted list of candidate
types printed by --list; which: index of the candidate within that type)"""
import os as _os, sys as _sys
if 'ROOT' not in globals():
    AN = _os.path.dirname(_os.path.abspath(__file__)); ROOT = _os.path.dirname(AN); _os.chdir(ROOT)
    __file__ = _os.path.join(ROOT, 'pcensus.py')
import sys, pickle, collections, numpy as np
name, p = sys.argv[1], int(sys.argv[2])
cd = pickle.load(open(f'logs/fixedlocus_cands_{name}_{p}.pkl', 'rb'))
types = sorted({c[0] for c in cd['cands']}, key=str)
if '--list' in sys.argv:
    for i, t in enumerate(types): print(i, t, sum(1 for c in cd['cands'] if c[0] == t))
    sys.exit()
ki, which = int(sys.argv[3]), int(sys.argv[4]) if len(sys.argv) > 4 else 0
sys.argv = ['x']
exec(open(_os.path.join(ROOT, 'pcensus.py')).read().split("if __name__ == '__main__':")[0])
ctx = setup(name, p); basis, d, M, DL = ctx['basis'], ctx['d'], ctx['M'], ctx['DL']
key_, Pq, Cb = [c for c in cd['cands'] if c[0] == types[ki]][which]
from sage.all__sagemath_singular import QQ, PolynomialRing, ideal, matrix
k = Cb.shape[0]
R = PolynomialRing(QQ, [f'y{i}' for i in range(k)], order='degrevlex'); ys = R.gens()
X = matrix(R, d, d, lambda i, j: Pq[i, j] + sum(int(Cb[t, i * d + j]) * ys[t] for t in range(k) if Cb[t, i * d + j]))
Id = matrix.identity(R, d)
def mul(x, y):
    out = [R(0)] * d
    for i in range(d):
        if x[i] == 0: continue
        for j in range(d):
            if y[j] == 0: continue
            for kk in np.nonzero(M[i, j])[0]: out[kk] += x[i] * y[j] * int(M[i, j, kk])
    return out
eqs = set(e for e in (X * X - X).list() if e != 0)
cols = [X.column(i) for i in range(d)]; qcols = [(Id - X).column(i) for i in range(d)]
for i in range(d):
    for j in range(d):
        h = mul(cols[i], cols[j]); hv = matrix(R, d, 1, h)
        eqs |= set(e for e in (X * hv - hv).list() if e != 0)
        q = mul(qcols[i], qcols[j]); qv = matrix(R, d, 1, q)
        eqs |= set(e for e in (X * qv).list() if e != 0)
for l in range(d):
    lhs = sum((X[m, l] * matrix(R, DL[m].tolist()) for m in range(d)), matrix(R, d, d))
    eqs |= set(e for e in (lhs - X * matrix(R, DL[l].tolist()) * X.transpose()).list() if e != 0)
I = ideal(R, list(eqs))
print(f"{name} F_{p}: candidate type {key_} #{which}; dim C_S = {k}; {len(eqs)} equations", flush=True)
G = I.groebner_basis(); I = ideal(R, G)
print(f"  Krull dimension of Folds ∩ C_S (all components): {I.dimension()}", flush=True)
comps = I.minimal_associated_primes()
origin = {y: 0 for y in ys}
through = [c for c in comps if all(g.subs(origin) == 0 for g in c.gens())]
print(f"  {len(comps)} components; {len(through)} through P; dims of those: {[c.dimension() for c in through]}", flush=True)
for c in through:
    lin = [g for g in c.groebner_basis() if g.degree() <= 1]
    nonlin = [g for g in c.groebner_basis() if g.degree() > 1]
    print(f"  component dim {c.dimension()}: {len(lin)} linear and {len(nonlin)} nonlinear generators; nonlinear: {nonlin[:6]}", flush=True)
pickle.dump(dict(key=key_, P=Pq, Cb=Cb, through=[[str(g) for g in c.gens()] for c in through]),
            open(f'logs/family_{name}_{p}_{ki}_{which}.pkl', 'wb'))
