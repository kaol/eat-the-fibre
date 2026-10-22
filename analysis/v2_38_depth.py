"""v2.38: Aut(kQ)-invariant row depths of folds, on the 19 families of 1⇉2⇉3 and on F1–F4 of 1⇉2⇇3.

Every φ ∈ Aut(kQ) satisfies ε∘φ = φ̄∘ε (ε : kQ → kQ/rad = k^{Q₀}, φ̄ a permutation of vertices) and φ(rad^n) = rad^n.
For φPφ⁻¹ this gives ε(φPφ⁻¹) = φ̄·εP·φ⁻¹. So the following are invariants, up to φ̄ acting on vertex labels:
  d_i      = max{n : ε_i P(rad^n) ≠ 0}  (−1 if row i of P vanishes)
  rkε_n    = rank εP|_{rad^n}
  im_n     = dim(im P ∩ rad^n),   ker_n = dim(ker P ∩ rad^n)
The script computes D over ℚ(t) (generic member) and at every t ∈ 𝔽_q, and asserts that D is constant along each
family (prediction 5a).
Usage: python3 analysis/v2_38_depth.py CERT.pkl PAIRS_IN_REPS_q.pkl
"""
import os, sys, pickle, collections, io, contextlib
import numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
cert_file, pin_file = sys.argv[1:3]

def rank_q(M, q):
    M = [[int(x) % q for x in row] for row in M]
    if not M or not M[0]: return 0
    return sp.Matrix(M).rank(iszerofunc=lambda x: x % q == 0) if False else _rank_mod(M, q)
def _rank_mod(M, q):
    M = [row[:] for row in M]; r = 0; ncol = len(M[0])
    for c in range(ncol):
        piv = next((i for i in range(r, len(M)) if M[i][c] % q), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]; inv = pow(M[r][c], q - 2, q)
        M[r] = [x * inv % q for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] % q:
                f = M[i][c]; M[i] = [(x - f * y) % q for x, y in zip(M[i], M[r])]
        r += 1
    return r
def rank_Q(M):
    return sp.Matrix(M).rank(simplify=True) if len(M) and len(M[0]) else 0

def D(P, nvert, radn, rank):
    """P: d×d (list of lists); radn[n] = list of basis indices of rad^n, n = 0.."""
    d = len(P); N = len(radn)
    def sub(rows, cols): return [[P[i][j] for j in cols] for i in rows]
    depths = []
    for i in range(nvert):
        dep = -1
        for n in range(N):
            if rank(sub([i], radn[n])) > 0: dep = n
        depths.append(dep)
    rP = rank(P)
    rk_eps = tuple(rank(sub(range(nvert), radn[n])) for n in range(N))
    im = tuple(rP - rank(sub([i for i in range(d) if i not in radn[n]], range(d))) if n else rP for n in range(N))
    ker = tuple(len(radn[n]) - rank(sub(range(d), radn[n])) for n in range(N))
    return tuple(depths), rk_eps, im, ker

# ---- 1⇉2⇉3: the 19 families
cert = pickle.load(open(cert_file, 'rb'))
Vvec = [[1] + [0] * 10, [0, 1] + [0] * 9, [0, 0, 1] + [0] * 8, [0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 1]]
t = sp.Symbol('t')
def famP(k):
    v = cert[k]; P = sp.zeros(11, 11)
    for j in v['nz']:
        c = [sp.sympify(x) for x in v['cols'][j]]
        for q_ in range(11): P[q_, j] = sp.expand(sum(c[a] * Vvec[a][q_] for a in range(6)))
    return P
rad_tp = [list(range(11)), list(range(3, 11)), list(range(7, 11))]
q = 7
print("== 1⇉2⇉3, the 19 families: D = (depths d_1..d_3, rank εP|rad^n, dim im∩rad^n, dim ker∩rad^n)")
res = {}
for k in sorted(cert):
    P = famP(k)
    gen = D(P.tolist(), 3, rad_tp, rank_Q)
    pts = {D(P.subs(t, tv).tolist(), 3, rad_tp, lambda M: rank_q(M, q)) for tv in range(q)}
    assert pts == {gen}, (k, gen, pts)
    res[k] = gen
    print(f"family {k:2d} ({cert[k]['tag']}): {gen}")
groups = collections.defaultdict(list)
for k in sorted(cert): groups[(cert[k]['tag'], res[k])].append(k)
print("distinct (signature, D):", len(groups), "of", len(cert))
for g, ks in groups.items():
    if len(ks) > 1: print("  shared:", ks, g)

# ---- 1⇉2⇇3: F1–F4 from the v2.36 inner-curve search
print("\n== 1⇉2⇇3 (pairs_in): families from analysis/family_curves.py")
saved = sys.argv; sys.argv = ['family_curves.py', pin_file]
g = {'__name__': 'fc', '__file__': os.path.join(HERE, 'family_curves.py')}
with contextlib.redirect_stdout(io.StringIO()):
    exec(open(os.path.join(HERE, 'family_curves.py')).read(), g)
sys.argv = saved
pq = g['p']; basis = g['basis']
nv = sum(1 for b in basis if b[0] == 'v')
radn = [list(range(len(basis)))]
n = 1
while True:
    cols = [i for i, b in enumerate(basis) if b[0] == 'p' and len(b[1]) >= n]
    if not cols: break
    radn.append(cols); n += 1
fams = [f for f in g['fams'] if len(f[0]) == pq]
for hit, (i, vec) in fams:
    Pz = g['lifted'](i); Inn, Innv = g['INN'][vec]
    Ps = (Inn * sp.Matrix(Pz.tolist()) * Innv).applyfunc(sp.expand)
    gen = D(Ps.tolist(), nv, radn, rank_Q)
    pts = {D(Ps.subs(g['t'], tv).tolist(), nv, radn, lambda M: rank_q(M, pq)) for tv in range(pq)}
    sig = g['sig_of_loc'][next(iter(hit))]
    print(f"rep {i} along {vec}, class {sig[0]}: generic D {gen}; constant on all t ∈ 𝔽_{pq}: {pts == {gen}}"
          + ("" if pts == {gen} else f" (values {sorted(pts)})"))
print("vertices of pairs_in in basis order:", [b for b in basis if b[0] == 'v'], "; Aut(Q) swaps the two sources")
