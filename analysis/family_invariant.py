"""Invariant of a family component (v2.26), exactly over Q.

Input: logs/family_QUIVER_p_KI_W.pkl from family.py (P, the commutant basis B of C_S, and the generators of the
component of Folds ∩ C_S through P). Steps:
 1. Parametrize the component: solve its linear generators, then each remaining nonlinear generator for a variable
    occurring linearly in it, leaving free parameters u_1..u_m. X(u) = P + Σ y_i(u) B_i.
 2. Verify symbolically that X(u) is a Δ-compatible fold for all u (idempotent, image and kernel subalgebras,
    ΔX = (X⊗X)Δ).
 3. Lie N = N_H(S): vector fields on the parameters, from [x, X(u)] = Σ_j (∂X/∂u_j) V_j(u).
 4. Rational invariants: f(u) with V f = 0 for all x ∈ Lie N. Since N·X(u) has codimension 1 in the component, the
    invariant field has transcendence degree 1. Polynomial and rational candidates are found from monomials up to a
    degree bound (linear algebra), then checked.
Usage: python3 analysis/family_invariant.py QUIVER p KI W"""
import os as _os, sys as _sys
if 'ROOT' not in globals():
    AN = _os.path.dirname(_os.path.abspath(__file__)); ROOT = _os.path.dirname(AN); _os.chdir(ROOT)
    __file__ = _os.path.join(ROOT, 'pcensus.py')
import sys, pickle, itertools, collections, numpy as np, sympy as sp
name, p, ki, w = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
fam = pickle.load(open(f'logs/family_{name}_{p}_{ki}_{w}.pkl', 'rb'))
sys.argv = ['x']
exec(open(_os.path.join(ROOT, 'pcensus.py')).read().split("if __name__ == '__main__':")[0])
ctx = setup(name, p); basis, arrows, d, M, DL = ctx['basis'], ctx['arrows'], ctx['d'], ctx['M'], ctx['DL']
Pq, Cb = fam['P'], fam['Cb']; k = Cb.shape[0]
ys = sp.symbols(f'y0:{k}')
gens_ = [sp.sympify(g, locals={f'y{i}': ys[i] for i in range(k)}) for g in fam['through'][0]]
# 1. parametrize: repeatedly substitute, and solve any generator for a variable it contains linearly with a constant
#    coefficient; then (rationally) for a variable occurring linearly with a polynomial coefficient.
sol = {}
def subst(e): return sp.expand(e.subs(sol))
def step(const_only):
    for g in gens_:
        g2 = subst(g)
        if g2 == 0: continue
        for y in sorted(g2.free_symbols, key=lambda s_: int(str(s_)[1:])):
            if sp.degree(g2, y) != 1: continue
            c = sp.expand(sp.diff(g2, y))
            if const_only and c.free_symbols: continue
            val = sp.simplify(-(g2 - c * y) / c)
            for a in list(sol): sol[a] = sp.simplify(sol[a].subs(y, val))
            sol[y] = val; return True
    return False
while step(True): pass
while step(False):
    while step(True): pass
left = [subst(g) for g in gens_ if subst(g) != 0]
assert not left, f"unsolved generators: {left}"
free = sorted(set(ys) - set(sol), key=lambda s: int(str(s)[1:]))
yv = [sp.simplify(sol.get(y, y)) for y in ys]
X = sp.Matrix(d, d, lambda i, j: Pq[i, j] + sum(int(Cb[t, i * d + j]) * yv[t] for t in range(k) if Cb[t, i * d + j]))
print(f"{name}: family component of type {fam['key']}: free parameters {free}")
# 2. verify fold equations symbolically
Id = sp.eye(d)
def mul(x, y):
    out = [0] * d
    for i in range(d):
        if x[i] == 0: continue
        for j in range(d):
            if y[j] == 0: continue
            for kk in np.nonzero(M[i, j])[0]: out[kk] += x[i] * y[j] * int(M[i, j, kk])
    return sp.Matrix(out)
ok = sp.simplify(X * X - X) == sp.zeros(d, d)
for i in range(d):
    for j in range(d):
        h = mul(X[:, i], X[:, j]); q = mul((Id - X)[:, i], (Id - X)[:, j])
        ok &= sp.simplify(X * h - h) == sp.zeros(d, 1) and sp.simplify(X * q) == sp.zeros(d, 1)
for l in range(d):
    lhs = sum((X[m, l] * sp.Matrix(DL[m].tolist()) for m in range(d)), sp.zeros(d, d))
    ok &= sp.simplify(lhs - X * sp.Matrix(DL[l].tolist()) * X.T) == sp.zeros(d, d)
print(f"  X(u) is a Δ-compatible fold for all parameter values: {ok}")
# 3. Lie H, S, N at the base point, vector fields
gens, E = lie_G(ctx)
cls = collections.defaultdict(list)
for a, st in enumerate(arrows): cls[st].append(a)
Ec = [np.diag([0 if b[0] == 'v' else sum(1 for x in b[1] if x in c) for b in basis]).astype(np.int64) for c in cls.values()]
LH = [sp.Matrix(x.tolist()) for x in gens + Ec]
P0 = sp.Matrix(Pq.tolist())
A = sp.Matrix([list(x * P0 - P0 * x) for x in LH]).T
S = [sum((v[i] * LH[i] for i in range(len(LH))), sp.zeros(d, d)) for v in A.nullspace()]
cs = sp.symbols(f'c0:{len(LH)}'); xs = sum((cs[i] * LH[i] for i in range(len(LH))), sp.zeros(d, d))
ms = sp.symbols(f'm0:{len(S)**2}'); neq = []
for a_, s_ in enumerate(S):
    T_ = xs * s_ - s_ * xs - sum((ms[a_ * len(S) + b_] * S[b_] for b_ in range(len(S))), sp.zeros(d, d))
    neq += [e for e in T_ if e != 0]
An, _ = sp.linear_eq_to_matrix(neq, list(cs) + list(ms))
Nv = [sum((v[i] * LH[i] for i in range(len(LH))), sp.zeros(d, d)) for v in An.nullspace()]
Jx = sp.Matrix([[sp.diff(X[i, j], u) for u in free] for i in range(d) for j in range(d)])     # (d², m)
V = []
for x in Nv:
    rhs = sp.Matrix([(x * X - X * x)[i, j] for i in range(d) for j in range(d)])
    vsol = sp.symbols(f'v0:{len(free)}')
    eqv = [e for e in (Jx * sp.Matrix(vsol) - rhs) if sp.simplify(e) != 0]
    sv = sp.solve(eqv, vsol, dict=True)
    if sv: V.append([sp.simplify(sv[0].get(v_, 0)) for v_ in vsol])
print(f"  Lie N: {len(Nv)} generators; vector fields on the parameters: {V}")
# 4. invariants: rational functions f = num/den of low degree with V(f) = 0
u = free
def invariants(maxdeg=2):
    mons = [sp.Mul(*[u[i] ** e[i] for i in range(len(u))]) for tot in range(maxdeg + 1)
            for e in itertools.product(range(tot + 1), repeat=len(u)) if sum(e) == tot]
    out = []
    for vf in V:
        pass
    # f = Σ a_m m / Σ b_m m ; V(f) = 0  <=>  den·V(num) − num·V(den) = 0 (bilinear) — search num, den among monomial ratios
    for m1, m2 in itertools.combinations(mons, 2):
        f = m1 / m2
        if all(sp.simplify(sum(vf[i] * sp.diff(f, u[i]) for i in range(len(u)))) == 0 for vf in V):
            out.append(sp.simplify(f))
    return out
inv = invariants()
print(f"  monomial-ratio invariants (degree ≤ 2): {inv[:6]}")
if not inv:
    # general: solve the first-order PDE system by characteristics through sympy's pdsolve is not available for systems;
    # fall back to polynomial invariants of degree ≤ 2 by linear algebra
    a = sp.symbols('a0:40')
    mons = [m for m in sorted({sp.Mul(*[u[i] ** e[i] for i in range(len(u))]) for e in itertools.product(range(3), repeat=len(u)) if sum(e) <= 2}, key=str)]
    f = sum(a[i] * mons[i] for i in range(len(mons)))
    eqs = []
    for vf in V:
        expr = sp.expand(sp.numer(sp.together(sum(vf[i] * sp.diff(f, u[i]) for i in range(len(u))))))
        eqs += sp.Poly(expr, *u).coeffs()
    s = sp.solve(eqs, a[:len(mons)], dict=True)
    print(f"  polynomial invariants (degree ≤ 2): {[sp.simplify(f.subs(si)) for si in s]}")
# 5. generic orbit rank of Lie H on the family (so a generic member has stabilizer of dim = dim S), the invariant in
#    shifted coordinates, and the family written in path names
rank_generic = sp.Matrix([list(x * X - X * x) for x in LH]).rank(simplify=True)
print(f"  generic rank of Lie H at X(u): {rank_generic}  (dim H = {len(LH)}, dim S = {len(S)})")
cand_t = None
for shift in itertools.product((-1, 0, 1), repeat=len(u)):
    z = [u[i] + shift[i] for i in range(len(u))]
    for e in itertools.product((-1, 0, 1), repeat=len(u)):
        if all(x == 0 for x in e): continue
        f = sp.Mul(*[z[i] ** e[i] for i in range(len(u))])
        if all(sp.simplify(sum(vf[i] * sp.diff(f, u[i]) for i in range(len(u)))) == 0 for vf in V):
            cand_t = f; break
    if cand_t is not None: break
print(f"  N-invariant (shifted monomial): t = {cand_t}")
names = []
cls_name = {}
letters = iter('abcdefgh')
for st, c in cls.items():
    L_ = next(letters)
    for j, a_ in enumerate(c): cls_name[a_] = f"{L_}{j+1}" if len(c) > 1 else L_
nm = [f"e{b[1]+1}" if b[0] == 'v' else "".join(cls_name[x] for x in b[1]) for b in basis]
print("  family X(u), columns P(b) = Σ X[i, b] i:")
for j in range(d):
    terms = [f"({sp.factor(X[i, j])})·{nm[i]}" for i in range(d) if X[i, j] != 0]
    if terms: print(f"    P({nm[j]}) = " + " + ".join(terms))
