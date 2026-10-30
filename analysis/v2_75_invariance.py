"""v2.75: H-invariance of φ on the family of v2_75_family (mult mode), symbolically.
S = the automorphism of kQ given by a generic element of H = ∏_classes (G_c × weighting): on class c, a_j ↦ s_c Σ_i
g^c_ij a_i with Σ_i g^c_ij = 1 (the census group), all g^c_ij (i < m−1) and s_c independent symbols. Checks
φ(S P S⁻¹) = φ(P) = λ_x m_x' / (λ_x' m_x) in the function field (cancel), for the given tails x, x'.
With INV_RANDOM=SEED: S at random rational g, s instead (a check, not a proof).
Usage: [INV_RANDOM=SEED] python3 analysis/v2_75_invariance.py QUIVER A0,A1,A2 X,X' C [OUT.json]"""
import os, sys, json, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sympy as sp
from v2_75_family import build
from v2_75_phi import phi
quiver = sys.argv[1]; cls = tuple(map(int, sys.argv[2].split(','))); x, x2 = map(int, sys.argv[3].split(',')); cmode = sys.argv[4]
out = sys.argv[5] if len(sys.argv) > 5 else None
ctx, R, Pr = build(quiver, cls, (x, x2), cmode, 'sigma0', 'mult')
d = ctx['d']; basis = ctx['basis']; arrows = ctx['arrows']; idx = {b: k for k, b in enumerate(basis)}
P = sp.Matrix(d, d, lambda r, c: Pr[r][c].as_expr())
classes = {}
for a, st in enumerate(arrows): classes.setdefault(st, []).append(a)
cl = list(classes.values()); cls_of = {a: k for k, arr in enumerate(cl) for a in arr}
G = []; Gi = []; s = sp.symbols(f's0:{len(cl)}')
for k, arr in enumerate(cl):
    m = len(arr)
    M = sp.Matrix(m, m, lambda i, j: sp.Symbol(f'g{k}_{i}{j}') if i < m - 1 else 0) if m > 1 else sp.Matrix([[1]])
    if m > 1:
        for j in range(m): M[m - 1, j] = 1 - sum(M[i, j] for i in range(m - 1))
    G.append(M); Gi.append(M.inv() if m > 1 else M)
def sigma(Gs, sc):
    S = sp.zeros(d, d)
    for col, b in enumerate(basis):
        if b[0] == 'v': S[col, col] = 1; continue
        choices = []
        for a in b[1]:
            k = cls_of[a]; arr = cl[k]; j = arr.index(a)
            choices.append([(arr[i], sc[k] * Gs[k][i, j]) for i in range(len(arr))])
        for combo in itertools.product(*choices):
            coef = sp.Integer(1)
            for _, cf in combo: coef *= cf
            S[idx[('p', tuple(a for a, _ in combo))], col] += coef
    return S
RND = os.environ.get('INV_RANDOM')                       # a seed: S at random rational g, s (a check, not a proof)
if RND:
    import random; rng = random.Random(int(RND))
    sub = {v: sp.Rational(rng.randint(-9, 9), rng.randint(1, 5)) for M in G for v in M.free_symbols}
    sub.update({v: sp.Rational(rng.choice([-1, 1]) * rng.randint(1, 9), rng.randint(1, 5)) for v in s})
    G = [M.subs(sub) for M in G]; assert all(M.det() != 0 for M in G); Gi = [M.inv() for M in G]; s = [sub[v] for v in s]
S = sigma(G, s); Si = sigma(Gi, [1 / v for v in s])
assert (S * Si).applyfunc(sp.cancel) == sp.eye(d)
P2 = (S * P * Si).applyfunc(sp.cancel)
lam = {y: sp.Symbol(f'l{y}') for y in (x, x2)}; mm = {y: sp.Symbol(f'm{y}') for y in (x, x2)}
expect = lam[x] * mm[x2] / (lam[x2] * mm[x])
f0, e0 = phi(P, basis, cl[cls_of[cls[0]]], x, x2)
f1, e1 = phi(P2, basis, cl[cls_of[cls[0]]], x, x2)
res = dict(random_seed=RND, quiver=quiver, cls=cls, tails=[x, x2], c=cmode, phi_P=str(f0), err_P=e0, phi_SPSinv=str(f1), err_SPSinv=e1,
           equal=(f0 is not None and f1 is not None and sp.cancel(f0 - f1) == 0),
           matches_formula=(f0 is not None and sp.cancel(f0 - expect) == 0), nsymbols=len((S * P).free_symbols))
print(json.dumps(res))
if out: json.dump(res, open(out, 'w'), indent=1)
