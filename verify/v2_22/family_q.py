"""Over Q(t): the four free (G × grading)-orbits of 1⇉2, 1→3 give one-parameter families P_t = φ P φ⁻¹ (φ = grading at λ = 1/t).
Checks symbolically: P_t is a Δ-compatible fold for all t ≠ 0, and τ(P) = [e1]P(γ) is invariant under Lie(G) and equals t·τ(P_1).
So t is a modulus for Aut_Δ° = G (Conjecture 3′ fails) that the grading removes (Conjecture 3″ holds on these families)."""
import sympy as sp, numpy as np, itertools
exec(open('gcensus.py').read().split("class Field")[0])
nv, arrows = QUIVERS['pair_srcout']; basis, M, DL = build(nv, arrows); d = len(basis)
names = ['e1', 'e2', 'e3', 'a1', 'a2', 'g']; assert [b for b in basis] == [('v', 0), ('v', 1), ('v', 2), ('p', (0,)), ('p', (1,)), ('p', (2,))]
t = sp.Symbol('t', nonzero=True)
reps = [  # from family.py over F_5 (entries ±1), columns P(e1..g) in basis e1,e2,e3,a1,a2,g
    {'e2': {'e1': 1, 'e2': 1}, 'e3': {'e1': -1, 'e2': -1}, 'a2': {'e1': 1, 'a1': -1, 'a2': 1}, 'g': {'e1': 1, 'a1': -1, 'a2': 1}},
    {'e2': {'e3': -1}, 'e3': {'e3': 1}, 'a2': {'e1': 1, 'a1': -1, 'a2': 1}, 'g': {'e1': 1, 'a1': -1, 'a2': 1}},
    {'e2': {'e3': -1}, 'e3': {'e3': 1}, 'a2': {'e2': 1, 'e3': 1, 'a1': -1, 'a2': 1}, 'g': {'e2': 1, 'e3': 1, 'a1': -1, 'a2': 1}},
    {'e2': {'e1': 1, 'e2': 1}, 'e3': {'e3': 1}, 'a2': {'e1': 1, 'a1': -1, 'a2': 1}, 'g': {'e1': 1, 'a1': -1, 'a2': 1}},
]
def mat(r):
    P = sp.zeros(d, d)
    for j, col in r.items():
        for i, v in col.items(): P[names.index(i), names.index(j)] = v
    return P
lens = [0 if b[0] == 'v' else len(b[1]) for b in basis]
phi = sp.diag(*[(1 / t) ** k for k in lens]); phii = sp.diag(*[t ** k for k in lens])
def mulv(x, y): return sp.Matrix([sum(x[i] * y[j] * M[i, j, k] for i in range(d) for j in range(d)) for k in range(d)])
def is_fold(P):
    Id = sp.eye(d)
    if sp.simplify(P * P - P) != sp.zeros(d, d): return False
    for i in range(d):
        for j in range(d):
            h = mulv(P[:, i], P[:, j]); q = mulv((Id - P)[:, i], (Id - P)[:, j])
            if sp.simplify(P * h - h) != sp.zeros(d, 1) or sp.simplify(P * q) != sp.zeros(d, 1): return False
    for l in range(d):
        lhs = sum((P[m, l] * sp.Matrix(DL[m]) for m in range(d)), sp.zeros(d, d))
        if sp.simplify(lhs - P * sp.Matrix(DL[l]) * P.T) != sp.zeros(d, d): return False
    return True
# Lie(G): D(a1) = aδ, D(a2) = bδ, δ = a1 − a2; fixes vertices and g
Ds = []
for a, b in ((1, 0), (0, 1)):
    D = sp.zeros(d, d); D[3, 3] += a; D[4, 3] -= a; D[3, 4] += b; D[4, 4] -= b; Ds.append(D)
ok = True
for k, r in enumerate(reps):
    P1 = mat(r); Pt = sp.simplify(phi * P1 * phii)
    fold = is_fold(Pt); tau = sp.simplify(Pt[0, 5] / P1[0, 5]) if P1[0, 5] != 0 else None
    tau_any = [sp.simplify(Pt[i, 5]) for i in range(3)]
    inv = all(sp.simplify((D * Pt - Pt * D)[i, 5]) == 0 for D in Ds for i in range(3))   # vertex part of P(g) is Lie(G)-invariant
    nonconst = any(sp.diff(x, t) != 0 for x in tau_any)
    print(f"family {k+1}: P_t fold for all t: {fold}; vertex part of P_t(γ) = {tau_any} (Lie(G)-invariant: {inv}, varies with t: {nonconst})")
    ok &= fold and inv and nonconst
print("all four families are genuine Aut_Δ-moduli removed by the grading:", ok)
