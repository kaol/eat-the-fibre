"""v2.33: every Theorem 2-4 family is one Aut(kQ)-orbit on its generic range; explicit witnesses, verified exactly.

For each family P_t (the kernel-checked matrices of lean/EatTheFibre/Gen/*.lean), an algebra automorphism
φ = Inn_u ∘ T (T an arrow torus, u = 1 + r) with φ P_t φ⁻¹ = P_t' is given in closed form, as rational functions of
(t, t'). It is checked symbolically: the identity φ P_t = P_t' φ holds in Q(t, t'), and φ is an algebra map (T is
multiplicative on paths, Inn_u is conjugation). Denominators vanish exactly at the degenerate parameters.
With δ_a = a1 − a2, δ_b = b1 − b2:
  Theorem 2 (1⇉2⇉3):  T = 1,                         u = 1 + (t'−t) δ_a δ_b
  Theorem 4 (1⇉2⇉3):  T: a_i ↦ (t−1)/(t'−1) a_i,      u = 1 − (t'−t)/(t'−1) δ_a δ_b              (t, t' ≠ 1)
  Theorem 3 F2, F4:   T = 1,                         u = 1 + (t'−t)(δ_b − δ_a)
  Theorem 3 F1:       T: a2 ↦ t/t' a2, b1 ↦ t/t' b1, b2 ↦ t(t+1)/(t'(t'+1)) b2,
                      u = 1 + (t'−t)/t' a2 + (t'−t)/(t'(t'+1)) b2                                (t, t' ∉ {0, −1})
  Theorem 3 F3:       F1 with the classes a and b exchanged in the roles of the arrows                (t, t' ∉ {0, −1})
Usage: python3 analysis/family_behaviour.py
"""
import os, re, sympy as sp
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def grab(src, name):
    i = src.index(f"def {name}"); j = src.index(":=", i) + 2
    k = src.find("\n\n", j); body = src[j:k] if k > 0 else src[j:]
    return eval(body.split("/--")[0].replace("![", "[").strip())

t, t2 = sp.symbols('t t2')
W = {  # family: (quiver file, torus scalings per arrow name, r as {path name: coefficient})
    'thm2P':   ('two_pairs', {}, {'a1b1': t2 - t, 'a1b2': -(t2 - t), 'a2b1': -(t2 - t), 'a2b2': t2 - t}),
    'thm4P':   ('two_pairs', {'a1': (t - 1) / (t2 - 1), 'a2': (t - 1) / (t2 - 1)},
                {k: -(t2 - t) / (t2 - 1) * v for k, v in {'a1b1': 1, 'a1b2': -1, 'a2b1': -1, 'a2b2': 1}.items()}),
    'thm3F2P': ('pairs_in', {}, {'a1': -(t2 - t), 'a2': t2 - t, 'b1': t2 - t, 'b2': -(t2 - t)}),
    'thm3F4P': ('pairs_in', {}, {'a1': -(t2 - t), 'a2': t2 - t, 'b1': t2 - t, 'b2': -(t2 - t)}),
    'thm3F1P': ('pairs_in', {'a2': t / t2, 'b1': t / t2, 'b2': t * (t + 1) / (t2 * (t2 + 1))},
                {'a2': (t2 - t) / t2, 'b2': (t2 - t) / (t2 * (t2 + 1))}),
    'thm3F3P': ('pairs_in', {'b2': t / t2, 'a1': t / t2, 'a2': t * (t + 1) / (t2 * (t2 + 1))},
                {'b2': (t2 - t) / t2, 'a2': (t2 - t) / (t2 * (t2 + 1))}),
}

def check(fam):
    q, tor, rr = W[fam]
    src = open(os.path.join(ROOT, 'lean', 'EatTheFibre', 'Gen', q + '.lean')).read()
    c = grab(src, 'c'); d = len(c); P = grab(src, fam)
    names = [x.split(': ')[1] for x in re.search(r"Basis \(index: path\): (.*)", src).group(1).split(', ')]
    Pt = lambda tv: sp.Matrix(d, d, lambda i, j: sum(a * tv**n for n, a in enumerate(P[i][j])))
    def mul(x, y): return [sp.expand(sum(c[i][j][m] * x[i] * y[j] for i in range(d) for j in range(d) if c[i][j][m])) for m in range(d)]
    E = lambda j: [1 if i == j else 0 for i in range(d)]
    # torus: diagonal, multiplicative on paths (names are concatenations of two-letter arrow names)
    tdiag = [sp.Integer(1) if n.startswith('e') else sp.prod([tor.get(n[k:k + 2], 1) for k in range(0, len(n), 2)]) for n in names]
    T = sp.diag(*tdiag)
    for i in range(d):                       # T is an algebra map: T(e_i e_j) = T e_i · T e_j
        for j in range(d):
            assert sp.simplify(sp.Matrix(mul(list(T[:, i]), list(T[:, j]))) - T * sp.Matrix(mul(E(i), E(j)))) == sp.zeros(d, 1)
    one = [1 if n.startswith('e') else 0 for n in names]
    r = [rr.get(n, 0) for n in names]
    u = [one[i] + r[i] for i in range(d)]
    uinv, pw = one[:], one[:]
    for _ in range(d): pw = mul(pw, [-x for x in r]); uinv = [sp.together(a + b) for a, b in zip(uinv, pw)]
    assert all(sp.simplify(x - y) == 0 for x, y in zip(mul(u, uinv), one))
    Inn = sp.Matrix([mul(mul(u, E(j)), uinv) for j in range(d)]).T
    Phi = Inn * T
    lhs = (Phi * Pt(t)).applyfunc(sp.simplify); rhs = (Pt(t2) * Phi).applyfunc(sp.simplify)
    ok = (lhs - rhs).applyfunc(sp.simplify) == sp.zeros(d, d)
    den = sorted({str(f) for e in list(Phi) for f in sp.factor_list(sp.denom(sp.together(e)))[1] for f in [f[0]]})
    return ok, den

if __name__ == '__main__':
    for fam in W:
        ok, den = check(fam)
        print(f"{fam:8} φ P_t = P_t' φ in Q(t, t'): {ok};  denominators of φ: {den or 'none'}")
