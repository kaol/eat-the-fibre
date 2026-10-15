"""Two-u consistency for identification folds, in exponential coordinates e_r = (1 - r a)^{-1}.
Elements of k[a] are kept as dicts r -> coeff meaning sum coeff * e_r (e_0 = 1). B(e_r, e_q) = C(r,q) e_r e_q.
Condition: Phi(x, z1, z2, y) = B(B(x,z1) z2, y) + B(x, z1 B(z2,y)) - B(x,z1) B(z2,y) depends only on z1*z2,
i.e. Phi(e_s, e_sig, e_tau, e_t) = Phi(e_s, 1, e_sig*e_tau, e_t)."""
import sympy as sp
s, sig, tau, t = sp.symbols('s sigma tau t')
Cf = sp.Function('C')
def emul(A, Bv):
    r = {}
    for p, cp in A.items():
        for q, cq in Bv.items():
            if p == q: raise ValueError("coincident exponents")
            for k, v in ((p, p / (p - q)), (q, -q / (p - q))):
                if k == 0 and False: pass
                r[k] = sp.together(r.get(k, 0) + cp * cq * v)
    return r
def mul(A, Bv):
    # handle e_0 = 1 exactly
    r = {}
    for p, cp in A.items():
        for q, cq in Bv.items():
            if p == 0: r[q] = r.get(q, 0) + cp * cq
            elif q == 0: r[p] = r.get(p, 0) + cp * cq
            else:
                for k, v in emul({p: 1}, {q: 1}).items(): r[k] = r.get(k, 0) + cp * cq * v
    return r
def Bil(A, Bv):
    r = {}
    for p, cp in A.items():
        for q, cq in Bv.items():
            for k, v in mul({p: 1}, {q: 1}).items(): r[k] = r.get(k, 0) + cp * cq * Cf(p, q) * v
    return r
def add(*Ds):
    r = {}
    for sgn, D in Ds:
        for k, v in D.items(): r[k] = r.get(k, 0) + sgn * v
    return r
def Phi(x, z1, z2, y):
    return add((1, Bil(mul(Bil(x, z1), z2), y)), (1, Bil(x, mul(z1, Bil(z2, y)))), (-1, mul(Bil(x, z1), Bil(z2, y))))
E = lambda r: {r: sp.Integer(1)}
lhs = Phi(E(s), E(sig), E(tau), E(t)); rhs = Phi(E(s), {0: 1}, mul(E(sig), E(tau)), E(t))
diff = add((1, lhs), (-1, rhs))
eqs = {k: sp.factor(sp.together(v)) for k, v in diff.items()}
for k, v in eqs.items(): print(f"coefficient of e_{k}:  {v}\n")
import pickle; pickle.dump(eqs, open('feq.pkl', 'wb'))
