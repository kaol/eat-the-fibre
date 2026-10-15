"""Exact rigidity modulo Delta-compatible automorphisms, chart by chart (standard affine charts are OPEN, so they
see local dimension):  on each chart pair (S_H, S_Q) the fold scheme is  I = X_H(y) + X_Q(z) + (t·det - 1).
Let k = Krull dim I.  If I + (k×k minors of the infinitesimal-orbit matrix) = (1), every fold in the chart has
orbit dim ≥ k ≥ local dim ≥ orbit dim, so local dim = orbit dim there.  Charts cover every fold."""
import sys, itertools, sympy as sp
from exactfold import QUIVERS, tensors, equations
def derivations(d, M, D):
    Dv = sp.Matrix(d, d, lambda i, j: sp.Symbol(f'D{i}_{j}')); eqs = []
    def mulv(u, v):
        out = [0] * d
        for i in range(d):
            for j in range(d):
                if u[i] != 0 and v[j] != 0:
                    for k, c in M[(i, j)]: out[k] += c * u[i] * v[j]
        return out
    I = sp.eye(d)
    for i, j in itertools.product(range(d), repeat=2):                        # D(b_i b_j) = D(b_i) b_j + b_i D(b_j)
        bij = mulv(list(I[:, i]), list(I[:, j]))
        lhs = Dv * sp.Matrix(bij); r1 = mulv(list(Dv[:, i]), list(I[:, j])); r2 = mulv(list(I[:, i]), list(Dv[:, j]))
        eqs += [lhs[k] - r1[k] - r2[k] for k in range(d)]
    for l in range(d):                                                         # Delta D b_l = (D⊗1 + 1⊗D) Delta b_l
        T = sp.zeros(d, d)
        for a, b, c in D[l]: T[a, b] += c
        L = sp.zeros(d, d)
        for m in range(d):
            for a, b, c in D[m]: L[a, b] += c * Dv[m, l]
        R = Dv * T + T * Dv.T
        eqs += list(L - R)
    syms = list(Dv); A, _ = sp.linear_eq_to_matrix(eqs, syms)
    return [sp.Matrix(d, d, list(v)) for v in A.nullspace()]
def chart_vars(d, S, pre):
    free = [(i, j) for i in range(len(S)) for j in range(d) if j not in S]
    xs = sp.symbols(f'{pre}0:{len(free)}') if free else ()
    B = sp.Matrix(len(S), d, lambda i, j: 1 if j == S[i] else 0)
    for (i, j), x in zip(free, xs): B[i, j] = x
    return xs, B, free
def chart_eqs(d, M, D, S, kind, pre):
    xs, B, free = chart_vars(d, S, pre)
    # reuse the cell equation builder with ALL non-pivot entries free (standard chart)
    import exactfold
    xs2, B2, eqs = exactfold.equations(d, M, D, S, free, kind)
    sub = dict(zip(xs2, xs)); return xs, B, [sp.expand(e.subs(sub)) for e in eqs]
def krull(G, xs):
    if list(G.exprs) == [1]: return -1
    lt = [sp.Poly(g, *xs).monoms(order='grevlex')[0] for g in G.exprs]
    for k in range(len(xs), -1, -1):
        for U in itertools.combinations(range(len(xs)), k):
            if all(any(m[i] > 0 for i in range(len(xs)) if i not in U) for m in lt): return k
def tangent(B, S, Dm, free):
    """infinitesimal action of D on the chart coordinates of row space of B"""
    V = B * Dm.T; V = V - V[:, list(S)] * B
    return [V[i, j] for (i, j) in free]
def run(name):
    nv, arrows = QUIVERS[name]; d, M, D = tensors(nv, arrows); g = derivations(d, M, D)
    t = sp.Symbol('t'); report = []; allok = True
    print(f"{name}: dim A = {d}, dim g = {len(g)} (exact)", flush=True)
    for r in range(d + 1):
        for SH in itertools.combinations(range(d), r):
            ys, BH, eH = chart_eqs(d, M, D, SH, 'H', 'y')
            if eH and list(sp.groebner(eH, *ys, order='grevlex').exprs) == [1]: continue
            for SQ in itertools.combinations(range(d), d - r):
                zs, BQ, eQ = chart_eqs(d, M, D, SQ, 'Q', 'z')
                det = BH.col_join(BQ).det() if d else sp.Integer(1)
                xs = list(ys) + list(zs) + [t]; I = eH + eQ + [sp.expand(t * det - 1)]
                G = sp.groebner(I, *xs, order='grevlex'); k = krull(G, xs)
                if k < 0: continue
                fH = [(i, j) for i in range(r) for j in range(d) if j not in SH]; fQ = [(i, j) for i in range(d - r) for j in range(d) if j not in SQ]
                O = sp.Matrix([tangent(BH, SH, Dm, fH) + tangent(BQ, SQ, Dm, fQ) for Dm in g]) if g else sp.zeros(0, len(ys) + len(zs))
                if k == 0: ok = True
                elif O.rows < k: ok = False
                else:
                    minors = [sp.expand(O.extract(list(R), list(C)).det()) for R in itertools.combinations(range(O.rows), k)
                              for C in itertools.combinations(range(O.cols), k)]
                    minors = [m for m in minors if m != 0]
                    ok = bool(minors) and list(sp.groebner(I + minors, *xs, order='grevlex').exprs) == [1]
                allok &= ok; report.append((r, SH, SQ, k, ok))
    dims = sorted({k for (_, _, _, k, _) in report})
    print(f"{name}: nonempty chart pairs = {len(report)}, Krull dims occurring = {dims}; "
          f"local dim = orbit dim on every chart: {allok}", flush=True)
    return allok
if __name__ == '__main__':
    for n in sys.argv[1:] or ['K2']: run(n)
