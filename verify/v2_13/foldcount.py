"""Exact fold counts for quivers whose X_H is finite: enumerate complementary pairs (H, Q), rebuild P and check the
original fold residual exactly; on positive-dimensional Q-pieces prove det(H; Q) ≡ 0 (no complement) by Groebner."""
import sys, itertools, sympy as sp
from exactfold import QUIVERS, tensors, cells, equations, solve_cell
def check_P(d, M, D, P):
    I = sp.eye(d); Pm = sp.Matrix(P)
    if Pm * Pm != Pm: return False
    def mul(u, v):
        out = sp.zeros(d, 1)
        for i in range(d):
            for j in range(d):
                if u[i] and v[j]:
                    for k, c in M[(i, j)]: out[k] += c * u[i] * v[j]
        return out
    for i, j in itertools.product(range(d), repeat=2):
        x, y = Pm[:, i], Pm[:, j]
        if sp.simplify(Pm * mul(x, y) - mul(x, y)) != sp.zeros(d, 1): return False
        if sp.simplify(Pm * mul(I[:, i] - x, I[:, j] - y)) != sp.zeros(d, 1): return False
    for j in range(d):            # Delta P = (P⊗P) Delta  on basis vector j
        lhs = sp.zeros(d, d)
        for l in range(d):
            if Pm[l, j]:
                for a, b, c in D[l]: lhs[a, b] += c * Pm[l, j]
        T = sp.zeros(d, d)
        for a, b, c in D[j]: T[a, b] += c
        if sp.simplify(lhs - Pm * T * Pm.T) != sp.zeros(d, d): return False
    return True
def run(name):
    nv, arrows = QUIVERS[name]; d, M, D = tensors(nv, arrows)
    pts = {'H': [], 'Q': []}; pieces = {'H': [], 'Q': []}
    for kind in 'HQ':
        for r in range(d + 1):
            for S, free in cells(d, r):
                xs, B, eqs = equations(d, M, D, S, free, kind); res = solve_cell(xs, eqs)
                if res[0] == 'finite':
                    for p in res[1]:
                        sub = dict(zip(xs, p)); pts[kind].append(sp.Matrix(B).subs(sub) if B else sp.zeros(0, d))
                elif res[0] == 'positive': pieces[kind].append((xs, B, eqs, res[1]))
    assert not pieces['H'], "X_H not finite: use the chart-pair engine"
    folds = []
    for Hm in pts['H']:
        for Qm in pts['Q']:
            if Hm.rows + Qm.rows != d: continue
            Bm = Hm.col_join(Qm).T                          # columns: basis of H then Q
            if sp.simplify(Bm.det()) == 0: continue
            r = Hm.rows; E = sp.diag(*([1] * r + [0] * (d - r))) if d else sp.zeros(0)
            folds.append(sp.simplify(Bm * E * Bm.inv()))
    ok = all(check_P(d, M, D, P) for P in folds)
    # positive-dimensional Q pieces: det(H;Q) vanishes identically for every H of complementary dimension
    no_comp = True
    t = sp.Symbol('t')
    for xs, B, eqs, k in pieces['Q']:
        for Hm in pts['H']:
            if Hm.rows + len(B) != d: continue
            det = sp.Matrix(Hm).col_join(sp.Matrix(B)).det()
            G = sp.groebner(eqs + [sp.expand(t * det - 1)], *xs, t, order='grevlex')
            no_comp &= (list(G.exprs) == [1])
    ranks = sorted(sp.Matrix(P).rank() for P in folds)
    print(f"{name}: |X_H| = {len(pts['H'])}, isolated Q-cell points = {len(pts['Q'])}, positive-dim Q pieces = {len(pieces['Q'])}; "
          f"folds = {len(folds)} (all satisfy the fold equations exactly: {ok}); positive-dim Q pieces admit no complement: {no_comp}; "
          f"rank distribution {dict((x, ranks.count(x)) for x in set(ranks))}", flush=True)
    return folds
if __name__ == '__main__':
    for n in sys.argv[1:] or ['A2']: run(n)
