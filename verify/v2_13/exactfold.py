"""Exact classification of folds on finite acyclic path algebras by decoupling:
   P is a fold  <=>  H = im P is a subalgebra and a subcoalgebra (Delta H ⊂ H⊗H),
                     Q = ker P is a subalgebra and a coideal (Delta Q ⊂ Q⊗A + A⊗Q),  and A = H ⊕ Q.
   So Folds = {(H, Q) complementary} ⊂ X_H × X_Q, an open subset of a product.  Each variety is cut out
   cell-by-cell on the Schubert-cell (RREF) decomposition of the Grassmannians, exactly over Q, with Groebner bases."""
import sys, itertools, sympy as sp
from quiver import build
QUIVERS = {'A2': (2, [(0, 1)]), 'A3': (3, [(0, 1), (1, 2)]), '1→2←3': (3, [(0, 1), (2, 1)]),
           'K2': (2, [(0, 1), (0, 1)]), 'K3': (2, [(0, 1), (0, 1), (0, 1)]), '1⇉2→3': (3, [(0, 1), (0, 1), (1, 2)])}
def tensors(nv, arrows):
    d, mult, Dl, allp = build(nv, arrows)
    M = {(i, j): [(k, int(mult[i, j, k])) for k in range(d) if mult[i, j, k]] for i in range(d) for j in range(d)}
    D = {l: [(a, b, int(Dl[l, a, b])) for a in range(d) for b in range(d) if Dl[l, a, b]] for l in range(d)}
    return d, M, D
def cells(d, r):
    for S in itertools.combinations(range(d), r):
        free = [(i, j) for i in range(r) for j in range(d) if j not in S and j > S[i]]
        yield S, free
def equations(d, M, D, S, free, kind):
    xs = sp.symbols(f'x0:{len(free)}') if free else ()
    r = len(S); B = [[sp.Integer(1 if j == S[i] else 0) for j in range(d)] for i in range(r)]
    for (i, j), x in zip(free, xs): B[i][j] = x
    def red(w):                      # kills V, injective on non-pivot coordinates
        return [w[j] - sum(w[S[i]] * B[i][j] for i in range(r)) for j in range(d)]
    def mul(u, v):
        out = [sp.Integer(0)] * d
        for i in range(d):
            if u[i] == 0: continue
            for j in range(d):
                if v[j] == 0: continue
                for k, c in M[(i, j)]: out[k] += c * u[i] * v[j]
        return out
    def delta(u):
        T = [[sp.Integer(0)] * d for _ in range(d)]
        for l in range(d):
            if u[l] == 0: continue
            for a, b, c in D[l]: T[a][b] += c * u[l]
        return T
    eqs = []
    for u, v in itertools.product(B, repeat=2): eqs += red(mul(u, v))           # subalgebra
    for u in B:
        T = delta(u)
        if kind == 'H':                                                         # Delta u ∈ V⊗V
            for b in range(d): eqs += red([T[a][b] for a in range(d)])
            for a in range(d): eqs += red(T[a])
        else:                                                                   # (red⊗red) Delta u = 0
            R1 = [red([T[a][b] for a in range(d)]) for b in range(d)]           # R1[b][a']
            for a2 in range(d): eqs += red([R1[b][a2] for b in range(d)])
    eqs = list({sp.expand(e) for e in eqs} - {0})
    return xs, B, eqs
def solve_cell(xs, eqs):
    """returns ('empty',), ('finite', [points]) or ('positive', dim) — dimension by maximal independent sets"""
    if not xs: return ('empty',) if eqs else ('finite', [()])
    if not eqs: return ('positive', len(xs))
    G = sp.groebner(eqs, *xs, order='grevlex')
    if list(G.exprs) == [1]: return ('empty',)
    if G.is_zero_dimensional:
        sols = sp.solve_poly_system(list(G.exprs), *xs)
        return ('finite', [tuple(s) for s in sols])
    # Krull dimension = size of a largest set of variables independent modulo the ideal
    lt = [sp.Poly(g, *xs).monoms(order='grevlex')[0] for g in G.exprs]   # leading monomials (grevlex)
    best = 0
    for k in range(len(xs), 0, -1):
        for U in itertools.combinations(range(len(xs)), k):
            if all(any(m[i] > 0 for i in range(len(xs)) if i not in U) for m in lt):
                return ('positive', k)                                      # Kredel–Weispfenning criterion on leading ideal
    return ('finite?', None)
def classify(name, kinds=('H', 'Q')):
    nv, arrows = QUIVERS[name]; d, M, D = tensors(nv, arrows); out = {}
    for kind in kinds:
        for r in range(0, d + 1):
            for S, free in cells(d, r):
                xs, B, eqs = equations(d, M, D, S, free, kind)
                res = solve_cell(xs, eqs)
                if res[0] == 'finite':
                    for p in res[1]:
                        sub = dict(zip(xs, p)); out.setdefault(kind, []).append(('pt', r, [[sp.nsimplify(sp.sympify(e).subs(sub)) for e in row] for row in B]))
                elif res[0] == 'positive':
                    out.setdefault(kind, []).append(('comp', r, S, res[1], len(xs)))
    return d, out
if __name__ == '__main__':
    for name in sys.argv[1:] or ['A2']:
        d, out = classify(name)
        for kind in ('H', 'Q'):
            pts = [o for o in out.get(kind, []) if o[0] == 'pt']; comps = [o for o in out.get(kind, []) if o[0] == 'comp']
            print(f"{name} X_{kind}: {len(pts)} isolated-cell points; positive-dim cell pieces: {[(c[1], c[2], c[3]) for c in comps]}", flush=True)
