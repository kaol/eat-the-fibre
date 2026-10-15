"""Exact: is the fold variety a finite union of orbits of the Delta-compatible automorphism group?
Chart pair (S_H, S_Q) (open affine charts of Gr × Gr), scheme I = X_H + X_Q + (t·det − 1).
finite_orbits(J): for every minimal prime P of J with k = dim P > 0, the infinitesimal-orbit matrix must have
rank k generically on P (some k-minor ∉ P: an orbit is open dense in V(P)); then recurse on P + (k-minors),
the G-stable closed locus where orbits are smaller.  Every stratum is then a finite union of orbits."""
import sys, itertools, sympy as sp
from sage.all__sagemath_singular import PolynomialRing, QQ, matrix
from exactfold import QUIVERS, tensors
from chartpair import derivations, chart_eqs, tangent
def finite_orbits(R, J, O, depth=0):
    fails = []
    for P in J.minimal_associated_primes():
        k = P.dimension() - (1 if False else 0)
        if k <= 0: continue
        if O.nrows() < k: return [('rank', k, 0)]
        mins = O.minors(k)
        nz = [m for m in mins if m not in P]
        if not nz: fails.append(('generic orbit rank < dim', k)); continue
        fails += finite_orbits(R, P + R.ideal(mins), O, depth + 1)
    return fails
def run(name, only_r=None):
    nv, arrows = QUIVERS[name]; d, M, D = tensors(nv, arrows); g = derivations(d, M, D)
    print(f"{name}: dim A = {d}, dim g = {len(g)}", flush=True)
    npairs = 0; bad = []; dims = set()
    for r in (range(d + 1) if only_r is None else only_r):
        for SH in itertools.combinations(range(d), r):
            ys, BH, eH = chart_eqs(d, M, D, SH, 'H', 'y')
            if eH and list(sp.groebner(eH, *ys, order='grevlex').exprs) == [1]: continue
            for SQ in itertools.combinations(range(d), d - r):
                zs, BQ, eQ = chart_eqs(d, M, D, SQ, 'Q', 'z')
                if eQ and zs and list(sp.groebner(eQ, *zs, order='grevlex').exprs) == [1]: continue
                names = [str(v) for v in ys] + [str(v) for v in zs] + ['t']
                R = PolynomialRing(QQ, len(names), names, order='degrevlex'); loc = dict(zip(names, R.gens()))
                det = BH.col_join(BQ).det() if d else sp.Integer(1)
                gens = [R(str(e).replace('**', '^')) if True else None for e in eH + eQ] + [R(str(sp.expand(sp.Symbol('t') * det - 1)).replace('**', '^'))]
                J = R.ideal(gens)
                if J.dimension() < 0: continue
                npairs += 1
                fH = [(i, j) for i in range(r) for j in range(d) if j not in SH]; fQ = [(i, j) for i in range(d - r) for j in range(d) if j not in SQ]
                rows = [[R(str(sp.expand(e)).replace('**', '^')) for e in tangent(BH, SH, Dm, fH) + tangent(BQ, SQ, Dm, fQ)] for Dm in g]
                O = matrix(R, rows) if rows else matrix(R, 0, len(names) - 1)
                dims.add(J.dimension())
                f = finite_orbits(R, J, O)
                if f: bad.append((r, SH, SQ, f))
    print(f"{name}: nonempty chart pairs = {npairs}, chart dims = {sorted(dims)}; fold variety is a finite union of G-orbits: {not bad}", flush=True)
    for b in bad[:5]: print("   FAIL", b, flush=True)
    return not bad
if __name__ == '__main__':
    for n in sys.argv[1:] or ['K2']: run(n)
