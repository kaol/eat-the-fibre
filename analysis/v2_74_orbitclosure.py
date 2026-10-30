"""v2.74: the closure of the orbit of rational points x of a component Z (char 0), by eliminating the orbit map of a
slice of H = G × T (all class weightings) through 1 transversal to the stabiliser (pivots of the differential).
Prints the extra generators of the closure ideal beyond I_Z (an orbit of codimension 1 in Z: one invariant level set).
Usage: python3 analysis/v2_74_orbitclosure.py CERT BLOCK COMP OUT [NPTS]"""
import os, sys, json, re
HERE = os.path.dirname(os.path.abspath(__file__))
cert, block, comp, out = sys.argv[1:5]; NP = int(sys.argv[5]) if len(sys.argv) > 5 else 2
sys.argv = ['v2_71_gcert.py', cert, block, os.path.join(os.path.dirname(out), 'empty.json'), out]
src = open(os.path.join(HERE, 'v2_71_gcert.py')).read(); src = src[:src.index("setup_primes = set()")]
G2 = {'__name__': 'v2_74o', '__file__': os.path.join(HERE, 'v2_71_gcert.py')}
exec(compile(src, 'v2_71_gcert', 'exec'), G2)
g = G2; G_ = g['G_']; sp = g['sp']; vs = g['vs']; d = g['d']; Pr = g['Pr']; NCL = g['NCL']; gvars = g['gvars']; cl = g['cl']; big = g['big']
C = [c for c in G_['comps'] if str(c['comp']) == comp][0]
gens = [x.replace('**', '^') for x in C['gens'] if x.strip() not in ('', '0')]
_, iv, _ = G_['info'](gens, 'oi')
Wint = [[1 if i == j else 0 for i in range(NCL)] for j in range(NCL)]
ts = sp.symbols(f't0:{NCL}'); us = sp.symbols(f'u0:{NCL}')
rat, _ = G_['points'](gens, iv, 4)
for pi, pt in enumerate(rat[:NP]):
    Px = sp.Matrix(d, d, lambda a, b: G_['evalr'](Pr[a][b], pt) if Pr[a][b] else 0)
    Sg = g['sigma'](g['Gm'], g['scal_of'](Wint, ts))
    at1 = {v: 0 for v in gvars}
    for c in big:
        for i in range(len(cl[c]) - 1): at1[sp.Symbol(f"g{c}_{i}{i}")] = 1
    at1.update({v: 1 for v in list(ts) + list(us)})
    hv0 = list(ts) + gvars
    Tm = sp.Matrix([[g['coord'](Dm * Px - Px * Dm, w) for w in vs] for Dm in [sp.diff(Sg, v).subs({v_: at1[v_] for v_ in hv0}) for v in hv0]]).T
    piv = Tm.rref()[1]
    hv1, yv, eqs = g['g1_graph']([hv0[k] for k in piv], Px, Sg, at1, ts, us, Wint)
    o = g['run_sing']([f"ring RQ = 0,({','.join(map(str, hv1 + yv))}),(dp({len(hv1)}),dp({len(yv)}));", "option(redSB);",
                       f"ideal I = {','.join(eqs)};", "ideal E = eliminate(I, " + '*'.join(map(str, hv1)) + ");",
                       f"ideal IZ = std(ideal({','.join(gens)})); ideal NE = reduce(std(E), IZ); int i;",
                       'for (i = 1; i <= ncols(NE); i++) { if (NE[i] != 0) { print("EXTRA " + string(NE[i])); } }',
                       'print("DIMS " + string(dim(std(E))) + " " + string(dim(IZ)));', 'print("END"); quit;'], f'oc{pi}')
    x = {v: str(pt[G_['gpos'][v]]) for v in vs if pt[G_['gpos'][v]] != 0}
    print(f"point {pi}: orbit dim {len(piv)}; {re.search(r'DIMS (.*)', o).group(1) if o and o not in ('TIMEOUT',) else o}; nonzero coords {x}")
    for m in re.finditer(r"^EXTRA (.*)$", o or '', re.M): print("   extra:", m.group(1)[:300])
