"""v2.74: the exact generic rank of the orbit matrix on components of a block (char 0): fraction-free Gaussian
elimination over R/I_Z (a domain), with 𝔥 = Der_Δ + the weightings W (default: all class weightings).
rank < dim Z at the generic point ⇒ every orbit in Z has dimension < dim Z ⇒ infinitely many orbits on Z.
Usage: python3 analysis/v2_74_genrank.py CERT BLOCK COMPS OUT"""
import os, sys, json, re
HERE = os.path.dirname(os.path.abspath(__file__))
cert, block, comps, out = sys.argv[1:5]
sys.argv = ['v2_67_orbits.py', cert, block, '0', comps, out + '_o', '--tmo', '6000']
src = open(os.path.join(HERE, 'v2_67_orbits.py')).read(); src = src[:src.index("stats = {'primes'")]
G_ = {'__name__': 'v2_74g', '__file__': os.path.join(HERE, 'v2_67_orbits.py')}
exec(compile(src, 'v2_67_orbits', 'exec'), G_)
sp = G_['sp']; ss = G_['ss']; NCL = G_['NCL']; vs = G_['vs']
Wb = [[sp.Integer(1 if i == j else 0) for i in range(NCL)] for j in range(NCL)]
act = G_['active'](Wb); MR = G_['combine'](act, G_['MB']); NH = len(MR); NC = len(MR[0])
for c in G_['comps']:
    gens = [g.replace('**', '^') for g in c['gens'] if g.strip() not in ('', '0')]
    lines = [f"ring R = 0,({','.join(vs)}),dp;", "option(redSB);", f"ideal Q = std(ideal({','.join(gens)}));",
             f"matrix M[{NH}][{NC}] = " + ",".join(ss(MR[a][b]) for a in range(NH) for b in range(NC)) + ";",
             "M = reduce(M, Q); int r = 0; int i; int j; int i2; int j2;",
             f"int nr = {NH}; int nc = {NC}; intvec usedr = 0:nr; intvec usedc = 0:nc; int found;",
             "for (j = 1; j <= nc; j++) { found = 0; for (i = 1; i <= nr && found == 0; i++) { if (usedr[i] == 0 && M[i, j] != 0) { found = i; } }",
             "  if (found > 0) { r++; usedr[found] = 1; usedc[j] = 1;",
             "    for (i2 = 1; i2 <= nr; i2++) { if (usedr[i2] == 0) { for (j2 = 1; j2 <= nc; j2++) { if (usedc[j2] == 0) {",
             "      M[i2, j2] = reduce(M[found, j] * M[i2, j2] - M[i2, j] * M[found, j2], Q); } } } } } }",
             'print("GENRANK " + string(r) + " DIM " + string(dim(Q)));', 'print("END"); quit;']
    o = G_['singular'](lines, f"gr{c['comp']}")
    m = re.search(r"GENRANK (\d+) DIM (-?\d+)", o or '')
    print(f"component {c['comp']}: dim {c['dim']}, generic orbit rank {m.group(1) if m else '?'} (Singular dim {m.group(2) if m else '?'}), 𝔥 = Der_Δ ({len(G_['DerD'])}) + {NCL} weightings", flush=True)
