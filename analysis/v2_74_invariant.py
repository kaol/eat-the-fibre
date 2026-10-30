"""v2.74: certify an H-invariant rational function φ = N/D on components of a block (char 0 and over ℤ): for every field
ξ of 𝔥 = Der_Δ + all class weightings (integral rows), ξ(N)·D − N·ξ(D) ∈ I_Z (exact reduction over ℚ, and over ℤ by a
strong basis: then the identity holds mod every p prime to its leading coefficients); D ∉ I_Z; φ non-constant on Z.
Usage: python3 analysis/v2_74_invariant.py CERT BLOCK COMPS N D OUT"""
import os, sys, json, re, math
HERE = os.path.dirname(os.path.abspath(__file__))
cert, block, comps, N, D, out = sys.argv[1:7]
sys.argv = ['v2_67_orbits.py', cert, block, '0', comps, out + '_o', '--tmo', '6000']
src = open(os.path.join(HERE, 'v2_67_orbits.py')).read(); src = src[:src.index("stats = {'primes'")]
G_ = {'__name__': 'v2_74i', '__file__': os.path.join(HERE, 'v2_67_orbits.py')}
exec(compile(src, 'v2_67_orbits', 'exec'), G_)
sp = G_['sp']; ss = G_['ss']; NCL = G_['NCL']; vs = G_['vs']; R = G_['R']
Wb = [[sp.Integer(1 if i == j else 0) for i in range(NCL)] for j in range(NCL)]
act = G_['active'](Wb); XI = G_['combine'](act, G_['XIB'])
def lcm(a, b): return a * b // math.gcd(a, b)
sc = []
for row in XI:
    den = 1
    for e in row:
        if e:
            for c in e.coeffs(): den = lcm(den, int(c.denominator))
    sc.append(den)
XI = [[e * R(sc[a]) if e else e for e in XI[a]] for a in range(len(XI))]
NH = len(XI); nv = len(vs)
for c in G_['comps']:
    gens = [g.replace('**', '^') for g in c['gens'] if g.strip() not in ('', '0')]
    head = [f"matrix XI[{NH}][{nv}] = " + ",".join(ss(XI[a][b]) for a in range(NH) for b in range(nv)) + ";",
            f"proc xi(poly h, int a) {{ poly g = 0; int b; for (b = 1; b <= {nv}; b++) {{ g = g + diff(h, var(b)) * XI[a, b]; }} return(g); }}",
            f"poly NN = {N}; poly DD = {D}; int a; ideal W;", f"for (a = 1; a <= {NH}; a++) {{ W = W + (xi(NN, a) * DD - NN * xi(DD, a)); }}"]
    o = G_['singular']([f"ring R = 0,({','.join(vs)}),dp;", "option(redSB);", f"ideal Q = std(ideal({','.join(gens)}));"] + head +
                       ['print("INV " + string(size(reduce(W, Q)) == 0));', 'print("DNOTIN " + string(reduce(DD, Q) != 0));',
                        'print("NNOTIN " + string(reduce(NN, Q) != 0));', 'print("END"); quit;'], f"iq{c['comp']}")
    oz = G_['singular']([f"ring R = integer,({','.join(vs)}),dp;", "option(redSB);", f"ideal Q = std(ideal({','.join(gens)}));",
                         'string s = "LC"; int i; for (i = 1; i <= size(Q); i++) { s = s + " " + string(leadcoef(Q[i])); } print(s);'] + head +
                        ['print("INVZ " + string(size(reduce(W, Q)) == 0));', 'print("END"); quit;'], f"iz{c['comp']}")
    lc = sorted({abs(int(x)) for x in re.search(r"LC(.*)", oz).group(1).split()} - {1}) if oz else None
    print(f"component {c['comp']}: φ = ({N})/({D}); ξ-invariant mod I_Z over ℚ: {re.search(r'INV (\d)', o).group(1) == '1'}; "
          f"D ∉ I_Z: {re.search(r'DNOTIN (\d)', o).group(1) == '1'}; N ∉ I_Z: {re.search(r'NNOTIN (\d)', o).group(1) == '1'}; "
          f"over ℤ: {re.search(r'INVZ (\d)', oz).group(1) == '1' if oz else None}, strong-basis leading coefficients ≠ ±1: {lc}", flush=True)
