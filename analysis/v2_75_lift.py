"""v2.75: integer certificates for the rigidity ideals of v2_75_rigid.py: Singular lift over ℤ of the target
(t − tp, or 1) from the entries E of h·P − P'·h; verified in sympy (Σ cof·E − target = 0).
Usage: python3 analysis/v2_75_lift.py OUTDIR"""
import re, subprocess, sys, json, sympy as sp
outd = sys.argv[1]; out = {}
for tag, target in (('src', 't-tp'), ('others', 't-tp'), ('src_vs_others', '1')):
    s = open(f'{outd}/rigid_{tag}.sing').read().split('\n')
    lines = [s[0].replace('ring r = 0', 'ring r = integer'), s[1], s[2], f"matrix T = lift(E, ideal({target}));",
             'int i; for (i = 1; i <= nrows(T); i++) { if (T[i,1] != 0) { print("COF " + string(i) + " " + string(T[i,1])); } }',
             'print("END"); quit;']
    open(f'{outd}/lift_{tag}.sing', 'w').write("\n".join(lines))
    o = subprocess.run(['Singular', '-q', '--no-rc', f'{outd}/lift_{tag}.sing'], capture_output=True, text=True).stdout
    cofs = {int(i): sp.sympify(c.replace('^', '**')) for i, c in re.findall(r'^COF (\d+) (.*)$', o, re.M)}
    E = [sp.sympify(e.replace('^', '**')) for e in re.search(r'ideal E = (.*);', s[2]).group(1).split(',')]
    res = sp.expand(sum(c * E[i - 1] for i, c in cofs.items()) - sp.sympify(target))
    integral = all(sp.Poly(c, *sorted(c.free_symbols, key=str)).coeffs() == [int(x) for x in sp.Poly(c, *sorted(c.free_symbols, key=str)).coeffs()] if c.free_symbols else c == int(c) for c in cofs.values())
    print(tag, 'target', target, {i: (str(c), str(E[i - 1])) for i, c in cofs.items()}, 'residual', res, 'integral', integral)
    out[tag] = dict(target=target, cofactors={str(i): str(c) for i, c in cofs.items()}, entries={str(i): str(E[i - 1]) for i in cofs},
                    residual=str(res), integral=integral)
json.dump(out, open(f'{outd}/certificates.json', 'w'), indent=1)
