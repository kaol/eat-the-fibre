"""v2.73: deduplicate the certified char-P nodes of a case-split block across its charts, in full coordinates mod P.
Phase 'export' (one chart): full coordinates mod P of the point of each certified node (coord(P_chart(x), name)).
Phase 'test' (one chart): for every exported point z of the other charts, and each certified node N of this chart of the
same dimension: z lies on N iff z is on the chart's graph (z_n = coord(P_chart(z|vs), n) mod P for every name) and the
generators of N (over 𝔽_P, chart coordinates) vanish at z|vs.
Phase 'count': union of mutual memberships (two certified nodes with the same closure: each point lies on the other's
closure, both being orbit closures of equal dimension) together with the fixed point leaves; prints the count.
Usage: python3 analysis/v2_73_dedupe_p.py export|test CERT BLOCK LEAF CHARPCERT.json CHARPTREE.json P OUTDIR
       python3 analysis/v2_73_dedupe_p.py count OUTDIR P"""
import os, sys, json, glob
HERE = os.path.dirname(os.path.abspath(__file__))
mode = sys.argv[1]
if mode == 'count':
    outd, P = sys.argv[2], int(sys.argv[3])
    pts = {}; tests = []
    for f in glob.glob(f'{outd}/dd_export_*_p{P}.json'): pts.update(json.load(open(f)))
    for f in glob.glob(f'{outd}/dd_test_*_p{P}.json'): tests += json.load(open(f))
    on = {(a, b) for a, b in tests}                         # point of node a lies on node b
    par = {k: k for k in pts}
    def find(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for a, b in on:
        if (b, a) in on and a in par and b in par: par[find(a)] = find(b)
    classes = {find(k) for k in pts}
    fixed = []
    for f in glob.glob(f'{outd}/leaf*.json'):
        L = json.load(open(f))
        if L['vs'] == 0 and L['ncomp'] and L['checks'].get('P_G_fixes') and L['checks'].get('P_W_fixes'): fixed.append(L)
    # a fixed point leaf is one more class unless it equals an exported point node (same coordinates mod P)
    import sympy as sp
    extra = 0
    for L in fixed:
        z = {n: int(sp.fraction(sp.Rational(v))[0]) * pow(int(sp.fraction(sp.Rational(v))[1]), -1, P) % P for n, v in L['point_full'].items()}
        if not any(pts[k]['dim'] == 0 and all(pts[k]['z'].get(n) == z[n] for n in z) for k in pts): extra += 1
    print(f"COUNT p = {P}: {len(classes)} certified orbit closures across charts ({len(pts)} certified nodes) + {extra} fixed point leaves = {len(classes) + extra}")
    json.dump(dict(count=len(classes) + extra, nodes=len(pts), classes=len(classes), fixed_extra=extra), open(f'{outd}/dedupe_p{P}.json', 'w'))
    sys.exit(0)
cert, block, leaf, cj, tj, P, outd = sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], int(sys.argv[7]), sys.argv[8]
os.environ['STAB_LEAF'] = leaf
sys.argv = ['v2_71_gcert.py', cert, block, tj, outd + f'/dd_{mode}_{leaf}']
src = open(os.path.join(HERE, 'v2_71_gcert.py')).read(); src = src[:src.index("setup_primes = set()")]
G2 = {'__name__': 'v2_73d', '__file__': os.path.join(HERE, 'v2_71_gcert.py')}
exec(compile(src, 'v2_71_gcert', 'exec'), G2)
G_ = G2['G_']; sp = G2['sp']; R = G2['R']; vs = G2['vs']; d = G2['d']; Pr = G2['Pr']; coord = G2['coord']; names = G_['S']['names']
red = lambda q: int(sp.fraction(sp.Rational(q))[0]) * pow(int(sp.fraction(sp.Rational(q))[1]), -1, P) % P
C = json.load(open(cj)); T = json.load(open(tj))
key_of = {nd['path']: k for k, nd in T['nodes'].items()}
tag = os.path.basename(cj).replace('.json', '')
def Pmat(vals):
    pt = [0] * len(R.gens)
    for v in vs: pt[G_['gpos'][v]] = sp.Integer(vals[v])
    return sp.Matrix(d, d, lambda a, b: red(G_['evalr'](Pr[a][b], pt)) if Pr[a][b] else 0)
if mode == 'export':
    out = {}
    for n in C['nodes']:
        if n['result'] == 'point certified':               # a 0-dimensional node: its point from the (linear) key mod P
            sol = {}
            for g in key_of[n['path']].split(','):
                e_ = sp.sympify(g.replace('^', '**')); v_ = list(e_.free_symbols)
                if len(v_) == 1: c1 = int(e_.coeff(v_[0])) % P; sol[str(v_[0])] = (-int(e_.subs(v_[0], 0)) * pow(c1, -1, P)) % P
            for v in vs: sol.setdefault(v, 0)
            n = dict(n, point=sol, dim=0)
        elif n['result'] != 'certified' or 'point' not in n: continue
        Px = Pmat(n['point'])
        out[f"{tag}:{n['path']}"] = dict(dim=n['dim'], z={nm: int(coord(Px, nm)) % P for nm in names})
    json.dump(out, open(f'{outd}/dd_export_{tag}_p{P}.json', 'w'))
    print(f"EXPORT {tag}: {len(out)} points")
else:
    pts = {}
    for f in glob.glob(f'{outd}/dd_export_*_p{P}.json'):
        if os.path.basename(f) != f'dd_export_{tag}_p{P}.json': pts.update(json.load(open(f)))
    mine = [n for n in C['nodes'] if n['result'] in ('certified', 'point certified')]
    found = []
    for k, z in pts.items():
        if not any(n['dim'] == z['dim'] for n in mine): continue
        vals = {v: z['z'][v] for v in vs}
        Px = Pmat(vals)
        if any((int(coord(Px, nm)) - z['z'][nm]) % P for nm in names): continue      # not on this chart's graph
        for n in mine:
            if n['dim'] != z['dim']: continue
            gens = [sp.sympify(g.replace('^', '**')) for g in key_of[n['path']].split(',') if g.strip() not in ('', '0')]
            if all(int(g.subs({sp.Symbol(v): vals[v] for v in vs})) % P == 0 for g in gens):
                found.append((k, f"{tag}:{n['path']}"))
    json.dump(found, open(f'{outd}/dd_test_{tag}_p{P}.json', 'w'))
    print(f"TEST {tag}: {len(found)} memberships")
