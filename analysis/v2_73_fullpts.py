"""v2.73: full coordinates of the certified points of a chart's tree (gcert JSON), via the chart's coordinate functions.
Usage: python3 analysis/v2_73_fullpts.py CERT BLOCK LEAF GC.json OUT.json"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
cert, block, leaf, gcj, outj = sys.argv[1:6]
os.environ['STAB_LEAF'] = leaf; os.environ['NO_POINT_EXIT'] = '1'
sys.argv = ['v2_71_gcert.py', cert, block, gcj, outj + '_x']
src = open(os.path.join(HERE, 'v2_71_gcert.py')).read(); src = src[:src.index("setup_primes = set()")]
G2 = {'__name__': 'v2_73f', '__file__': os.path.join(HERE, 'v2_71_gcert.py')}
exec(compile(src, 'v2_71_gcert', 'exec'), G2)
G_ = G2['G_']; sp = G2['sp']; R = G2['R']; vs = G2['vs']; d = G2['d']; Pr = G2['Pr']; coord = G2['coord']; names = G_['S']['names']
gc = json.load(open(gcj)); out = {}
for n in gc['nodes']:
    pts = []
    for p in n.get('points', []):
        if 'x' not in p: pts.append(None); continue
        pt = [0] * len(R.gens)
        for v, x in zip(vs, p['x']): pt[G_['gpos'][v]] = sp.Rational(x)
        Px = sp.Matrix(d, d, lambda a, b: G_['evalr'](Pr[a][b], pt) if Pr[a][b] else 0)
        pts.append({nm: str(coord(Px, nm)) for nm in names})
    out[n['path']] = pts
json.dump(out, open(outj, 'w'))
print(f"FULLPTS leaf {leaf}: {sum(len(v) for v in out.values())} points")
