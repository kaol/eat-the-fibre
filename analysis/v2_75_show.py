"""v2.75: show rational points of a component Z of a vertex block as folds: P in the path basis, a basis of Im P and
of Ker P (reduced row echelon form, labelled by paths), and the values of given coordinate functions.
With SHOW_PHI='A0;X;X2;FORMULA': also the intrinsic φ (v2_75_phi) and FORMULA in stratum coordinates.
Usage: [SHOW_PHI=...] python3 analysis/v2_75_show.py CERT BLOCK COMP OUTPREFIX [NPTS]"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
cert, block, comp, out = sys.argv[1:5]; NP = int(sys.argv[5]) if len(sys.argv) > 5 else 2
PHI = os.environ.get('SHOW_PHI', '').split(';') if os.environ.get('SHOW_PHI') else None   # A0;X;X2;FORMULA
if PHI: PHI = [int(PHI[0]), int(PHI[1]), int(PHI[2]), PHI[3]]
os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
emp = os.path.join(os.path.dirname(out) or '.', 'empty.json')
if not os.path.exists(emp): open(emp, 'w').write('{"nodes": {}}')
sys.argv = ['v2_71_gcert.py', cert, block, emp, out]
src = open(os.path.join(HERE, 'v2_71_gcert.py')).read(); src = src[:src.index("setup_primes = set()")]
g = {'__name__': 'v2_75s', '__file__': os.path.join(HERE, 'v2_71_gcert.py')}
exec(compile(src, 'v2_71_gcert', 'exec'), g)
G_ = g['G_']; sp = g['sp']; vs = g['vs']; d = g['d']; Pr = g['Pr']
basis = G_['S']['ctx']['basis']; arrows = G_['S']['ctx']['arrows']
def lab(b):
    if b[0] == 'v': return f"e{b[1]}"
    return '·'.join(f"a{i}" for i in b[1])
L = [lab(b) for b in basis]
C = [c for c in G_['comps'] if str(c['comp']) == comp][0]
gens = [x.replace('**', '^') for x in C['gens'] if x.strip() not in ('', '0')]
_, iv, _ = G_['info'](gens, 'sh')
rat, _ = G_['points'](gens, iv, max(NP, 2))
print("arrows:", {f"a{i}": a for i, a in enumerate(arrows)})
def vecstr(v): return ' + '.join(f"({c})·{L[i]}" if c != 1 else L[i] for i, c in enumerate(v) if c != 0) or '0'
for pi, pt in enumerate(rat[:NP]):
    Px = sp.Matrix(d, d, lambda a, b: G_['evalr'](Pr[a][b], pt) if Pr[a][b] else 0)
    assert Px * Px == Px
    print(f"--- point {pi}: rank P = {Px.rank()}")
    for c in range(d):
        col = Px[:, c]
        if any(col[r] != (1 if r == c else 0) for r in range(d)): print(f"  P({L[c]}) = {vecstr(list(col))}")
    im = Px.T.rref()[0]; ke = (sp.eye(d) - Px).T.rref()[0]
    print("  Im P :", ' | '.join(vecstr(list(im[r, :])) for r in range(im.rows) if any(im[r, :])))
    print("  Ker P:", ' | '.join(vecstr(list(ke[r, :])) for r in range(ke.rows) if any(ke[r, :])))
    if PHI:                                             # v2.75: the intrinsic φ against a stratum-coordinate formula
        from v2_75_phi import phi
        cls_arrows = [a for a, st in enumerate(arrows) if st == arrows[PHI[0]]]
        f_, err = phi(Px, basis, cls_arrows, PHI[1], PHI[2])
        val = sp.sympify(PHI[3], locals={v: g['coord'](Px, v) for v in vs})
        print(f"  φ intrinsic = {f_} {err or ''}; formula {PHI[3]} = {val}; equal: {f_ is not None and sp.simplify(f_ - val) == 0}")
