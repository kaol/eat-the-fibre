"""v2.44: for the line families in a lines log (P(t) = P₀ + t·N), the merged orbit and (b, e) class of each point
t ∈ 𝔽_q, and the overlap between lines. Usage: python3 analysis/v2_44_lineclass.py REPS_p.pkl REP_INDEX LINES_LOG REPS_q.pkl"""
import os, sys, io, re, ast, contextlib
import numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
rf, ri, logf, qf = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
def load_fc(f):
    sys.argv = ['family_curves.py', f]
    g = {'__name__': 'fc', '__file__': os.path.join(HERE, 'family_curves.py')}
    with contextlib.redirect_stdout(io.StringIO()):
        exec(open(os.path.join(HERE, 'family_curves.py')).read().split("fams = []; nlift")[0], g)
    return g
g = load_fc(rf); gq = load_fc(qf); q, d = gq['p'], g['d']
sys.argv = ['fibre_rigidity.py']
fr = {'__name__': 'fr', '__file__': os.path.join(HERE, 'fibre_rigidity.py')}
exec(open(os.path.join(HERE, 'fibre_rigidity.py')).read().split("\ndef main():")[0], fr)
C = fr['setup'](gq['r']['name'], q)
P0 = np.array(g['lifted'](ri), dtype=np.int64)
Ns = [ast.literal_eval(m.group(1)) for m in re.finditer(r"N \(nonzero entries\): (\{.*\})", open(logf).read())]
sets = []
for k, Nd in enumerate(Ns):
    N = np.zeros((d, d), np.int64)
    for (i, j), v in Nd.items(): N[i, j] = int(sp.Rational(v))
    pts = []
    for s in range(q):
        P = (P0 + s * N) % q
        l = gq['key2loc'][gq['key'](P)]; a = fr['at'](P, C, q)
        pts.append((s, l, (a['b'], a['e'])))
    sets.append({l for _, l, _ in pts})
    print(f"line #{k} over F_{q}: " + ", ".join(f"t={s}: orbit {l} {be}" for s, l, be in pts), flush=True)
for a in range(len(sets)):
    for b in range(a + 1, len(sets)):
        print(f"lines #{a}, #{b} share orbits: {sorted(sets[a] & sets[b])}")
