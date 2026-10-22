"""v2.44: union of the merged orbits met by all line families in driver logs, per prime, broken down by signature
group and (b, e) class. Usage: python3 analysis/v2_44_coverage.py REPS_p.pkl LOG... -- REPS_q.pkl..."""
import os, sys, io, re, ast, contextlib, collections
import numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
args = sys.argv[1:]; k = args.index('--'); rf, logs, files = args[0], args[1:k], args[k + 1:]
def load_fc(f):
    sys.argv = ['family_curves.py', f]
    g = {'__name__': 'fc', '__file__': os.path.join(HERE, 'family_curves.py')}
    with contextlib.redirect_stdout(io.StringIO()):
        exec(open(os.path.join(HERE, 'family_curves.py')).read().split("fams = []; nlift")[0], g)
    return g
g = load_fc(rf); d = g['d']
sys.argv = ['fibre_rigidity.py']
fr = {'__name__': 'fr', '__file__': os.path.join(HERE, 'fibre_rigidity.py')}
exec(open(os.path.join(HERE, 'fibre_rigidity.py')).read().split("\ndef main():")[0], fr)
lines = []
for lf in logs:
    for m in re.finditer(r"from orbit \d+ rep (\d+), .*?N = (\{.*\})", open(lf).read()):
        N = np.full((d, d), sp.Integer(0), dtype=object)
        for (i, j), v in ast.literal_eval(m.group(2)).items(): N[i, j] = sp.Rational(v)
        lines.append((np.array(g['lifted'](int(m.group(1))), dtype=object), N))
print(f"{len(lines)} lines", flush=True)
for f in files:
    gq = load_fc(f); q = gq['p']; C = fr['setup'](gq['r']['name'], q); sig = gq['sig_of_loc']
    U = set()
    for P0, N in lines:
        Nq = np.vectorize(lambda z: int(z.p * pow(int(z.q), -1, q)) % q)(N).astype(np.int64)
        for s in range(q): U.add(gq['key2loc'][gq['key']((P0.astype(np.int64) + s * Nq) % q)])
    rep = {}
    for i, l in enumerate(gq['r']['loc']):
        if l in U and l not in rep: rep[l] = i
    cnt = collections.Counter()
    for l, i in rep.items():
        a = fr['at'](gq['r']['greps'][i], C, q)
        s = sig[l]; grp = ('L' if "'Uki'), (7" in str(s) else 'R' if "'Uik'), (7" in str(s) else 'other')
        cnt[(grp, (a['b'], a['e']))] += 1
    print(f"q={q}: lines meet {len(U)} orbits: {dict(sorted(cnt.items()))}", flush=True)
