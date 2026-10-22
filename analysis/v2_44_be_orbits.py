"""v2.44: v2.34's (b, e) per merged orbit (one representative each) in given signature groups, at several primes, and
at sample points of the inner-curve families from a family_curves log.
Usage: python3 analysis/v2_44_be_orbits.py CURVES_LOG REPS_p_of_log.pkl SIG_SUBSTRING REPS_q.pkl...
"""
import os, sys, re, ast, io, contextlib, collections, pickle, json
import numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
logf, rpf, tg = sys.argv[1:4]; files = sys.argv[4:]
def load_fc(f):
    sys.argv = ['family_curves.py', f]
    g = {'__name__': 'fc', '__file__': os.path.join(HERE, 'family_curves.py')}
    with contextlib.redirect_stdout(io.StringIO()):
        exec(open(os.path.join(HERE, 'family_curves.py')).read().split("fams = []; nlift")[0], g)
    return g
sys.argv = ['fibre_rigidity.py']
fr = {'__name__': 'fr', '__file__': os.path.join(HERE, 'fibre_rigidity.py')}
exec(open(os.path.join(HERE, 'fibre_rigidity.py')).read().split("\ndef main():")[0], fr)

# (8) sanity: points of the inner-curve families that meet q orbits (from the log), over the log's prime
gp = load_fc(rpf); p0 = gp['p']; Cp = fr['setup'](gp['r']['name'], p0)
fams = [(int(m.group(1)), int(m.group(2)), ast.literal_eval(m.group(3))) for m in
        re.finditer(r"family: (\d+) local orbits, from rep (\d+) along (\([^)]*\))", open(logf).read())]
hist = collections.Counter()
for n, i, vec in fams:
    if n != p0: continue
    Pz = gp['lifted'](i); Inn, Innv = gp['INN'][vec]
    Ps = (Inn * sp.Matrix(Pz.tolist()) * Innv).applyfunc(sp.expand)
    for tv in (1, 2):
        a = fr['at'](np.array(Ps.subs(gp['t'], tv).tolist(), dtype=np.int64) % p0, Cp, p0)
        hist[(a['b'], a['e'])] += 1
print(f"(8) inner-family points over F_{p0} (t = 1, 2): (b, e) histogram {dict(hist)}", flush=True)

# (9) per-orbit (b, e) in the target groups
res = {}
for f in files:
    g = load_fc(f); r, sig, p = g['r'], g['sig_of_loc'], g['p']
    C = fr['setup'](r['name'], p)
    seen = {}
    for i, l in enumerate(r['loc']):
        if l in seen or tg not in str(sig[l]): continue
        a = fr['at'](r['greps'][i], C, p); seen[l] = (a['b'], a['e'])
    res[p] = collections.Counter(seen.values())
    print(f"q={p}: {len(seen)} orbits; per (b, e): {dict(sorted(res[p].items()))}", flush=True)
ps = sorted(res)
print("growing (b, e) classes:", {k: [res[p][k] for p in ps] for k in sorted(set().union(*res.values()))
                                  if len({res[p][k] for p in ps}) > 1})
