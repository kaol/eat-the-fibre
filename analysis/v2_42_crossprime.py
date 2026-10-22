"""v2.42: closed forms found by analysis/family_curves.py over 𝔽_p (from its log), lifted to ℤ, reduced mod another
prime q and located in that prime's census (merged orbit labels).
Usage: python3 analysis/v2_42_crossprime.py CURVES_LOG REPS_p.pkl REPS_q.pkl
"""
import os, sys, re, ast, pickle, io, contextlib
import numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
logf, rpf, rqf = sys.argv[1:4]
fams = [(int(m.group(2)), ast.literal_eval(m.group(3))) for m in
        re.finditer(r"family: (\d+) local orbits, from rep (\d+) along (\([^)]*\))", open(logf).read())]
rp = pickle.load(open(rpf, 'rb')); p0 = rp['p']
sys.argv = ['family_curves.py', rqf]
g = {'__name__': 'fc', '__file__': os.path.join(HERE, 'family_curves.py')}
src = open(os.path.join(HERE, 'family_curves.py')).read().split("fams = []; nlift")[0]
with contextlib.redirect_stdout(io.StringIO()):
    exec(src, g)
q, key, key2loc, INN, t, sig_of_loc = g['p'], g['key'], g['key2loc'], g['INN'], g['t'], g['sig_of_loc']
allhit = []
for i, vec in fams:
    Pz = ((rp['greps'][i].astype(np.int64) % p0 + p0 // 2) % p0) - p0 // 2
    Inn, Innv = INN[vec]
    Ps = (Inn * sp.Matrix(Pz.tolist()) * Innv).applyfunc(sp.expand)
    hit = set()
    for sv in range(q):
        Pn = np.array(Ps.subs(t, sv).tolist(), dtype=np.int64) % q
        l = key2loc.get(key(Pn))
        assert l is not None, "family member is not a census fold"
        hit.add(l)
    allhit.append(hit)
    print(f"rep {i} along {vec}: meets {len(hit)} orbits over 𝔽_{q}; classes {sorted({sig_of_loc[l][0] for l in hit})}")
inter = sum(len(a & b) for k_, a in enumerate(allhit) for b in allhit[k_ + 1:])
n = len(set(g['r']['loc']))
print(f"{len(allhit)} families over 𝔽_{q}: pairwise shared orbits {inter}; orbits met {len(set().union(*allhit))}; "
      f"census orbits {n}; remaining {n - len(set().union(*allhit))}")
