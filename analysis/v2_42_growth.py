"""v2.42: where is the growth not explained by the families found? Counts of merged orbits per v2.35 signature, at two
primes, with the orbits met by the given genuine families removed.
Usage: python3 analysis/v2_42_growth.py CURVES_LOG REPS_p.pkl REPS_q.pkl REPS_q2.pkl  (families = curves meeting q orbits at every prime)"""
import os, sys, re, ast, pickle, io, contextlib, collections
import numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
logf, rpf = sys.argv[1:3]; others = sys.argv[3:]
fams = [(int(m.group(2)), ast.literal_eval(m.group(3))) for m in
        re.finditer(r"family: (\d+) local orbits, from rep (\d+) along (\([^)]*\))", open(logf).read())]
rp = pickle.load(open(rpf, 'rb')); p0 = rp['p']
res = {}
for rqf in others:
    sys.argv = ['family_curves.py', rqf]
    g = {'__name__': 'fc', '__file__': os.path.join(HERE, 'family_curves.py')}
    with contextlib.redirect_stdout(io.StringIO()):
        exec(open(os.path.join(HERE, 'family_curves.py')).read().split("fams = []; nlift")[0], g)
    q = g['p']; hits = []
    for i, vec in fams:
        Pz = ((rp['greps'][i].astype(np.int64) % p0 + p0 // 2) % p0) - p0 // 2
        Inn, Innv = g['INN'][vec]
        Ps = (Inn * sp.Matrix(Pz.tolist()) * Innv).applyfunc(sp.expand)
        hits.append({g['key2loc'][g['key'](np.array(Ps.subs(g['t'], sv).tolist(), dtype=np.int64) % q)] for sv in range(q)})
    res[q] = (g['sig_of_loc'], hits)
qs = sorted(res)
genuine = [k for k in range(len(fams)) if all(len(res[q][1][k]) == q for q in qs)]
print("genuine families (meet q orbits at every prime):", len(genuine), "of", len(fams))
cnt = {}
for q in qs:
    sig, hits = res[q]
    gone = set().union(*(hits[k] for k in genuine))
    cnt[q] = collections.Counter(str(s) for l, s in sig.items() if l not in gone)
    print(f"q={q}: remaining orbits {sum(cnt[q].values())}")
q0, q1 = qs[0], qs[-1]
for s in sorted(set(cnt[q0]) | set(cnt[q1])):
    if cnt[q0][s] != cnt[q1][s]: print(f"  grows: {s}: {cnt[q0][s]} -> {cnt[q1][s]}")
