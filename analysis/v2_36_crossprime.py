"""v2.36: the 19 closed forms found over 𝔽₇, reduced mod another prime, located in that prime's census.
Usage: python3 analysis/v2_36_crossprime.py NF.pkl REPS_two_pairs_7.pkl REPS_two_pairs_q.pkl"""
import os, sys, pickle, collections
import numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
nf, r7f, rqf = sys.argv[1:4]
sys.argv = ['x', rqf]
src = open(os.path.join(HERE, 'family_normal_forms.py')).read().split("growing = {}")[0]
src = src.replace("_NF_FILE = sys.argv[1]; sys.argv = ['x']", "_NF_FILE = sys.argv[1]; sys.argv = ['x']; __file__ = os.path.join(HERE, 'family_normal_forms.py')")
exec(src)                                   # census of the target prime: key(), key2loc, sig_of_loc, INN, s_
r7 = pickle.load(open(r7f, 'rb')); fams = pickle.load(open(nf, 'rb'))
allhit = []
for tag in sorted(fams):
    for locs, (i, xn, _) in sorted(fams[tag].items()):
        Pz = ((r7['greps'][i].astype(np.int64) % 7 + 3) % 7) - 3
        Inn, Innv = INN[xn]
        Ps = (Inn * sp.Matrix(Pz.tolist()) * Innv).applyfunc(sp.expand)
        hit = set()
        for sv in range(p):
            Pn = np.array(Ps.subs(s_, sv).tolist(), dtype=np.int64) % p
            l = key2loc.get(key(Pn))
            assert l is not None, "family member is not a census fold"
            hit.add(l)
        sigs = {sig_of_loc[l][0] for l in hit}
        allhit.append(hit)
        print(f"{tag[:40]:40} x={xn:5} meets {len(hit)} local orbits over 𝔽_{p}; classes {sigs}")
inter = sum(len(a & b) for k, a in enumerate(allhit) for b in allhit[k + 1:])
print(f"pairwise shared local orbits among the {len(allhit)} families over 𝔽_{p}: {inter}")
