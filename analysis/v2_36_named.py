"""v2.36: which of the 19 families are Theorems 2 and 4 (located through the census over 𝔽₇).
Usage: python3 analysis/v2_36_named.py NF.pkl REPS_two_pairs_7.pkl"""
import os, sys, pickle, re, numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
nf, r7f = sys.argv[1:3]; sys.argv = ['x', r7f]
src = open(os.path.join(HERE, 'family_normal_forms.py')).read().split("growing = {}")[0]
exec(src)
fams = pickle.load(open(nf, 'rb'))
famlist = [(tag, set(locs)) for tag in sorted(fams) for locs in sorted(fams[tag])]
gsrc = open(os.path.join(ROOT, 'lean', 'EatTheFibre', 'Gen', 'two_pairs.lean')).read()
def grab(nm):
    i = gsrc.index(f"def {nm}"); j = gsrc.index(":=", i) + 2; k = gsrc.find("\n\n", j)
    return eval(gsrc[j:k].split("/--")[0].replace("![", "[").strip())
for nm in ('thm2P', 'thm4P'):
    Pp = grab(nm); hits = set()
    for tv in range(p):
        Pn = np.array([[sum(a * tv ** n for n, a in enumerate(Pp[i][j])) % p for j in range(d)] for i in range(d)], np.int64)
        hits.add(key2loc.get(key(Pn)))
    owners = [k + 1 for k, (tag, f) in enumerate(famlist) if f & hits]
    print(f"{nm}: meets {len(hits - {None})} local orbits over 𝔽_{p}; shared with family {owners}; "
          f"orbits inside that family: {[len(hits & famlist[o - 1][1]) for o in owners]}")
