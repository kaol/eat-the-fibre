"""v2.36: print the closed forms P_t of the 19 families (nonzero columns, readable).
Usage: python3 analysis/v2_36_show.py NF.pkl REPS_two_pairs_7.pkl [family numbers]"""
import os, sys, pickle, numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
nf, r7f = sys.argv[1:3]; which = set(map(int, sys.argv[3:]))
sys.argv = ['x', r7f]
src = open(os.path.join(HERE, 'family_normal_forms.py')).read().split("# ---- H = G × T² as matrices")[0]
exec(src)
src2 = open(os.path.join(HERE, 'family_normal_forms.py')).read()
exec(src2[src2.index("# ---- inner curves"):src2.index("def is_fold_family")])
fams = pickle.load(open(nf, 'rb')); t = sp.Symbol('t'); k = 0
for tag in sorted(fams):
    for locs, (i, xn, _) in sorted(fams[tag].items()):
        k += 1
        if which and k not in which: continue
        Pz = ((r['greps'][i].astype(np.int64) % 7 + 3) % 7) - 3
        Inn, Innv = INN[xn]
        Ps = (Inn * sp.Matrix(Pz.tolist()) * Innv).applyfunc(sp.expand).subs(s_, t)
        print(f"family {k}: {tag}, curve along {xn}")
        for j in range(d):
            col = Ps[:, j]
            if all(c == 0 for c in col): continue
            terms = [f"({sp.factor(c)})·{names[q]}" if c not in (1, -1) else ('' if c == 1 else '−') + names[q]
                     for q, c in enumerate(col) if c != 0]
            print(f"   P({names[j]}) = " + " + ".join(terms))
