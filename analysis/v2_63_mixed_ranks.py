"""v2.63: ranks of the mixed-kind folds in nine leaf extensions (census representatives, q = 3)."""
import os, sys, collections
pairs = [("K2","pair_after"),("K2","pair_srcout"),("K2","pair_tgtin"),("pair_after","pair_after2"),("pair_after","pair_tgtfork"),("pair_srcout","pair_srcboth"),("K3","K3_after"),("K3","K3_srcout"),("pairs_in_short","pairs_in_short_tail")]
HERE = '/home/claude/eat-the-fibre/analysis'; sys.path.insert(0, os.path.dirname(HERE)); sys.path.insert(0, HERE)
import numpy as np
src = open(os.path.join(HERE, 'v2_60_leaf_pairs.py')).read().split("key0 = {}")[0]
for Qn, Qp in pairs:
    sys.argv = ['x', Qn, Qp, '3', '/home/claude/eat-the-fibre/logs/v2_60/reps']
    G = {'__file__': os.path.join(HERE, 'v2_60_leaf_pairs.py'), '__name__': 'x'}
    exec(src, G)
    R1, Qx, Jx, c1, q = G['R1'], G['Qx'], G['Jx'], G['c1'], G['q']
    ranks = collections.defaultdict(set); allr = collections.Counter()
    for P, lab in zip(R1['greps'], R1['loc']):
        P = P.astype(np.int64) % q
        J = not P[np.ix_(Qx, Jx)].any(); Qs = not P[np.ix_(Jx, Qx)].any()
        if J or Qs: continue
        ranks[lab].add(c1['F'].rank(P))
    print(Qn, '→', Qp, ': neither classes', len(ranks), '; ranks', dict(collections.Counter(frozenset(v) for v in ranks.values())))
