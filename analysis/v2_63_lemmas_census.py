"""v2.63: Lemmas A and B (λ = 0 one-sided folds of a sink-leaf extension) checked on census representatives, q = 3."""
import os, sys, collections, itertools
pairs = [("K2","pair_after"),("K2","pair_srcout"),("pair_after","pair_after2"),("pair_after","pair_tgtfork"),("K3","K3_after"),("K3","K3_srcout"),("pairs_in_short","pairs_in_short_tail")]
HERE = '/home/claude/eat-the-fibre/analysis'; sys.path.insert(0, os.path.dirname(HERE)); sys.path.insert(0, HERE)
import numpy as np
src = open(os.path.join(HERE, 'v2_60_leaf_pairs.py')).read().split("key0 = {}")[0]
for Qn, Qp in pairs:
    sys.argv = ['x', Qn, Qp, '3', '/home/claude/eat-the-fibre/logs/v2_60/reps']
    G = {'__file__': os.path.join(HERE, 'v2_60_leaf_pairs.py'), '__name__': 'x'}
    exec(src, G)
    R1, Qx, Jx, c1, c0, q, ev = G['R1'], G['Qx'], G['Jx'], G['c1'], G['c0'], G['q'], G['ev']
    M0 = np.asarray(c0['M']); d0 = len(Qx); F = c1['F']
    K0 = np.asarray(c0['Kmat']) % q
    res = collections.Counter()
    for P, lab in zip(R1['greps'], R1['loc']):
        P = P.astype(np.int64) % q
        J = not P[np.ix_(Qx, Jx)].any(); Qs = not P[np.ix_(Jx, Qx)].any()
        if J == Qs: continue                                  # one-sided only
        lam = int(P[ev, ev])
        if lam != 0: continue
        Pb = P[np.ix_(Qx, Qx)]
        PJ0 = not P[np.ix_(Jx, Jx)].any()
        if J:   # Lemma A: xi = e_v row on kQ; xi∘P̄ = xi; xi multiplicative on im P̄
            xi = P[ev, Qx]
            inv = not ((xi @ Pb - xi) % q).any()
            mult = True
            for i, j in itertools.product(range(d0), repeat=2):
                a, b = Pb[:, i], Pb[:, j]
                ab = np.einsum('i,j,ijk->k', a, b, M0) % q
                if (xi @ ab - (xi @ a) * (xi @ b)) % q: mult = False; break
            res[('J-kind λ=0', 'P_J = 0' if PJ0 else 'P_J ≠ 0', 'ξ∘P̄ = ξ' if inv else 'ξ∘P̄ ≠ ξ', 'ξ multiplicative on im P̄' if mult else 'not multiplicative')] += 1
        else:   # Lemma B: y = P(e_v) ∈ kQ; −y idempotent, y ∈ im P̄ ∩ K
            y = P[Qx, ev]
            idem = not ((np.einsum('i,j,ijk->k', y, y, M0) + y) % q).any()        # y² = −y
            inim = not ((Pb @ y - y) % q).any()
            inK = F.rank(np.vstack([K0.T, y])) == F.rank(K0.T)
            YinIm = not ((Pb @ P[np.ix_(Qx, Jx)] - P[np.ix_(Qx, Jx)]) % q).any()
            res[('Q-kind λ=0', 'P_J = 0' if PJ0 else 'P_J ≠ 0', '−y idempotent' if idem else 'not', 'y ∈ im P̄' if inim else 'y ∉ im', 'y ∈ K' if inK else 'y ∉ K', 'Y(J) ⊆ im P̄' if YinIm else 'not')] += 1
    print(Qn, '→', Qp)
    for k, v in sorted(res.items()): print('   ', k, v)
