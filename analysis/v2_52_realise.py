"""v2.52: one integer fold of k(1⇉2⇇3 + 1→3) for each of its vertex blocks, from the census checkpoints.

For each vertex block found in the checkpoint of QUIVER over F_p, scan the folds with that block. Lift each
symmetrically to an integer matrix (standard basis: P = X·B⁻¹) and keep the first lift that is a fold over ℚ. The fold
check is exact (sympy integers): P² = P, image and kernel closed under products, and Σ_m P[m,l] DL[m] = P DL[l] Pᵀ for
every l. Writes JSON {block: P}.
Usage: python3 analysis/v2_52_realise.py QUIVER P OUT.json"""
import os, sys, json, pickle, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
import numpy as np
import pcensus as pc
name, p, outf = sys.argv[1], int(sys.argv[2]), sys.argv[3]
ctx = pc.setup(name, p)
d, basis, kdim = ctx['d'], ctx['basis'], ctx['kdim']
M = np.asarray(ctx['M'], dtype=object); DL = np.asarray(ctx['DL'], dtype=object)
Binv = np.asarray(ctx['Binv'], dtype=np.int64) % p
vid = [i for i, b in enumerate(basis) if b[0] == 'v']
def is_fold(P):
    """exact over ℤ (object arrays)"""
    I = np.eye(d, dtype=object)
    if ((P.dot(P) - P) != 0).any(): return False
    Q = I - P
    prod = lambda x, y: np.einsum('i,j,ijk->k', x, y, M)
    for a in range(d):
        for b in range(d):
            u = prod(P[:, a], P[:, b])
            if ((P.dot(u) - u) != 0).any(): return False
            w = prod(Q[:, a], Q[:, b])
            if (P.dot(w) != 0).any(): return False
    for l in range(d):
        lhs = sum(P[m, l] * DL[m] for m in range(d))
        if ((lhs - P.dot(DL[l]).dot(P.T)) != 0).any(): return False
    return True
ck = pickle.load(open(os.path.join(os.path.dirname(HERE), 'checkpoints', f'{name}_{p}.pkl'), 'rb'))
found, tried = {}, collections.Counter()
for r, fib in ck['fibres'].items():
    for X in fib:
        Pm = (np.asarray(X, dtype=np.int64) @ Binv) % p
        Ps = (Pm + p // 2) % p - p // 2
        V = tuple(map(tuple, Ps[np.ix_(vid, vid)].tolist()))
        if V in found: continue
        tried[V] += 1
        if tried[V] > 2000: continue
        Po = Ps.astype(object)
        if is_fold(Po): found[V] = Ps.tolist()
blocks = set(tried)
print(f"{name} F_{p}: {len(blocks)} vertex blocks; integer folds found for {len(found)}; "
      f"missing: {sorted(blocks - set(found))}")
json.dump({repr(V): P for V, P in sorted(found.items())}, open(outf, 'w'))
