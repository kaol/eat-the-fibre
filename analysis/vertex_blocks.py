"""Vertex-block accounting of orbit growth (v2.27).

For a quiver with converging classes A: u → v, B: w → v (u, w, v = the first such pair found by theorem3.family), and
for each prime given: H-orbits of folds (H = G × the class weightings preserving the fold set) grouped by vertex block
B(P) = π₀P|_{kQ₀} (an H-invariant, by the vertex-block lemma). Prints
  - the histogram {orbits per block: number of blocks}; growth confined to single-orbit blocks shows up as identical
    histograms except at 1;
  - the single-orbit blocks, classified by their restriction to {u, v, w}: one of the families F1–F4 of Theorem 3 at a
    generic parameter, or 'other';
  - the orbit total on all other blocks.
Usage: python3 analysis/vertex_blocks.py QUIVER p [p ...]   (needs logs/folds_QUIVER_p.pkl)"""
import os as _os, sys as _sys, runpy, collections
import numpy as np
ARGS = list(_sys.argv)
AN = _os.path.dirname(_os.path.abspath(__file__))
T3 = runpy.run_path(_os.path.join(AN, 'theorem3.py'), run_name='t3')
name, primes = ARGS[1], [int(x) for x in ARGS[2:]]
basis = T3['setup'](name, 3)['basis']
_, _, _, _, fams = T3['family'](name)
_, u, v, w, _ = fams[0][0]
vi = {b[1]: i for i, b in enumerate(basis) if b[0] == 'v'}
iu, iv, iw = vi[u - 1], vi[v - 1], vi[w - 1]
nv = len(vi)

def classify(B, p):
    """B[i] = P(e_i) restricted to vertices (row i), entries mod p"""
    r1 = B[iu][[iu, iv, iw]] % p; r2 = B[iv] % p; r3 = B[iw][[iu, iv, iw]] % p
    if r2.any(): return 'other'
    if r1[0] == 1 and r1[2] == 0 and ((r3 + r1) % p == 0).all() and r1[1] not in (0, 1): return 'F1'
    if r3[2] == 1 and r3[0] == 0 and ((r1 + r3) % p == 0).all() and r3[1] not in (0, 1): return 'F3'
    if r1[1] == 0 and (r1[2] - (r1[0] - 1)) % p == 0 and ((r3 + r1) % p == 0).all() and r1[0] not in (0, 1): return 'F2'
    if (r1[0], r1[2], r3[0], r3[2]) == (1, 0, 0, 1) and (r1[1] + r3[1] - 1) % p == 0 and r1[1] not in (0, 1): return 'F4'
    return 'other'

for p in primes:
    r = T3['vertex_blocks'](name, p)
    pv = r['_per_vb']
    print(f"{name} F_{p}: N = {r['N']}, H-orbits = {r['orbits']}, vertex blocks = {r['vertex_blocks']}, F1–F4 members in census: {r['families_in_census']}")
    print(f"  orbits per block: {r['orbits_per_block_hist']}")
    cls = collections.Counter(); others = []
    for k, n in pv.items():
        if n != 1: continue
        B = np.frombuffer(k, dtype=np.int8).reshape(nv, nv).astype(int).T        # row i = P(e_i)
        c = classify(B, p); cls[c] += 1
        if c == 'other': others.append(np.where(B > p // 2, B - p, B).tolist())
    print(f"  single-orbit blocks: {sum(cls.values())}, by restriction to (u, v, w) = ({u}, {v}, {w}): {dict(sorted(cls.items()))}")
    for o in others: print(f"    other: {o}")
    print(f"  orbits on all other blocks: {r['orbits'] - sum(cls.values())}", flush=True)
