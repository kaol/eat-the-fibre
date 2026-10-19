"""Vertex blocks of local-torus orbits, from a --save-reps dump (v2.28).

For each local-torus orbit (a component of G-orbits joined by local class-torus moves), take the vertex block
B(P) = π₀P|_{kQ₀} of its G-orbit representatives. Every move fixes the vertices, so B is constant on a component;
this is asserted. Prints the number of distinct blocks and the histogram {orbits per block: number of blocks}.
Usage: python3 analysis/loc_blocks.py logs/reps_QUIVER_p.pkl [more dumps ...]"""
import sys, pickle, collections
import numpy as np

def blocks(path):
    D = pickle.load(open(path, 'rb'))
    p, basis, greps, loc = D['p'], D['basis'], D['greps'], D['loc']
    vert = [i for i, b in enumerate(basis) if b[0] == 'v']
    comp_block = {}
    for i, P in enumerate(greps):
        k = (P.astype(np.int64)[np.ix_(vert, vert)] % p).astype(np.int8).tobytes()
        assert comp_block.setdefault(loc[i], k) == k, "vertex block not constant on a local-torus orbit"
    per_block = collections.Counter(comp_block.values())
    return D, comp_block, per_block

if __name__ == '__main__':
    for path in sys.argv[1:]:
        D, cb, pb = blocks(path)
        print(f"{D['name']} F_{D['p']}: {len(cb)} local-torus orbits, {len(pb)} vertex blocks; "
              f"orbits per block: {dict(sorted(collections.Counter(pb.values()).items()))}", flush=True)
