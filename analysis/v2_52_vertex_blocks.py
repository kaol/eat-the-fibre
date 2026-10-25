"""v2.52: the vertex blocks of all folds in a census checkpoint.

The vertex block of a fold P is the matrix V[i, j] = coefficient of e_i in P(e_j), over the vertices. The vertices are
the first basis vectors of K (pcensus.setup), and the fibres store X = P·B, so V = X[vertex rows, vertex columns] for
every fold of every fibre. Entries are lifted symmetrically to integers in (−p/2, p/2]. The vertex block is constant
along grading-torus orbits, since φ_l is the identity on vertices.
Usage: python3 analysis/v2_52_vertex_blocks.py QUIVER P [P ...]"""
import os, sys, pickle, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
import numpy as np
import pcensus as pc
name, ps = sys.argv[1], [int(x) for x in sys.argv[2:]]
sets = {}
for p in ps:
    ctx = pc.setup(name, p)
    vid = [i for i, b in enumerate(ctx['basis']) if b[0] == 'v']
    kv = [j for j in range(ctx['kdim']) if np.count_nonzero(np.asarray(ctx['Kb'][j]) % p) == 1 and
          int(np.nonzero(np.asarray(ctx['Kb'][j]) % p)[0][0]) in vid]
    assert len(kv) == len(vid), "vertices are not K basis vectors"
    order = [int(np.nonzero(np.asarray(ctx['Kb'][j]) % p)[0][0]) for j in kv]
    ck = pickle.load(open(os.path.join(os.path.dirname(HERE), 'checkpoints', f'{name}_{p}.pkl'), 'rb'))
    assert ck['nreps'] == len(ck['fibres'])
    cnt = collections.Counter()
    for r, fib in ck['fibres'].items():
        for X in fib:
            X = np.asarray(X, dtype=np.int64) % p
            V = X[np.ix_(order, kv)]                      # rows: vertex coefficients; columns: P(e_v)
            V = (V + p // 2) % p - p // 2
            cnt[tuple(map(tuple, V.tolist()))] += 1
    sets[p] = cnt
    ents = sorted({x for V in cnt for row in V for x in row})
    print(f"{name} F_{p}: {sum(cnt.values())} fibre folds (torus reps), {len(cnt)} distinct vertex blocks; entries {ents}",
          flush=True)
common = set.intersection(*[set(c) for c in sets.values()])
print(f"vertex blocks common to all primes: {len(common)}; per prime: {[len(c) for c in sets.values()]}")
if len(ps) > 1:
    for p in ps:
        extra = set(sets[p]) - common
        print(f"  F_{p}: {len(extra)} blocks not in every prime" + (f", e.g. {sorted(extra)[:3]}" if extra else ""))
if len(common) <= 60:
    for V in sorted(common): print("   ", V)
