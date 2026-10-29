"""v2.53: the census of Q₀ = 1⇉2⇇3 + 1→3 split by vertex block.
For each prime: from the saved G-orbit representatives (standard basis), the vertex block of each representative, the
number of folds with each block (the sum of G-orbit sizes), and the number of local-torus orbits with each block.
Usage: python3 -I analysis/v2_53_block_census.py REPS_q1.pkl REPS_q2.pkl ..."""
import sys, pickle, math, collections
import numpy as np
class SafeU(pickle.Unpickler):                        # v2.67: numpy arrays only (some pickles were uploaded)
    OK = {('numpy._core.numeric', '_frombuffer'), ('numpy', 'dtype'), ('numpy._core.multiarray', '_reconstruct'), ('numpy', 'ndarray'),
          ('numpy.core.numeric', '_frombuffer'), ('numpy.core.multiarray', '_reconstruct')}
    def find_class(self, mod, name):
        if (mod, name) in self.OK: return super().find_class(mod, name)
        raise pickle.UnpicklingError(f"blocked {mod}.{name}")
rows = {}
for f in sys.argv[1:]:
    R = SafeU(open(f, 'rb')).load(); p = R['p']; basis = R['basis']
    vid = [i for i, b in enumerate(basis) if b[0] == 'v']
    N = collections.Counter(); loc = collections.defaultdict(set); gorb = collections.Counter()
    for P, sz, l in zip(R['greps'], R['gsizes'], R['loc']):
        P = (np.asarray(P, dtype=np.int64) % p + p // 2) % p - p // 2
        V = tuple(map(tuple, P[np.ix_(vid, vid)].tolist()))
        N[V] += sz; loc[V].add(l); gorb[V] += 1
    rows[p] = (N, {V: len(s) for V, s in loc.items()}, gorb)
ps = sorted(rows)
blocks = sorted(set().union(*[set(rows[p][0]) for p in ps]))
fmt = lambda V: " ".join("".join(str(x) for x in r) for r in V)
print(f"{'vertex block':13} | " + " | ".join(f"F_{p}: folds, G-orb, loc" for p in ps) + " | deg est | loc constant")
tot = collections.Counter()
for V in blocks:
    cells = []; locs = []
    for p in ps:
        N, L, Go = rows[p]; cells.append(f"{N.get(V, 0):>8} {Go.get(V, 0):>4} {L.get(V, 0):>4}"); locs.append(L.get(V, 0)); tot[p] += L.get(V, 0)
    a, b = (ps[-2], ps[-1]) if len(ps) > 1 else (ps[0], ps[0])
    deg = None if a == b else round(math.log(rows[b][0][V] / rows[a][0][V]) / math.log(b / a)) if rows[a][0].get(V) and rows[b][0].get(V) else None
    print(f"{fmt(V):13} | " + " | ".join(cells) + f" | {deg} | {len(set(locs)) == 1}")
print("local-torus orbits per prime:", dict(tot))
