"""Conjecture 4 on random acyclic quivers without parallel arrows: linear certificate (nullity 0 at every length)
and arrow-rule fold counts.  Arrows i→j only for i<j (acyclic), each pair at most once (no parallel arrows)."""
import sys, random
from lin import system
from extend import vertex_folds
random.seed(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
NQ, NV = int(sys.argv[1]), int(sys.argv[2])
bad = 0; done = 0
for trial in range(NQ):
    nv = NV; pairs = [(i, j) for i in range(nv) for j in range(i + 1, nv)]
    arrows = sorted(random.sample(pairs, random.randint(1, len(pairs))))
    perm = list(range(nv)); random.shuffle(perm); arrows = [(perm[i], perm[j]) for i, j in arrows]
    A = set(arrows); cnt = 0; worst = 0
    for P0 in vertex_folds(nv):
        if not all(not [a for a in range(nv) if P0[a, s]] or not [b for b in range(nv) if P0[b, t]] or
                   all((a, b) in A for a in range(nv) if P0[a, s] for b in range(nv) if P0[b, t]) for s, t in arrows): continue
        cnt += 1
        for m in range(1, nv):
            k = system(nv, arrows, P0, m)
            if k is not None: worst = max(worst, k)
    bad += worst > 0; done += 1
    print(f"{arrows}: {cnt} folds, max nullity {worst}", flush=True)
print(f"{done} random quivers on {NV} vertices: Conjecture 4 holds on all: {bad == 0}")
