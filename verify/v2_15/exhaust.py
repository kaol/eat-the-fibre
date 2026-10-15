import sys
from dags import dags
from lin import system
from extend import vertex_folds
import podge
n = int(sys.argv[1]); allc = True; alld = True; cnt = 0
for arrows in dags(n):
    if not arrows: continue
    A = set(arrows); worst = 0; dok = True; nf = 0
    for P0 in vertex_folds(n):
        if not all(not [a for a in range(n) if P0[a, s]] or not [b for b in range(n) if P0[b, t]] or
                   all((a, b) in A for a in range(n) if P0[a, s] for b in range(n) if P0[b, t]) for s, t in arrows): continue
        nf += 1
        for m in range(1, n):
            k = system(n, arrows, P0, m)
            if k is not None: worst = max(worst, k)
        if 'D' in sys.argv[2:]:
            for m, (kind, ok, tot) in podge.analyse(n, arrows, P0).items(): dok &= ok
    allc &= worst == 0; alld &= dok; cnt += 1
    print(f"{arrows}: {nf} folds, nullity {worst}, spanning lemma {dok}", flush=True)
print(f"all {cnt} DAGs on {n} vertices: Conjecture 4 certified: {allc}; spanning lemma (D1, D2): {alld}")
