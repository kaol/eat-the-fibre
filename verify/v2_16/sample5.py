import random, sys
from dags import dags
import proofcheck
from extend import vertex_folds
random.seed(int(sys.argv[2])); qs = random.sample([a for a in dags(5) if a], int(sys.argv[1]))
allok = [True] * 4; nf = 0; nt = 0
for arrows in qs:
    A = set(arrows)
    for P0 in vertex_folds(5):
        if not all(not [a for a in range(5) if P0[a, s]] or not [b for b in range(5) if P0[b, t]] or
                   all((a, b) in A for a in range(5) if P0[a, s] for b in range(5) if P0[b, t]) for s, t in arrows): continue
        r = proofcheck.check(5, arrows, P0); nf += 1; nt += r[4]
        for i in range(4): allok[i] &= r[i]
    print(arrows, flush=True)
print(f"{len(qs)} random 5-vertex DAGs, {nf} vertex folds, {nt} run terms: P1 {allok[0]}, P2 {allok[1]}, P3 {allok[2]}, P4 {allok[3]}")
