"""Predictor for non-uniqueness on loopless quivers: some path p (length ≥ 1, all vertices with S_v ≠ ∅) has
S_start ∩ S_end ≠ ∅ (in acyclic quivers this is impossible: it would close a walk).  Compare with the cyclic sweep:
the predictor should flag exactly the vertex folds with positive linear nullity."""
import sys, itertools
from extend import vertex_folds
nv, L = int(sys.argv[1]), int(sys.argv[2])
exec(open('cyclic.py').read().split("allok = True")[0].split("L = int")[0])   # reuse has_cycle (imports only)
def has_cycle(n, E):
    adj = {i: [b for a, b in E if a == i] for i in range(n)}; color = [0] * n
    def dfs(u):
        color[u] = 1
        for w in adj[u]:
            if color[w] == 1 or (color[w] == 0 and dfs(w)): return True
        color[u] = 2; return False
    return any(color[i] == 0 and dfs(i) for i in range(n))
cand = [(i, j) for i in range(nv) for j in range(nv) if i != j]; seen = set(); qs = []
for mask in range(1, 1 << len(cand)):
    E = [cand[k] for k in range(len(cand)) if mask >> k & 1]
    if not has_cycle(nv, E): continue
    canon = min(tuple(sorted((p[a], p[b]) for a, b in E)) for p in itertools.permutations(range(nv)))
    if canon not in seen: seen.add(canon); qs.append(list(canon))
flag = []
for arrows in qs:
    A = set(arrows)
    for P0 in vertex_folds(nv):
        if not all(not [a for a in range(nv) if P0[a, s]] or not [b for b in range(nv) if P0[b, t]] or
                   all((a, b) in A for a in range(nv) if P0[a, s] for b in range(nv) if P0[b, t]) for s, t in arrows): continue
        S = [{a for a in range(nv) if P0[a, v]} for v in range(nv)]
        paths = [[v] for v in range(nv) if S[v]]; hit = False
        for _ in range(L):
            new = []
            for p in paths:
                for (a, b) in arrows:
                    if a == p[-1] and S[b]: new.append(p + [b])
            for p in new:
                if p[0] != p[-1] and S[p[0]] & S[p[-1]]: hit = True
            paths = new
        if hit: flag.append((arrows, [list(P0.column(v)) for v in range(nv)]))
print(f"predictor flags {len(flag)} vertex folds:")
for f in flag: print("  ", f)
