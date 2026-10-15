"""Loopless quivers WITH cycles, no parallel arrows: linear uniqueness certificate (lengths 1..L) for every arrow-rule
vertex fold; positive nullities are then re-examined with the exact truncated fibre (which includes Δ at length L+1)."""
import sys, itertools
import lin
from loops_build import tbuild
from extend import vertex_folds
L = int(sys.argv[2]); nv = int(sys.argv[1])
lin.build = lambda nv_, arrows_: tbuild(nv_, arrows_, L)
def has_cycle(n, E):
    adj = {i: [b for a, b in E if a == i] for i in range(n)}
    color = [0] * n
    def dfs(u):
        color[u] = 1
        for w in adj[u]:
            if color[w] == 1 or (color[w] == 0 and dfs(w)): return True
        color[u] = 2; return False
    return any(color[i] == 0 and dfs(i) for i in range(n))
cand = [(i, j) for i in range(nv) for j in range(nv) if i != j]
seen = set(); qs = []
for mask in range(1, 1 << len(cand)):
    E = [cand[k] for k in range(len(cand)) if mask >> k & 1]
    if not has_cycle(nv, E): continue
    canon = min(tuple(sorted((p[a], p[b]) for a, b in E)) for p in itertools.permutations(range(nv)))
    if canon not in seen: seen.add(canon); qs.append(list(canon))
allok = True; suspicious = []
for arrows in qs:
    A = set(arrows); worst = 0; nf = 0
    for P0 in vertex_folds(nv):
        if not all(not [a for a in range(nv) if P0[a, s]] or not [b for b in range(nv) if P0[b, t]] or
                   all((a, b) in A for a in range(nv) if P0[a, s] for b in range(nv) if P0[b, t]) for s, t in arrows): continue
        nf += 1
        nul = [lin.system(nv, arrows, P0, m) for m in range(1, L + 1)]
        if any(nul): suspicious.append((arrows, [list(P0.column(v)) for v in range(nv)], nul))
        worst = max([worst] + [x for x in nul if x])
    print(f"{arrows}: {nf} vertex folds, max linear nullity {worst}", flush=True)
print(f"{len(qs)} loopless cyclic quivers on {nv} vertices; vertex folds with positive linear nullity: {len(suspicious)}")
for s in suspicious: print("  SUSPICIOUS", s)
