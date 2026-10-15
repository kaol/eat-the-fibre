"""Random loopless cyclic quivers on 4 vertices: for every arrow-rule vertex fold compare the predictor
(∃ N-path x→…→y, x ≠ y, S_x ∩ S_y ≠ ∅) with the linear certificate (positive nullity at some length ≤ L)."""
import sys, random, itertools
import lin
from loops_build import tbuild
from extend import vertex_folds
nv, NQ, L, seed = 4, int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); random.seed(seed)
lin.build = lambda a, b: tbuild(a, b, L)
def has_cycle(n, E):
    adj = {i: [b for a, b in E if a == i] for i in range(n)}; color = [0] * n
    def dfs(u):
        color[u] = 1
        for w in adj[u]:
            if color[w] == 1 or (color[w] == 0 and dfs(w)): return True
        color[u] = 2; return False
    return any(color[i] == 0 and dfs(i) for i in range(n))
cand = [(i, j) for i in range(nv) for j in range(nv) if i != j]
qs = []
while len(qs) < NQ:
    E = sorted(random.sample(cand, random.randint(3, 6)))
    if has_cycle(nv, E) and E not in qs: qs.append(E)
agree = True; pos = 0; flagged = 0
for arrows in qs:
    A = set(arrows)
    for P0 in vertex_folds(nv):
        if not all(not [a for a in range(nv) if P0[a, s]] or not [b for b in range(nv) if P0[b, t]] or
                   all((a, b) in A for a in range(nv) if P0[a, s] for b in range(nv) if P0[b, t]) for s, t in arrows): continue
        S = [{a for a in range(nv) if P0[a, v]} for v in range(nv)]
        paths = [[v] for v in range(nv) if S[v]]; hit = False
        for _ in range(L):
            paths = [p + [b] for p in paths for (a, b) in arrows if a == p[-1] and S[b]]
            hit |= any(p[0] != p[-1] and S[p[0]] & S[p[-1]] for p in paths)
        nul = any(lin.system(nv, arrows, P0, m) for m in range(1, L + 1))
        pos += nul; flagged += hit
        if nul and not hit: agree = False; print("  linear nullity without predictor:", arrows, [list(P0.column(v)) for v in range(nv)], flush=True)
        if hit and not nul: print("  predictor without linear nullity:", arrows, [list(P0.column(v)) for v in range(nv)], flush=True)
    print(arrows, flush=True)
print(f"{NQ} random cyclic 4-vertex quivers: positive linear nullity {pos}, predictor flags {flagged}; every positive nullity is predicted: {agree}")
