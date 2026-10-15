"""Confirm the 4-vertex cyclic sample: for every (H)-violating vertex fold, compute the exact truncated fibre dimension
(L given) and the candidate count Σ over violating head–block pairs (x, y) of |S_c|; for every non-violating fold
the linear certificate already proves uniqueness."""
import sys, random, itertools
import cfibre, loops_build
from hcheck import violates_H, podge_blocks
from extend import vertex_folds
nv, NQ, L, seed = 4, int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); random.seed(seed)
def has_cycle(n, E):
    adj = {i: [b for a, b in E if a == i] for i in range(n)}; color = [0] * n
    def dfs(u):
        color[u] = 1
        for w in adj[u]:
            if color[w] == 1 or (color[w] == 0 and dfs(w)): return True
        color[u] = 2; return False
    return any(color[i] == 0 and dfs(i) for i in range(n))
cand = [(i, j) for i in range(nv) for j in range(nv) if i != j]; qs = []
while len(qs) < NQ:
    E = sorted(random.sample(cand, random.randint(3, 6)))
    if has_cycle(nv, E) and E not in qs: qs.append(E)
def pairs(arrows, P0, Lp):
    S = [{a for a in range(nv) if P0[a, v]} for v in range(nv)]
    blocks = podge_blocks(nv, P0); found = set()
    paths = [[v] for v in range(nv) if S[v]]
    for _ in range(Lp):
        paths = [p + [b] for p in paths for (a, b) in arrows if a == p[-1] and S[b]]
        for p in paths:
            for B, heads in blocks:
                for c in heads:
                    if (p[0] == c and p[-1] in B) or (p[-1] == c and p[0] in B): found.add((p[0], p[-1], c))
    return sum(len(S[c]) for (_, _, c) in found), found
rows = []
for arrows in qs:
    A = set(arrows)
    for P0 in vertex_folds(nv):
        if not all(not [a for a in range(nv) if P0[a, s]] or not [b for b in range(nv) if P0[b, t]] or
                   all((a, b) in A for a in range(nv) if P0[a, s] for b in range(nv) if P0[b, t]) for s, t in arrows): continue
        if not violates_H(nv, arrows, P0, L + 1): continue
        R, J, n = cfibre.fibre(nv, arrows, P0, L)
        dim = J.dimension(); cnt, pr = pairs(arrows, P0, L + 1)
        rows.append((arrows, [list(P0.column(v)) for v in range(nv)], dim, cnt, sorted(pr)))
        print(f"{arrows} P0={rows[-1][1]}: fibre dim {dim}, pair count {cnt}, pairs {sorted(pr)}", flush=True)
pos = sum(1 for r in rows if r[2] > 0)
print(f"(H)-violating vertex folds: {len(rows)}; with positive fibre: {pos}; pair-count formula matches: {sum(1 for r in rows if r[2] == r[3])}/{len(rows)}")
