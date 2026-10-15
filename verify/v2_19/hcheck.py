"""Hypothesis (H): no path of length ≥ 1 through vertices with S_v ≠ ∅ joins a Podge head c to a vertex of its own block C
(pairs (c,d), (d,c), (c,c)).  Under (H) the uniqueness proof of v2.16 applies verbatim (loopless, no parallel arrows).
Compare (H)-violations with the predictor and with the truth on all cyclic 3-vertex quivers."""
import sys, itertools
from sage.all__sagemath_singular import QQ, matrix, vector, identity_matrix
from extend import vertex_folds
def podge_blocks(nv, P0):
    Q0 = (identity_matrix(QQ, nv) - matrix(QQ, P0)).column_space(); blocks = []
    for S in itertools.chain.from_iterable(itertools.combinations(range(nv), r) for r in range(1, nv + 1)):
        v = vector(QQ, [1 if i in S else 0 for i in range(nv)])
        if v in Q0 and not any(set(S) > B for B in blocks): blocks.append(set(S))
    hodge = {v for v in range(nv) if P0[v, v] != 0 and all(P0[a, v] >= 0 for a in range(nv))}
    return [(B, [c for c in B if c not in hodge]) for B in blocks]
def violates_H(nv, arrows, P0, L):
    S = [{a for a in range(nv) if P0[a, v]} for v in range(nv)]
    paths = [[v] for v in range(nv) if S[v]]; viol = False
    blocks = podge_blocks(nv, P0)
    for _ in range(L):
        paths = [p + [b] for p in paths for (a, b) in arrows if a == p[-1] and S[b]]
        for p in paths:
            for B, heads in blocks:
                for c in heads:
                    if (p[0] == c and p[-1] in B) or (p[-1] == c and p[0] in B): viol = True
    return viol
def predictor(nv, arrows, P0, L):
    S = [{a for a in range(nv) if P0[a, v]} for v in range(nv)]
    paths = [[v] for v in range(nv) if S[v]]; hit = False
    for _ in range(L):
        paths = [p + [b] for p in paths for (a, b) in arrows if a == p[-1] and S[b]]
        hit |= any(p[0] != p[-1] and S[p[0]] & S[p[-1]] for p in paths)
    return hit
if __name__ == '__main__':
    nv, L = 3, 4
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
    hv = []; pr = []
    for arrows in qs:
        A = set(arrows)
        for P0 in vertex_folds(nv):
            if not all(not [a for a in range(nv) if P0[a, s]] or not [b for b in range(nv) if P0[b, t]] or
                       all((a, b) in A for a in range(nv) if P0[a, s] for b in range(nv) if P0[b, t]) for s, t in arrows): continue
            key = (tuple(arrows), tuple(tuple(P0.column(v)) for v in range(nv)))
            if violates_H(nv, arrows, P0, L): hv.append(key)
            if predictor(nv, arrows, P0, L): pr.append(key)
    print(f"cyclic 3-vertex quivers: (H) violated by {len(hv)} vertex folds; predictor flags {len(pr)}; same set: {set(hv) == set(pr)}")
