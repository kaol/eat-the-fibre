"""Lightweight arrow-rule count directly on Cayley forests (parent arrays; roots 'H' and 'P'), using the explicit
description S_v = subtree(v) (Hodge v) or subtree(v)∖{v} (Podge v).  No matrices: independent of extend.vertex_folds."""
import sys, itertools
def count(n, arrows):
    E = set(arrows); total = 0
    choices = list(range(n)) + ['H', 'P']
    for par in itertools.product(choices, repeat=n):
        side = [None] * n; ok = True
        for v in range(n):                         # walk to root, detect cycles
            seen = set(); u = v; k = 0
            while isinstance(u, int):
                if u in seen: ok = False; break
                seen.add(u); u = par[u]; k += 1
            if not ok: break
            # depth k ≥ 1 from root u: depth-1 nodes take the root's side, sides alternate
            side[v] = u if k % 2 == 1 else ('P' if u == 'H' else 'H')
        if not ok: continue
        ch = [[w for w in range(n) if par[w] == v] for v in range(n)]
        sub = [None] * n
        def subtree(v):
            if sub[v] is None:
                s = {v}
                for w in ch[v]: s |= subtree(w)
                sub[v] = s
            return sub[v]
        S = [subtree(v) if side[v] == 'H' else subtree(v) - {v} for v in range(n)]
        if all(not S[s] or not S[t] or all((a, b) in E for a in S[s] for b in S[t]) for s, t in arrows): total += 1
    return total
if __name__ == '__main__':
    fam = sys.argv[1]; N = int(sys.argv[2])
    f = {'star': lambda n: [(0, i) for i in range(1, n)], 'linear': lambda n: [(i, i + 1) for i in range(n - 1)],
         'zigzag': lambda n: [(i, i + 1) if i % 2 == 0 else (i + 1, i) for i in range(n - 1)]}[fam]
    print(fam, [count(n, f(n)) for n in range(1, N + 1)], flush=True)
