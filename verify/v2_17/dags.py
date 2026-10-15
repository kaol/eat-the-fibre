"""All DAGs on n vertices up to isomorphism (no parallel arrows), canonical form by brute-force relabelling."""
import itertools
def dags(n):
    seen = set(); out = []
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for mask in range(1 << len(pairs)):
        E = [pairs[k] for k in range(len(pairs)) if mask >> k & 1]
        canon = min(tuple(sorted((p[a], p[b]) for a, b in E)) for p in itertools.permutations(range(n)))
        if canon not in seen: seen.add(canon); out.append(list(canon))
    return out
if __name__ == '__main__':
    import sys; n = int(sys.argv[1]); print(n, len(dags(n)))
