"""Fold counts (= arrow-rule Cayley forests, by the theorem) for further quiver families."""
import sys
from arrowrule import count
fams = {'linear A_n': lambda n: [(i, i + 1) for i in range(n - 1)],
        'zigzag A_n': lambda n: [(i, i + 1) if i % 2 == 0 else (i + 1, i) for i in range(n - 1)],
        'source star': lambda n: [(0, i) for i in range(1, n)],
        'sink star': lambda n: [(i, 0) for i in range(1, n)]}
N = int(sys.argv[1])
for name, f in fams.items():
    print(name, [count(n, f(n)) for n in range(1, N + 1)], flush=True)
