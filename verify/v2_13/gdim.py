"""dim of the Lie algebra of Delta-compatible derivations (exact) vs the parallel-arrow count  Σ_{i,j} m_ij (m_ij − 1)."""
import sys
from collections import Counter
from exactfold import tensors
from chartpair import derivations
Qs = {'A2': (2, [(0, 1)]), 'A3': (3, [(0, 1), (1, 2)]), '1→2←3': (3, [(0, 1), (2, 1)]), 'K2': (2, [(0, 1)] * 2),
      'K3': (2, [(0, 1)] * 3), 'K4': (2, [(0, 1)] * 4), '1⇉2→3': (3, [(0, 1), (0, 1), (1, 2)]), '1⇉2⇉3': (3, [(0, 1), (0, 1), (1, 2), (1, 2)]),
      'diamond': (4, [(0, 1), (0, 2), (1, 3), (2, 3)]), '1→2, 1→3, 2→3': (3, [(0, 1), (0, 2), (1, 2)]), 'K2+K3 (1⇉2, 2⇛3)': (3, [(0, 1)] * 2 + [(1, 2)] * 3)}
if __name__ == '__main__':
    ok = True
    for name, (nv, arrows) in Qs.items():
        d, M, D = tensors(nv, arrows); g = len(derivations(d, M, D))
        pred = sum(m * (m - 1) for m in Counter(arrows).values()); ok &= (g == pred)
        print(f"{name}: dim A = {d}, dim g = {g}, Σ m(m−1) = {pred}", flush=True)
    print("dim g = Σ m_ij(m_ij − 1) on all examples:", ok)
