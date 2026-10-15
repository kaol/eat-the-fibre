"""Arrow rule: a vertex fold P0 (Cayley forest) passes iff for every arrow s→t, supp P0(e_s) × supp P0(e_t) ⊂ arrows."""
import sys
from extend import vertex_folds
QS = {'A1': (1, []), 'A2': (2, [(0, 1)]), 'A3': (3, [(0, 1), (1, 2)]), '1→2←3': (3, [(0, 1), (2, 1)]), '1←2→3': (3, [(1, 0), (1, 2)]),
      'A4': (4, [(0, 1), (1, 2), (2, 3)]), 'A4 alt': (4, [(0, 1), (2, 1), (2, 3)]), 'D4': (4, [(0, 3), (1, 3), (2, 3)]),
      '1→2→3, 1→3': (3, [(0, 1), (1, 2), (0, 2)]), 'diamond': (4, [(0, 1), (0, 2), (1, 3), (2, 3)]),
      'A5': (5, [(0, 1), (1, 2), (2, 3), (3, 4)]), 'A4 →→←': (4, [(0, 1), (1, 2), (3, 2)]), 'A4 ←←→': (4, [(1, 0), (2, 1), (2, 3)]), 'A4 →←←': (4, [(0, 1), (2, 1), (3, 2)]), 'swapfail': (5, [(0, 2), (1, 3), (2, 4), (3, 4)]), 'A6': (6, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)])}
def count(nv, arrows):
    A = set(arrows); c = 0
    for P0 in vertex_folds(nv):
        ok = True
        for (s, t) in arrows:
            S = [a for a in range(nv) if P0[a, s] != 0]; T = [b for b in range(nv) if P0[b, t] != 0]
            if S and T and not all((a, b) in A for a in S for b in T): ok = False; break
        c += ok
    return c
if __name__ == '__main__':
    for name in sys.argv[1:] or list(QS): print(name, count(*QS[name]), flush=True)
