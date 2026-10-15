"""Representatives of the free (G × grading)-orbits of pair_srcout (1⇉2, 1→3) over F_p, printed in path names."""
import sys, collections, numpy as np
sys.argv = ['x', 'pair_srcout', sys.argv[1] if len(sys.argv) > 1 else '5']
exec(open("graded.py").read().split("print(f\"{name} F_{p}")[0])
names = ['e1', 'e2', 'e3', 'a1', 'a2', 'g']
def show(P):
    out = []
    for j in range(d):
        terms = [f"{'' if P[i, j] == 1 else ('-' if P[i, j] == p - 1 else int(P[i, j]))}{names[i]}" for i in range(d) if P[i, j]]
        out.append(f"P({names[j]}) = {' + '.join(terms) if terms else '0'}")
    return "; ".join(out)
seen = set(); big = []
for P in folds:
    if key(P) in seen: continue
    o = set()
    for g, gi in G:
        Q = g @ P @ gi % p
        for e, ei in Gr: o.add(key(e @ Q @ ei))
    seen |= o
    if len(o) == p * (p - 1) ** 2: big.append(P)
for P in big: print(show(P))
