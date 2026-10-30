"""v2.77: the path family of v2_76_pathfamily with a set of vertices K killed: every basis path touching a vertex of K
(as a vertex, start or end of an arrow) is sent to 0, all else as in the family. Checks the fold equations identically
and the invariance of φ (generic symbolic S).
Usage: python3 analysis/v2_77_killsink.py QUIVER A0,A1,A2 PI1;PI2 K1,K2,... [OUT.json]"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import numpy as np, sympy as sp
import fold_ideal as fi
from v2_76_pathfamily import build
quiver = sys.argv[1]; cls = tuple(map(int, sys.argv[2].split(',')))
paths = [tuple(map(int, p.split('.'))) for p in sys.argv[3].split(';')]; K = set(map(int, sys.argv[4].split(',')))
out = sys.argv[5] if len(sys.argv) > 5 else None
ctx, R, P = build(quiver, cls, paths)
d = ctx['d']; basis = ctx['basis']; arrows = ctx['arrows']
def touches(b):
    if b[0] == 'v': return b[1] in K
    return any(arrows[a][0] in K or arrows[a][1] in K for a in b[1])
dead = [k for k, b in enumerate(basis) if touches(b)]
for c in dead:
    for r in range(d): P[r][c] = R(0)
for r in dead:
    for c in range(d): assert P[r][c] == 0 or c in dead, (basis[r], basis[c])
M = np.asarray(ctx['M'], dtype=np.int64); DL = np.asarray(ctx['DL'], dtype=np.int64)
nzM = [(a, b, k, int(M[a, b, k])) for a in range(d) for b in range(d) for k in range(d) if M[a, b, k]]
res = dict(quiver=quiver, paths=paths, killed=sorted(K), dead=[str(basis[k]) for k in dead])
for part in ('idem', 'delta', 'closure'): res[part] = len(fi._equations(P, R, d, nzM, DL, (part,)))
print(json.dumps(res))
if out: json.dump(res, open(out, 'w'), indent=1)
