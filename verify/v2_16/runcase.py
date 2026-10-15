"""Hand-built case with a run of two consecutive Podge vertices: x→c→c′→y, x→d, d→d′, d′→y (d child of c, d′ child of c′).
Build P0 from the intended Cayley forest, check the arrow rule, run the proof check and the linear certificate."""
import itertools
from sage.all__sagemath_singular import QQ, matrix
from proofcheck import check
from lin import system
x, c, d, cp, dp, y = range(6); nv = 6
arrows = [(x, c), (x, d), (c, cp), (d, dp), (cp, y), (dp, y)]
par = {x: 'R1', y: 'R1', c: 'R2', d: c, cp: 'R2', dp: cp}
def side(v):
    k = 0; u = v
    while isinstance(par[u], int): u = par[u]; k += 1
    root = par[u]; base = 'H' if root == 'R1' else 'P'
    return base if k % 2 == 0 else ('P' if base == 'H' else 'H')
blocksH = [[v] + [w for w in range(nv) if par[w] == v] for v in range(nv) if side(v) == 'H']
blocksQ = [[v] + [w for w in range(nv) if par[w] == v] for v in range(nv) if side(v) == 'P']
B = matrix(QQ, [[1 if v in blk else 0 for blk in blocksH + blocksQ] for v in range(nv)])
P0 = B * matrix.diagonal(QQ, [1] * len(blocksH) + [0] * len(blocksQ)) * B.inverse()
A = set(arrows)
ar = all(not [a for a in range(nv) if P0[a, s]] or not [b for b in range(nv) if P0[b, t]] or
         all((a, b) in A for a in range(nv) if P0[a, s] for b in range(nv) if P0[b, t]) for s, t in arrows)
print("sides:", {v: side(v) for v in range(nv)}, " arrow rule:", ar)
r = check(nv, arrows, P0)
print("proof check P1..P4:", r[:4], " run terms:", r[4])
print("linear certificate nullities by length:", [system(nv, arrows, P0, m) for m in (1, 2, 3)])
