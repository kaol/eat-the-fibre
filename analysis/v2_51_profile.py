"""v2.51: per-column frontier profile of the torch fibre engine, one fibre at a time (diagnosis of job 14's OOM).

For each sampled source fibre r (a torus rep whose fibre the engine computes, after the G × T symmetry plan) this runs the
engine's own `level` and `quotient` breadth-first, column by column, and records the frontier size after each column
(after the N-quotient, as the engine stores it). A fibre is abandoned when its frontier exceeds CAP rows or it uses more
than TCAP seconds; the profile up to that point is kept. Breadth-first is not the engine's order (depth-first) but it
gives the per-column widths that the engine's stack must hold for one parent batch at each column.

Output: one JSON line per fibre: r, rank PK, nnz PK, profile [(column, frontier)], leaves (complete N-orbit reps), time,
status (done | cap | time).
Usage: python3 analysis/v2_51_profile.py QUIVER P OUT.jsonl [--sample N] [--seed S] [--only r1,r2,...] [--cap ROWS] [--tcap SEC]
       env KCACHE=path for the kernel-stage cache (default checkpoints/QUIVER_P.kernel.npz)"""
import os, sys, time, json, argparse, random, io, contextlib
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
import numpy as np, torch
import pcensus as pc, fibre_torch as ft
ap = argparse.ArgumentParser(); ap.add_argument('quiver'); ap.add_argument('p', type=int); ap.add_argument('out')
ap.add_argument('--sample', type=int, default=50); ap.add_argument('--seed', type=int, default=51)
ap.add_argument('--only'); ap.add_argument('--cap', type=float, default=2e6); ap.add_argument('--tcap', type=float, default=120)
a = ap.parse_args()
t0 = time.time()
ctx = pc.setup(a.quiver, a.p); pc.C.update(ctx)
phis = pc.grading(ctx)
cache = os.environ.get('KCACHE', os.path.join(ROOT, 'checkpoints', f'{a.quiver}_{a.p}.kernel.npz'))
Kf, reps = pc.kernel_stage(ctx, phis, 'auto', None, cache)
src, trans, ctx['N'] = pc.symmetry_plan(ctx, reps, phis)
d, kdim, p = ctx['d'], ctx['kdim'], ctx['p']
print(f"{a.quiver} F_{p}: d={d}, kdim={kdim}, {len(Kf)} kernel folds, {len(reps)} torus reps, {len(src)} sources, "
      f"|N| = {len(ctx['N'][0])}; setup {time.time() - t0:.0f}s", flush=True)
dev = torch.device('cpu'); P = ft.Plan(ctx, dev)
Kmat = np.asarray(ctx['Kmat'])
if a.only: todo = [int(x) for x in a.only.split(',')]
else:
    random.seed(a.seed); todo = sorted(random.sample(src, min(a.sample, len(src))))
out = open(a.out, 'a')
for r in todo:
    PK = np.asarray(reps[r][0]) % p
    X0 = np.zeros((1, d, d), dtype=np.int8); X0[0, :, :kdim] = Kmat @ PK % p
    rid, X = torch.zeros(1, dtype=torch.long), torch.as_tensor(X0)
    prof = []; status = 'done'; ts = time.time(); leaves = 0
    for c in range(kdim, d):
        nb = P.parent_batch(c); kr, kx = [], []; n = 0
        for s in range(0, len(rid), nb):
            for r2, x2 in ft.level(P, rid[s:s + nb], X[s:s + nb], c):
                if P.sym and c + 1 < d and len(r2): r2, x2 = ft.quotient(P, r2, x2, c + 1)
                kr.append(r2); kx.append(x2); n += len(r2)
                if n > a.cap: status = 'cap'; break
                if time.time() - ts > a.tcap: status = 'time'; break
            if status != 'done': break
        prof.append((c + 1, n))
        if status != 'done' or n == 0: break
        rid, X = torch.cat(kr), torch.cat(kx)
    if status == 'done' and prof and prof[-1][0] == d: leaves = prof[-1][1]
    rec = dict(r=r, rankPK=int(np.linalg.matrix_rank(PK.astype(float))), nnzPK=int((PK != 0).sum()), profile=prof,
               leaves=leaves, t=round(time.time() - ts, 2), status=status)
    out.write(json.dumps(rec) + "\n"); out.flush()
    print(json.dumps(rec), flush=True)
