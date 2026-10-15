"""Orbit dimension versus orbit size for every (G × T)-orbit of folds over F_p (T = class-weight torus, all weightings;
only for quivers where every fold admits all of them, i.e. local_torus_sizes = {(p−1)^classes: all}).

For each orbit: its size |H(F_p)|/|Stab(F_p)| (H = G × T), and the rank r of the Lie-algebra action at the representative
(mod p: rank{[D, P]} over D ∈ Lie G plus the class generators E_c), which is the orbit dimension when the rank is not
accidentally lowered mod p. dim Stab = dim H − r. A stabilizer of dimension 0 with |Stab(F_p)| > 1 is a finite group:
its F_q-point count can depend on q mod n (e.g. μ₄ ⊂ F_q^× iff 4 | q − 1), so the number of F_q-orbits it produces can
change with q without any continuous family.
Usage: python3 analysis/stabtype.py QUIVER p   (needs logs/folds_QUIVER_p.pkl)"""
import os as _os, sys as _sys
if 'ROOT' not in globals():
    AN = _os.path.dirname(_os.path.abspath(__file__)); ROOT = _os.path.dirname(AN); _os.chdir(ROOT)
    __file__ = _os.path.join(ROOT, 'pcensus.py')
import sys, pickle, itertools, collections, numpy as np
name, p = sys.argv[1], int(sys.argv[2]); sys.argv = ['x']
exec(open(_os.path.join(ROOT, 'pcensus.py')).read().split("if __name__ == '__main__':")[0])
ctx = setup(name, p); C.update(ctx); F = ctx['F']; basis, arrows, d = ctx['basis'], ctx['arrows'], ctx['d']
folds = [f.astype(np.int64) for f in pickle.load(open(f'logs/folds_{name}_{p}.pkl', 'rb'))['folds']]
G = build_G(ctx)
cls = collections.defaultdict(list)
for a, st in enumerate(arrows): cls[st].append(a)
cl = list(cls.values())
T = []
for w in itertools.product(range(1, p), repeat=len(cl)):
    wa = {a: w[k] for k, c in enumerate(cl) for a in c}
    dg = [1 if b[0] == 'v' else int(np.prod([wa[x] for x in b[1]])) % p for b in basis]
    T.append((np.diag(dg), np.diag([pow(x, p - 2, p) for x in dg])))
H = [(t @ g % p, gi @ ti % p) for g, gi in G for t, ti in T]
Hg = (np.array([h for h, _ in H]), np.array([hi for _, hi in H]))
lab, reps, sizes = orbit_label(folds, Hg, p)
gens, E = lie_G(ctx)
Ec = [np.diag([0 if b[0] == 'v' else sum(1 for x in b[1] if x in c) for b in basis]).astype(np.int64) for c in cl]
dimH = len(gens) + len(cl)
tab = collections.Counter()
for P, sz in zip(reps, sizes):
    r = F.rank(np.array([((D @ P - P @ D) % p).ravel() for D in gens + Ec]))
    tab[(sz, r)] += 1
print(f"{name} F_{p}: |H| = {len(H)}, dim H = {dimH}, {len(reps)} H-orbits")
print("  (orbit size, orbit dimension r) : orbits    [dim Stab = dim H − r; |Stab(F_p)| = |H|/size]")
for (sz, r), n in sorted(tab.items()):
    print(f"  size {sz:6d}  r = {r}  dim Stab = {dimH - r}  |Stab(F_p)| = {len(H) // sz:5d} : {n}")
