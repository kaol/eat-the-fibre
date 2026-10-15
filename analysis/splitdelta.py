"""Is every Δ-compatible fold also compatible with each partial deletion Δ_c (delete only arrows of class c) separately?
If so, conjugation by any class-wise weighting φ_w (Δφ_w = (φ_w⊗φ_w)Σ_c w_c Δ_c) maps folds to folds: a torus of rank
#classes acts, larger than the single grading. Classes: parallel classes (and single arrows). Also tested per arrow."""
import os as _os, sys as _sys
if 'ROOT' not in globals():
    AN = _os.path.dirname(_os.path.abspath(__file__)); ROOT = _os.path.dirname(AN); _os.chdir(ROOT)
    __file__ = _os.path.join(ROOT, 'pcensus.py')   # pcensus.py locates gcensus.py relative to its own __file__
import sys, pickle, collections, numpy as np
q, p = sys.argv[1], int(sys.argv[2]); sys.argv = ['x']
exec(open(_os.path.join(ROOT, 'pcensus.py')).read().split("if __name__ == '__main__':")[0])
ctx = setup(q, p); basis, arrows, d = ctx['basis'], ctx['arrows'], ctx['d']
idx = {b: i for i, b in enumerate(basis)}
folds = [f.astype(np.int64) for f in pickle.load(open(f'logs/folds_{q}_{p}.pkl', 'rb'))['folds']]
def partial_delta(sel):
    DLs = np.zeros((d, d, d), dtype=np.int64)
    for l, b in enumerate(basis):
        if b[0] == 'v': continue
        s = b[1]
        for k, a in enumerate(s):
            if a not in sel: continue
            pre = idx[('p', s[:k])] if k > 0 else idx[('v', arrows[a][0])]
            suf = idx[('p', s[k + 1:])] if k + 1 < len(s) else idx[('v', arrows[a][1])]
            DLs[l, pre, suf] += 1
    return DLs
def ok(P, DLs): return all(((np.einsum('m,mab->ab', P[:, l], DLs) - P @ DLs[l] @ P.T) % p == 0).all() for l in range(d))
cls = collections.defaultdict(list)
for a, st in enumerate(arrows): cls[st].append(a)
for name_, sel in [(f"class {st}", set(c)) for st, c in cls.items()] + [(f"arrow {a}", {a}) for a in range(len(arrows))]:
    DLs = partial_delta(sel); n = sum(ok(P, DLs) for P in folds)
    print(f"{q}: {name_}: {n}/{len(folds)} folds compatible", flush=True)
