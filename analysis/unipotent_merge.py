"""Orbits of Aut_Δ° where it is larger than G (v2.28).

On a quiver with consecutive classes A: u→v, B: v→w and an arrow c: u→w, Lie(Aut_Δ) has one more generator than Lie G:
the derivation c ↦ δ_aδ_b (theorem3.lie_aut_delta). Its group U = {φ_s : c ↦ c + s·δ_aδ_b} consists of algebra
automorphisms with Δφ_s = (φ_s⊗φ_s)Δ (δ_aδ_b ∈ ker Δ, and c composes with no arrow here), normalized by G and by the
class weightings. This script takes a --save-reps dump (G-orbit representatives and local-torus components), merges
the components along φ_s for s ∈ F_p, and reports the orbit count under G, U and the local torus together.
Usage: python3 analysis/unipotent_merge.py QUIVER p   (needs logs/reps_QUIVER_p.pkl)"""
import os as _os, sys as _sys, pickle, collections
import numpy as np
ARGS = list(_sys.argv)
AN = _os.path.dirname(_os.path.abspath(__file__)); ROOT = _os.path.dirname(AN); _os.chdir(ROOT)
_sys.argv = ['x']; __file__ = _os.path.join(ROOT, 'pcensus.py')
exec(open(__file__).read().split("if __name__ == '__main" + "__':")[0])
name, p = ARGS[1], int(ARGS[2])
ctx = setup(name, p); C.update(ctx); basis, arrows, d = ctx['basis'], ctx['arrows'], ctx['d']
idx = {b: i for i, b in enumerate(basis)}
cls = collections.defaultdict(list)
for a, st in enumerate(arrows): cls[st].append(a)
pairs = [(A, B, cls[(u, w)]) for (u, v), A in cls.items() for (v2, w), B in cls.items()
         if v2 == v and len(A) >= 2 and len(B) >= 2 and (u, w) in cls]
assert pairs, "no consecutive classes with a shortcut arrow"
D = pickle.load(open(f'logs/reps_{name}_{p}.pkl', 'rb')); greps, loc = D['greps'].astype(np.int64), D['loc']
key = lambda A_: (A_ % p).astype(np.int8).tobytes()
gk = {key(R_): i for i, R_ in enumerate(greps)}
G = build_G(ctx); Gm = np.array([g for g, _ in G]); Gi = np.array([gi for _, gi in G])
par = {c: c for c in set(loc)}
def find(x_):
    while par[x_] != x_: par[x_] = par[par[x_]]; x_ = par[x_]
    return x_
for A, B, Cs in pairs:
    a1, a2, b1, b2 = A[0], A[1], B[0], B[1]
    dadb = np.zeros(d, np.int64)
    for (x_, y_), sg in (((a1, b1), 1), ((a1, b2), -1), ((a2, b1), -1), ((a2, b2), 1)): dadb[idx[('p', (x_, y_))]] += sg
    for c in Cs:
        for s in range(1, p):
            phi = np.eye(d, dtype=np.int64); phi[:, idx[('p', (c,))]] += s * dadb
            phii = np.eye(d, dtype=np.int64); phii[:, idx[('p', (c,))]] -= s * dadb
            for i, P in enumerate(greps):
                Q = phi @ P @ phii % p
                orb = (Gm @ Q @ Gi) % p                                      # its G-orbit; find a stored representative
                j = next(gk[key(o)] for o in orb if key(o) in gk)
                a_, b_ = find(loc[i]), find(loc[j])
                if a_ != b_: par[a_] = b_
n_loc = len(set(loc)); n_u = len({find(c) for c in set(loc)})
print(f"{name} F_{p}: local-torus orbits {n_loc}; with U = {{c ↦ c + s·δ_aδ_b}} as well: {n_u}")
