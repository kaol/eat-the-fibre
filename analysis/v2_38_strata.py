"""v2.38: stratify the local-torus orbits of 1⇉2⇉3 by the (G × T²)-invariant τ and locate the 19 families.

τ(P) = (vertex block εPι mod q, and for each vertex i and block X ∈ {A, B, AB}: [ε_i P|_X ≠ 0], [P|_X ι e_i ≠ 0]).
Invariance: ε∘h = ε, h∘ι = ι, h block-diagonal and invertible on A, B, AB for h ∈ G × T².
Usage: python3 analysis/v2_38_strata.py CERT.pkl REPS_q.pkl [REPS_q' ...]
"""
import sys, pickle, collections, json
import numpy as np, sympy as sp
cert = pickle.load(open(sys.argv[1], 'rb'))
Vvec = [[1] + [0] * 10, [0, 1] + [0] * 9, [0, 0, 1] + [0] * 8, [0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 1]]
BL = {'A': [3, 4], 'B': [5, 6], 'AB': [7, 8, 9, 10]}
def tau(P, q):
    P = np.asarray(P, np.int64) % q
    vb = tuple(int(x) for x in ((P[:3, :3] + q // 2) % q - q // 2).flatten())   # symmetric residues: comparable across q
    rows = tuple(int(P[i, BL[X]].any()) for i in range(3) for X in BL)
    cols = tuple(int(P[BL[X], i].any()) for i in range(3) for X in BL)
    return (vb, rows, cols)
tsym = sp.Symbol('t')
def famP(k, tv):
    v = cert[k]; P = np.zeros((11, 11), np.int64)
    for j in v['nz']:
        c = [int(sp.sympify(x).subs(tsym, tv)) for x in v['cols'][j]]
        for qq in range(11): P[qq, j] = sum(c[a] * Vvec[a][qq] for a in range(6))
    return P
out = {}
for f in sys.argv[2:]:
    r = pickle.load(open(f, 'rb')); q = r['p']
    assert [str(b) for b in r['basis']][:3] and len(r['basis']) == 11
    seen = {}
    for P, l in zip(r['greps'], r['loc']):
        tp = tau(P, q)
        if l in seen: assert seen[l] == tp, "τ not constant on a local orbit"
        seen[l] = tp
    strata = collections.Counter(seen.values())
    famin = collections.defaultdict(collections.Counter)
    for k in cert:
        for tv in range(q): famin[tau(famP(k, tv), q)][k] += 1
    out[q] = (strata, famin)
    print(f"q={q}: local orbits {len(seen)}, strata {len(strata)}, family points outside census strata "
          f"{sum(sum(c.values()) for s, c in famin.items() if s not in strata)}")
qs = sorted(out)
if len(qs) == 2:
    q0, q1 = qs; S0, F0 = out[q0]; S1, F1 = out[q1]
    keys = sorted(set(S0) | set(S1), key=str)
    print("strata only at one prime:", sum((s in S0) != (s in S1) for s in keys))
    bad = []; famstrata = 0
    for s in keys:
        nf = len(F1.get(s, {}))
        d_ = S1.get(s, 0) - S0.get(s, 0)
        if nf: famstrata += 1
        if d_ != (q1 - q0) * nf: bad.append((s, S0.get(s, 0), S1.get(s, 0), dict(F0.get(s, {})), dict(F1.get(s, {}))))
    print("strata containing family points:", famstrata, "; prediction-4 (as written) violations:", len(bad))
    for b in bad[:12]: print("  ", b)
    ref = [(s, S0.get(s, 0), S1.get(s, 0), dict(F0.get(s, {})), dict(F1.get(s, {}))) for s in keys
           if S1.get(s, 0) - S0.get(s, 0) != sum(F1.get(s, {}).values()) - sum(F0.get(s, {}).values())]
    print("refined check (growth of each stratum = growth of its family points) violations:", len(ref))
    for b in ref: print("  ", b)
    nonfam0 = {s: S0.get(s, 0) - sum(F0.get(s, {}).values()) for s in keys}
    print("non-family orbits: total", sum(nonfam0.values()), "in", sum(v > 0 for v in nonfam0.values()), "strata;",
          "strata mixing family and non-family orbits:", sum(1 for s in keys if nonfam0[s] > 0 and s in F1))
    print("max stratum size", max(S1.values()), "; singleton strata", sum(v == 1 for v in S1.values()))
    print("family strata detail:")
    for s in keys:
        if s in F1: print("  ", S0.get(s, 0), S1.get(s, 0), dict(F0.get(s, {})), dict(F1[s]))
