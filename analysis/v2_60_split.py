"""v2.60 (Program 3): the split folds of Q₁ = Q₀ + ℓ (2→4) over the folds of Q₀, on census data.

A split fold is P = P̄ ⊕ P_J (P(kQ₀) ⊆ kQ₀ and P(J) ⊆ J). From ΔP = (P⊗P)Δ (by hand, doc §2):
  P_J(e₄) = λ e₄ with λ ∈ {0, 1}, and P_J(pℓ) = φ(p)ℓ + μ(p) e₄ for p ∈ B (paths of Q₀ ending at 2), with
  φ(p) = λ·P̄(p) + Σ_{Δp = Σ pre ⊗ suf} μ(suf)·P̄(pre)       (*)
This checks (*) on every split G-orbit representative, records (λ, μ = 0?) per census class, and checks that every
Q₀ class has the trivial lift P̄ ⊕ 0 among the split classes over it.
Usage: python3 analysis/v2_60_split.py q"""
import os, sys, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE)); sys.path.insert(0, HERE)
import numpy as np
src = open(os.path.join(HERE, 'v2_60_leaf_census.py')).read().split("stats = collections.Counter()")[0]
exec(src)
idx1 = {b: i for i, b in enumerate(b1)}; idx0 = {b: i for i, b in enumerate(b0)}
e4 = idx1[('v', 3)]
B = [b for b in b0 if (b[0] == 'v' and b[1] == 1) or (b[0] == 'p' and c0['arrows'][b[1][-1]][1] == 1)]   # paths ending at vertex 2
Bl = [idx1[('p', ((b[1],) if False else ()) + (5,))] if b[0] == 'v' else idx1[('p', b[1] + (5,))] for b in B]
B0 = [idx0[b] for b in B]
ok = 0; bad = 0; per = collections.defaultdict(set); trivial = collections.defaultdict(bool)
for P, lab in zip(R1['greps'], R1['loc']):
    P = P.astype(np.int64) % q
    if P[np.ix_(Qx, Jx)].any() or P[np.ix_(Jx, Qx)].any(): continue
    Pb = P[np.ix_(Qx, Qx)]
    lam = int(P[e4, e4]); col_e4_ok = not P[[i for i in Jx if i != e4], e4].any()
    mu = {B0[k]: int(P[e4, Bl[k]]) for k in range(len(B))}
    good = col_e4_ok and lam in (0, 1)
    for k, b in enumerate(B):
        phi = np.zeros(d0, dtype=np.int64)
        for k2, bl2 in enumerate(Bl): phi[B0[k2]] = P[bl2, Bl[k]]
        rhs = lam * Pb[:, B0[k]]
        Dl = np.asarray(DL0[B0[k]])
        for a in range(d0):
            for bb in range(d0):
                if Dl[a, bb] and bb in mu: rhs = rhs + Dl[a, bb] * mu[bb] * Pb[:, a]
        if ((rhs - phi) % q).any(): good = False
    ok += good; bad += not good
    per[lab].add((lam, any(mu.values())))
    if lam == 0 and not any(mu.values()) and not P[np.ix_(Jx, Jx)].any(): trivial[lab] = True
key = lambda lab: per[lab]
types = collections.Counter(frozenset(v) for v in per.values())
print(f"q = {q}: split G-orbit representatives satisfying (*): {ok}, violating: {bad}")
print(f"split classes {len(per)}; by (λ, μ ≠ 0) of their G-orbits: " + ", ".join(f"{sorted(k)}: {v}" for k, v in sorted(types.items(), key=lambda kv: -kv[1])))
print(f"classes containing a trivial lift P̄ ⊕ 0: {sum(trivial.values())}")
# L4: λ = 1, μ = 0 split lifts exist exactly over the Q₀ classes with P̄(B) ⊆ B (B = paths ending at vertex 2)
notB = [i for i in range(d0) if i not in B0]
pres = {}
for P, lab in zip(R0['greps'], R0['loc']):
    P = P.astype(np.int64) % q
    pres.setdefault(lab, set()).add(not P[np.ix_(notB, B0)].any())
nB = sum(1 for v in pres.values() if v == {True}); mixed = sum(1 for v in pres.values() if len(v) > 1)
lam1mu0 = sum(1 for v in per.values() if v == {(1, False)})
print(f"Q₀ classes with P̄(B) ⊆ B: {nB} (classes mixed: {mixed}); split classes with λ = 1, μ = 0: {lam1mu0}")
