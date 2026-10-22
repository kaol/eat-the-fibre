"""v2.42: family search along shears σ_w : c ↦ c + t·w (w ∈ ⟨a_i b_j⟩), optionally composed with Inn(1 + t·x),
x ∈ {δ_a, δ_b, δ_aδ_b}, on 1⇉2⇉3 + 1→3. σ_w is an algebra automorphism (c composes with nothing; w ∈ e₁ kQ e₃) but
not Δ-compatible unless w ∈ k·δ_aδ_b. Curves P_t = φ_t P φ_t⁻¹ are tested as folds over ℚ[t] and located in the
merged censuses at two primes.
Usage: python3 analysis/v2_42_shear_search.py REPS_p.pkl REPS_q.pkl [SIGNATURE_SUBSTRING...]"""
import os, sys, io, contextlib, itertools
import numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
rpf, rqf = sys.argv[1:3]; targets = sys.argv[3:]
def load(f):
    sys.argv = ['family_curves.py', f]
    g = {'__name__': 'fc', '__file__': os.path.join(HERE, 'family_curves.py')}
    with contextlib.redirect_stdout(io.StringIO()):
        exec(open(os.path.join(HERE, 'family_curves.py')).read().split("fams = []; nlift")[0], g)
    return g
gp, gq = load(rpf), load(rqf)
p, q, d, basis, t = gp['p'], gq['p'], gp['d'], gp['basis'], gp['t']
idx = {b: i for i, b in enumerate(basis)}
arrows = gp['arrows']
cidx = [i for i, b in enumerate(basis) if b[0] == 'p' and len(b[1]) == 1 and arrows[b[1][0]] == (0, 2)]
assert len(cidx) == 1; ci = cidx[0]
ab = [i for i, b in enumerate(basis) if b[0] == 'p' and len(b[1]) == 2]
def shear(w, s):
    S = sp.eye(d)
    for j, wj in zip(ab, w): S[j, ci] += s * wj
    return S
r = gp['r']; sig = gp['sig_of_loc']
tlocs = sorted(l for l, s in sig.items() if not targets or any(tg in str(s) for tg in targets))
SH, NSH = (int(x) for x in os.environ.get('SHARD', '0/1').split('/'))
tlocs = tlocs[SH::NSH]
print(f"orbits searched at q={p}: {len(tlocs)}", flush=True)
import time; T0 = time.time()
ws = [w for w in itertools.product(range(-1, 2), repeat=4) if any(w) and next(c for c in w if c) > 0
      and not (w[0] == -w[1] == -w[2] == w[3])]                      # exclude multiples of δ_aδ_b (census symmetry)
inns = [None] + list(gp['INN'].items())
BIG = 1000003
DLn = np.array(gp['DL'], dtype=object)
def numfold(Pn):
    """Δ-compatibility of a numeric matrix mod BIG (image/kernel closure hold along automorphism orbits)"""
    Pn = np.array(Pn, dtype=object) % BIG
    for l in range(d):
        lhs = sum(Pn[m, l] * DLn[m] for m in range(d))
        if ((lhs - Pn.dot(DLn[l]).dot(Pn.T)) % BIG).any(): return False
    return True
tsamp = [3, 17, 101, 9999]
INNnum = {vec: [(np.array(Inn.subs(t, tv).tolist(), dtype=object) % BIG, np.array(Innv.subs(t, tv).tolist(), dtype=object) % BIG)
                for tv in tsamp] for vec, (Inn, Innv) in gp['INN'].items()}
def shear_num(w, s_):
    S_ = np.eye(d, dtype=object)
    for j, wj in zip(ab, w): S_[j, ci] = (S_[j, ci] + s_ * wj) % BIG
    return S_
def num_curve(w, vec, P):
    Pn = np.array(P.tolist(), dtype=object) % BIG
    for k_, tv in enumerate(tsamp):
        A = shear_num(w, tv); B = shear_num(w, -tv)
        if vec is not None:
            I_, Iv = INNnum[vec][k_]; A = I_.dot(A) % BIG; B = B.dot(Iv) % BIG
        if not numfold(A.dot(Pn).dot(B) % BIG): return False
    return True
found = []
for l in tlocs:
    i = next((i for i, li in enumerate(r['loc']) if li == l and gp['lifted'](i) is not None), None)
    if i is None: continue
    P = sp.Matrix(gp['lifted'](i).tolist())
    print(f"  orbit {l} rep {i} ({time.time() - T0:.0f}s)", flush=True)
    for w in ws:
        for inn in inns:
            if not num_curve(w, inn[0] if inn else None, P): continue
            Phi = shear(w, t); Phiv = shear(w, -t)
            if inn is not None:
                (vec, (Inn, Innv)) = inn
                Phi = Inn * Phi; Phiv = Phiv * Innv
            Ps = (Phi * P * Phiv).applyfunc(sp.expand)
            if Ps == P or not gp['is_fold_family'](Ps): continue
            hp = {gp['key2loc'].get(gp['key'](np.array(Ps.subs(t, s).tolist(), dtype=np.int64) % p)) for s in range(p)}
            if len(hp) < 3 or any(hp & f[0] for f in found): continue
            hq = {gq['key2loc'].get(gq['key'](np.array(Ps.subs(t, s).tolist(), dtype=np.int64) % q)) for s in range(q)}
            found.append((hp, hq, (i, w, inn[0] if inn else None)))
            print(f"curve from rep {i} (orbit {l}, class {sig[l][0]}): w = {w}, inner {inn[0] if inn else None}; "
                  f"meets {len(hp)} orbits at q={p}, {len(hq)} at q={q}", flush=True)
gen = [f for f in found if len(f[0]) == p and len(f[1]) == q]
print(f"genuine (q orbits at both primes): {len(gen)}; curves kept: {len(found)}")
