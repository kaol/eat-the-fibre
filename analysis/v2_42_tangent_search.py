"""v2.42: targeted family search in given signature groups, along inner directions x ∈ rad found by a tangent computation.

For each local orbit in the target groups and each of its G-orbit representatives that lifts to an integer fold
(symmetric residues):
  L = {x ∈ rad : d/dt|₀ Inn(1 + t·x) P Inn⁻¹ satisfies the linearised Δ-compatibility}  (exact, over ℚ);
then, for small integer combinations of a basis of L, the full curve P_t = Inn(1 + t·x) P Inn(1 + t·x)⁻¹ is tested as a
fold over ℚ[t] and its 𝔽_p- and 𝔽_q-points are located in the merged censuses.
Usage: python3 analysis/v2_42_tangent_search.py REPS_p.pkl REPS_q.pkl SIGNATURE_SUBSTRING...
"""
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
p, q = gp['p'], gq['p']
M, DL, d, basis = gp['M'], gp['DL'], gp['d'], gp['basis']
rad = [i for i, b in enumerate(basis) if b[0] == 'p']
one = [1 if b[0] == 'v' else 0 for b in basis]
t = sp.Symbol('t')
def smul(x, y): return [sp.expand(sum(int(M[i, j, m]) * x[i] * y[j] for i in range(d) for j in range(d) if M[i, j, m])) for m in range(d)]
E = lambda j: [1 if i == j else 0 for i in range(d)]
def admat(x):   # y ↦ x y − y x
    return sp.Matrix([[a - b for a, b in zip(smul(x, E(j)), smul(E(j), x))] for j in range(d)]).T
def innmat(x):
    u = [one[i] + t * x[i] for i in range(d)]; uinv = one[:]; pw = one[:]
    for _ in range(d):
        pw = smul(pw, [-t * y for y in x])
        if all(z == 0 for z in pw): break
        uinv = [sp.expand(a + b) for a, b in zip(uinv, pw)]
    return (sp.Matrix([smul(smul(u, E(j)), uinv) for j in range(d)]).T,)
DLs = [sp.Matrix(DL[m].tolist()) for m in range(d)]
def lin_delta(P, X):
    return [sum((X[m, l] * DLs[m] for m in range(d)), sp.zeros(d, d)) - X * DLs[l] * P.T - P * DLs[l] * X.T for l in range(d)]
def is_fold_family(Ps):
    for l in range(d):
        lhs = sum((Ps[m, l] * DLs[m] for m in range(d)), sp.zeros(d, d))
        if (lhs - Ps * DLs[l] * Ps.T).applyfunc(sp.expand) != sp.zeros(d, d): return False
    return True
r = gp['r']; sig = gp['sig_of_loc']
tlocs = [l for l, s in sig.items() if any(tg in str(s) for tg in targets)]
print(f"target orbits at q={p}: {len(tlocs)}", flush=True)
xs = sp.symbols(f'x0:{len(rad)}')
found = []
for l in sorted(tlocs):
    reps = [i for i, li in enumerate(r['loc']) if li == l]
    for i in reps:
        Pz = gp['lifted'](i)
        if Pz is None: continue
        P = sp.Matrix(Pz.tolist())
        xvec = [0] * d
        for k_, j in enumerate(rad): xvec[j] = xs[k_]
        A = admat(xvec); X = A * P - P * A
        eqs = [e for Mx in lin_delta(P, X) for e in Mx if e != 0]
        Lsol = sp.linsolve(eqs, xs) if eqs else None
        sol = list(Lsol)[0] if Lsol else xs
        free = sorted(set().union(*[sp.sympify(c_).free_symbols for c_ in sol]), key=str)
        Lbasis = [[sp.sympify(c_).subs({f: (1 if f == f0 else 0) for f in free}) for c_ in sol] for f0 in free]
        # drop directions whose curve is constant
        tried = 0
        for coeffs in itertools.product(range(-1, 2), repeat=len(Lbasis)):
            if not any(coeffs) or next(c_ for c_ in coeffs if c_) < 0: continue
            xv = [0] * d
            for c_, bvec in zip(coeffs, Lbasis):
                for k_, j in enumerate(rad): xv[j] += c_ * bvec[k_]
            if any(not sp.sympify(v).is_integer for v in xv): continue
            Inn = innmat(xv)[0]; Innv = innmat([-v for v in xv])[0]
            Ps = (Inn * P * Innv).applyfunc(sp.expand)
            if Ps == P or not is_fold_family(Ps): continue
            tried += 1
            hp = {gp['key2loc'].get(gp['key'](np.array(Ps.subs(t, s).tolist(), dtype=np.int64) % p)) for s in range(p)}
            hq = {gq['key2loc'].get(gq['key'](np.array(Ps.subs(t, s).tolist(), dtype=np.int64) % q)) for s in range(q)}
            if len(hp) >= 3 and not any(hp & f for f, _, _ in found):
                found.append((hp, hq, (i, xv)))
                print(f"curve from rep {i} (orbit {l}): x = {[(basis[j], xv[j]) for j in rad if xv[j]]}; "
                      f"meets {len(hp)} orbits at q={p}, {len(hq)} at q={q}; classes {sorted({str(sig[z][0]) for z in hp})}", flush=True)
        break   # one liftable representative per orbit
gen = [f for f in found if len(f[0]) == p and len(f[1]) == q]
print(f"genuine (q orbits at both primes): {len(gen)}; total curves kept: {len(found)}")
