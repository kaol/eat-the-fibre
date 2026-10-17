"""Theorem 3 (v2.27): two parallel classes with a common target (or a common source) give a one-parameter family of
folds whose invariant is a vertex coefficient.

Let A: u → v and B: w → v be parallel classes with |A|, |B| ≥ 2 and u ≠ w (converging), with no arrow between u and w
and no loop at u or w. δ_a = a₁ − a₂, δ_b = b₁ − b₂. For s ∈ k put z_s = −s·e_v + s·δ_a + δ_b and
    P_s(e_u) = e_u + δ_a + z_s,  P_s(e_w) = −P_s(e_u),  P_s(a₂) = z_s,  P_s(b₂) = −(1 + s)·z_s,  P_s = 0 on every other path.
Diverging classes (A: v → u, B: v → w) take the same formula (the fold conditions are invariant under A ↦ A^op).
Invariant: the e_v-coefficient of P(e_u) is −s. Every h ∈ G × T fixes each vertex and maps paths of positive length into
their span, so (hPh⁻¹)(e_u) ≡ P(e_u) modulo arrows: the whole vertex block P₀₀ is an invariant.
Checks:
 1. P_s is a Δ-compatible fold for all s (symbolic over Q(s)), on quivers satisfying the hypothesis and on controls.
 2. Lie(Aut_Δ): every derivation D with ΔD = (D⊗1 + 1⊗D)Δ kills the vertices (so the identity component of Aut_Δ, not
    only G × T, fixes the vertex block), and dim Lie(Aut_Δ) = dim Lie(G × T).
 3. Over F_p (folds files): every P_s is a census fold; the vertex block takes finitely many values except for the
    parametric ones; orbit counts split by vertex block.
Usage: python3 analysis/theorem3.py"""
import os as _os, sys as _sys
if 'ROOT' not in globals():
    AN = _os.path.dirname(_os.path.abspath(__file__)); ROOT = _os.path.dirname(AN); _os.chdir(ROOT)
    __file__ = _os.path.join(ROOT, 'pcensus.py')
import sys, collections, itertools, pickle, numpy as np, sympy as sp
sys.argv = ['x']
exec(open(_os.path.join(ROOT, 'pcensus.py')).read().split("if __name__ == '__main__':")[0])
s = sp.Symbol('s')

c_, x_ = sp.symbols('c x')

def members(kind_P, eu, ev, ew, pa1, pa2, pb1, pb2, d):
    """the three families on converging (or diverging) classes A ∋ a₁, a₂ (u–v) and B ∋ b₁, b₂ (w–v); P = 0 on every
    basis path not listed. Vertex blocks (the invariants):
      F1(s): P(e_u) = e_u − s·e_v + …,                P(e_w) = −P(e_u)
      F2(c): P(e_u) = c·e_u + (c − 1)·e_w + …,         P(e_w) = −P(e_u)
      F4(x): P(e_u) = e_u + x·e_v + …,                 P(e_w) = e_w + (1 − x)·e_v + …
    (F3 is F1 with u and w exchanged.)"""
    def vec(u=0, v=0, w=0, da=0, db=0):
        z = sp.zeros(d, 1); z[eu] += u; z[ev] += v; z[ew] += w; z[pa1] += da; z[pa2] -= da; z[pb1] += db; z[pb2] -= db
        return z
    P = sp.zeros(d, d)
    if kind_P == 'F1':
        z = vec(0, -s, 0, s, 1); P[:, eu] = vec(1, -s, 0, s + 1, 1); P[:, ew] = -P[:, eu]; P[:, pa2] = z; P[:, pb2] = -(1 + s) * z
    elif kind_P == 'F2':
        c = c_; z = vec(1, 0, 1, c - 1, -c); P[:, eu] = vec(c, 0, c - 1, c * (c - 1), -c * (c - 1)); P[:, ew] = -P[:, eu]
        P[:, pa2] = z; P[:, pb2] = z
    elif kind_P == 'F4':
        x = x_; z = vec(0, 1, 0, -x, x - 1); P[:, eu] = vec(1, x, 0, x * (1 - x), -x * (1 - x))
        P[:, ew] = vec(0, 1 - x, 1, -x * (1 - x), x * (1 - x)); P[:, pa2] = z; P[:, pb2] = z
    return P

def family(name, kinds=('F1', 'F2', 'F4')):
    """all converging/diverging pairs of classes (A, B) with |A|, |B| ≥ 2 sharing an endpoint v, u ≠ w; the families"""
    nv, arrows = QUIVERS[name]; basis, M, DL = build(nv, arrows); d = len(basis); idx = {b: i for i, b in enumerate(basis)}
    cls = collections.defaultdict(list)
    for a, st in enumerate(arrows): cls[st].append(a)
    out = []
    for (A_st, A), (B_st, B) in itertools.permutations(cls.items(), 2):
        if len(A) < 2 or len(B) < 2: continue
        if A_st[1] == B_st[1] and A_st[0] != B_st[0]: kind, u, v, w = 'in', A_st[0], A_st[1], B_st[0]
        elif A_st[0] == B_st[0] and A_st[1] != B_st[1]: kind, u, v, w = 'out', A_st[1], A_st[0], B_st[1]
        else: continue
        hyp = not any((x, y) in ((u, w), (w, u), (u, u), (w, w)) for x, y in arrows)
        ev, eu, ew = idx[('v', v)], idx[('v', u)], idx[('v', w)]
        pa1, pa2, pb1, pb2 = (idx[('p', (x,))] for x in (A[0], A[1], B[0], B[1]))
        for kp in kinds:
            out.append(((kind, u + 1, v + 1, w + 1, kp), hyp, members(kp, eu, ev, ew, pa1, pa2, pb1, pb2, d)))
    return basis, M, DL, d, out

def is_fold(P, M, DL, d):
    Id = sp.eye(d)
    def mul(xv, yv):
        o = [0] * d
        for i in range(d):
            if xv[i] == 0: continue
            for j in range(d):
                if yv[j] == 0: continue
                for kk in np.nonzero(M[i, j])[0]: o[kk] += xv[i] * yv[j] * int(M[i, j, kk])
        return sp.Matrix(o)
    ok = sp.expand(P * P - P) == sp.zeros(d, d)
    for i in range(d):
        for j in range(d):
            h = mul(P[:, i], P[:, j]); q = mul((Id - P)[:, i], (Id - P)[:, j])
            ok &= sp.expand(P * h - h) == sp.zeros(d, 1) and sp.expand(P * q) == sp.zeros(d, 1)
    for l in range(d):
        lhs = sum((P[m, l] * sp.Matrix(DL[m].tolist()) for m in range(d)), sp.zeros(d, d))
        ok &= sp.expand(lhs - P * sp.Matrix(DL[l].tolist()) * P.T) == sp.zeros(d, d)
    return bool(ok)

def lie_aut_delta(name):
    """dim of {D : D derivation of kQ, ΔD = (D⊗1 + 1⊗D)Δ}, and whether all such D kill the vertices (exact, over Q)"""
    nv, arrows = QUIVERS[name]; basis, M, DL = build(nv, arrows); d = len(basis)
    Ds = sp.symbols(f'D0:{d * d}'); D = sp.Matrix(d, d, Ds)
    eqs = []
    for i in range(d):
        for j in range(d):
            prod = sp.Matrix(M[i, j].tolist())                       # b_i b_j
            lhs = D * prod
            rhs = sum((D[k, i] * sp.Matrix(M[k, j].tolist()) + D[k, j] * sp.Matrix(M[i, k].tolist()) for k in range(d)), sp.zeros(d, 1))
            eqs += list(lhs - rhs)
    for l in range(d):                                              # Δ(D b_l) = (D⊗1 + 1⊗D) Δ(b_l);  Δ(b)_ab = DL[b][a, b]
        lhs = sum((D[m, l] * sp.Matrix(DL[m].tolist()) for m in range(d)), sp.zeros(d, d))
        rhs = D * sp.Matrix(DL[l].tolist()) + sp.Matrix(DL[l].tolist()) * D.T
        eqs += list(lhs - rhs)
    A_, _ = sp.linear_eq_to_matrix([e for e in eqs if e != 0], Ds)
    N = A_.nullspace()
    vert = [i for i, b in enumerate(basis) if b[0] == 'v']
    kills = all(all(v_[r * d + c] == 0 for r in range(d) for c in vert) for v_ in N)   # D[:, vertex] = 0
    cls = collections.Counter(arrows)
    dimG = sum(m * (m - 1) for m in cls.values() if m > 1)        # Lie G: column-sum-0 matrices on each class
    return len(N), dimG, kills

def vertex_blocks(name, p, orbits=True):
    """over F_p: P_s ∈ census; distinct vertex blocks P₀₀; H-orbits (H = G × the torus elements preserving the fold
    set) grouped by vertex block"""
    ctx = setup(name, p); C.update(ctx); basis, arrows, d = ctx['basis'], ctx['arrows'], ctx['d']
    folds = [f.astype(np.int64) % p for f in pickle.load(open(f'logs/folds_{name}_{p}.pkl', 'rb'))['folds']]
    key = lambda A: (A % p).astype(np.int8).tobytes()
    FS = {key(f) for f in folds}
    vert = [i for i, b in enumerate(basis) if b[0] == 'v']
    _, _, _, _, fams = family(name)
    mem = {}
    fams = [f_ for f_ in fams if f_[0][1] == 1 or f_[0][4] == 'F1'][:4]          # F1, F2, F4 at u = 1; F1 at u = 3 (= F3)
    for k_, h_, Ps in fams:
        par = list(Ps.free_symbols)
        mem[k_[4] if k_[1] == 1 else 'F3'] = [np.array([[int(sp.sympify(x).subs({q_: sv for q_ in par})) % p for x in row] for row in Ps.tolist()], dtype=np.int64)
                      for sv in range(p)]
    in_census = {kp: all(key(P) in FS for P in L) for kp, L in mem.items()}
    VB = collections.Counter(key(f[np.ix_(vert, vert)]) for f in folds)
    out = dict(p=p, N=len(folds), families_in_census=in_census, vertex_blocks=len(VB))
    if orbits:
        G = build_G(ctx); cls = collections.defaultdict(list)
        for a, st in enumerate(arrows): cls[st].append(a)
        cl = list(cls.values()); T = []
        for w in itertools.product(range(1, p), repeat=len(cl)):
            wa = {a: w[k] for k, c in enumerate(cl) for a in c}
            dg = [1 if b[0] == 'v' else int(np.prod([wa[x] for x in b[1]])) % p for b in basis]
            t_, ti = np.diag(dg), np.diag([pow(x, p - 2, p) for x in dg])
            if all(key(t_ @ f @ ti) in FS for f in folds): T.append((t_, ti))
        H = [(t_ @ g % p, gi @ ti % p) for g, gi in G for t_, ti in T]
        lab, reps, sizes = orbit_label(folds, (np.array([h for h, _ in H]), np.array([hi for _, hi in H])), p)
        per_vb = collections.Counter(key(r[np.ix_(vert, vert)]) for r in reps)
        out.update(orbits=len(reps), family_orbits={kp: len({lab[key(P)] for P in L}) for kp, L in mem.items()},
                   family_orbit_sizes={kp: [sizes[lab[key(P)]] for P in L] for kp, L in mem.items()},
                   orbits_per_block_hist=dict(sorted(collections.Counter(per_vb.values()).items())))
        single = {k for k, v in per_vb.items() if v == 1}
        generic = {kp: {key(P[np.ix_(vert, vert)]) for P in L if per_vb[key(P[np.ix_(vert, vert)])] == 1} for kp, L in mem.items()}
        out['generic_members_per_family'] = {kp: len(v) for kp, v in generic.items()}
        out['single_orbit_blocks'] = len(single)
        out['single_orbit_blocks_are_generic_family_members'] = single == set().union(*generic.values())
        out['other_blocks_orbit_total'] = len(reps) - len(single)
        out['_per_vb'] = per_vb; out['_reps'] = reps; out['_sizes'] = sizes
    out['_VB'] = VB
    return out

if __name__ == '__main__':
    print("1. P_s is a Δ-compatible fold for all s (symbolic):")
    for q in ('pairs_in', 'pairs_out', 'pairs_in_tail', 'pairs_in_extra', 'pairs_in_pre', 'pairs_in_side',
              'triple_pair_in', 'pairs_in_short', 'two_pairs', 'two_pairs_sep'):
        _, M, DL, d, fams = family(q)
        r = [(k, h, is_fold(P, M, DL, d)) for k, h, P in fams]
        if not r: print(f"   {q:15s}: no converging or diverging pair of classes"); continue
        byk = collections.defaultdict(list)
        for k, h, ok in r: byk[(k[0], k[1], k[2], k[3], h)].append(f"{k[4]} {ok}")
        print(f"   {q:15s}: " + "; ".join(f"{k[0]} u={k[1]}, w={k[3]}, v={k[2]}, hypothesis {k[4]}: " + ", ".join(v_)
                                         for k, v_ in byk.items()))
    print("2. Lie(Aut_Δ) (Δ-derivations): dimension vs dim Lie(G); all kill the vertices (the class torus is not in")
    print("   Aut_Δ: it rescales Δ, and moves folds to folds only through the class-weight lemma):")
    for q in ('pairs_in', 'pairs_out', 'two_pairs', 'pairs_in_tail'):
        n, g, k = lie_aut_delta(q); print(f"   {q:15s}: dim {n}, dim Lie(G) {g}, vertices fixed: {k}")

    print("3. Over F_p (pairs_in): families in the census, vertex blocks, orbits (H = G × fold-preserving torus):")
    for p in (3, 5):
        r = vertex_blocks('pairs_in', p)
        print(f"   F_{p}: " + ", ".join(f"{k} = {v}" for k, v in r.items() if not k.startswith('_')), flush=True)
    r = vertex_blocks('pairs_in', 7, orbits=False)
    print(f"   F_7: " + ", ".join(f"{k} = {v}" for k, v in r.items() if not k.startswith('_')), flush=True)
