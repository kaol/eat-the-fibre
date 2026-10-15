"""Unique extension as linear algebra.  Conicity: conjugation by λ^{length} preserves Δ-compatible folds (Δ has degree −1)
and fixes P₀, acting on the unknowns by a_x ↦ λ^{−|x|} a_x.  By induction on length the fibre equations at length m are
LINEAR in a_m.  Test: do the vertex-product equations alone force a_m = 0?
  (I)  a(Tz) = (1−P₀) a(z)
  (HL) a(h·Tz) = (1−P₀)(h·a(z)),  (HR) a(Tz·h) = (1−P₀)(a(z)·h)          h ∈ H₀ = im P₀
  (QL) a(q·(z−Tz)) = P₀(q·a(z)),   (QR) a((z−Tz)·q) = P₀(a(z)·q)           q ∈ Q₀ = ker P₀"""
import sys
from sage.all__sagemath_singular import QQ, matrix, vector, identity_matrix
from quiver import build
from extend import vertex_folds
from arrowrule import QS
def system(nv, arrows, P0, m, use=('I', 'HL', 'HR', 'QL', 'QR', 'HH', 'QQ')):
    d, mult, Dl, allp = build(nv, arrows)
    plen = lambda p: 0 if p[0] != 'a' else len(p) - 1
    Am = [i for i, p in enumerate(allp) if plen(p) == m]
    if not Am: return None
    src = lambda p: p[0] if p[0] != 'a' else arrows[p[1]][0]
    tgt = lambda p: p[0] if p[0] != 'a' else arrows[p[-1]][1]
    vseq = lambda p: [src(p)] + [arrows[k][1] for k in p[1:]] if p[0] == 'a' else [p[0]]
    index = {tuple(vseq(allp[i])): j for j, i in enumerate(Am)}
    N = len(Am)
    # T on A_m in the A_m basis: P0^{⊗(m+1)} on vertex sequences
    T = matrix(QQ, N, N)
    import itertools
    for j, i in enumerate(Am):
        vs = vseq(allp[i])
        supp = [[a for a in range(nv) if P0[a, v] != 0] for v in vs]
        for seq in itertools.product(*supp):
            c = 1
            for a, v in zip(seq, vs): c *= P0[a, v]
            T[index[tuple(seq)], j] += c
    S = [src(allp[i]) for i in Am]; Tg = [tgt(allp[i]) for i in Am]
    L = lambda f: matrix.diagonal(QQ, [f[S[j]] for j in range(N)])      # left mult by a vertex function on A_m
    R = lambda f: matrix.diagonal(QQ, [f[Tg[j]] for j in range(N)])
    I0 = identity_matrix(QQ, nv); Pm = matrix(QQ, P0)
    H0 = Pm.column_space().basis(); Q0 = (I0 - Pm).column_space().basis()
    # unknown a = n × N matrix; vectorise column-major: a(z_j) = column j. Equation "a∘M − K∘a = 0" ↦ (Mᵀ ⊗ I_n − I_N ⊗ K)
    rows = []
    def eq(M, K):                                    # a·M − K·a  with a: n×N
        rows.append(M.transpose().tensor_product(identity_matrix(QQ, nv)) - identity_matrix(QQ, N).tensor_product(K))
    if 'I' in use: eq(T, I0 - Pm)
    if 'ZQ' in use: eq(identity_matrix(QQ, N) - T, matrix(QQ, nv, nv))     # restrict a to the Hodge part: a∘(1−T) = 0
    if 'ZH' in use: eq(T, matrix(QQ, nv, nv))                              # restrict a to the Podge part: a∘T = 0
    for h in H0:
        D = matrix.diagonal(QQ, list(h))
        if 'HL' in use: eq(L(h) * T, (I0 - Pm) * D)
        if 'HR' in use: eq(T * R(h), (I0 - Pm) * D)
    for q in Q0:
        D = matrix.diagonal(QQ, list(q))
        if 'QL' in use: eq(L(q) * (identity_matrix(QQ, N) - T), Pm * D)
        if 'QR' in use: eq((identity_matrix(QQ, N) - T) * R(q), Pm * D)
    if 'HH' in use or 'QQ' in use:                   # a(Tx·Ty) = 0 and a((x−Tx)(y−Ty)) = 0 for |x|, |y| ≥ 1, |x|+|y| = m
        Tfull = {}
        def Tvec(i):                                 # T on a path index i (any length) as a dict over path indices
            vs = vseq(allp[i]); supp = [[a for a in range(nv) if P0[a, v] != 0] for v in vs]; out = {}
            byseq = {tuple(vseq(p)): k for k, p in enumerate(allp)}
            for seq in itertools.product(*supp):
                c = 1
                for a, v in zip(seq, vs): c *= P0[a, v]
                out[byseq[tuple(seq)]] = out.get(byseq[tuple(seq)], 0) + c
            return out
        def prod(u, v):
            out = {}
            for i, ci in u.items():
                for j, cj in v.items():
                    for k in range(d):
                        if mult[i, j, k]: out[k] = out.get(k, 0) + ci * cj * int(mult[i, j, k])
            return out
        for i, p in enumerate(allp):
            for j, q in enumerate(allp):
                if plen(p) >= 1 and plen(q) >= 1 and plen(p) + plen(q) == m and tgt(p) == src(q):
                    tx, ty = Tvec(i), Tvec(j)
                    for kind in ('HH', 'QQ'):
                        if kind not in use: continue
                        if kind == 'HH': w = prod(tx, ty)
                        else:
                            ux = dict((k, -c) for k, c in tx.items()); ux[i] = ux.get(i, 0) + 1
                            uy = dict((k, -c) for k, c in ty.items()); uy[j] = uy.get(j, 0) + 1
                            w = prod(ux, uy)
                        vec = vector(QQ, N)
                        for k, c in w.items():
                            if plen(allp[k]) == m: vec[Am.index(k)] += c
                        rows.append(vec.row().tensor_product(identity_matrix(QQ, nv)))   # a(w) = 0
    M = matrix(QQ, 0, nv * N)
    for r in rows: M = M.stack(r)
    return nv * N - M.rank()
if __name__ == '__main__':
    use = tuple(sys.argv[2].split(',')) if len(sys.argv) > 2 else ('I', 'HL', 'HR', 'QL', 'QR', 'HH', 'QQ')
    for name in sys.argv[1].split(';'):
        nv, arrows = QS[name]; A = set(arrows); worst = 0; cnt = 0
        for P0 in vertex_folds(nv):
            if not all(not [a for a in range(nv) if P0[a, s]] or not [b for b in range(nv) if P0[b, t]] or
                       all((a, b) in A for a in range(nv) if P0[a, s] for b in range(nv) if P0[b, t]) for s, t in arrows): continue
            cnt += 1
            for m in range(1, nv):
                k = system(nv, arrows, P0, m, use)
                if k is not None: worst = max(worst, k)
        print(f"{name}: {cnt} arrow-rule vertex folds; max nullity of the vertex-product system over all lengths, using {use}: {worst}", flush=True)
