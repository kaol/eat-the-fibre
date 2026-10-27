"""v2.61 (Program 3, item 1): the μ-cones of the split λ = 1 lifts, Q₀ → Q₁ = Q₀ + ℓ (2→4), over 𝔽_q.

For each G-orbit representative P̄ of Folds(𝔽_q Q₀) (census dump), and every μ ∈ B* (B = paths of Q₀ ending at 2,
dim 7), form P = P̄ ⊕ P_J with P_J(e₄) = e₄ and P_J(pℓ) = φ(p)ℓ + μ(p)e₄, φ from the normal form (*) of v2.60
(φ(p) = P̄(p) + Σ μ(suf)·P̄(pre)); discard μ for which some φ(p) ∉ B; test whether P is a fold of kQ₁ (vectorised).
M(P̄) = {μ : P is a fold}. Reports per census class of Q₀ with a μ ≠ 0: |M(P̄)|, whether M(P̄) is a linear subspace
(closed under addition, contains 0), its dimension, and whether it is spanned by coordinate functionals p*.
Usage: python3 analysis/v2_61_mu.py q"""
import os, sys, pickle, itertools, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
import numpy as np
import pcensus as pc
q = int(sys.argv[1])
class SafeU(pickle.Unpickler):
    OK = {('numpy._core.numeric', '_frombuffer'), ('numpy.core.numeric', '_frombuffer'), ('numpy', 'dtype'),
          ('numpy._core.multiarray', '_reconstruct'), ('numpy.core.multiarray', '_reconstruct'), ('numpy', 'ndarray')}
    def find_class(self, mod, name):
        if (mod, name) not in self.OK: raise pickle.UnpicklingError(f"refused global {mod}.{name}")
        return super().find_class(mod, name)
root = os.path.dirname(HERE)
R0 = SafeU(open(os.path.join(root, 'logs', f'reps_pairs_in_short_{q}.pkl'), 'rb')).load()
c1, c0 = pc.setup('pairs_in_short_tail', q), pc.setup('pairs_in_short', q)
b1, b0 = c1['basis'], c0['basis']; d1, d0 = len(b1), len(b0)
idx1 = {b: i for i, b in enumerate(b1)}; idx0 = {b: i for i, b in enumerate(b0)}
M1 = np.asarray(c1['M']); DL1 = np.asarray(c1['DL']); DL0 = np.asarray(c0['DL'])
e4 = idx1[('v', 3)]
B = [b for b in b0 if (b[0] == 'v' and b[1] == 1) or (b[0] == 'p' and c0['arrows'][b[1][-1]][1] == 1)]
B0 = [idx0[b] for b in B]; Bpos = {i: k for k, i in enumerate(B0)}
Bl = [idx1[('p', (5,))] if b[0] == 'v' else idx1[('p', b[1] + (5,))] for b in B]
Qx = [idx1[b] for b in b0]
nb = len(B)
name = lambda b: (f"e{b[1]+1}" if b[0] == 'v' else "".join(["a", "a'", "b", "b'", "c", "ℓ"][x] for x in b[1]))
def fold_mask(Ps):
    """Ps: (N, d, d) over 𝔽_q → boolean (N,) fold test for kQ₁"""
    d = d1; P = Ps % q
    ok = ~(((np.einsum('nij,njk->nik', P, P) - P) % q).any(axis=(1, 2)))       # staged: cheap tests first
    sel = np.nonzero(ok)[0]
    if len(sel):
        Q_ = P[sel]
        lhs = np.einsum('nml,mab->nlab', Q_, DL1) % q
        rhs = np.einsum('nax,lxy,nby->nlab', Q_, DL1, Q_, optimize=True) % q
        ok[sel] &= ~(((lhs - rhs) % q).any(axis=(1, 2, 3)))
    sel = np.nonzero(ok)[0]
    if len(sel):
        Q_ = P[sel]; I = np.eye(d, dtype=np.int64)
        for A, img in ((Q_, True), ((I - Q_) % q, False)):
            prod = np.einsum('nia,njb,ijk->nkab', A, A, M1, optimize=True) % q
            Pu = np.einsum('nlk,nkab->nlab', Q_, prod, optimize=True) % q
            bad = ((Pu - prod) % q if img else Pu).any(axis=(1, 2, 3))
            ok[sel] &= ~bad
    return ok
mus = np.array(list(itertools.product(range(q), repeat=nb)), dtype=np.int64)          # all μ ∈ 𝔽_q^7
out = collections.defaultdict(list); nreps = 0; nonempty_reps = 0
seen_lab = set()
for Pb, lab in zip(R0['greps'], R0['loc']):
    if lab in seen_lab and os.environ.get('ONE_PER_CLASS') == '1': continue
    seen_lab.add(lab)
    Pb = Pb.astype(np.int64) % q; nreps += 1
    # φ for every μ: φ(p) = P̄p + Σ_{(a,b)} DL0[p][a,b] μ(b) P̄[:, a]   (b ∈ B)
    base = Pb[:, B0]                                                    # (d0, nb): P̄ p
    corr = np.zeros((len(mus), d0, nb), dtype=np.int64)
    for k, pi in enumerate(B0):
        Dl = DL0[pi]
        for a, bb in zip(*np.nonzero(Dl)):
            if bb in Bpos: corr[:, :, k] += int(Dl[a, bb]) * np.outer(mus[:, Bpos[bb]], Pb[:, a])
    phi = (base[None] + corr) % q                                       # (N, d0, nb)
    notB = [i for i in range(d0) if i not in Bpos]
    inB = ~phi[:, notB, :].any(axis=(1, 2))
    cand = np.nonzero(inB)[0]
    if not len(cand): continue
    Ps = np.zeros((len(cand), d1, d1), dtype=np.int64)
    Ps[:, np.ix_(Qx, Qx)[0], np.ix_(Qx, Qx)[1]] = Pb[None]
    Ps[:, e4, e4] = 1
    for k in range(nb):
        for k2 in range(nb): Ps[:, Bl[k2], Bl[k]] = phi[cand, B0[k2], k]
        Ps[:, e4, Bl[k]] = mus[cand, k]
    good = []
    for s in range(0, len(cand), 2000):
        m = fold_mask(Ps[s:s + 2000]); good += [int(cand[s + i]) for i in np.nonzero(m)[0]]
    Mset = {tuple(mus[i]) for i in good}
    if not (Mset - {tuple([0] * nb)}): continue                     # only the identity lift (or nothing)
    out[lab].append((Pb, Mset)); nonempty_reps += 1
print(f"q = {q}: {nreps} G-orbit representatives of Q₀ ({nonempty_reps} with a μ ≠ 0 lift); {len(out)} census classes with a μ ≠ 0 lift")
summ = collections.Counter()
for lab, lst in sorted(out.items()):
  for Pb, Mset in lst:
      zero = tuple([0] * nb) in Mset
      m0 = np.array(next(iter(Mset)))
      shifted = {tuple((np.array(x) - m0) % q) for x in Mset}
      affine = all(tuple((np.array(x) + np.array(y)) % q) in shifted for x in shifted for y in shifted)
      closed = all(tuple((np.array(x) + np.array(y)) % q) in Mset for x in Mset for y in Mset)
      vecs = np.array(sorted(Mset))
      rank = c0['F'].rank(vecs) if len(vecs) else 0
      coord = all(sum(1 for t in x if t) <= 1 for x in Mset) or closed and all(
          any(all(t == 0 for j, t in enumerate(x) if j != k) for x in Mset if x[k]) for k in range(nb) if any(x[k] for x in Mset))
      supp = sorted({name(B[k]) for x in Mset for k in range(nb) if x[k]})
      summ[(zero, closed, affine, rank, coord)] += 1
      print(f"  class {lab}: |M| = {len(Mset)}, 0 ∈ M: {zero}, linear: {closed}, affine: {affine}, dim {rank}, coordinate span: {coord}, support {supp}")
print("summary (0 ∈ M, linear, affine, dim of span, coordinate span):", dict(summ))
