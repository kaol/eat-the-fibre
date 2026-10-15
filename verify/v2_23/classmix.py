"""1⇉2, 1→3 over F_3: split compatibility (with the pair's class deletion Δ_{α}) against (r_G, r_H) for every fold.
Shows the family folds (r_H > r_G) are all split compatible, and the split failures are exactly the folds mixing the two classes
beyond ker Δ: σ(α-part of P(γ)) ≠ 0 (e.g. P(γ) = α₁ + α₂), or γ occurring in some P(α_i)."""
import sys, pickle, collections, numpy as np
sys.argv = ['x']
exec(open('pcensus.py').read().split("if __name__ == '__main__':")[0])
q, p = 'pair_srcout', 3
ctx = setup(q, p); C.update(ctx); basis, arrows, d = ctx['basis'], ctx['arrows'], ctx['d']
folds = [f.astype(np.int64) for f in pickle.load(open(f'logs/folds_{q}_{p}.pkl', 'rb'))['folds']]
idx = {b: i for i, b in enumerate(basis)}
D01 = np.zeros((d, d, d), dtype=np.int64)
for l, b in enumerate(basis):
    if b[0] == 'v': continue
    for k, a in enumerate(b[1]):
        if a not in (0, 1): continue
        pre = idx[('p', b[1][:k])] if k > 0 else idx[('v', arrows[a][0])]
        suf = idx[('p', b[1][k + 1:])] if k + 1 < len(b[1]) else idx[('v', arrows[a][1])]
        D01[l, pre, suf] += 1
ok = lambda P: all(((np.einsum('m,mab->ab', P[:, l], D01) - P @ D01[l] @ P.T) % p == 0).all() for l in range(d))
gens, E = lie_G(ctx); F = ctx['F']
def ranks(P):
    rows = [((D @ P - P @ D) % p).ravel() for D in gens]
    return F.rank(np.array(rows)), F.rank(np.array(rows + [((E @ P - P @ E) % p).ravel()]))
g_col = idx[('p', (2,))]; a_rows = [idx[('p', (0,))], idx[('p', (1,))]]
tab = collections.Counter()
for P in folds:
    rG, rH = ranks(P); s = ok(P)
    mixing = bool((P[a_rows, g_col].sum() % p) or P[g_col, a_rows].any())   # σ-part of P(γ) on α's ≠ 0, or γ in some P(α_i)
    tab[('split ok' if s else 'split fails', 'family' if rH > rG else 'no family', 'mixes beyond ker Δ' if mixing else 'no mixing')] += 1
print(sorted(tab.items()))
fam_split_ok = all(k[0] == 'split ok' for k in tab if k[1] == 'family')
fails_mixing = all(k[2] == 'mixes beyond ker Δ' for k in tab if k[0] == 'split fails')
exact = all((k[0] == "split fails") == (k[2] == "mixes beyond ker Δ") for k in tab)
print(f"family folds all split compatible: {fam_split_ok}; split failure ⇔ mixing beyond ker Δ: {exact}")
