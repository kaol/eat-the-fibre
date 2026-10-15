"""Show transversal (G × grading)-orbit representatives over Q for a quiver: folds whose H-orbit carries a 1-parameter
family of G-orbits. Prints P in path names, marking entries whose (output length − input length) < 0: those entries are
rescaled by the grading (λ^{Δlength}) and are the candidates for the modulus."""
import sys, pickle, collections
import numpy as np, sympy as sp
name, p, nshow = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 3
sys.argv = ['x']
exec(open('pcensus.py').read().split("if __name__ == '__main__':")[0])
ctx = setup(name, p); C.update(ctx)
folds = [f.astype(np.int64) for f in pickle.load(open(f'logs/folds_{name}_{p}.pkl', 'rb'))['folds']]
basis, arrows, d, M, DL = ctx['basis'], ctx['arrows'], ctx['d'], ctx['M'], ctx['DL']
aname = {}
cls = collections.defaultdict(list)
for ai, st in enumerate(arrows): cls[st].append(ai)
letters = iter('abcdefgh')
for st, c in cls.items():
    L = next(letters)
    for k, ai in enumerate(c): aname[ai] = f"{L}{k+1}" if len(c) > 1 else L
nm = [f"e{b[1]+1}" if b[0] == 'v' else "".join(aname[x] for x in b[1]) for b in basis]
lens = [0 if b[0] == 'v' else len(b[1]) for b in basis]
gens, E = lie_G(ctx)
def qfold(P):
    Pq = np.where(P == p - 1, -1, P).astype(np.int64); I = np.eye(d, dtype=np.int64)
    if (Pq @ Pq != Pq).any(): return None
    for i in range(d):
        for j in range(d):
            h = np.einsum('i,j,ijk->k', Pq[:, i], Pq[:, j], M); q = np.einsum('i,j,ijk->k', (I - Pq)[:, i], (I - Pq)[:, j], M)
            if (Pq @ h != h).any() or (Pq @ q).any(): return None
    for l in range(d):
        if (np.einsum('m,mab->ab', Pq[:, l], DL) != Pq @ DL[l] @ Pq.T).any(): return None
    return Pq
shown = 0; seen_kind = set()
for P in folds:
    Pq = qfold(P)
    if Pq is None: continue
    rows = [(D @ Pq - Pq @ D).ravel().tolist() for D in gens]
    rG = sp.Matrix(rows).rank(); rH = sp.Matrix(rows + [(E @ Pq - Pq @ E).ravel().tolist()]).rank()
    if rH <= rG: continue
    drops = tuple(sorted((nm[i], nm[j]) for i in range(d) for j in range(d) if Pq[i, j] and lens[i] < lens[j]))
    if drops in seen_kind: continue
    seen_kind.add(drops); shown += 1
    print(f"--- transversal fold (r_G = {rG}, r_H = {rH}); length-lowering entries (out ← in): {drops}")
    for j in range(d):
        terms = [f"{'' if Pq[i, j] == 1 else ('-' if Pq[i, j] == -1 else Pq[i, j])}{nm[i]}" for i in range(d) if Pq[i, j]]
        if terms: print(f"   P({nm[j]}) = {' + '.join(terms)}")
    if shown >= nshow: break
