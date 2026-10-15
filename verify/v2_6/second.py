import sys, time, itertools
from fold import *
from tangent import spec, gvec
def analyse(kind, Lw, delta, Llow, show=True):
    S = spec(kind, Lw + 2 * delta + 2); F = Fold(S, Lw, delta)
    eqs = equations(F, Lw); J = linear_part(eqs)
    K, E = kernel(J, len(F.names))
    low = [j for j, key in enumerate(F.names) if key[0] == 'c' and F.wlev(key[1]) <= Llow]
    G = [gvec(F, h) for h in range(len(S.B)) if S.lev[h] <= delta]
    # basis of tangent mod G (on low coords): start echelon with G projections, add K vectors
    Ech = Echelon(); reps = []
    for g in G: Ech.add({j: g[j] for j in low if j in g})
    for v in K:
        if Ech.add({j: v[j] for j in low if j in v}): reps.append(v)
    print(f"{kind} Lw={Lw} δ={delta}: kernel={len(K)}, G={len(G)}, extras on low coords={len(reps)}")
    def describe(v):
        out = []
        for j in sorted(low):
            if j in v and v[j]:
                _, w, o = F.names[j]
                val = v[j] if v[j] < MOD // 2 else v[j] - MOD
                out.append(f"c{tuple(S.side[i][0]+str(S.lev[i]) for i in w)}→{S.side[o][0]}{S.lev[o]}:{val}")
        return out
    if show:
        for r in reps: print("   ", describe(r)[:6])
    return S, F, eqs, E, K, G, reps, low
def quadrics(F, eqs, taus):
    """obstruction quadrics in parameters t_a for tau = sum t_a taus[a]; eta ranges over all vars."""
    T = len(taus); nv = len(F.names)
    pairs = [(a, b) for a in range(T) for b in range(a, T)]
    aux = {p: nv + i for i, p in enumerate(pairs)}
    Ech = Echelon()
    for e in eqs:
        row = {k[0]: v for k, v in e.items() if len(k) == 1}
        quad = [(k, v) for k, v in e.items() if len(k) == 2]
        if quad:
            for (a, b) in pairs:
                s = 0
                for (i, j), cf in quad:
                    ta_i, tb_j = taus[a].get(i, 0), taus[b].get(j, 0)
                    s += cf * ta_i * tb_j
                    if a != b: s += cf * taus[b].get(i, 0) * taus[a].get(j, 0)
                s %= MOD
                if s: row[aux[(a, b)]] = s
        if row: Ech.add(row)
    Q = [r for c, r in Ech.piv.items() if c >= nv]
    return Q, pairs, nv
if __name__ == '__main__':
    kind = sys.argv[1]; Lw, delta, Llow = map(int, sys.argv[2:5])
    S, F, eqs, E, K, G, reps, low = analyse(kind, Lw, delta, Llow)
    taus = G + reps
    Q, pairs, nv = quadrics(F, eqs, taus)
    print(f"  {len(Q)} independent obstruction quadrics on span(G ∪ extras), T={len(taus)} (first {len(G)} = conjugation)")
    def fmt(r):
        terms = []
        for c, v in sorted(r.items()):
            a, b = pairs[c - nv]; vv = v if v < MOD // 2 else v - MOD
            terms.append(f"{vv}·t{a}t{b}")
        return " + ".join(terms)
    for r in Q[:12]: print("   ", fmt(r))
