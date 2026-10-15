import sys, time, sympy as sp
from fold import *
from tangent import spec, gvec
from second import quadrics
def ratrec(a, m=MOD):
    a %= m; r0, r1, s0, s1 = m, a, 0, 1
    while r1 * r1 * 2 > m: q = r0 // r1; r0, r1, s0, s1 = r1, r0 - q * r1, s1, s0 - q * s1
    return sp.Rational(r1, s1)
def wtype(S, w):
    x, y = w; return S.side[x] + S.side[y]       # 'HH','HQ' (left edge: Q-run right of a),'QH' (right edge),'QQ'
def build(kind, Lw, delta):
    S = spec(kind, Lw + 2 * delta + 4); F = Fold(S, Lw, 2 * delta)
    eqs = equations(F, Lw); J = linear_part(eqs)
    big = [j for j, k in enumerate(F.names) if S.lev[k[2]] > F.wlev(k[1]) + delta]
    return S, F, eqs, J, big
def sub_tangent(S, F, J, big, kill_types):
    kill = [j for j, k in enumerate(F.names) if k[0] == 'c' and wtype(S, k[1]) in kill_types]
    K, _ = kernel(J + [{j: 1} for j in big + kill], len(F.names)); return K
def modG(S, F, vecs, G, Llow):
    low = [j for j, k in enumerate(F.names) if k[0] == 'c' and F.wlev(k[1]) <= Llow]
    E = Echelon(); [E.add({j: g[j] for j in low if j in g}) for g in G]
    reps = [v for v in vecs if E.add({j: v[j] for j in low if j in v})]
    return reps
def qpolys(F, eqs, taus, names):
    Q, pairs, nv = quadrics(F, eqs, taus)
    t = sp.symbols(names)
    return [sp.expand(sum(ratrec(v) * t[pairs[c - nv][0]] * t[pairs[c - nv][1]] for c, v in r.items())) for r in Q]
if __name__ == '__main__':
    kind = sys.argv[1]; Lw, delta, Llow = map(int, sys.argv[2:5]); t0 = time.time()
    S, F, eqs, J, big = build(kind, Lw, delta)
    G = [gvec(F, h) for h in range(len(S.B)) if S.lev[h] <= delta]
    full = sub_tangent(S, F, J, big, set())
    TL = sub_tangent(S, F, J, big, {'QH'})    # vanish on right-edge words
    TR = sub_tangent(S, F, J, big, {'HQ'})    # vanish on left-edge words
    TH = sub_tangent(S, F, J, big, {'QH', 'HQ', 'QQ'})  # only Hodge-side corrections
    rF, rL, rR, rH = (len(modG(S, F, X, G, Llow)) for X in (full, TL, TR, TH))
    print(f"{kind} Lw={Lw} δ={delta}: extras mod G: all={rF}  left-only={rL}  right-only={rR}  hodge-only={rH}   [{time.time()-t0:.0f}s]")
    for nm, X in (('left', TL), ('right', TR), ('hodge', TH)):
        reps = modG(S, F, X, G, Llow)
        if not reps: continue
        P = qpolys(F, eqs, reps, [f'{nm[0]}{i}' for i in range(len(reps))])
        print(f"   {nm}-only family: {len(reps)} directions, obstruction quadrics restricted to it: {len(P)}", P[:3])
    repsL = modG(S, F, TL, G, Llow); repsR = modG(S, F, TR, G, Llow)
    P = qpolys(F, eqs, repsL + repsR, [f'l{i}' for i in range(len(repsL))] + [f'r{i}' for i in range(len(repsR))])
    print(f"   left ⊕ right: {len(P)} quadrics; all of mixed type l·r: {all(all(any(str(s)[0]=='l' for s in m.free_symbols) and any(str(s)[0]=='r' for s in m.free_symbols) for m in sp.Add.make_args(p)) for p in P)}")
    for p in P[:6]: print("     ", p)
