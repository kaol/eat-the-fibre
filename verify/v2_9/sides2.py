import sys, time, sympy as sp
import fold as FD
from fold import *
from tangent import spec, gvec
from second import quadrics
from sides import ratrec, wtype, modG, qpolys
def build(kind, Leq, Lu, delta, D):
    S = spec(kind, Lu + 2 * D + 4); F = Fold(S, Lu, D)
    eqs = equations(F, Leq); assert FD.SKIPS[0] == 0, FD.SKIPS[0]
    J = linear_part(eqs)
    big = [j for j, k in enumerate(F.names) if S.lev[k[2]] > F.wlev(k[1]) + delta]
    return S, F, eqs, J, big
def sub(S, F, J, big, kill_types):
    kill = [j for j, k in enumerate(F.names) if k[0] == 'c' and wtype(S, k[1]) in kill_types]
    return kernel(J + [{j: 1} for j in big + kill], len(F.names))[0]
if __name__ == '__main__':
    kind = sys.argv[1]; Leq, Lu, delta, D, Llow = map(int, sys.argv[2:7]); t0 = time.time()
    S, F, eqs, J, big = build(kind, Leq, Lu, delta, D)
    G = [gvec(F, h) for h in range(len(S.B)) if S.lev[h] <= delta]
    fam = {nm: modG(S, F, sub(S, F, J, big, kt), G, Llow) for nm, kt in
           (('all', set()), ('left', {'QH'}), ('right', {'HQ'}), ('hodge', {'QH', 'HQ', 'QQ'}), ('edges', {'HH', 'QQ'}))}
    print(f"{kind} Leq={Leq} Lu={Lu} δ={delta} D={D}: extras mod G: " + ", ".join(f"{k}={len(v)}" for k, v in fam.items()) + f"  [{time.time()-t0:.0f}s]")
    for nm in ('left', 'right', 'hodge'):
        if fam[nm]:
            P = qpolys(F, eqs, fam[nm], [f'{nm[0]}{i}' for i in range(len(fam[nm]))])
            print(f"   {nm}-only: {len(fam[nm])} dirs, {len(P)} quadrics", P[:4])
    L, R = fam['left'], fam['right']
    P = qpolys(F, eqs, L + R, [f'l{i}' for i in range(len(L))] + [f'r{i}' for i in range(len(R))])
    print(f"   left ⊕ right: {len(P)} quadrics")
    for p in P[:8]: print("     ", sp.factor(p))
    Pall = qpolys(F, eqs, G + fam['all'], [f'g{i}' for i in range(len(G))] + [f't{i}' for i in range(len(fam['all']))])
    print(f"   whole tangent (G + extras): {len(Pall)} quadrics; G-variables appear: {any(str(s)[0]=='g' for p in Pall for s in p.free_symbols)}")
