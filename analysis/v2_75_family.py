"""v2.75: the one-parameter families on a class of ≥ 3 arrows a0, a1, a2: i→j with two tails b, b' out of j.

P = the projection onto ⟨paths⟩ ∖ {a0, a0·b, a0·b'} along ⟨a0 − a1, k_b, k_b'⟩,
    k_x = a0·x − a1·x − λ_x c − m_x ℓ₀,
so P(a0) = a1, P(a0·x) = a1·x + λ_x c + m_x ℓ₀ (x = b, b'), and P fixes every other path.
  c   'src'   e_i;  'others'  Σ_{v ≠ i} e_v;  'local'  e_j + e_t(b) + e_t(b')
  ℓ₀  'sigma0'  a1 − a2 (in ker σ);  'a1'  a1 (not in ker σ)
Checks the fold equations (fold_ideal._equations: P² = P, Δ-compatibility, image and kernel closed under products)
identically in ℚ[λ, m, λ', m'] and, if some fail, prints a Gröbner basis of the ideal they generate.
  c 'tail': e_t(x) for the tail x.
Usage: [FAMILY_MODE=naive|mult] python3 analysis/v2_75_family.py QUIVER A0,A1,A2 B,B' C ELL [OUT.json]
(arrow indices of pcensus.QUIVERS)"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE)); sys.path.insert(0, HERE)
import numpy as np, sympy as sp
from sympy import QQ
from sympy.polys.rings import ring
import pcensus as pc, fold_ideal as fi

def build(quiver, cls, tails, cmode, ell, mode='naive'):
    """mode 'naive': only a0, a0·b, a0·b' move (b, b' = tails).
    mode 'mult': every path through a0 moves: P(u·a0·r) = u·a1·r if r = e_j; for r = x·q (x any arrow out of j, with
    its own parameters λ_x, m_x): u·(a1·x + λ_x c_x + m_x ℓ₀)·q, where an empty u or q is not multiplied in.
    In 'mult' mode `tails` is ignored: every arrow out of j gets parameters."""
    ctx = pc.setup(quiver, 7); d = ctx['d']; basis = ctx['basis']; idx = {b: i for i, b in enumerate(basis)}
    arrows = ctx['arrows']; a0, a1, a2 = cls
    i, j = arrows[a0]
    assert arrows[a1] == arrows[a2] == (i, j)
    outs = [x for x, (s_, t_) in enumerate(arrows) if s_ == j] if mode == 'mult' else list(tails)
    assert all(arrows[x][0] == j for x in outs)
    R, *gens = ring(','.join(f'l{x},m{x}' for x in outs), QQ)
    par = {x: (gens[2 * k], gens[2 * k + 1]) for k, x in enumerate(outs)}
    nv = sum(1 for x in basis if x[0] == 'v')
    M = np.asarray(ctx['M'], dtype=np.int64)
    nzM = [(a, b, k, int(M[a, b, k])) for a in range(d) for b in range(d) for k in range(d) if M[a, b, k]]
    def mul(x, y):
        out = [R(0)] * d
        for a, b, k, c in nzM:
            if x[a] and y[b]: out[k] += c * x[a] * y[b]
        return out
    unit = lambda k: [R(1) if r == k else R(0) for r in range(d)]
    def cvec(x):
        if cmode == 'src': vs_ = [i]
        elif cmode == 'tail': vs_ = [arrows[x][1]]
        elif cmode == 'others': vs_ = [v for v in range(nv) if v != i]
        elif cmode == 'local': vs_ = [j] + [arrows[y][1] for y in outs]
        else: raise ValueError(cmode)
        return {idx[('v', v)]: 1 for v in vs_}
    lvec = {idx[('p', (a1,))]: 1, idx[('p', (a2,))]: -1} if ell == 'sigma0' else {idx[('p', (a1,))]: 1}
    def E(x):                                           # P(a0·x)
        l, m = par[x]; v = unit(idx[('p', (a1, x))])
        for r, c in cvec(x).items(): v[r] += l * c
        for r, c in lvec.items(): v[r] += m * c
        return v
    P = [[R(1) if r == c else R(0) for c in range(d)] for r in range(d)]
    for col, bb in enumerate(basis):
        if bb[0] != 'p' or a0 not in bb[1]: continue
        k = bb[1].index(a0); u, r = bb[1][:k], bb[1][k + 1:]
        if mode == 'naive' and (u or (r and (len(r) > 1 or r[0] not in outs))): continue
        if not r: v = unit(idx[('p', u + (a1,))])
        else:
            x, q = r[0], r[1:]
            v = E(x)
            if u: v = mul(unit(idx[('p', u)]), v)
            if q: v = mul(v, unit(idx[('p', q)]))
        for rr in range(d): P[rr][col] = v[rr]
    return ctx, R, P

if __name__ == '__main__':
    quiver = sys.argv[1]; cls = tuple(map(int, sys.argv[2].split(','))); tails = tuple(map(int, sys.argv[3].split(',')))
    cmode, ell = sys.argv[4], sys.argv[5]; out = sys.argv[6] if len(sys.argv) > 6 else None
    mode = os.environ.get('FAMILY_MODE', 'naive')
    ctx, R, P = build(quiver, cls, tails, cmode, ell, mode)
    d = ctx['d']; M = np.asarray(ctx['M'], dtype=np.int64); DL = np.asarray(ctx['DL'], dtype=np.int64)
    nzM = [(a, b, k, int(M[a, b, k])) for a in range(d) for b in range(d) for k in range(d) if M[a, b, k]]
    res = dict(quiver=quiver, mode=mode, cls=cls, tails=tails, c=cmode, ell=ell, d=d, params=[str(g) for g in R.gens])
    for part in ('idem', 'delta', 'closure'):
        eqs = fi._equations(P, R, d, nzM, DL, (part,))
        res[part] = len(eqs)
    eqs = fi._equations(P, R, d, nzM, DL)
    if eqs:
        G = sp.groebner([fi.to_expr(e) for e in eqs], *[sp.Symbol(str(g)) for g in R.gens], order='grevlex')
        res['groebner'] = [str(g) for g in G.exprs]
    print(json.dumps(res))
    if out: json.dump(res, open(out, 'w'), indent=1)
