"""v2.61 (Program 3, item 1): the split λ = 1 μ-fibre of Q₀ → Q₁ over the generic point of each component of Q₀.

For each component Z of Folds(ℚQ₀) (certificates logs/v2_58/certs_q0.jsonl): solve the prime's generators for the
dependent coordinates in terms of an independent set t (Z must be a graph over it; otherwise reported and skipped), so
P̄ = P̄(t) has entries in ℚ(t). Build P(μ) = P̄(t) ⊕ P_J(λ = 1, μ) over ℚ(t)[m₀, …, m₆] (normal form (*), v2.60) and let
Singular (ring with parameters t) generate the fold equations of P(μ) and the conditions φ(B) ⊆ B, and compute the
reduced Gröbner basis over ℚ(t). Reported per component: empty / affine of dimension k (all basis elements of degree
≤ 1) / non-affine, and whether μ = 0 (the identity lift) is a solution.
This is the generic fibre only: special points of Z may have other fibres.
Usage: python3 analysis/v2_61_mu_generic.py OUT_PREFIX [--only BLOCK]"""
import os, sys, json, subprocess, collections, argparse
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE)); sys.path.insert(0, HERE)
import numpy as np, sympy as sp
import fold_ideal as fi, pcensus as pc
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--only'); A = ap.parse_args()
root = os.path.dirname(HERE)
recs = [json.loads(l) for l in open(os.path.join(root, 'logs', 'v2_58', 'certs_q0.jsonl'))]
groups = collections.OrderedDict()
for r in recs:
    g = groups.setdefault(r['block'], {'stratum': None, 'comps': []})
    if r['kind'] == 'stratum': g['stratum'] = r
    else: g['comps'].append(r)
c1 = pc.setup('pairs_in_short_tail', 3); b1 = c1['basis']; d1 = len(b1)
M1 = np.asarray(c1['M']); DL1 = np.asarray(c1['DL'])
idx1 = {b: i for i, b in enumerate(b1)}; e4 = idx1[('v', 3)]
ss = lambda e: str(e).replace('**', '^')
def fold_code():
    """Singular code: ideal F = fold equations of the matrix PP (d1 × d1) in the current ring."""
    nzM = [(i + 1, j + 1, k + 1, int(M1[i, j, k])) for i in range(d1) for j in range(d1) for k in range(d1) if M1[i, j, k]]
    nzD = [(l + 1, a + 1, b + 1, int(DL1[l, a, b])) for l in range(d1) for a in range(d1) for b in range(d1) if DL1[l, a, b]]
    return "\n".join([
        f"intmat MI[{len(nzM)}][4] = " + ",".join(f"{a},{b},{k},{c}" for a, b, k, c in nzM) + ";",
        f"intmat DI[{len(nzD)}][4] = " + ",".join(f"{l},{a},{b},{c}" for l, a, b, c in nzD) + ";",
        f"matrix QQm[{d1}][{d1}]; int ci; for (ci = 1; ci <= {d1}; ci++) {{ QQm[ci,ci] = 1; }} QQm = QQm - PP;",
        "ideal F = ideal(PP * PP - PP); int ca; int cb; int ct; int cl;",
        f"for (ca = 1; ca <= {d1}; ca++) {{ for (cb = 1; cb <= {d1}; cb++) {{",
        f"  matrix u[{d1}][1]; matrix w[{d1}][1];",
        f"  for (ct = 1; ct <= {len(nzM)}; ct++) {{",
        "    u[MI[ct,3],1] = u[MI[ct,3],1] + MI[ct,4] * PP[MI[ct,1],ca] * PP[MI[ct,2],cb];",
        "    w[MI[ct,3],1] = w[MI[ct,3],1] + MI[ct,4] * QQm[MI[ct,1],ca] * QQm[MI[ct,2],cb]; }",
        "  F = F + ideal(PP * u - u) + ideal(PP * w); kill u; kill w; } }",
        # Δ-compatibility: Σ_m PP[m,l] DL[m] = PP DL[l] PPᵀ, with DL[l] as d×d matrices
        f"list DLm; for (cl = 1; cl <= {d1}; cl++) {{ matrix Dm[{d1}][{d1}]; DLm[cl] = Dm; kill Dm; }}",
        f"for (ct = 1; ct <= {len(nzD)}; ct++) {{ matrix Dm = DLm[DI[ct,1]]; Dm[DI[ct,2],DI[ct,3]] = Dm[DI[ct,2],DI[ct,3]] + DI[ct,4]; DLm[DI[ct,1]] = Dm; kill Dm; }}",
        f"for (cl = 1; cl <= {d1}; cl++) {{ matrix L[{d1}][{d1}]; int cm; for (cm = 1; cm <= {d1}; cm++) {{ L = L + PP[cm,cl] * DLm[cm]; }} kill cm;",
        "  matrix Dl = DLm[cl]; F = F + ideal(L - PP * Dl * transpose(PP)); kill L; kill Dl; }",
        "F = simplify(F, 2 + 4 + 8);"])
FOLD = fold_code()
out_rows = []; tally = collections.Counter()
for block, g in groups.items():
    if A.only and block != A.only: continue
    st = g['stratum']; vs = st['vs']
    V = [[-1 if ch == 'm' else int(ch) for ch in row] for row in block.split(',')]
    S = fi.stratum('pairs_in_short', V); R = S['R']; c0 = S['ctx']; b0 = c0['basis']; d0 = c0['d']
    DL0 = np.asarray(c0['DL'])
    gl = {str(x): x for x in R.gens}
    parse = lambda s_: R(sp.sympify(s_, locals={k: sp.Symbol(k) for k in gl}))
    Pr = [row[:] for row in S['P']]
    for x, v in st['subs']:
        xv, vv = gl[x], parse(v); Pr = [[e.compose(xv, vv) for e in row] for row in Pr]
    Sv = sp.symbols(vs) if vs else []
    PrE = sp.Matrix([[e.as_expr() for e in row] for row in Pr])
    idx0 = {b: i for i, b in enumerate(b0)}
    B = [b for b in b0 if (b[0] == 'v' and b[1] == 1) or (b[0] == 'p' and c0['arrows'][b[1][-1]][1] == 1)]
    B0 = [idx0[b] for b in B]; Bpos = {i: k for k, i in enumerate(B0)}; nb = len(B)
    Bl = [idx1[('p', (5,))] if b[0] == 'v' else idx1[('p', b[1] + (5,))] for b in B]
    Qx = [idx1[b] for b in b0]
    for c in g['comps']:
        free = [Sv[j] for j in range(len(vs)) if c['indep'] and c['indep'][j] == 1] if vs else []
        gens = [sp.sympify(x, locals={str(v): v for v in Sv}) for x in c['gens'] if x.strip() not in ('0', '')]
        dep = [v for v in Sv if v not in free]
        if gens:
            sol = sp.solve(gens, dep, dict=True)
            if len(sol) != 1: out_rows.append((block, c['comp'], c['dim'], f'not a graph ({len(sol)} solutions)')); tally['skipped'] += 1; continue
            sol = sol[0]
        else: sol = {}
        par = {v: sol.get(v, v) for v in Sv}
        Pt = PrE.subs(par).applyfunc(sp.cancel)
        tv = sorted(set().union(*[e.free_symbols for e in Pt]) | set(), key=str)
        mu = sp.symbols(f'm0:{nb}')
        # φ(p) = P̄p + Σ μ(suf) P̄(pre)
        phi = sp.Matrix(d0, nb, lambda r, k: Pt[r, B0[k]])
        for k, pi in enumerate(B0):
            for a, bb in zip(*np.nonzero(DL0[pi])):
                if bb in Bpos:
                    for r in range(d0): phi[r, k] += int(DL0[pi][a, bb]) * mu[Bpos[bb]] * Pt[r, a]
        lin = [phi[r, k] for r in range(d0) if r not in Bpos for k in range(nb) if phi[r, k] != 0]
        P = sp.zeros(d1, d1)
        for i, qi in enumerate(Qx):
            for j, qj in enumerate(Qx): P[qi, qj] = Pt[i, j]
        P[e4, e4] = 1
        for k in range(nb):
            for k2 in range(nb): P[Bl[k2], Bl[k]] = phi[B0[k2], k]
            P[e4, Bl[k]] = mu[k]
        pars = ','.join(map(str, tv)) or 'xdummy'
        src = [f"ring r = (0,{pars}),({','.join(map(str, mu))}),dp; option(redSB);",
               f"matrix PP[{d1}][{d1}] = " + ",".join(ss(sp.together(P[i, j])) for i in range(d1) for j in range(d1)) + ";",
               FOLD, f"ideal LIN = {','.join(ss(sp.together(e)) for e in lin) or '0'};",
               "ideal G = std(F + LIN); int i; int md = 0; for (i = 1; i <= size(G); i++) { if (deg(G[i]) > md) { md = deg(G[i]); } }",
               "ideal Z0 = subst(F + LIN, " + ", ".join(f"{m}, 0" for m in mu) + "); Z0 = simplify(Z0, 2);",
               'print("DIM " + string(dim(G)) + " MAXDEG " + string(md) + " ZEROLIFT " + string(size(Z0) == 0)); quit;']
        fn = f"{A.out}_{block.replace(',', '')}_{c['comp']}.sing"; open(fn, 'w').write("\n".join(src))
        try:
            rr = subprocess.run(["Singular", "-q", "--no-rc", fn], capture_output=True, text=True, timeout=1200); res = rr.stdout + rr.stderr
        except subprocess.TimeoutExpired: res = ''
        if "DIM" not in res or "? " in res: out_rows.append((block, c['comp'], c['dim'], 'Singular failed: ' + res[-200:])); tally['failed'] += 1; continue
        w = res.split("DIM ")[1].split(); dm, md, z = int(w[0]), int(w[2]), w[4] == '1'
        kind = 'empty' if dm == -1 else (f'affine dim {dm}' if md <= 1 else f'NON-AFFINE dim {dm} degree {md}')
        tally[(kind, z)] += 1
        out_rows.append((block, c['comp'], c['dim'], kind + (', identity lift' if z else '')))
        print(f"{block} component {c['comp']} (dim {c['dim']}): generic μ-fibre {kind}{', identity lift present' if z else ''}", flush=True)
print("TALLY", dict(tally))
