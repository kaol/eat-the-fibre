"""v2.75: rigidity of the one-parameter normal forms on 1⇛2, 2→3, 2→4 under the census symmetries G × T³, as an ideal
membership: the entries of h·P_t − P_t'·h (h = the substitution by g ∈ G (3×3, column sums 1) scaled by the class
weight wa, the tails scaled by wb, wc) generate an ideal containing u·(t − t') for an explicit unit-type u.
P_t (c = 'src'): λ = m = λ' = 1, m' = t; P⁻_t (c = 'others'): the same with vertex part e2 + e3 + e4.
Also P_t against P⁻_t': the ideal contains a unit-type element (never conjugate).
Writes Singular scripts and runs them (lift gives the certificate).
Usage: python3 analysis/v2_75_rigid.py OUTDIR"""
import os, sys, json, itertools, subprocess, re
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sympy as sp
from v2_75_family import build
outd = sys.argv[1]; os.makedirs(outd, exist_ok=True)
def fam(cmode, tv):
    ctx, R, Pr = build('K3_tgtfork', (0, 1, 2), (3, 4), cmode, 'sigma0', 'mult')
    sub = {sp.Symbol('l3'): 1, sp.Symbol('m3'): 1, sp.Symbol('l4'): 1, sp.Symbol('m4'): tv}
    return ctx, sp.Matrix(ctx['d'], ctx['d'], lambda r, c: Pr[r][c].as_expr().subs(sub))
t, tp = sp.symbols('t tp')
ctx, Pt = fam('src', t); _, Ptp = fam('src', tp); _, Mt = fam('others', t); _, Mtp = fam('others', tp)
d = ctx['d']; basis = ctx['basis']; idx = {b: k for k, b in enumerate(basis)}
g = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'g{i}{j}') if i < 2 else 0)
for j in range(3): g[2, j] = 1 - g[0, j] - g[1, j]
wa, wb, wc = sp.symbols('wa wb wc')
cl = [[0, 1, 2], [3], [4]]; Gs = [g, sp.Matrix([[1]]), sp.Matrix([[1]])]; sc = [wa, wb, wc]
cls_of = {a: k for k, arr in enumerate(cl) for a in arr}
h = sp.zeros(d, d)
for col, b in enumerate(basis):
    if b[0] == 'v': h[col, col] = 1; continue
    choices = [[(cl[cls_of[a]][i], sc[cls_of[a]] * Gs[cls_of[a]][i, cl[cls_of[a]].index(a)]) for i in range(len(cl[cls_of[a]]))] for a in b[1]]
    for combo in itertools.product(*choices):
        cf = sp.Integer(1)
        for _, x in combo: cf *= x
        h[idx[('p', tuple(a for a, _ in combo))], col] += cf
json.dump(dict(h=[[str(sp.expand(h[r, c])) for c in range(d)] for r in range(d)]), open(f'{outd}/subst.json', 'w'))
vars_ = 'g00,g01,g02,g10,g11,g12,wa,wb,wc,t,tp'
res = {}
for tag, A, B in (('src', Pt, Ptp), ('others', Mt, Mtp), ('src_vs_others', Pt, Mtp)):
    E = [sp.expand(e) for e in (h * A - B * h) if sp.expand(e) != 0]
    lines = [f"ring r = 0,({vars_}),dp;", "option(redSB);", f"ideal E = {','.join(str(e).replace('**', '^') for e in E)};",
             "ideal S = std(E);", 'print("STD"); print(S);',
             "ideal Q = sat(E, ideal(wa*wb*wc))[1];", "Q = sat(Q, ideal(g00*g11 - g01*g10 + g00*g12 - g02*g10 + g01*g12 - g02*g11 + 0))[1];",
             'print("SAT"); print(std(Q));', 'print("END"); quit;']
    fn = f'{outd}/rigid_{tag}.sing'; open(fn, 'w').write("\n".join(lines))
    o = subprocess.run(['Singular', '-q', '--no-rc', fn], capture_output=True, text=True, timeout=3000).stdout
    res[tag] = dict(n_entries=len(E), std=o.split('STD')[1].split('SAT')[0].strip()[:2000], sat=o.split('SAT')[1].split('END')[0].strip()[:2000])
    print(tag, len(E), 'entries;\n  std:', res[tag]['std'][:600].replace('\n', ' '), '\n  saturated:', res[tag]['sat'][:600].replace('\n', ' '))
json.dump(res, open(f'{outd}/rigid.json', 'w'), indent=1)
