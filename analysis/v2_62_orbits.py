"""v2.62: explicit representatives of the census orbits of Q₀ in characteristic 0, from the v2.58 certificates.

For a component Z of a stratum and its group H_Z (Lie algebra 𝔥_Z = Der_Δ ⊕ E_{W_Z}), every H_Z-orbit O ⊆ Z is
- the open orbit (dim O = dim Z), or
- of dimension r < dim Z, and then O ⊆ Z_r ⊆ V(J_{Z,r}) (gens + the certified (r+1)-minors), which has dimension ≤ r, so
  the closure of O is an r-dimensional irreducible component C of V(J_{Z,r}), and O = {P ∈ C : rank [𝔥_Z, P] = r}
  (rank = orbit dimension in characteristic 0; such points form one open dense orbit of C).
So: for each Z, a rational point of Z of rank dim Z, and for each r and each r-dimensional component C of V(J_{Z,r}), a
rational point of C of rank r, if one is found. Components C with no point of rank r (all sampled points of lower rank)
are not orbit closures. Components without a rational point found are reported (their orbit may still exist; a point
over a number field would be needed).
Output: JSON lines (block, comp, r, point {coordinate: value}, rank) and a summary.
Usage: python3 analysis/v2_62_orbits.py OUT.jsonl [--tries N]"""
import os, sys, json, subprocess, collections, argparse, random, re
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE)); sys.path.insert(0, HERE)
import sympy as sp
import fold_ideal as fi, singular_nohilb
singular_nohilb.install()
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--tries', type=int, default=60); ap.add_argument('--only')
A = ap.parse_args(); rng = random.Random(62)
root = os.path.dirname(HERE)
recs = [json.loads(l) for l in open(os.path.join(root, 'logs', 'v2_58', 'certs_q0.jsonl'))]
groups = collections.OrderedDict()
for r in recs:
    g = groups.setdefault(r['block'], {'stratum': None, 'comps': []})
    if r['kind'] == 'stratum': g['stratum'] = r
    else: g['comps'].append(r)
out = open(A.out, 'w'); summ = collections.Counter()
def singular_components(vs, gens, r):
    src = [f'LIB "primdec.lib"; ring R = 0,({",".join(vs)}),dp;', f"ideal J = {','.join(g.replace('**', '^') for g in gens) or '0'};",
           "list L = minAssGTZ(J); int i;",
           f'for (i = 1; i <= size(L); i++) {{ ideal G = std(L[i]); if (dim(G) == {r}) {{ print("COMP " + string(G)); intvec iv = indepSet(G); print("INDEP " + string(iv)); kill iv; }} kill G; }}',
           'print("END"); quit;']
    fn = A.out + '.sing'; open(fn, 'w').write("\n".join(src))
    res = subprocess.run(["Singular", "-q", "--no-rc", fn], capture_output=True, text=True, timeout=3000).stdout
    assert "END" in res and "? " not in res, res[-300:]
    comps = re.findall(r"COMP (.*)\nINDEP (.*)", res)
    return [(c.split(','), [int(x) for x in iv.split(',')]) for c, iv in comps]
def points(Sv, gl, iv, k):
    pts = []
    for t in range(A.tries):
        vals = {Sv[j]: sp.Integer(rng.randint(-4, 4)) for j in range(len(Sv)) if iv[j] == 1}
        g2 = [g for g in (sp.expand(g.subs(vals)) for g in gl) if g != 0]
        unk = sorted(set().union(*[g.free_symbols for g in g2]), key=str) if g2 else []
        try: sols = sp.solve(g2, unk, dict=True) if g2 else [{}]
        except Exception: sols = []
        for so in sols:
            if len(so) < len(unk) or not all(x.is_rational for x in so.values()): continue
            pt = dict(vals); pt.update(so)
            for v in Sv: pt.setdefault(v, sp.Integer(rng.randint(-4, 4)))
            pts.append(pt); break
        if len(pts) >= k: break
    return pts
for block, g in groups.items():
    if A.only and block != A.only: continue
    st = g['stratum']; vs = st['vs']
    if not vs:
        for c in g['comps']:
            out.write(json.dumps(dict(block=block, comp=c['comp'], r=0, point={}, rank=0)) + "\n"); summ['point strata'] += 1
        continue
    V = [[-1 if ch == 'm' else int(ch) for ch in row] for row in block.split(',')]
    S = fi.stratum('pairs_in_short', V); R = S['R']; d = S['ctx']['d']
    gl_ = {str(x): x for x in R.gens}
    parse = lambda s_: R(sp.sympify(s_, locals={k: sp.Symbol(k) for k in gl_}))
    Pr = [row[:] for row in S['P']]
    for x, v in st['subs']:
        xv, vv = gl_[x], parse(v); Pr = [[e.compose(xv, vv) for e in row] for row in Pr]
    Sv = sp.symbols(vs); loc = {str(v): v for v in Sv}
    PrE = sp.Matrix([[e.as_expr() for e in row] for row in Pr])
    DerD = [sp.Matrix(d, d, [sp.Rational(x) for x in D]) for D in st['DerD']]
    Ec = [sp.diag(*e) for e in st['Ec']]
    for c in g['comps']:
        if c['dim'] <= 0:
            pt = {}
            gl = [sp.sympify(x, locals=loc) for x in c['gens']]
            so = sp.solve(gl, Sv, dict=True)
            if so and all(v.is_rational for v in so[0].values()):
                out.write(json.dumps(dict(block=block, comp=c['comp'], r=0, point={str(k): str(v) for k, v in so[0].items()}, rank=0)) + "\n"); summ['orbits (points)'] += 1
            else: summ['point components without a rational point'] += 1
            continue
        W = [[sp.Rational(x) for x in w] for w in c['W']]
        h = DerD + [sum((w[k] * Ec[k] for k in range(len(Ec))), sp.zeros(d, d)) for w in W]
        Mfun = sp.lambdify(Sv, sp.Matrix([list(D * PrE - PrE * D) for D in h]), 'sympy')
        rank_at = lambda pt: sp.Matrix(Mfun(*[pt[v] for v in Sv])).rank()
        targets = [(c['dim'], [(c['gens'], c['indep'])])]
        for cut in c['cuts']:
            comps = singular_components(vs, c['gens'] + [m['det'] for m in cut['minors']], cut['r'])
            targets.append((cut['r'], comps))
        for r, comps in targets:
            for gens, iv in comps:
                gl = [sp.sympify(x, locals=loc) for x in gens if x.strip() not in ('0', '')]
                pts = points(Sv, gl, iv, 6)
                if not pts: summ['components without a rational point found'] += 1; print(f"{block} comp {c['comp']} r {r}: no rational point", flush=True); continue
                hit = [pt for pt in pts if rank_at(pt) == r]
                if hit:
                    out.write(json.dumps(dict(block=block, comp=c['comp'], r=r, point={str(k): str(v) for k, v in hit[0].items()}, rank=r)) + "\n"); out.flush()
                    summ['orbits'] += 1
                else: summ['components that are not orbit closures (no point of rank r)'] += 1
    print(f"{block}: {dict(summ)}", flush=True)
print("SUMMARY", dict(summ))
