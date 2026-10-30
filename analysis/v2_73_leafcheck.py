"""v2.73: per-leaf ℤ-checks for a case-split block (v2.57 case split, leaf k of CERT).
  B  the leaf's points (rem = 0, ∏ nz ≠ 0, in the chart vs) lie in its leaf prime mod p: for each integral generator ρ of
     the leaf prime, the ideal (rem) + (1 − T·∏nz·ρ) over ℤ contains an integer N_ρ; candidates: primes of N_ρ.
     (An empty leaf, ncomp = 0: 1 ∈ (rem) + (1 − T·∏nz), integer N.)
  C  if the leaf prime is not one of the H-stable roots ROOTS (full-coordinate keys): a root R containing it over ℚ, and
     over ℤ: every generator of R's key, pulled back to the chart (each full coordinate is coord(P_leaf, name), a
     polynomial in vs), reduces to 0 modulo a strong ℤ-basis of the leaf prime; candidates: its leading coefficients
     (and the denominators of the chart's coordinate functions).
  P  a point leaf (vs = []) fixed by H in char 0: G fixes it and every compatible weighting move fixes it, as identities
     over ℤ (in g, and in λ over a ℤ-basis of the compatible weightings); candidates: Smith primes of the defect matrix.
Usage: python3 analysis/v2_73_leafcheck.py CERT BLOCK LEAF ROOTS.json OUT"""
import os, sys, json, re, math, subprocess
HERE = os.path.dirname(os.path.abspath(__file__))
cert, block, leaf, rootsf, out = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
os.environ['STAB_LEAF'] = leaf; os.environ['NO_POINT_EXIT'] = '1'
sys.argv = ['v2_71_gcert.py', cert, block, rootsf, out]
src = open(os.path.join(HERE, 'v2_71_gcert.py')).read(); src = src[:src.index("setup_primes = set()")]
G2 = {'__name__': 'v2_73l', '__file__': os.path.join(HERE, 'v2_71_gcert.py')}
exec(compile(src, 'v2_71_gcert', 'exec'), G2)
G_ = G2['G_']; sp = G2['sp']; R = G2['R']; vs = G2['vs']; st = G2['st']; d = G2['d']; Pr = G2['Pr']; coord = G2['coord']
primes_of = G2['primes_of']; smith_primes = G2['smith_primes']; sigma = G2['sigma']; Gm = G2['Gm']; cl = G2['cl']
defects = G2['defects']; nzd = G2['nzd']; NCL = G2['NCL']; S = G_['S']
ROOTS = json.load(open(rootsf))['roots']                    # full-coordinate keys of the H-stable roots
names = S['names']
comp = None
for l in open(cert):
    if f'"block": "{block}"' not in l: continue
    r = json.loads(l)
    if r['kind'] == 'component' and r['block'] == block and r['leaf'] == int(leaf): comp = r
ssy = lambda e: str(e).replace('**', '^')
def ints(gl):
    out_ = []
    for g in gl:
        e = sp.sympify(str(g).replace('^', '**'))
        if e == 0: continue
        num, _ = sp.fraction(sp.together(sp.expand(e))); num = sp.expand(num)
        c = sp.Poly(num, *[sp.Symbol(v) for v in vs]).content() if vs else abs(num)
        out_.append(sp.expand(num / c) if c not in (0, 1) else num)
    return out_
def run(lines, tag):
    fn = f"{out}_{tag}.sing"; open(fn, 'w').write("\n".join(lines))
    try: o = subprocess.run(["Singular", "-q", "--no-rc", fn], capture_output=True, text=True, timeout=3000)
    except subprocess.TimeoutExpired: return None
    o = o.stdout + o.stderr
    return None if ("? error" in o or "END" not in o) else o
res = dict(leaf=int(leaf), vs=len(vs), ncomp=st['ncomp'], branch_consts=st.get('branch_consts'), cand={}, checks={})
def addc(tag, ps):
    for q in ps: res['cand'].setdefault(str(q), []).append(tag)
for a, b in st.get('branch_consts') or []: addc('A', primes_of(a) | primes_of(b))
# the chart's coordinate functions: full coordinate n = coord(P_leaf(vs), n)
CF = {n: sp.expand(coord(sp.Matrix(d, d, lambda a, b: Pr[a][b].as_expr() if Pr[a][b] else 0), n)) for n in names}
for e in CF.values():
    for c in sp.Poly(e, *[sp.Symbol(v) for v in vs]).coeffs() if vs else [e]:
        addc('chart', primes_of(sp.fraction(sp.Rational(c))[1]))
nzp = sp.Integer(1)
for f in st['nz']: nzp *= sp.sympify(f.replace('^', '**'))
rem = ints(st['rem'])
P_L = ints(comp['gens']) if comp and comp.get('gens') else []
RING = f"ring RZ = integer,({','.join(vs + ['TT'])}),dp;" if vs else "ring RZ = integer,(TT),dp;"
# B: leaf points lie in the leaf prime
if st['ncomp'] == 0 or comp is None: targets = [sp.Integer(1)]
else: targets = P_L
lines = [RING, "option(redSB);", f"ideal REM = {','.join(ssy(e) for e in rem) or '0'};"]
for j, rho in enumerate(targets):
    lines.append(f"ideal S{j} = std(REM + ideal(1 - TT*({ssy(sp.expand(nzp))})*({ssy(rho)}))); if (deg(S{j}[1]) == 0) {{ print(\"RAB {j} \" + string(S{j}[1])); }} else {{ print(\"RABFAIL {j}\"); }}")
lines.append('print("END"); quit;')
o = run(lines, 'B') if targets else 'END'
if o is None: res['checks']['B'] = 'singular error/timeout'
else:
    rab = [int(x) for x in re.findall(r"^RAB \d+ (-?\d+)", o, re.M)]; fail = o.count("RABFAIL")
    res['checks']['B'] = 'ok' if not fail else f'{fail} generators not in the radical'
    for N in rab: addc('B', primes_of(N))
# C: a leaf prime that is not an H-stable root lies in one
full = lambda gl: [sp.sympify(g.replace('^', '**')) for g in gl.split(',') if g.strip() not in ('', '0')]
def pull(key):
    sub = {sp.Symbol(n): CF[n] for n in names}
    return [sp.expand(g.subs(sub, simultaneous=True)) for g in full(key)]
mykey = None
if comp is not None and st['ncomp'] and vs:
    # this leaf prime's own full key: its generators plus the graph of the chart
    pass
res['root_of_leaf'] = None
if st['ncomp'] and comp is not None:
    for ri, key in enumerate(ROOTS):
        pb = [g for g in pull(key) if g != 0]
        if not pb: res['root_of_leaf'] = ri; res['checks']['C'] = 'ok (identically)'; break
        lines = [RING.replace('integer', '0'), "option(redSB);", f"ideal PL = std(ideal({','.join(ssy(e) for e in P_L) or '0'}));",
                 f"ideal PB = {','.join(ssy(e) for e in pb)};", 'print("INQ " + string(size(reduce(PB, PL)) == 0));', 'print("END"); quit;']
        o = run(lines, f'Cq{ri}')
        if o is None or re.search(r"INQ (\d)", o).group(1) != '1': continue
        lines = [RING, "option(redSB);", G2['LCPROC'], f"ideal PL = std(ideal({','.join(ssy(e) for e in P_L) or '0'}));", 'print("LC" + lcs(PL));',
                 f"ideal PB = {','.join(ssy(sp.fraction(sp.together(e))[0]) for e in pb)};", 'print("INZ " + string(size(reduce(PB, PL)) == 0));', 'print("END"); quit;']
        o = run(lines, f'Cz{ri}')
        if o is None: continue
        for c in re.search(r"LC(.*)", o).group(1).split(): addc('C', primes_of(c))
        res['root_of_leaf'] = ri; res['checks']['C'] = 'ok' if re.search(r"INZ (\d)", o).group(1) == '1' else 'not exact over ℤ'
        break
    if res['root_of_leaf'] is None: res['checks']['C'] = 'no containing root found'
# P: a point leaf: fixed by G and by the compatible weighting moves, identically over ℤ
if not vs and st['ncomp']:
    Px = sp.Matrix(d, d, lambda a, b: sp.Rational(str(Pr[a][b])) if Pr[a][b] else 0)
    for x in Px: addc('Pden', primes_of(sp.fraction(x)[1]))
    Sg = sigma(Gm, lambda c: sp.Integer(1))
    res['checks']['P_G_fixes'] = all(sp.expand(x) == 0 for x in (Sg * Px - Px * Sg))
    Dx = sp.Matrix([[sp.Rational(str(defects[c][k])) if defects[c][k] else 0 for c in range(NCL)] for k in nzd]) if nzd else sp.zeros(0, NCL)
    wp, wr = smith_primes(Dx); addc('PW', wp)
    ker = Dx.nullspace() if Dx.rows else [sp.Matrix([1 if i == j else 0 for i in range(NCL)]) for j in range(NCL)]
    lam = sp.symbols(f'lam0:{len(ker)}') if ker else ()
    s_ = [sum((lam[k] * ker[k][c] for k in range(len(ker))), sp.Integer(0)) for c in range(NCL)]
    Sw = sigma([sp.eye(len(arr)) for arr in cl], lambda c: s_[c])
    res['checks']['P_W_fixes'] = all(sp.expand(x) == 0 for x in (Sw * Px - Px * Sw))
    res['point_full'] = {n: str(CF[n]) for n in names}
json.dump(res, open(out + '.json', 'w'), indent=1)
print(f"LEAF {leaf}: vs {len(vs)}, ncomp {st['ncomp']}, checks {res['checks']}, root {res['root_of_leaf']}, candidates {res['cand']}")
