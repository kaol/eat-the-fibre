"""v2.71 (from v2_70_gcert.py; separation per point, target primes ≥ 3): per-node ℤ-certificates for Conjecture G (Lemma L of v2.69 at a prime p), for one block's orbit tree.

For every dense node Y (group H_Y = ∏_c G_c × T_{W(Y)}, W(Y) a partition span) and each of NPTS rational points x of
the dense orbit (the node's minor ≠ 0 at x), the candidate bad primes of:
  den  denominators of x and of P_x (so x_p, P_{x_p} are defined and the cleared equations reduce correctly);
  W    the defect matrix at x (rows: defect entries, columns: classes): its kernel is W(Y) over ℚ (checked); primes of
       its Smith invariants (rank mod p = rank over ℚ ⇒ the compatible weightings at x_p are W(Y) mod p);
  S    the stabiliser ideal in ℤ[g, zz, t, u]: affine-linear in g (checked: Stab ∩ G is open in an affine space, so
       connected, in every characteristic); a strong ℤ-basis for the block order (g, zz) > (t, u), whose
       leading coefficients' primes are candidates (then its (t, u)-part reduces to the ideal of the image torus,
       which is the split subtorus of a saturated lattice: connected); the Jacobian at 1 has rank
       #vars − (dim H − dim Y) over ℚ (checked), primes of its Smith invariants are candidates (tangent space of
       Stab_p at 1 has the char 0 dimension ⇒ Stab_p smooth of dimension dim H − dim Y);
  G1   the graph ideal of h ↦ h·x in ℤ[g, z_c, t, u, y]: a strong ℤ-basis for the block order h > y; its y-part
       equals I_Y over ℚ (checked); primes of its leading coefficients are candidates (then Y_p is the closure of
       H_p·x_p: geometrically irreducible, with x_p in its dense orbit).
A prime is a candidate for Y only if it is one at every point. Writes OUT.json.
Usage: python3 analysis/v2_70_gcert.py CERT BLOCK TREE.json OUT [NPTS]"""
import os, sys, json, re, time, itertools, math, subprocess, zlib
HERE = os.path.dirname(os.path.abspath(__file__))
cert, block, tree, out = sys.argv[1:5]; NPTS = int(sys.argv[5]) if len(sys.argv) > 5 else 3; t0 = time.time()
TMO = int(os.environ.get('GCERT_TMO', '1200')); MAXPTS = int(os.environ.get('GCERT_MAXPTS', '40')); SMALLP = {int(x) for x in os.environ.get('GCERT_EXEMPT', '2').split(',') if x}
sys.argv = ['v2_67_orbits.py', cert, block, os.environ.get('STAB_LEAF', '0'), 'all', out + '_o', '--tmo', str(TMO)] + (['--leafmode'] if os.environ.get('STAB_LEAF') else [])
src = open(os.path.join(HERE, 'v2_67_orbits.py')).read(); src = src[:src.index("stats = {'primes'")]
G_ = {'__name__': 'v2_70', '__file__': os.path.join(HERE, 'v2_67_orbits.py')}
exec(compile(src, 'v2_67_orbits', 'exec'), G_)
sp = G_['sp']; R = G_['R']; vs = G_['vs']; st = G_['st']; ctx = G_['ctx']; d = ctx['d']; basis = ctx['basis']; Pr = G_['Pr']
defects = G_['defects']; nzd = G_['nzd']; NCL = G_['NCL']; compl = G_['compl']
Kl = G_['Kl']; Kmat = G_['Kmat']
T = json.load(open(tree))
cl = st['classes']; arrow_class = {a: c for c, arr in enumerate(cl) for a in arr}
idx = {b: i for i, b in enumerate(basis)}
from sympy.matrices.normalforms import smith_normal_form
def primes_of(n):
    n = abs(int(n)); return set(sp.factorint(n)) if n > 1 else set()
def den_primes(q): return primes_of(sp.fraction(sp.Rational(q))[1])
def smith_primes(M):
    """primes p with rank_p(M) < rank_ℚ(M), M a rational matrix: clear row denominators, Smith form over ℤ"""
    if M.rows == 0 or M.cols == 0: return set(), 0
    rows = []
    for i in range(M.rows):
        den = 1
        for x in M.row(i): den = den * int(sp.fraction(x)[1]) // math.gcd(den, int(sp.fraction(x)[1]))
        rows.append([int(x * den) for x in M.row(i)])
    Z = sp.Matrix(rows); r = Z.rank()
    if r == 0: return set(), 0
    S = smith_normal_form(Z, domain=sp.ZZ); out = set()
    for i in range(r): out |= primes_of(S[i, i])
    return out, r
# group coordinates (as v2_69_stab.py)
gvars = []; Gm = []
for c, arr in enumerate(cl):
    m = len(arr); M = [[None] * m for _ in range(m)]
    for j in range(m):
        for i in range(m - 1):
            v = sp.Symbol(f"g{c}_{i}{j}"); gvars.append(v); M[i][j] = v
        M[m - 1][j] = 1 - sum(M[i][j] for i in range(m - 1))
    Gm.append(sp.Matrix(M))
dimG = len(gvars); assert dimG == len(st['DerD'])
big = [c for c in range(len(cl)) if len(cl[c]) > 1]
zc = {c: sp.Symbol(f"z{c}") for c in big}
Ginv = [(zc[c] * Gm[c].adjugate()).applyfunc(sp.expand) if c in zc else Gm[c] for c in range(len(cl))]
def sigma(Gs, scal):
    """d × d matrix of the automorphism a_j ↦ scal(c) Σ_i Gs[c][i, j] a_i on arrows of class c"""
    S = sp.zeros(d, d)
    for col, b in enumerate(basis):
        if b[0] == 'v': S[col, col] = 1; continue
        choices = []
        for a in b[1]:
            c = arrow_class[a]; arr = cl[c]; j = arr.index(a)
            choices.append([(arr[i], scal(c) * Gs[c][i, j]) for i in range(len(arr))])
        for combo in itertools.product(*choices):
            path = ('p', tuple(x for x, _ in combo)); coef = sp.Integer(1)
            for _, cf in combo: coef *= cf
            S[idx[path], col] += coef
    return S
def scal_of(Wint, ts):
    def s(c):
        e = sp.Integer(1)
        for k in range(len(Wint)): e *= ts[k] ** Wint[k][c]
        return e
    return s
def coord(C, v):
    """the stratum coordinate v of a d × d matrix C (as G_['xi'], over sympy)"""
    m = re.fullmatch(r'y(\d+)_(\d+)', v)
    if m: return C[int(m.group(1)), compl[int(m.group(2))]]
    m = re.fullmatch(r'k(\d)(\d)', v); i, j = int(m.group(1)), int(m.group(2))
    return sp.expand((Kl[i, :] * C * Kmat[:, j])[0, 0])
ssy = lambda e: str(e).replace('**', '^')
def run_sing(lines, tag):
    fn = f"{out}_{tag}.sing"; open(fn, 'w').write("\n".join(lines))
    try: r = subprocess.run(["Singular", "-q", "--no-rc", fn], capture_output=True, text=True, timeout=TMO)
    except subprocess.TimeoutExpired: return 'TIMEOUT'
    o = r.stdout + r.stderr
    return None if ("? error" in o or "END" not in o) else o
def clear(e):
    num, den = sp.fraction(sp.together(sp.expand(e))); return sp.expand(num)
LCPROC = 'proc lcs(ideal G) { string s = ""; int i; for (i = 1; i <= size(G); i++) { if (G[i] != 0) { s = s + " " + string(leadcoef(G[i])); } } return(s); }'
def g1_graph(keep, Px, Sg, at1, ts, us, Wint):
    """the graph ideal of the slice map (coordinates keep of (g, t), the rest at 1): Σ → Y, s ↦ s·x"""
    r_ = len(ts); hv0 = list(ts) + gvars
    fix = {v: at1[v] for v in hv0 if v not in keep}
    keep_t = [k_ for k_ in range(r_) if ts[k_] in keep]
    fix.update({us[k_]: 1 for k_ in range(r_) if k_ not in keep_t})
    zk = [c for c in big if any(v in keep for v in Gm[c].free_symbols)]
    fix.update({zc[c]: 1 for c in big if c not in zk})
    Si = sigma(Ginv, scal_of(Wint, us))
    Cx = (Sg.subs(fix) * Px * Si.subs(fix))
    yv = [sp.Symbol(v) for v in vs]
    gr = [clear(yv[j_] - coord(Cx, v)) for j_, v in enumerate(vs)]
    hv1 = list(keep) + [zc[c] for c in zk] + [us[k_] for k_ in keep_t]
    rel1 = [sp.expand((zc[c] * Gm[c].det()).subs(fix) - 1) for c in zk] + [ts[k_] * us[k_] - 1 for k_ in keep_t]
    return hv1, yv, [ssy(q) for q in gr + rel1]
def pivots_mod(Tm, q):
    """pivot columns of a rational matrix over 𝔽_q (denominators prime to q)"""
    M = [[int(sp.fraction(x)[0]) * pow(int(sp.fraction(x)[1]), -1, q) % q for x in Tm.row(i)] for i in range(Tm.rows)]
    piv = []; r = 0
    for c in range(Tm.cols):
        k = next((i for i in range(r, len(M)) if M[i][c] % q), None)
        if k is None: continue
        M[r], M[k] = M[k], M[r]; inv = pow(M[r][c], -1, q); M[r] = [x * inv % q for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]: f_ = M[i][c]; M[i] = [(a - f_ * b) % q for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
    return piv
setup_primes = set()
for row in Pr:
    for e in row:
        if e:
            for c in e.coeffs(): setup_primes |= primes_of(c.denominator)
for M_ in (Kl, Kmat):
    for x in M_: setup_primes |= den_primes(x)
res = []; tally = {'nodes': 0, 'certified_points': 0, 'checks_failed': 0, 'timeouts': 0}
# separation, per point: for every other dense node Y' of the same dimension, the primes dividing gcd_q q(x) over the
# integral generators q of I_Y' (x_p ∉ Y'_p unless p divides it). Added to the point's candidates, so a prime survives
# only if no single point certifies every claim at once.
IYR = {}
for nd_ in T['nodes'].values():
    if nd_.get('result') == 'dense orbit':
        IYR[nd_['path']] = (nd_['dim'], [G_['parse'](g) for g in nd_['gens'] if g.strip() not in ('', '0')])
def sep_primes(path, dim_, pt):
    out_ = set(); bad_ = []
    for p2, (d2, gl) in IYR.items():
        if p2 == path or d2 != dim_: continue
        g_ = 0
        for q in gl:
            val = G_['evalr'](q, pt); g_ = math.gcd(g_, abs(int(sp.fraction(val)[0])))
            if g_ == 1: break
        if g_ == 0: bad_.append(p2)
        else: out_ |= primes_of(g_)
    return out_, bad_
pnodes = []
def point_node(nd):
    """a 0-dimensional node: one rational point x, fixed by G (char 0, hence by a polynomial identity in g over
    ℤ_(p)). Census moves at x_p: the weightings s ∈ W_{x_p}(𝔽_p) with s_c ≠ 0, W_{x_p} = ker D(x_p). Candidates: den(x),
    Smith primes of D(x) (then W_{x_p} = W_x ⊗ 𝔽_p, W_x = ker D(x) ∩ ℤ^cl saturated); checked: G fixes x, and
    φ_{Σ λ_k w_k} commutes with P_x identically in λ over ℤ (so every census weighting move fixes x_p)."""
    Sv = [sp.Symbol(v) for v in vs]
    sol = sp.solve([sp.sympify(g.replace('^', '**')) for g in nd['gens']], Sv, dict=True)
    if len(sol) != 1 or any(v not in sol[0] for v in Sv): return dict(path=nd['path'], dim=0, candidates=None, checks=[('rational_point', False)])
    pt = [0] * len(R.gens)
    for v in Sv: pt[G_['gpos'][str(v)]] = sp.Rational(sol[0][v])
    cand = set(setup_primes); src = {}
    def addc(tag, ps):
        for q in ps: cand.add(q); src.setdefault(q, set()).add(tag)
    for x in pt: addc('den', den_primes(x))
    Px = sp.Matrix(d, d, lambda a_, b_: G_['evalr'](Pr[a_][b_], pt) if Pr[a_][b_] else 0)
    for x in Px: addc('den', den_primes(x))
    checks = []
    Sg = sigma(Gm, lambda c: sp.Integer(1))
    checks.append(('G_fixes', (Sg * Px - Px * Sg).applyfunc(sp.expand).is_zero_matrix))
    Dx = sp.Matrix([[G_['evalr'](defects[c][k], pt) if defects[c][k] else 0 for c in range(NCL)] for k in nzd]) if nzd else sp.zeros(0, NCL)
    wp, wr = smith_primes(Dx); addc('W', wp)
    ker = Dx.nullspace() if Dx.rows else [sp.Matrix([1 if i == j else 0 for i in range(NCL)]) for j in range(NCL)]
    lam = sp.symbols(f'lam0:{len(ker)}') if ker else ()
    s_ = [sum((lam[k] * ker[k][c] for k in range(len(ker))), sp.Integer(0)) for c in range(NCL)]
    Sw = sigma([sp.eye(len(arr)) for arr in cl], lambda c: s_[c])
    checks.append(('W_fixes', (Sw * Px - Px * Sw).applyfunc(sp.expand).is_zero_matrix))
    ok = all(v for _, v in checks)
    return dict(path=nd['path'], dim=0, W_x=len(ker), candidates=sorted(cand) if ok else None, src={str(q): sorted(v) for q, v in src.items()},
                checks=checks, _x=[pt[G_['gpos'][v]] for v in vs])
for key, nd in T['nodes'].items():
    if nd.get('result') == 'points':
        pn = point_node(nd); pnodes.append(pn)
        print(f"  {nd['path']}: point, W_x {pn.get('W_x')}: candidates {pn.get('src')}{'' if pn['candidates'] is not None else ' CHECKS ' + str(pn['checks'])}", flush=True)
        continue
    if nd.get('result') != 'dense orbit': continue
    tally['nodes'] += 1; n_ = tally['nodes']
    gens = nd['gens']; e = nd['dim']
    Wb = [[sp.Rational(x) for x in w] for w in nd.get('weights') or []]
    Wm = sp.Matrix(Wb) if Wb else sp.zeros(0, NCL)
    Wint = [list(r) for r in (Wm.rref()[0].tolist() if Wb else [])]       # partition spans: rref rows = block indicators
    assert all(all(x in (0, 1) for x in r) for r in Wint) and (not Wint or sum(sum(r) for r in Wint) == NCL), ('W not a partition span', Wint)
    r_ = len(Wint); ts = sp.symbols(f't0:{r_}') if r_ else (); us = sp.symbols(f'u0:{r_}') if r_ else ()
    f = G_['parse'](nd['minor']) if nd.get('minor') else None
    _, iv, _ = G_['info'](gens, f'i{n_}')
    node = dict(path=nd['path'], dim=e, W=r_, points=[])
    Iy = [clear(sp.sympify(g.replace('^', '**'))) for g in gens if g.strip() not in ('', '0')]
    seen = []; want = set()
    def gen_points():
        for attempt in range(400):
            rat, _ = G_['points'](gens, iv, 2)
            for p_ in rat:
                if p_ in seen: continue
                seen.append(p_)
                if want and any(q in den_primes(x) for x in p_ for q in want): continue   # targeted: q-integral points only
                if f is None or G_['evalr'](f, p_) != 0: yield p_
    pgen = gen_points(); pts = []; pi = -1; g1in = {}; sdata = {}
    while True:                                             # adaptive: more points while a candidate ≥ 5 survives
        goodp = [p for p in node['points'] if all(v for _, v in p['checks'])]
        if len(node['points']) >= MAXPTS or node.get('G1_timeout'): break
        if len(goodp) >= NPTS:
            resid = set.intersection(*[set(p['cand']) for p in goodp]) - SMALLP
            want = {q for q in resid if not any(set(p['src'].get(str(q), [])) <= {'G1lc', 'Slc'} for p in goodp)}
            if not want: break                              # every residual prime is directly checkable at some point
        pt = next(pgen, None)
        if pt is None: break
        pts.append(pt); pi += 1
        cand = set(setup_primes); checks = []; src = {}
        def addc(tag, ps):
            for q in ps: cand.add(q); src.setdefault(q, set()).add(tag)
        for x in pt: addc('den', den_primes(x))
        Px = sp.Matrix(d, d, lambda a, b: G_['evalr'](Pr[a][b], pt) if Pr[a][b] else 0)
        for x in Px: addc('den', den_primes(x))
        xv = [pt[G_['gpos'][v]] for v in vs]
        checks.append(('coords', all(sp.simplify(coord(Px, v) - xv[j]) == 0 for j, v in enumerate(vs))))
        # W: defect matrix at x
        Dx = sp.Matrix([[G_['evalr'](defects[c][k], pt) if defects[c][k] else 0 for c in range(NCL)] for k in nzd]) if nzd else sp.zeros(0, NCL)
        wp, wr = smith_primes(Dx); addc('W', wp)
        checks.append(('W_kernel', NCL - wr == r_ and all((Dx * sp.Matrix(w)).is_zero_matrix for w in Wint) if Dx.rows else r_ == NCL))
        # S: stabiliser
        Sg = sigma(Gm, scal_of(Wint, ts))
        eqs = [clear(q) for q in (Sg * Px - Px * Sg) if q != 0]
        eqs = list(dict.fromkeys(q for q in eqs if q != 0))
        checks.append(('S_affine_g', all(sp.Poly(q, *gvars).total_degree() <= 1 for q in eqs) if gvars else True))
        zz = sp.Symbol('zz'); dets = sp.Integer(1)
        for c in big: dets *= Gm[c].det()
        rel = [sp.expand(zz * dets - 1)] + [ts[i] * us[i] - 1 for i in range(r_)]
        allS = eqs + rel; hv = gvars + [zz]; tv = list(ts) + list(us)
        sdata[pi] = (hv, tv, [ssy(q) for q in eqs + rel])
        at1 = {v: 0 for v in gvars}
        for c in big:
            for i in range(len(cl[c]) - 1): at1[sp.Symbol(f"g{c}_{i}{i}")] = 1
        at1.update({zz: 1}); at1.update({v: 1 for v in tv})
        Jm = sp.Matrix([[sp.diff(q, v).subs(at1) for v in hv + tv] for q in allS])
        jp, jr = smith_primes(Jm); addc('Sjac', jp)
        dimH = dimG + r_
        checks.append(('S_tangent', len(hv + tv) - jr == dimH - e))
        if hv + tv:
            o = run_sing([f"ring RZ = integer,({','.join(map(str, hv + tv))}),(dp({len(hv)}),dp({len(tv) or 1}));" if tv else
                          f"ring RZ = integer,({','.join(map(str, hv))}),dp;", "option(redSB);", LCPROC,
                          f"ideal I = {','.join(ssy(q) for q in allS)};", "ideal G = std(I);", 'print("LC" + lcs(G));',
                          'print("END"); quit;'], f's{n_}_{pi}')
            if o == 'TIMEOUT' or o is None: tally['timeouts'] += 1; checks.append(('S_zbasis', False)); o = ''
            else:
                for c in re.search(r"LC(.*)", o).group(1).split(): addc('Slc', primes_of(c))
                checks.append(('S_zbasis', True))
        # G1: orbit closure, on a slice Σ through 1 transversal to Stab: the e coordinates of (g, t) at which the
        # differential of h ↦ h·x has full rank. Over ℚ: the reduced Gröbner basis G of the graph ideal (block order
        # Σ-side > y). If G is p-integral (p ∤ leading coefficients of its primitive integral elements), the inputs reduce
        # to 0 mod G by p-integral steps, so Γ_p ⊆ (Ḡ), Ḡ is a Gröbner basis, and ker φ_p* = Γ_p ∩ 𝔽_p[y] ⊆ (Ē) = I_Y ⊗ 𝔽_p;
        # the other inclusion holds since I_Y vanishes on the orbit. So Y_p = closure(φ_p(Σ_p)): irreducible, with x_p in
        # its dense orbit.
        gsort = sorted(gvars, key=lambda v: (len(cl[int(str(v)[1:].split('_')[0])]), str(v)))   # cheapest slice first
        if pi % 3 == 1: gsort = gsort[::-1]                 # vary the slice between points (any slice is a certificate)
        elif pi % 3 == 2: gsort = sorted(gsort, key=lambda v: zlib.crc32(f'{v},{pi}'.encode()))
        hv0 = list(ts) + gsort
        at1s = {v: at1[v] for v in hv0}
        Tm = sp.Matrix([[coord(Dm * Px - Px * Dm, w) for w in vs] for Dm in [sp.diff(Sg, v).subs(at1s) for v in hv0]]).T
        piv = Tm.rref()[1]
        checks.append(('G1_rank', len(piv) == e))
        hv1, yv, g1eqs = g1_graph([hv0[k_] for k_ in piv], Px, Sg, at1, ts, us, Wint)
        g1in[pi] = (hv0, Tm, Px, Sg, at1, ts, us, Wint)
        if hv1:
            o = run_sing([f"ring RQ = 0,({','.join(map(str, hv1 + yv))}),(dp({len(hv1)}),dp({len(yv)}));", "option(redSB);",
                          f"ideal I = {','.join(g1eqs)};", "ideal G = std(I);",
                          'int i; string s = "LC"; for (i = 1; i <= size(G); i++) { s = s + " " + string(leadcoef(cleardenom(G[i]))); } print(s);',
                          f"int j; int s0; intvec ev; ideal E; for (i = 1; i <= size(G); i++) {{ ev = leadexp(G[i]); s0 = 0; for (j = 1; j <= {len(hv1)}; j++) {{ s0 = s0 + ev[j]; }} if (s0 == 0) {{ E = E + G[i]; }} }}",
                          "ideal SE = std(E);", f"ideal IY = {','.join(ssy(q) for q in Iy) or '0'}; ideal SY = std(IY);",
                          'print("EQ " + string(size(reduce(IY, SE)) == 0 && size(reduce(E, SY)) == 0));', 'print("NG " + string(size(G)));', 'print("END"); quit;'], f'g{n_}_{pi}')
            if o == 'TIMEOUT' or o is None: tally['timeouts'] += 1; checks.append(('G1', False)); node['G1_timeout'] = True
            else:
                for c in re.search(r"LC(.*)", o).group(1).split(): addc('G1lc', primes_of(c))
                checks.append(('G1', re.search(r"EQ (\d)", o).group(1) == '1'))
        sp_, sb_ = sep_primes(nd['path'], e, pt); addc('sep', sp_)
        checks.append(('separated', not sb_))
        ok = all(v for _, v in checks)
        if not ok: tally['checks_failed'] += 1
        else: tally['certified_points'] += 1
        node['points'].append(dict(x=[str(pt[G_['gpos'][v]]) for v in vs], cand=sorted(cand), src={str(q): sorted(v) for q, v in src.items()}, checks=checks))
    good = [set(p['cand']) for p in node['points'] if all(v for _, v in p['checks'])]
    node['candidates'] = sorted(set.intersection(*good)) if good else None
    # residual primes not exempt that come from G1 alone at some point: check G1 directly in characteristic p there
    # (std over 𝔽_p of the graph ideal; its elimination ideal must equal I_Y ⊗ 𝔽_p; I_Y has a strong ℤ-basis with unit
    # leading coefficients by C1, so its integral generators reduce to I_Y ⊗ 𝔽_p)
    def direct_G1(q, pi_):
        hv0_, Tm_, Px_, Sg_, at1_, ts_, us_, Wint_ = g1in[pi_]
        pv = pivots_mod(Tm_, q)
        if len(pv) != e: return False                       # the orbit map has lower rank mod q at this point
        hv1_, yv_, eq_ = g1_graph([hv0_[k_] for k_ in pv], Px_, Sg_, at1_, ts_, us_, Wint_)
        o = run_sing([f"ring RP = {q},({','.join(map(str, hv1_ + yv_))}),(dp({len(hv1_)}),dp({len(yv_)}));", "option(redSB);",
                      f"ideal I = {','.join(eq_)};", "ideal G = std(I);",
                      f"int i; int j; int s0; intvec ev; ideal E; for (i = 1; i <= size(G); i++) {{ ev = leadexp(G[i]); s0 = 0; for (j = 1; j <= {len(hv1_)}; j++) {{ s0 = s0 + ev[j]; }} if (s0 == 0) {{ E = E + G[i]; }} }}",
                      "ideal SE = std(E);", f"ideal IY = {','.join(ssy(t_) for t_ in Iy) or '0'}; ideal SY = std(IY);",
                      'print("EQ " + string(size(reduce(IY, SE)) == 0 && size(reduce(E, SY)) == 0));', 'print("END"); quit;'], f'gp{n_}_{pi_}_{q}')
        return o not in (None, 'TIMEOUT') and re.search(r"EQ (\d)", o).group(1) == '1'
    def direct_S(q, pi_):
        """Stab over 𝔽_q: one minimal prime (the identity component, defined over 𝔽_q) of dimension dim H − dim Y"""
        hv_, tv_, eq_ = sdata[pi_]
        o = run_sing(['LIB "primdec.lib";', f"ring RP = {q},({','.join(map(str, hv_ + tv_))}),dp;", f"ideal I = {','.join(eq_)};",
                      "list L = minAssGTZ(I); int i;", 'print("NCOMP " + string(size(L)));',
                      'for (i = 1; i <= size(L); i++) { print("CDIM " + string(dim(std(L[i])))); }', 'print("END"); quit;'], f'sp{n_}_{pi_}_{q}')
        if o in (None, 'TIMEOUT'): return False
        cd_ = [int(x) for x in re.findall(r"CDIM (-?\d+)", o)]
        return int(re.search(r"NCOMP (\d+)", o).group(1)) == 1 and cd_ == [dimG + r_ - e]
    # residual primes ≥ 5 that come only from G1 / the stabiliser's ℤ-basis at some point: check those directly in
    # characteristic q at that point (C1 gives I_Y a strong ℤ-basis with unit leading coefficients, so its integral
    # generators reduce to I_Y ⊗ 𝔽_q; the stabiliser's tangent space at 1 is certified by Sjac)
    for q in [q for q in (node['candidates'] or []) if q not in SMALLP]:
        for pi_, p_ in enumerate(node['points']):
            srcq = set(p_['src'].get(str(q), []))
            if not all(v for _, v in p_['checks']) or not srcq or not srcq <= {'G1lc', 'Slc'} or pi_ not in g1in: continue
            if ('G1lc' in srcq and not direct_G1(q, pi_)) or ('Slc' in srcq and not direct_S(q, pi_)): continue
            node['candidates'] = [c for c in node['candidates'] if c != q]
            node.setdefault('direct', []).append([q, pi_, sorted(srcq)]); break
    res.append(node)
    print(f"  {nd['path']}: dim {e}, W {r_}: candidates {[p['src'] for p in node['points']]} -> {node['candidates']}"
          f"{'' if all(all(v for _, v in p['checks']) for p in node['points']) else ' CHECKS ' + str([p['checks'] for p in node['points']])}; t = {time.time() - t0:.0f}s", flush=True)
for a_ in pnodes:                                           # distinct points stay distinct mod p
    for b_ in pnodes:
        if a_ is b_ or '_x' not in a_ or '_x' not in b_: continue
        g_ = 0
        for u_, v_ in zip(a_['_x'], b_['_x']): g_ = math.gcd(g_, abs(int(sp.fraction(sp.Rational(u_ - v_))[0])))
        sep = primes_of(g_) if g_ else {'INSEP'}
        if sep:
            a_.setdefault('sep', {})[b_['path']] = sorted(map(str, sep))
            if a_['candidates'] is not None: a_['candidates'] = sorted(set(map(str, a_['candidates'])) | set(map(str, sep)))
for a_ in pnodes: a_.pop('_x', None)
tally['points'] = len(pnodes); tally['points_failed'] = sum(1 for a_ in pnodes if a_['candidates'] is None)
res += pnodes
allc = set()
for n in res:
    if n['candidates'] is None: allc.add('UNCERTIFIED')
    else: allc |= set(map(str, n['candidates']))
json.dump(dict(block=block, tally=tally, setup=sorted(setup_primes), candidates=sorted(map(str, allc)), nodes=res), open(out + '.json', 'w'))
print(f"RESULT {block}: {tally}; setup primes {sorted(setup_primes)}; candidates {sorted(map(str, allc))}; t = {time.time() - t0:.0f}s", flush=True)
