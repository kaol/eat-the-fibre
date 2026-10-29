"""v2.69: stabilisers at the dense orbits of an orbit tree (v2.67), and the local weightings.

For each positive-dimensional node Y of BLOCK's tree:
  * a rational point x of Y in the dense orbit (sampled; the node's minor is non-zero at x, so rank = dim Y);
  * H_Y = ∏_classes {g ∈ GL_m : column sums 1} × (torus of W(Y)); σ_{g,t}: a_j ↦ s_c Σ_i g_ij a_i on the arrows of
    class c (s_c = ∏_k t_k^{w_k,c}), extended multiplicatively to paths;
  * Stab(x) = {σ : σ P_x = P_x σ} (determinants and t inverted), its minimal primes over ℚ and their dimensions.
    The identity component is defined over ℚ with a smooth rational point, so Stab(x) is connected iff there is
    exactly one minimal prime. Checked: its dimension is dim H_Y − dim Y.
  * the weightings compatible on Y by the defect condition alone, against W(Y) (which also asks for tangency).
Usage: python3 analysis/v2_69_stab.py CERT BLOCK TREE.json OUT"""
import os, sys, json, re, time, itertools, math
HERE = os.path.dirname(os.path.abspath(__file__))
cert, block, tree, out = sys.argv[1:5]; t0 = time.time()
sys.argv = ['v2_67_orbits.py', cert, block, os.environ.get('STAB_LEAF', '0'), 'all', out + '_o', '--tmo', '3000'] + (['--leafmode'] if os.environ.get('STAB_LEAF') else [])
src = open(os.path.join(HERE, 'v2_67_orbits.py')).read(); src = src[:src.index("stats = {'primes'")]
G_ = {'__name__': 'v2_69', '__file__': os.path.join(HERE, 'v2_67_orbits.py')}
exec(compile(src, 'v2_67_orbits', 'exec'), G_)
sp = G_['sp']; R = G_['R']; vs = G_['vs']; st = G_['st']; ctx = G_['ctx']; d = ctx['d']; basis = ctx['basis']; Pr = G_['Pr']
T = json.load(open(tree))
cl = st['classes']; arrow_class = {a: c for c, arr in enumerate(cl) for a in arr}
idx = {b: i for i, b in enumerate(basis)}
# group coordinates: for a class of size m, an m × m matrix g with column sums 1 (last row eliminated)
gvars = []; Gm = []
for c, arr in enumerate(cl):
    m = len(arr); M = [[None] * m for _ in range(m)]
    for j in range(m):
        for i in range(m - 1):
            v = sp.Symbol(f"g{c}_{i}{j}"); gvars.append(v); M[i][j] = v
        M[m - 1][j] = 1 - sum(M[i][j] for i in range(m - 1))
    Gm.append(sp.Matrix(M))
dimG = len(gvars)
assert dimG == len(st['DerD']), (dimG, len(st['DerD']))      # Der_Δ is the Lie algebra of these substitution groups
def sigma(Wint):
    """the d × d matrix of σ_{g,t} (columns: images of basis vectors), polynomial in g, t, u = 1/t"""
    r = len(Wint); ts = sp.symbols(f't0:{r}') if r else (); us = sp.symbols(f'u0:{r}') if r else ()
    def s(c):
        e = sp.Integer(1)
        for k in range(r):
            w = Wint[k][c]; e *= ts[k] ** w if w >= 0 else us[k] ** (-w)
        return e
    S = sp.zeros(d, d)
    for col, b in enumerate(basis):
        if b[0] == 'v': S[col, col] = 1; continue
        choices = []
        for a in b[1]:
            c = arrow_class[a]; arr = cl[c]; j = arr.index(a)
            choices.append([(arr[i], s(c) * Gm[c][i, j]) for i in range(len(arr))])
        for combo in itertools.product(*choices):
            path = ('p', tuple(x for x, _ in combo)); coef = sp.Integer(1)
            for _, cf in combo: coef *= cf
            S[idx[path], col] += coef
    return S, ts, us
def integral_lattice(Wb):
    out_ = []
    for w in Wb:
        den = 1
        for x in w: den = den * int(sp.fraction(x)[1]) // math.gcd(den, int(sp.fraction(x)[1]))
        wi = [int(x * den) for x in w]; g = 0
        for x in wi: g = math.gcd(g, abs(x))
        out_.append([x // g for x in wi] if g > 1 else wi)
    return out_
res = []; tally = {'nodes': 0, 'connected': 0, 'disconnected': 0, 'dimension_mismatch': 0, 'no_point': 0, 'W_defect_larger': 0}
for k, nd in T['nodes'].items():
    if nd.get('result') != 'dense orbit': continue
    tally['nodes'] += 1
    gens = nd['gens']; dim_ = nd['dim']
    Wb = [[sp.Rational(x) for x in w] for w in nd.get('weights') or []]
    G_['local_weights'](gens, f'w{tally["nodes"]}'); wdef = G_['LASTW'].get('defect_only')
    if wdef is not None and wdef > len(Wb): tally['W_defect_larger'] += 1
    f = G_['parse'](nd['minor']) if nd.get('minor') else None
    _, iv, _ = G_['info'](gens, f'i{tally["nodes"]}')
    pt = None
    for attempt in range(4):
        rat, _ = G_['points'](gens, iv, 4)
        for p_ in rat:
            if f is None or G_['evalr'](f, p_) != 0: pt = p_; break
        if pt: break
    if pt is None: tally['no_point'] += 1; res.append(dict(path=nd['path'], dim=dim_, result='no rational point in the dense orbit')); continue
    Px = sp.Matrix(d, d, lambda a, b: G_['evalr'](Pr[a][b], pt) if Pr[a][b] else 0)
    Wint = integral_lattice(Wb)
    S, ts, us = sigma(Wint)
    eqs = [e for e in (S * Px - Px * S) if e != 0]
    z = sp.Symbol('zz')
    dets = sp.Integer(1)
    for M in Gm:
        if M.shape[0] > 1: dets *= M.det()
    # denominators cleared (Singular reads x^2/31 as x^(2/31))
    eqs = [sp.expand(sp.fraction(sp.together(sp.expand(e)))[0]) for e in eqs] + [sp.expand(z * dets - 1)] + [ts[i] * us[i] - 1 for i in range(len(ts))]
    eqs = list(dict.fromkeys(eqs))
    allv = gvars + list(ts) + list(us) + [z]
    ssy = lambda e: str(e).replace('**', '^')
    o = G_['singular']([f'LIB "primdec.lib"; ring RS = 0,({",".join(map(str, allv))}),dp;', "option(redSB);",
                        f"ideal I = {','.join(ssy(e) for e in eqs)};", "list L = minAssGTZ(I); int i;",
                        'print("NCOMP " + string(size(L)));', 'for (i = 1; i <= size(L); i++) { print("CDIM " + string(dim(std(L[i])))); }',
                        'print("END"); quit;'], f's{tally["nodes"]}')
    if o is None: res.append(dict(path=nd['path'], dim=dim_, result='singular error')); continue
    ncomp = int(re.search(r"NCOMP (\d+)", o).group(1)); cd = [int(x) for x in re.findall(r"CDIM (-?\d+)", o)]
    dimH = dimG + len(Wint); expect = dimH - dim_
    ok_dim = (max(cd) == expect) if cd else False
    if not ok_dim: tally['dimension_mismatch'] += 1
    tally['connected' if ncomp == 1 else 'disconnected'] += 1
    res.append(dict(path=nd['path'], dim=dim_, dimH=dimH, stab_components=ncomp, stab_dims=cd, expected=expect, W=len(Wb), W_defect=wdef))
    if ncomp != 1 or not ok_dim:
        print(f"  {nd['path']}: dim {dim_}, dim H {dimH}: stabiliser components {ncomp}, dims {cd} (expected {expect})", flush=True)
json.dump(dict(block=block, tally=tally, nodes=res), open(out + '.json', 'w'))
print(f"RESULT {block}: {tally}; t = {time.time() - t0:.0f}s", flush=True)
