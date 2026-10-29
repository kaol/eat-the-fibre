"""v2.71: Lemma L directly in characteristic p, on a char-p orbit tree (v2_68_charp.py output), for the blocks whose
char 0 tree does not specialise at p.

The char-p tree (roots: minimal primes of the stratum over 𝔽_p; children: minimal primes of closures under the fields
of the fixed group H_fix = G × T_{W_root}) gives, by Lemma C at p, that every H_fix-orbit of S_p is a top orbit of a
node. For each node Z of dimension e > 0 this script finds an 𝔽_p-point x of Z and checks, over 𝔽_p:
  W   W_x = ker D(x) (the census weightings at x) is a partition span; H_x = G × T_{W_x} ⊇ H_fix;
  rk  the orbit matrix of H_x at x has rank e (pivots mod p), so dim H_x·x = e;
  G1  the elimination ideal of the graph of a slice Σ ⊂ H_x through 1 (pivot coordinates) equals I_Z over 𝔽_p: Z is
      the closure of H_x·x, geometrically irreducible; as H_fix ◁ H_x and every H_fix-orbit in H_x·x is a translate of
      the open H_fix-orbit, H_x·x = O_Z (one census class at most, as moves at points of O_Z lie in H_x);
  S   the stabiliser of x in H_x over 𝔽_p: affine-linear in g, one minimal prime of dimension dim H_x − e, tangent
      space at 1 of that dimension (Jacobian rank mod p): smooth and connected, so Lang gives exactly one class;
  sep x ∉ Z' for every other node Z' of dimension e (some generator of I_Z' non-zero at x).
Or, for a node Z that is not H-stable at p (possible in characteristic p, where the closures are only ξ-stable): a
point where the orbit matrix has rank > e. The generic rank on Z is then > e, so no orbit of dimension e lies in Z (a
dense orbit O ⊆ Z would meet the open set of rank > e, but rank ≤ dim O at every point of O): Z contributes no class,
and its lower orbits lie in its children by Lemma C.
For a node of dimension 0: a single 𝔽_p-point, fixed by G and by every weighting move compatible at it.
Then the census count of the block over every 𝔽_{p^k} equals the number of certified nodes.
Usage: python3 analysis/v2_71_charp_cert.py CERT BLOCK P CHARP_TREE.json OUT"""
import os, sys, json, re, time, random
HERE = os.path.dirname(os.path.abspath(__file__))
cert, block, P, tree, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5]
sys.argv = ['v2_71_gcert.py', cert, block, tree, out]
src = open(os.path.join(HERE, 'v2_71_gcert.py')).read(); src = src[:src.index("setup_primes = set()")]
G2 = {'__name__': 'v2_71c', '__file__': os.path.join(HERE, 'v2_71_gcert.py')}
exec(compile(src, 'v2_71_gcert', 'exec'), G2)
globals().update({k: G2[k] for k in ('G_', 'sp', 'R', 'vs', 'st', 'd', 'Pr', 'defects', 'nzd', 'NCL', 'cl', 'Gm', 'gvars', 'big', 'zc',
                                      'sigma', 'scal_of', 'coord', 'run_sing', 'clear', 'g1_graph', 'pivots_mod', 'ssy', 'dimG', 'T', 't0')})
rng = random.Random(71 + P)
red = lambda q: int(sp.fraction(sp.Rational(q))[0]) * pow(int(sp.fraction(sp.Rational(q))[1]), -1, P) % P
def kernel_mod(rows, ncols):
    """basis of {w ∈ 𝔽_p^ncols : rows·w = 0}"""
    M = [[x % P for x in r] for r in rows]; piv = []; r = 0
    for c in range(ncols):
        k = next((i for i in range(r, len(M)) if M[i][c]), None)
        if k is None: continue
        M[r], M[k] = M[k], M[r]; inv = pow(M[r][c], -1, P); M[r] = [x * inv % P for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]: f_ = M[i][c]; M[i] = [(a - f_ * b) % P for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
    free = [c for c in range(ncols) if c not in piv]; basis = []
    for fc in free:
        w = [0] * ncols; w[fc] = 1
        for i, pc in enumerate(piv): w[pc] = (-M[i][fc]) % P
        basis.append(w)
    return basis
def rank_mod(rows):
    return len(rows[0]) - len(kernel_mod(rows, len(rows[0]))) if rows and rows[0] else 0
def partition_rows(basis):
    """the block indicators if span(basis) is a partition span over 𝔽_p, else None"""
    if not basis: return None
    rows = [list(r) for r in basis]; ncol = len(rows[0])
    # reduced row echelon form mod p
    M = [r[:] for r in rows]; piv = []; r = 0
    for c in range(ncol):
        k = next((i for i in range(r, len(M)) if M[i][c] % P), None)
        if k is None: continue
        M[r], M[k] = M[k], M[r]; inv = pow(M[r][c], -1, P); M[r] = [x * inv % P for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] % P: f_ = M[i][c]; M[i] = [(a - f_ * b) % P for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
    M = M[:r]
    if not all(all(x in (0, 1) for x in row) for row in M): return None
    if any(sum(row[c] for row in M) != 1 for c in range(ncol)): return None
    return M
def find_points(gens, tag, want=4, rounds=4):
    """𝔽_p-points of V(gens): random values on an independent set, then the linear minimal primes"""
    pts = []
    for rd in range(rounds):
        tries = []
        for t in range(25): tries.append([rng.randrange(P) for _ in vs])
        lines = ['LIB "primdec.lib";', f"ring Rp = {P},({','.join(vs)}),dp;", "option(redSB);", f"ideal Z = std(ideal({gens}));",
                 'intvec iv = indepSet(Z); print("IV " + string(iv)); int i; int j; ideal S; list L;']
        for t, vals in enumerate(tries):
            lines.append(f"ideal A{t} = {','.join(f'{v}-({vals[j]})' for j, v in enumerate(vs))};")
            lines.append(f"S = 0; for (j = 1; j <= {len(vs)}; j++) {{ if (iv[j] == 1) {{ S = S + A{t}[j]; }} }} S = std(Z + S);"
                         f' if (dim(S) == 0) {{ L = minAssGTZ(S); for (i = 1; i <= size(L); i++) {{ if (vdim(std(L[i])) == 1) {{ print("PT " + string(std(L[i]))); }} }} }}')
        lines.append('print("END"); quit;')
        o = run_sing(lines, f'{tag}_pts{rd}')
        if o in (None, 'TIMEOUT'): continue
        for m in re.finditer(r"^PT (.*)$", o, re.M):
            sol = {}
            for g in m.group(1).split(','):
                e_ = sp.sympify(g.replace('^', '**')); v = [s for s in e_.free_symbols]
                if len(v) != 1: break
                c1 = e_.coeff(v[0]); c0 = e_.subs(v[0], 0)
                sol[str(v[0])] = int(-c0 * pow(int(c1), -1, P)) % P
            if len(sol) == len(vs) and sol not in pts: pts.append(sol)
        if len(pts) >= want: break
    return pts
NODES = T['nodes']; byd = {}
ONLY = os.environ.get('CHARP_ONLY')
for k, nd in NODES.items(): byd.setdefault(nd['dim'], []).append(k)
def gens_eval(k, sol):
    out_ = []
    for g in k.split(','):
        e_ = sp.sympify(g.replace('^', '**')); out_.append(int(e_.subs({sp.Symbol(v): sol[v] for v in vs})) % P)
    return out_
res = []; tally = dict(nodes=0, certified=0, no_top_orbit=0, failed=0, points=0)
for k, nd in NODES.items():
    e = nd['dim']; n_ = len(res) + 1; rec = dict(path=nd['path'], dim=e)
    if ONLY and nd['path'] not in ONLY.split(','): continue
    if e <= 0:
        # a 0-dimensional node over 𝔽_p: its geometric points are conjugate. It contributes one class per geometric point
        # (all defined over 𝔽_p iff vdim = 1) if the points are fixed by G and by every weighting compatible at them;
        # if G moves them, they lie in a larger orbit and the node contributes nothing (no orbit is contained in it).
        tally['points'] += 1
        Sg = sigma(Gm, lambda c: sp.Integer(1))
        PxS = sp.Matrix(d, d, lambda a, b: Pr[a][b].as_expr() if Pr[a][b] else 0)
        com = [sp.expand(x) for x in (Sg * PxS - PxS * Sg) if x != 0]
        o = run_sing([f"ring RV = {P},({','.join(vs)}),dp;", f"ideal QV = std(ideal({k}));", 'print("VDIM " + string(vdim(QV)));',
                      f"ring RP = {P},({','.join(vs + [str(g) for g in gvars])}),dp;", f"ideal Q = std(ideal({k}));",
                      f"ideal C = {','.join(ssy(x) for x in com) or '0'};",
                      'print("GFIX " + string(size(reduce(C, Q)) == 0));', 'print("END"); quit;'], f'pt{n_}')
        vd = int(re.search(r"VDIM (-?\d+)", o).group(1)); gfix = re.search(r"GFIX (\d)", o).group(1) == '1'
        rec.update(vdim=vd, G_fixes=gfix)
        if not gfix: rec['result'] = 'no orbit (G moves the point)'; tally['no_top_orbit'] += 1
        elif vd != 1: rec['result'] = 'FAILED (fixed, not rational)'; tally['failed'] += 1
        else:
            sol = find_points(k, f'p{n_}', want=1, rounds=1)[0]; pt = [0] * len(R.gens)
            for v in vs: pt[G_['gpos'][v]] = sp.Integer(sol[v])
            Px = sp.Matrix(d, d, lambda a, b: red(G_['evalr'](Pr[a][b], pt)) if Pr[a][b] else 0)
            Dx = [[red(G_['evalr'](defects[c][kk], pt)) if defects[c][kk] else 0 for c in range(NCL)] for kk in nzd]
            ker = kernel_mod(Dx, NCL) if Dx else [[1 if i == j else 0 for i in range(NCL)] for j in range(NCL)]
            lam = sp.symbols(f'lam0:{len(ker)}') if ker else ()
            s_ = [sum((lam[j] * ker[j][c] for j in range(len(ker))), sp.Integer(0)) for c in range(NCL)]
            Sw = sigma([sp.eye(len(arr)) for arr in cl], lambda c: s_[c])
            wfix = all((sp.Poly(x, *lam, modulus=P).is_zero if lam else x % P == 0) for x in (Sw * Px - Px * Sw).applyfunc(sp.expand))
            if wfix: rec['result'] = 'point certified'; tally['certified'] += 1
            else: rec['result'] = 'no orbit (a weighting move moves the point)'; tally['no_top_orbit'] += 1
        res.append(rec); print(f"  {nd['path']}: point: {rec['result']}", flush=True); continue
    tally['nodes'] += 1
    rec['tried'] = []; seen_ = []
    def stream():
        for batch in range(int(os.environ.get('CHARP_BATCHES', '8'))):     # fresh points until one certifies
            for sol in find_points(k, f'n{n_}_{batch}', want=4, rounds=4):
                if sol not in seen_: seen_.append(sol); yield sol
    for sol in stream():
        pt = [0] * len(R.gens)
        for v in vs: pt[G_['gpos'][v]] = sp.Integer(sol[v])
        Px = sp.Matrix(d, d, lambda a, b: red(G_['evalr'](Pr[a][b], pt)) if Pr[a][b] else 0)
        chk = {}
        chk['coords'] = all((coord(Px, v) - sol[v]) % P == 0 for v in vs) if all(True for _ in vs) else True
        Dx = [[red(G_['evalr'](defects[c][kk], pt)) if defects[c][kk] else 0 for c in range(NCL)] for kk in nzd]
        Wx = partition_rows(kernel_mod(Dx, NCL) if Dx else [[1 if i == j else 0 for i in range(NCL)] for j in range(NCL)])
        chk['W_partition'] = Wx is not None
        if Wx is None: rec['tried'].append(chk); continue
        r_ = len(Wx); ts = sp.symbols(f't0:{r_}'); us = sp.symbols(f'u0:{r_}')
        Sg = sigma(Gm, scal_of(Wx, ts))
        at1 = {v: 0 for v in gvars}
        for c in big:
            for i in range(len(cl[c]) - 1): at1[sp.Symbol(f"g{c}_{i}{i}")] = 1
        zz = sp.Symbol('zz'); at1.update({zz: 1}); at1.update({v: 1 for v in list(ts) + list(us)})
        hv0 = list(ts) + sorted(gvars, key=lambda v: (len(cl[int(str(v)[1:].split('_')[0])]), str(v)))
        Tm = sp.Matrix([[coord(Dm * Px - Px * Dm, w) for w in vs] for Dm in [sp.diff(Sg, v).subs({v_: at1[v_] for v_ in hv0}) for v in hv0]]).T
        pv = pivots_mod(Tm, P); chk['rank'] = len(pv) == e; chk['rank_value'] = len(pv)
        if len(pv) > e:                                     # generic rank on Z > dim Z: no orbit of dimension e lies in Z
            rec['tried'].append(chk); rec.update(result='no top orbit', overrank_point=sol, rank=len(pv)); break
        if not chk['rank']: rec['tried'].append(chk); continue
        hv1, yv, g1eqs = g1_graph([hv0[j] for j in pv], Px, Sg, at1, ts, us, Wx)
        o = run_sing([f"ring RP = {P},({','.join(map(str, hv1 + yv))}),(dp({len(hv1)}),dp({len(yv)}));", "option(redSB);",
                      f"ideal I = {','.join(g1eqs)};", "ideal G = std(I);",
                      f"int i; int j; int s0; intvec ev; ideal E; for (i = 1; i <= size(G); i++) {{ ev = leadexp(G[i]); s0 = 0; for (j = 1; j <= {len(hv1)}; j++) {{ s0 = s0 + ev[j]; }} if (s0 == 0) {{ E = E + G[i]; }} }}",
                      f"ideal SE = std(E); ideal IZ = {k}; ideal SZ = std(IZ);",
                      'print("EQ " + string(size(reduce(IZ, SE)) == 0 && size(reduce(E, SZ)) == 0));', 'print("END"); quit;'], f'g{n_}')
        chk['G1'] = o not in (None, 'TIMEOUT') and re.search(r"EQ (\d)", o).group(1) == '1'
        eqs = list(dict.fromkeys(q for q in (clear(q) for q in (Sg * Px - Px * Sg) if q != 0) if q != 0))
        chk['S_affine_g'] = all(sp.Poly(q, *gvars).total_degree() <= 1 for q in eqs) if gvars else True
        dets = sp.Integer(1)
        for c in big: dets *= Gm[c].det()
        allS = eqs + [sp.expand(zz * dets - 1)] + [ts[i] * us[i] - 1 for i in range(r_)]
        hv = gvars + [zz] + list(ts) + list(us)
        Jrows = [[int(sp.diff(q, v).subs(at1)) % P for v in hv] for q in allS]
        dimH = dimG + r_
        chk['S_tangent'] = len(hv) - rank_mod(Jrows) == dimH - e
        o = run_sing(['LIB "primdec.lib";', f"ring RP = {P},({','.join(map(str, hv))}),dp;", f"ideal I = {','.join(ssy(q) for q in allS)};",
                      "list L = minAssGTZ(I); int i;", 'print("NCOMP " + string(size(L)));',
                      'for (i = 1; i <= size(L); i++) { print("CDIM " + string(dim(std(L[i])))); }', 'print("END"); quit;'], f's{n_}')
        chk['S_connected'] = o not in (None, 'TIMEOUT') and int(re.search(r"NCOMP (\d+)", o).group(1)) == 1 and \
            [int(x) for x in re.findall(r"CDIM (-?\d+)", o)] == [dimH - e]
        chk['separated'] = all(any(gens_eval(k2, sol)) for k2 in byd[e] if k2 != k)
        rec['tried'].append(chk)
        if all(v for k_, v in chk.items() if k_ != 'rank_value'): rec.update(result='certified', W_x=r_, point=sol); break
    if rec.get('result') not in ('certified', 'no top orbit'): rec['result'] = 'FAILED'
    tally[{'certified': 'certified', 'no top orbit': 'no_top_orbit'}.get(rec['result'], 'failed')] += 1
    res.append(rec)
    print(f"  {nd['path']}: dim {e}: {rec['result']} ({len(rec['tried'])} points tried{'' if rec['result'] == 'certified' else ': ' + str(rec['tried'][-2:])}); t = {time.time() - t0:.0f}s", flush=True)
json.dump(dict(block=block, p=P, tally=tally, nodes=res), open(out + '.json', 'w'))
print(f"RESULT {block} at p = {P}: {tally}; census count over every 𝔽_(p^k) = {tally['certified'] if tally['failed'] == 0 and not ONLY else 'NOT ESTABLISHED' if tally['failed'] else 'partial run'}; t = {time.time() - t0:.0f}s", flush=True)
