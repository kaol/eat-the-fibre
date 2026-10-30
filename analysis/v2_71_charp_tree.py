"""v2.71 (from v2_68_charp.py: one group per root, its weightings computed over 𝔽_p on that root): orbit stratification directly in characteristic p ("contained orbits"), for the primes where the char 0
tree does not specialise.

Over 𝔽̄_p (p > 0) the Lie algebra no longer controls the group, so this version never uses H-stability below the
roots. Its invariant: finitely many H-orbits are *contained* in each node Z. With e = dim Z and an e × e minor f of
the orbit matrix with f ∉ I_Z:
  * orbits contained in Z of dimension e are dense in components of Z: at most #components of them;
  * an orbit O ⊆ Z of dimension < e has rank < e at each point (rank ≤ orbit dimension in any characteristic), so
    O ⊆ V(f); O is H-stable and every field ξ̃ is the tangent of a curve in the smooth group H_p, so every element
    of the closure J of I_Z + (f) under the fields vanishes on O; O is irreducible, so it lies in one minimal prime
    of J: recurse.
Roots: the minimal primes of the stratum ideal over 𝔽_p (computed here, so no coverage claim is needed). A root Z is
stable under H_Z = Der_Δ + (weightings compatible on Z): Der_Δ is connected and preserves the stratum, and the
torus of a compatible weighting preserves the closed set of w-compatible folds, of which Z is a component. So all
orbits through points of Z are contained in Z, and the classes meeting the stratum are finitely many if every node
gets a minor. A node without one (generic rank < dim) is reported: the method is inconclusive there.
Minors are found at closed points of Z ∩ (random 𝔽_p-linear section), by fraction-free elimination modulo the point.
Usage: python3 analysis/v2_68_charp.py CERT BLOCK P TREE.json OUT"""
import os, sys, json, re, time, random, math
HERE = os.path.dirname(os.path.abspath(__file__))
cert, block, P, tree, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5]
t0 = time.time(); rng = random.Random(68 + P)
sys.argv = ['v2_67_orbits.py', cert, block, os.environ.get('STAB_LEAF', '0'), 'all', out + '_o', '--tmo', os.environ.get('CHARP_TMO', '600')] + (['--leafmode'] if os.environ.get('STAB_LEAF') else [])   # v2.73: leaf charts
src = open(os.path.join(HERE, 'v2_67_orbits.py')).read(); src = src[:src.index("stats = {'primes'")]
G_ = {'__name__': 'v2_68p', '__file__': os.path.join(HERE, 'v2_67_orbits.py')}
exec(compile(src, 'v2_67_orbits', 'exec'), G_)
sp = G_['sp']; R = G_['R']; vs = G_['vs']; nv = len(vs); ss = G_['ss']; st = G_['st']
T = json.load(open(tree))
def lcm(a, b): return a * b // math.gcd(a, b)
def row_scale(row):
    den = 1
    for e in row:
        if e:
            for c in e.coeffs(): den = lcm(den, int(c.denominator))
    return den
def kernel_mod(rows, ncols):
    M = [[x % P for x in r] for r in rows]; piv = []; r = 0
    for c in range(ncols):
        k = next((i for i in range(r, len(M)) if M[i][c]), None)
        if k is None: continue
        M[r], M[k] = M[k], M[r]; inv = pow(M[r][c], -1, P); M[r] = [x * inv % P for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]: f_ = M[i][c]; M[i] = [(a - f_ * b) % P for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
    out_ = []
    for fc in [c for c in range(ncols) if c not in piv]:
        w = [0] * ncols; w[fc] = 1
        for i, pc in enumerate(piv): w[pc] = (-M[i][fc]) % P
        out_.append(w)
    return out_
def partition_rows(basis):
    if not basis: return None
    M = [r[:] for r in basis]; ncol = len(M[0]); r = 0
    for c in range(ncol):
        k = next((i for i in range(r, len(M)) if M[i][c] % P), None)
        if k is None: continue
        M[r], M[k] = M[k], M[r]; inv = pow(M[r][c], -1, P); M[r] = [x * inv % P for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] % P: f_ = M[i][c]; M[i] = [(a - f_ * b) % P for a, b in zip(M[i], M[r])]
        r += 1
    M = M[:r]
    if not all(all(x in (0, 1) for x in row) for row in M) or any(sum(row[c] for row in M) != 1 for c in range(ncol)): return None
    return M
NCL_ = G_['NCL']; nzd_ = G_['nzd']; defects_ = G_['defects']
def root_weightings(gens, tag):
    """the weightings w ∈ 𝔽_p^classes with Σ_c w_c (Δ_c-defect) ∈ I_Z (computed mod p), as block indicators"""
    lines = [f"ring Rp = {P},({','.join(vs)}),dp;", "option(redSB);", f"ideal GZ = std(ideal({gens}));"]
    for c in range(NCL_):
        lines.append(f"ideal D{c} = {','.join(ss(defects_[c][k]) for k in nzd_) or '0'};")
        lines.append(f'ideal N{c} = reduce(D{c}, GZ); int i{c}; for (i{c} = 1; i{c} <= ncols(N{c}); i{c}++) {{ print("NF {c} " + string(i{c}) + " " + string(N{c}[i{c}])); }}')
    lines.append('print("END"); quit;')
    o = G_['singular'](lines, tag)
    eqs = {}
    for m in re.finditer(r"^NF (\d+) (\d+) (.*)$", o, re.M):
        c, k, poly = int(m.group(1)), int(m.group(2)), m.group(3).strip()
        if poly == '0': continue
        for mon, co in G_['parse'](poly).terms():
            v = eqs.setdefault((k, mon), [0] * NCL_); v[c] = (v[c] + int(co.numerator) * pow(int(co.denominator), -1, P)) % P
    ker = kernel_mod(list(eqs.values()), NCL_) if eqs else [[1 if i == j else 0 for i in range(NCL_)] for j in range(NCL_)]
    return partition_rows(ker)
GROUPS = {}
def group(Wrows):
    key = json.dumps(Wrows)
    if key in GROUPS: return GROUPS[key]
    Wb = [[sp.Integer(x) for x in w] for w in Wrows]
    act = G_['active'](Wb); XI = G_['combine'](act, G_['XIB']); MR = G_['combine'](act, G_['MB'])
    sc = [lcm(row_scale(XI[a]), row_scale(MR[a])) for a in range(len(act))]
    XI = [[e * R(sc[a]) if e else e for e in XI[a]] for a in range(len(act))]
    MR = [[e * R(sc[a]) if e else e for e in MR[a]] for a in range(len(act))]
    NH = len(XI); NC = len(MR[0])
    nterm = lambda e: len(e.terms()) if e else 0
    corder = sorted(range(NC), key=lambda j: (sum(nterm(MR[i][j]) for i in range(NH)), j))
    HEAD = [f'LIB "primdec.lib";', f"ring Rp = {P},({','.join(vs)}),dp;", "option(redSB);",
            f"matrix XI[{NH}][{nv}] = " + ",".join(ss(XI[a][b]) for a in range(NH) for b in range(nv)) + ";",
            f"proc xi(poly h, int a) {{ poly g = 0; int b; for (b = 1; b <= {nv}; b++) {{ g = g + diff(h, var(b)) * XI[a, b]; }} return(g); }}",
            f"matrix MM[{NH}][{NC}] = " + ",".join(ss(MR[a][corder[b]]) for a in range(NH) for b in range(NC)) + ";"]
    GROUPS[key] = dict(HEAD=HEAD, NH=NH, NC=NC, MR=MR, corder=corder, W=Wrows); return GROUPS[key]
def run(lines, tag): return G_['singular'](lines, tag)
def roots_at_p():
    HEAD = [f'LIB "primdec.lib";', f"ring Rp = {P},({','.join(vs)}),dp;", "option(redSB);"]
    rem = [r.replace('**', '^') for r in (json.loads(os.environ['CHARP_ROOTGENS']) if os.environ.get('CHARP_ROOTGENS') else st['rem'])] or ['0']   # v2.73: an H-stable root's leaf prime
    o = run(HEAD + [f"ideal I = {','.join(rem)};", "list L = minAssGTZ(I); int i;",
                    'for (i = 1; i <= size(L); i++) { ideal Gi = std(L[i]); print("COMP " + string(dim(Gi))); print("GENS " + string(Gi)); kill Gi; }',
                    'print("END"); quit;'], 'roots')
    return [(int(b.split()[0]), b.split("GENS ")[1].split("\n")[0]) for b in o.split("COMP ")[1:]]
def find_minor(gens, e, tag, g_, tries=6):
    HEAD, NH, NC = g_['HEAD'], g_['NH'], g_['NC']
    """pivot rows/columns of an e × e minor non-vanishing at a closed point of Z ∩ (random linear section)"""
    # first: the generic point itself (fraction-free elimination modulo the prime: exact in the domain R/I_Z)
    o = run(HEAD + [f"ideal Q = std(ideal({gens}));", "matrix M = reduce(MM, Q); int r = 0; int i; int j; int i2; int j2; intvec pr; intvec pc;",
                    f"int nr = {NH}; int nc = {NC}; intvec usedr = 0:nr; intvec usedc = 0:nc; int found;",
                    f"for (j = 1; j <= nc && r < {e}; j++) {{ found = 0; for (i = 1; i <= nr && found == 0; i++) {{ if (usedr[i] == 0 && M[i, j] != 0) {{ found = i; }} }}",
                    "  if (found > 0) { r++; usedr[found] = 1; usedc[j] = 1; pr[r] = found; pc[r] = j;",
                    "    for (i2 = 1; i2 <= nr; i2++) { if (usedr[i2] == 0) { for (j2 = 1; j2 <= nc; j2++) { if (usedc[j2] == 0) {",
                    "      M[i2, j2] = reduce(M[found, j] * M[i2, j2] - M[i2, j] * M[found, j2], Q); } } } } } }",
                    'print("RANK " + string(r)); print("PR " + string(pr)); print("PC " + string(pc));', 'print("END"); quit;'], f'{tag}_gen')
    if o is not None:
        r = int(re.search(r"RANK (\d+)", o).group(1))
        if r >= e:
            return ([int(x) for x in re.search(r"PR (.*)", o).group(1).split(',')][:e],
                    [int(x) for x in re.search(r"PC (.*)", o).group(1).split(',')][:e], r)
        print(f"    generic rank {r} < dim {e}", flush=True)
    for t in range(tries):
        lin = [' + '.join(f"{rng.randrange(P)}*{v}" for v in vs) + f" + {rng.randrange(P)}" for _ in range(e)]
        o = run(HEAD + [f"ideal Z = {gens};", f"ideal S = std(Z + ideal({','.join(lin)}));",
                        'if (dim(S) != 0) { print("SKIP"); print("END"); quit; }', "list L = minAssGTZ(S);",
                        "ideal Q = std(L[1]); matrix M = reduce(MM, Q); int r = 0; int i; int j; int i2; int j2; intvec pr; intvec pc;",
                        f"int nr = {NH}; int nc = {NC}; intvec usedr = 0:nr; intvec usedc = 0:nc; int found;",
                        "for (j = 1; j <= nc; j++) { found = 0; for (i = 1; i <= nr && found == 0; i++) { if (usedr[i] == 0 && M[i, j] != 0) { found = i; } }",
                        "  if (found > 0) { r++; usedr[found] = 1; usedc[j] = 1; pr[r] = found; pc[r] = j;",
                        "    for (i2 = 1; i2 <= nr; i2++) { if (usedr[i2] == 0) { for (j2 = 1; j2 <= nc; j2++) { if (usedc[j2] == 0) {",
                        "      M[i2, j2] = reduce(M[found, j] * M[i2, j2] - M[i2, j] * M[found, j2], Q); } } } } } }",
                        'print("RANK " + string(r)); print("PR " + string(pr)); print("PC " + string(pc));', 'print("END"); quit;'], f'{tag}_pt{t}')
        if o is None or "SKIP" in o: continue
        r = int(re.search(r"RANK (\d+)", o).group(1))
        if r < e: continue
        pr = [int(x) for x in re.search(r"PR (.*)", o).group(1).split(',')][:e]
        pc = [int(x) for x in re.search(r"PC (.*)", o).group(1).split(',')][:e]
        return pr, pc, r
    return None
CAND = {}
def char0_pivots(e, g_):
    MR, corder = g_['MR'], g_['corder']
    """pivot rows/columns (MM positions) of e × e minors non-vanishing at rational points of the char 0 nodes of
    dimension e: candidates to try first (each is then tested exactly in characteristic p)"""
    if (e, json.dumps(g_['W'])) in CAND: return CAND[(e, json.dumps(g_['W']))]
    out = []
    for k, nd0 in T['nodes'].items():
        if nd0.get('dim') != e or nd0.get('result') != 'dense orbit': continue
        try:
            d0, iv0, _ = G_['info'](nd0['gens'], f'c0_{len(out)}')
            rat, _ = G_['points'](nd0['gens'], iv0, 2)
        except Exception: continue
        for pt in rat:
            Mv = sp.Matrix([[G_['evalr'](x, pt) for x in row] for row in MR])
            if Mv.rank() < e: continue
            _, pcs = Mv[:, corder].rref(); pcs = list(pcs)[:e]          # positions in MM's column order
            _, prs = Mv[:, [corder[j] for j in pcs]].T.rref(); prs = [i + 1 for i in list(prs)[:e]]
            cand = (tuple(prs), tuple(j + 1 for j in pcs))
            if cand not in out: out.append(cand)
        if len(out) >= 12: break
    CAND[(e, json.dumps(g_['W']))] = out; return out
def try_candidates(gens, e, tag, g_):
    HEAD = g_['HEAD']
    cands = char0_pivots(e, g_)
    if not cands: return None
    lines = HEAD + [f"ideal GZ = std(ideal({gens})); int c = 0; poly f;"]
    for n, (pr, pc) in enumerate(cands):
        lines.append(f"if (c == 0) {{ matrix SM{n}[{e}][{e}] = " + ",".join(f"MM[{i},{j}]" for i in pr for j in pc) +
                     f"; f = det(SM{n}); if (reduce(f, GZ) != 0) {{ c = {n + 1}; }} }}")
    lines += ['print("CAND " + string(c));', 'print("END"); quit;']
    o = run(lines, tag)
    if o is None: return None
    c = int(re.search(r"CAND (\d+)", o).group(1))
    return (list(cands[c - 1][0]), list(cands[c - 1][1]), -1) if c else None
NODES = {}; stats = {'nodes': 0, 'minor': 0, 'points': 0, 'inconclusive': 0, 'error': 0}
def node(gens, dim_, path, g_=None, root=False):
    if gens in NODES: return gens
    nd = NODES[gens] = dict(path=path, dim=dim_); stats['nodes'] += 1; n_ = stats['nodes']
    if dim_ <= 0: stats['points'] += 1; nd['result'] = 'points'; return gens
    if root:   # the root's own group: the weightings compatible on it over 𝔽_p (a partition span)
        Wr = root_weightings(gens, f'rw{n_}')
        if Wr is None: stats['inconclusive'] += 1; nd['result'] = 'root weightings not a partition span'; return gens
        g_ = group(Wr); nd['weights'] = Wr
    HEAD, NH = g_['HEAD'], g_['NH']
    piv = try_candidates(gens, dim_, f'tc{n_}', g_) or find_minor(gens, dim_, f'n{n_}', g_)
    if piv is None:
        stats['inconclusive'] += 1; nd['result'] = 'no minor of size dim found'; print(f"  {path}: dim {dim_}: NO MINOR (inconclusive)", flush=True); return gens
    pr, pc, r = piv
    o = run(HEAD + [f"ideal I = {gens}; ideal G = std(I); int a; int i; poly g;",
                    "matrix SM[" + str(dim_) + "][" + str(dim_) + "] = " + ",".join(f"MM[{i},{j}]" for i in pr for j in pc) + ";",
                    "poly f = det(SM);", 'if (reduce(f, G) == 0) { print("FININ"); print("END"); quit; }',
                    "G = std(G + f); ideal NEW = f; int rounds = 0;",
                    "while (size(NEW) > 0) { rounds++; ideal NX = 0; for (i = 1; i <= size(NEW); i++) { for (a = 1; a <= " + str(NH) +
                    "; a++) { g = reduce(xi(NEW[i], a), G); if (g != 0) { NX = NX + g; G = std(G + g); } } } NEW = simplify(NX, 2); kill NX; }",
                    'print("CLOSED " + string(rounds) + " DIM " + string(dim(G)));', "list L = minAssGTZ(G);",
                    'for (i = 1; i <= size(L); i++) { ideal Gi = std(L[i]); print("COMP " + string(dim(Gi))); print("GENS " + string(Gi)); kill Gi; }',
                    'print("END"); quit;'], f'c{n_}')
    if o is None or "FININ" in o:
        stats['error'] += 1; nd['result'] = 'singular error/timeout' if o is None else 'minor in I_Z'; return gens
    kids = [(int(b.split()[0]), b.split("GENS ")[1].split("\n")[0]) for b in o.split("COMP ")[1:]]
    kids = [k for k in kids if k[0] >= 0]
    assert all(k[0] < dim_ for k in kids)
    stats['minor'] += 1; nd.update(result='minor', rank_at_point=r, children=[])
    print(f"  {path}: dim {dim_}: minor found (rank {r} at a closed point); boundary {[k[0] for k in kids]}; t = {time.time() - t0:.0f}s", flush=True)
    nd['children'] = [node(g, dm, f"{path}.{k + 1}", g_) for k, (dm, g) in enumerate(kids)]
    return gens
roots = roots_at_p()
print(f"{block} at p = {P}: {len(roots)} components of the stratum, dims {[d for d, _ in roots]}; t = {time.time() - t0:.0f}s", flush=True)
for k, (dm, g) in enumerate(roots): node(g, dm, f"Z{k + 1}", root=True)
ok = stats['inconclusive'] == 0 and stats['error'] == 0
json.dump(dict(block=block, p=P, roots=[g for _, g in roots], nodes=NODES, stats=stats), open(out + '.json', 'w'))
print(f"RESULT p = {P}: {'finitely many classes' if ok else 'NOT ESTABLISHED'}: {stats}; t = {time.time() - t0:.0f}s", flush=True)
