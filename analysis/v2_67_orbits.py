"""v2.67: orbit stratification by invariant boundaries (char 0), for one component of a vertex-block stratum.

Let H be the connected group with Lie algebra 𝔥 = Der_Δ + (class weightings W_Z of the component), acting on the
stratum by conjugation; its vector fields in the stratum coordinates vs are ξ_D(x) = the x-coordinate of [D, P(vs)]
(x = y_{r,j}: entry (r, c_j) of [D,P]; x = k_{ij}: entry (i,j) of Kl·[D,P]·Kmat).
For an H-stable prime X:
  * e = generic orbit dimension. An e × e minor f of the orbit matrix (rows D ∈ 𝔥, columns the entries of [D,P]) with
    f ∉ I_X certifies e ≥ rank; e ≤ dim X always. If the sampled points never give rank dim X, X is reported as
    "no dense orbit" (infinitely many orbits on X; not certified here).
  * If e = dim X, X is one dense orbit O plus points of orbit dimension < e, and those lie in the H-core of V(f) ∩ X,
    which is contained in V(J) with J the closure of I_X + (f) under all ξ_D (Taylor expansion of f along exp(tD)).
    J is ξ-stable, so its minimal primes are H-stable (char 0, H connected); recurse on them.
Finitely many H-orbits on X iff every prime met has a dense orbit. Checks done in Singular: ξ_D(I_X) ⊆ I_X for the
input prime (tangency: validates the vector fields), f ∉ I_X, the closure loop to stability, dims of the children.
Writes OUT.json (the tree: gens, dim, minor, children) and prints a summary.
Usage: python3 analysis/v2_67_orbits.py CERT BLOCK LEAF COMP OUT [--seed S] [--tmo SECONDS]"""
import os, sys, json, lzma, re, random, subprocess, time, argparse
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import numpy as np, sympy as sp
import fold_ideal as fi, singular_nohilb
singular_nohilb.install()
ap = argparse.ArgumentParser(); ap.add_argument('cert'); ap.add_argument('block'); ap.add_argument('leaf', type=int)
ap.add_argument('comp'); ap.add_argument('out'); ap.add_argument('--seed', type=int, default=67)
ap.add_argument('--tmo', type=int, default=20000)
ap.add_argument('--leafmode', action='store_true', help='case-split leaf: skip roots that are not H-stable; record full-coordinate keys')
A = ap.parse_args(); rng = random.Random(A.seed); t0 = time.time()
op = lzma.open if A.cert.endswith('.xz') else open
st = None; comps = []
for l in op(A.cert, 'rt'):
    if f'"block": "{A.block}"' not in l: continue
    r = json.loads(l)
    if r['block'] != A.block or r['leaf'] != A.leaf: continue
    if r['kind'] == 'stratum': st = r
    elif A.comp == 'all' or r.get('comp') in {int(x) for x in A.comp.split(',')}: comps.append(r)
assert st and comps, "stratum/component not found"
if not st['vs']:                                            # a point stratum: one fold (or none)
    print(f"RESULT finitely many orbits: point stratum; geometric H-orbits {st['ncomp']}; t = 0s"); sys.exit(0)
Q = st['quiver']; V = [[-1 if ch == 'm' else int(ch) for ch in row] for row in A.block.split(',')]
S = fi.stratum(Q, V); R = S['R']; ctx = S['ctx']; d = ctx['d']; kdim = ctx['kdim']; compl = S['complement']
gl = {str(g): g for g in R.gens}; loc = {k: sp.Symbol(k) for k in gl}
parse = lambda s_: R(sp.sympify(s_.replace('^', '**'), locals=loc))
Pr = [row[:] for row in S['P']]
for x, v in st['subs']:
    xv, vv = gl[x], parse(v); Pr = [[e.compose(xv, vv) for e in row] for row in Pr]
vs = st['vs']; nv = len(vs)
# 𝔥 at a prime Y: Der_Δ plus the class weightings w with Σ_c w_c (Δ_c-defect) ≡ 0 on Y whose vector fields are
# tangent to Y (the local torus of v2.23 at the generic point of Y); computed per prime, so it grows down the tree
DerD = [sp.Matrix(d, d, [sp.Rational(x) for x in D]) for D in st['DerD']]
Ec = [sp.diag(*e) for e in st['Ec']]
NCL = len(Ec); ND = len(DerD); BASE = DerD + Ec
def rc(q): return sp.QQ(int(sp.fraction(q)[0]), int(sp.fraction(q)[1]))
def comm(D):
    Dr = [[rc(D[a, b]) for b in range(d)] for a in range(d)]
    return [[sum((Pr[k][b] * Dr[a][k] for k in range(d) if Dr[a][k]), R(0)) - sum((Pr[a][k] * Dr[k][b] for k in range(d) if Dr[k][b]), R(0))
             for b in range(d)] for a in range(d)]
CB = [comm(D) for D in BASE]
p7 = 7
Kmat = sp.Matrix(((np.asarray(ctx['Kmat']) + p7 // 2) % p7 - p7 // 2).tolist())
Bm = sp.Matrix.hstack(Kmat, sp.Matrix([[1 if i == j else 0 for j in compl] for i in range(d)]))
Kl = Bm.inv()[:kdim, :]
KlQ = [[rc(Kl[i, j]) for j in range(d)] for i in range(kdim)]; KmQ = [[rc(Kmat[i, j]) for j in range(kdim)] for i in range(d)]
def xi(Cd, x):
    m = re.fullmatch(r'y(\d+)_(\d+)', x)
    if m: return Cd[int(m.group(1))][compl[int(m.group(2))]]
    m = re.fullmatch(r'k(\d)(\d)', x); i, j = int(m.group(1)), int(m.group(2))
    return sum((Cd[a][b] * KlQ[i][a] * KmQ[b][j] for a in range(d) if KlQ[i][a] for b in range(d) if KmQ[b][j] and Cd[a][b]), R(0))
XIB = [[xi(Cd, x) for x in vs] for Cd in CB]
ELIM = [x for x in S['names'] if x not in vs] if A.leafmode else []
XCE = [xi(Pr, x) for x in ELIM]                            # the eliminated coordinates as functions of vs
XFB = [[xi(Cd, x) for x in ELIM] for Cd in CB]             # the true field on them
ss = lambda e: str(e).replace('**', '^') if e else '0'
cols = [(a, b) for a in range(d) for b in range(d) if any(Cd[a][b] for Cd in CB)]
MB = [[Cd[a][b] for a, b in cols] for Cd in CB]
nterm = lambda e: len(e.terms()) if e else 0
# class Δ-defects (as in v2_55_pipeline.py): Σ_m P[m,l] Δ_c[m] − P Δ_c[l] Pᵀ
arrows = ctx['arrows']; basis = ctx['basis']; idx = {b_: i for i, b_ in enumerate(basis)}
cl = st['classes']
def DLclass(c):
    D_ = np.zeros((d, d, d), dtype=np.int64)
    for l, b_ in enumerate(basis):
        if b_[0] == 'v': continue
        for k, x in enumerate(b_[1]):
            if x not in c: continue
            pre = idx[('p', b_[1][:k])] if k > 0 else idx[('v', arrows[x][0])]
            suf = idx[('p', b_[1][k + 1:])] if k + 1 < len(b_[1]) else idx[('v', arrows[x][1])]
            D_[l, pre, suf] += 1
    return D_
def defect(Dl):
    out = []
    for l in range(d):
        PD = [[sum((Pr[a][x] * int(Dl[l][x, y]) for x in range(d) if Dl[l][x, y]), R(0)) for y in range(d)] for a in range(d)]
        for a in range(d):
            for b in range(d):
                lhs = sum((Pr[m][l] * int(Dl[m][a, b]) for m in range(d) if Dl[m][a, b]), R(0))
                out.append(lhs - sum((PD[a][y] * Pr[b][y] for y in range(d) if PD[a][y]), R(0)))
    return out
defects = [defect(DLclass(c)) for c in cl]
nzd = [k for k in range(len(defects[0])) if any(defects[c][k] for c in range(NCL))]
DFSTR = ','.join(ss(defects[c][k]) for c in range(NCL) for k in nzd) or '0'
XESTR = f"matrix XE[{NCL}][{nv}] = " + ",".join(ss(XIB[ND + c][b]) for c in range(NCL) for b in range(nv)) + ";"
def active(Wb):
    """coefficient rows over BASE of the active fields: Der_Δ, then Σ_c w_c E_c for w in Wb"""
    return [[1 if j == i else 0 for j in range(ND + NCL)] for i in range(ND)] + [[0] * ND + [rc(x) for x in w] for w in Wb]
def combine(rows, B):
    out = []
    for co in rows:
        out.append([sum((rc(sp.Rational(str(co[j]))) * B[j][t] for j in range(len(co)) if co[j] and B[j][t]), R(0)) for t in range(len(B[0]))])
    return out
gpos = {str(g): i for i, g in enumerate(R.gens)}
from sympy.polys.matrices import DomainMatrix as _DM
KD = R.to_domain()
print(f"{Q} V = {A.block} leaf {A.leaf} comps {A.comp}: {nv} coordinates, dim Der_Δ = {ND}, classes {NCL}, "
      f"{len(cols)} orbit-matrix columns; t = {time.time() - t0:.0f}s", flush=True)
from sympy.polys.matrices import DomainMatrix as _DM
KD = R.to_domain()
gpos = {str(g): i for i, g in enumerate(R.gens)}
def singular(lines, tag):
    fn = f"{A.out}_{tag}.sing"; open(fn, 'w').write("\n".join(lines))
    try: r = subprocess.run(["Singular", "-q", "--no-rc", fn], capture_output=True, text=True, timeout=A.tmo)
    except subprocess.TimeoutExpired: return None
    o = r.stdout + r.stderr
    return None if ("? error" in o or "END" not in o) else o
RING = [f'LIB "primdec.lib";', f"ring R = 0,({','.join(vs)}),dp;", "option(redSB);"]
def head(XI):
    return RING + [f"matrix XI[{len(XI)}][{nv}] = " + ",".join(ss(XI[a][b]) for a in range(len(XI)) for b in range(nv)) + ";",
                   f"proc xi(poly h, int a) {{ poly g = 0; int b; for (b = 1; b <= {nv}; b++) {{ g = g + diff(h, var(b)) * XI[a, b]; }} return(g); }}"]
def evalr(e, pt):
    out = sp.Rational(0)
    for mon, c in e.items():
        t = sp.Rational(int(c.numerator), int(c.denominator))
        for i, k in enumerate(mon):
            if k: t *= pt[i] ** k
        out += t
    return out
def evalc(e, pt):
    out = 0j
    for mon, c in e.items():
        t = complex(float(c))
        for i, k in enumerate(mon):
            if k: t *= pt[i] ** k
        out += t
    return out
def points(gens, iv, k=3):
    Sv = [sp.Symbol(v) for v in vs]; gs = [sp.sympify(g.replace('^', '**'), locals=loc) for g in gens]
    rat, num = [], []
    for t in range(16 * k):
        vals = {Sv[j]: sp.Integer(rng.randint(-60, 60)) for j in range(nv) if iv[j] == 1}
        g2 = [g for g in (sp.expand(g.subs(vals)) for g in gs) if g != 0]
        unk = sorted(set().union(*[g.free_symbols for g in g2]), key=str) if g2 else []
        try: sols = sp.solve(g2, unk, dict=True) if g2 else [{}]
        except Exception: sols = []
        for so in sols:
            if len(so) < len(unk): continue
            pt = dict(vals); pt.update(so)
            for v in Sv: pt.setdefault(v, sp.Integer(rng.randint(-60, 60)))
            full = [0] * len(R.gens)
            if all(x.is_rational for x in pt.values()):
                for s_, v in pt.items(): full[gpos[str(s_)]] = sp.Rational(v)
                rat.append(full)
            else:
                try:
                    for s_, v in pt.items(): full[gpos[str(s_)]] = complex(sp.N(v, 30))
                    num.append(full)
                except Exception: pass
            break
        if len(rat) >= k: break
    return rat, num
def choose_minor(gens, iv, e, Mrows):
    """an e × e minor of the orbit matrix not vanishing at a sampled point of V(gens) (None if none found)"""
    nr = len(Mrows)
    corder = sorted(range(len(cols)), key=lambda j: (sum(nterm(Mrows[i][j]) for i in range(nr)), j))
    rorder = sorted(range(nr), key=lambda i: (sum(nterm(x) for x in Mrows[i]), i))
    if e > nr: return None, None
    for attempt in range(4):                               # rank can drop at special points: sample again
        f_, how_ = _choose(gens, iv, e, Mrows, corder, rorder)
        if f_ is not None: return f_, how_
    return None, None
def _choose(gens, iv, e, Mrows, corder, rorder):
    rat, num = points(gens, iv)
    for pt in rat:
        Mv = sp.Matrix([[evalr(x, pt) for x in row] for row in Mrows])
        if Mv.rank() < e: continue
        _, pc = Mv[:, corder].rref(); pc = [corder[j] for j in list(pc)[:e]]
        _, pr = Mv[rorder, pc].T.rref(); pr = [rorder[i] for i in list(pr)[:e]]
        return _DM([[Mrows[i][j] for j in pc] for i in pr], (e, e), KD).det(), 'rational'
    for pt in num:
        Mn = np.array([[evalc(x, pt) for x in row] for row in Mrows])[np.ix_(rorder, corder)]
        tol = 1e-7 * max(1.0, np.abs(Mn).max()); rk = lambda B: np.linalg.matrix_rank(B, tol=tol) if B.size else 0
        if rk(Mn) < e: continue
        pc = []
        for j in range(Mn.shape[1]):
            if rk(Mn[:, pc + [j]]) > len(pc): pc.append(j)
            if len(pc) == e: break
        pr = []
        for i in range(Mn.shape[0]):
            if rk(Mn[pr + [i]][:, pc]) > len(pr): pr.append(i)
            if len(pr) == e: break
        pr, pc = [rorder[i] for i in pr], [corder[j] for j in pc]
        return _DM([[Mrows[i][j] for j in pc] for i in pr], (e, e), KD).det(), 'numeric'
    return None, None
def info(gens, tag):
    o = singular(RING + [f"ideal G = std(ideal({','.join(gens) or '0'}));", 'print("DIM " + string(dim(G)));',
                         'print("VDIM " + string(vdim(G)));', 'intvec iv = indepSet(G); print("INDEP " + string(iv));', 'print("END"); quit;'], tag)
    return (int(re.search(r"DIM (-?\d+)", o).group(1)), [int(x) for x in re.search(r"INDEP (.*)", o).group(1).split(',')],
            int(re.search(r"VDIM (-?\d+)", o).group(1)))
def local_weights(gens, tag):
    """basis of {w : Σ_c w_c defect_c ∈ I_Y and Σ_c w_c ξ_{E_c}(I_Y) ⊆ I_Y}"""
    lines = RING + [XESTR, f"ideal I = {','.join(gens) or '0'}; ideal G = std(I);", f"ideal DF = {DFSTR};",
                    'ideal NFD = reduce(DF, G); int c; int i; int b; poly g;',
                    'for (c = 1; c <= ncols(DF); c++) { print("D " + string(c) + " " + string(NFD[c])); }',
                    f'for (c = 1; c <= {NCL}; c++) {{ for (i = 1; i <= size(I); i++) {{ g = 0; for (b = 1; b <= {nv}; b++) '
                    '{ g = g + diff(I[i], var(b)) * XE[c, b]; } print("T " + string(c) + " " + string(i) + " " + string(reduce(g, G))); } }',
                    'print("END"); quit;']
    o = singular(lines, tag)
    nk = len(nzd); eqs = {}
    def add(key, c, poly):
        if poly.strip() == '0': return
        for mon, co in parse(poly).terms():
            eqs.setdefault((key, mon), [sp.Rational(0)] * NCL)[c] += sp.Rational(int(co.numerator), int(co.denominator))
    for m in re.finditer(r"^D (\d+) (.*)$", o, re.M) if nk else []:
        j = int(m.group(1)) - 1; add(('d', j % nk), j // nk, m.group(2))
    for m in re.finditer(r"^T (\d+) (\d+) (.*)$", o, re.M):
        add(('t', int(m.group(2))), int(m.group(1)) - 1, m.group(3))
    rows = list(eqs.values())
    drows = [v for k, v in eqs.items() if k[0][0] == 'd']
    LASTW['defect_only'] = NCL - (sp.Matrix(drows).rank() if drows else 0)
    if not rows: return [[sp.Integer(1) if i == j else sp.Integer(0) for i in range(NCL)] for j in range(NCL)]
    Nm = sp.Matrix(rows).nullspace()
    return [[v[i] for i in range(NCL)] for v in Nm]
LASTW = {}

stats = {'primes': 0, 'dense': 0, 'points': 0, 'nodense': 0, 'unresolved': 0}
NODES = {}                                                 # canonical reduced Gröbner basis → node (primes are shared)
def stratify(gens, path, check_tangent=False):
    key = ','.join(gens)
    if key in NODES: return key
    n_ = stats['primes'] = stats['primes'] + 1
    dim_, iv, vd = info(gens, f'i{n_}')
    node = NODES[key] = dict(path=path, dim=dim_, gens=gens)
    if dim_ <= 0:
        stats['points'] += 1; node.update(result='points', geometric=vd); return key
    Wb = local_weights(gens, f'w{n_}')
    act = active(Wb); XI = combine(act, XIB); Mrows = combine(act, MB)
    node['weights'] = [[str(x) for x in w] for w in Wb]; node['weights_defect_only'] = LASTW.get('defect_only')
    f, how = choose_minor(gens, iv, dim_, Mrows)
    if f is None:
        stats['nodense'] += 1; node['result'] = 'no dense orbit found at sampled points'
        print(f"  {path}: dim {dim_}: NO DENSE ORBIT at sampled points (dim 𝔥_Y = {len(act)})", flush=True); return key
    NH = len(XI)
    lines = head(XI) + [f"ideal I = {','.join(gens) or '0'};", "ideal G = std(I); int a; int i; poly g;",
                        f"poly f = {ss(f)};", 'if (reduce(f, G) == 0) { print("FININ"); print("END"); quit; }']
    if check_tangent:
        lines += ['for (i = 1; i <= size(I); i++) { for (a = 1; a <= ' + str(NH) + '; a++) { if (reduce(xi(I[i], a), G) != 0) '
                  '{ print("NOTTANGENT " + string(i) + " " + string(a)); print("END"); quit; } } }']
        if ELIM:   # a leaf's graph need not be H-stable: the field on the eliminated coordinates must match too
            XF = combine(act, XFB)
            lines += [f"matrix XC[1][{len(ELIM)}] = " + ",".join(ss(e) for e in XCE) + ";",
                      f"matrix XF[{NH}][{len(ELIM)}] = " + ",".join(ss(XF[a][e]) for a in range(NH) for e in range(len(ELIM))) + ";",
                      'int e; for (a = 1; a <= ' + str(NH) + '; a++) { for (e = 1; e <= ' + str(len(ELIM)) + '; e++) { '
                      'if (reduce(XF[a, e] - xi(XC[1, e], a), G) != 0) { print("NOTTANGENT graph " + string(e) + " " + string(a)); print("END"); quit; } } }']
    lines += ["G = std(G + f); ideal NEW = f; int rounds = 0;",
              "while (size(NEW) > 0) { rounds++; ideal NX = 0; for (i = 1; i <= size(NEW); i++) { for (a = 1; a <= " + str(NH) +
              "; a++) { g = reduce(xi(NEW[i], a), G); if (g != 0) { NX = NX + g; G = std(G + g); } } } NEW = simplify(NX, 2); kill NX; }",
              'print("CLOSED " + string(rounds) + " DIM " + string(dim(G)));',
              "list L = minAssGTZ(G);",
              'for (i = 1; i <= size(L); i++) { ideal Gi = std(L[i]); print("COMP " + string(dim(Gi))); print("GENS " + string(Gi)); kill Gi; }',
              'print("END"); quit;']
    o = singular(lines, f'n{n_}')
    if o is None:
        stats['unresolved'] += 1; node['result'] = 'singular timeout/error'
        print(f"  {path}: dim {dim_}: SINGULAR timeout/error; t = {time.time() - t0:.0f}s", flush=True); return key
    assert "FININ" not in o, "chosen minor lies in I_X"
    if "NOTTANGENT" in o and A.leafmode:                   # a leaf prime that is not H-stable is not a component
        stats['primes'] -= 1; stats['dense'] -= 0; node['result'] = 'not H-stable (not a component of the stratum)'
        print(f"  {path}: not H-stable: skipped", flush=True); return key
    assert "NOTTANGENT" not in o, "vector fields not tangent to the input prime: " + o[o.index("NOTTANGENT"):][:40]
    m = re.search(r"CLOSED (\d+) DIM (-?\d+)", o)
    kids = [(int(b.split()[0]), b.split("GENS ")[1].split("\n")[0].split(',')) for b in o.split("COMP ")[1:]]
    kids = [k for k in kids if k[0] >= 0]                 # the unit ideal (empty boundary) has dim -1
    assert all(k[0] < dim_ for k in kids), "boundary not of lower dimension"
    stats['dense'] += 1
    node.update(result='dense orbit', minor=ss(f), minor_point=how, geometric=1 if how == 'rational' else None,
                closure_rounds=int(m.group(1)), boundary_dim=int(m.group(2)), h_dim=NH)
    print(f"  {path}: dim {dim_}: dense orbit (dim 𝔥_Y {NH}, minor {len(ss(f))} chars, {how} point); boundary dim {m.group(2)} after "
          f"{m.group(1)} closure rounds, {len(kids)} primes {[k[0] for k in kids]}; t = {time.time() - t0:.0f}s", flush=True)
    node['children'] = [stratify(g, f"{path}.{k + 1}") for k, (dm, g) in enumerate(kids)]
    return key

roots = []
for c in comps:
    gens = [g.replace('**', '^') for g in c['gens'] if g.strip() not in ('', '0')]
    if c['dim'] <= 0 and not gens:                          # a point stratum: one fold
        roots.append(None); stats['points'] += 1; stats['primes'] += 1; continue
    roots.append(stratify(gens, f"X{c['comp']}", check_tangent=True))
    nd = NODES[roots[-1]]
    assert nd['dim'] == c['dim'], (nd['dim'], c['dim'])
    if c.get('W') is not None and nd.get('weights') is not None and len(nd['weights']) != len(c['W']):
        print(f"  note: X{c['comp']} local weightings {len(nd['weights'])} vs pipeline W_Z {len(c['W'])}", flush=True)
for k_, nd in NODES.items():                               # dense orbit seen only at algebraic points: count the
    if nd.get('result') == 'dense orbit' and nd.get('geometric') is None:   # geometric components (absPrimdecGTZ)
        o = singular(RING + [f"ideal I = {k_};", "def SA = absPrimdecGTZ(I); setring SA; int i;",
                             'for (i = 1; i <= size(absolute_primes); i++) { print("CONJ " + string(absolute_primes[i][2])); }',
                             'print("END"); quit;'], 'abs' + str(abs(hash(k_)) % 10**8))
        if o: nd['geometric'] = sum(int(x) for x in re.findall(r"CONJ (\d+)", o)); nd['geometric_by'] = 'absPrimdecGTZ'
if A.leafmode:                                             # canonical key in all coordinates (graph of the leaf's substitutions)
    rel = [f"{x}-({v.replace('**', '^')})" for x, v in st['subs']]
    for k_, nd in NODES.items():
        if nd.get('result') not in ('dense orbit', 'points'): continue
        o = singular([f"ring RF = 0,({','.join(S['names'])}),dp;", "option(redSB);", f"ideal I = {','.join(([k_] if k_ else []) + rel) or '0'};",
                      'ideal G = std(I); print("FULL " + string(G));', 'print("END"); quit;'], 'full' + str(len(k_)) + '_' + str(abs(hash(k_)) % 10**8))
        nd['fullkey'] = re.search(r"FULL (.*)", o).group(1)
NODES_OK = {k_: nd for k_, nd in NODES.items() if not str(nd.get('result', '')).startswith('not H-stable')}
ok = stats['nodense'] == 0 and stats['unresolved'] == 0
geo = [n.get('geometric') for n in NODES_OK.values()]
count = sum(geo) + roots.count(None) if all(g is not None for g in geo) else None
json.dump(dict(block=A.block, leaf=A.leaf, comps=[c['comp'] for c in comps], quiver=Q, roots=roots, nodes=NODES,
               stats=stats, orbits=count), open(A.out + '.json', 'w'))
print(f"RESULT {'finitely many orbits' if ok else 'NOT ESTABLISHED'}: {stats}; geometric H-orbits {count}; t = {time.time() - t0:.0f}s", flush=True)
