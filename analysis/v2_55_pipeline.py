"""v2.55: the v2.53–v2.54 pipeline for any quiver, built on fold_ideal.py (fast sparse polynomials).

For a vertex block V of QUIVER:
1. the stratum's fold ideal, presolved (fold_ideal.stratum, fold_ideal.presolve);
2. its components over ℚ (Singular minAssGTZ), with dimension, independent set and the normal forms of the class
   Δ-defects modulo each prime;
3. per component Z: W_Z (weightings compatible on all of Z), 𝔥_Z = Der_Δ ⊕ E_{W_Z};
4. the orbit-dimension criterion (v2.54): for r = dim Z − 1, …, 0, an upper bound dim Z_r ≤ r by point-guided
   (r+1)-minors of M(P) = ([h, P])_{h ∈ 𝔥_Z}. r = dim Z − 1 failing with "rank ≤ r at every sampled point" means no
   open orbit; lower r failing means a boundary locus of dimension > r with orbits of dimension ≤ r (a modulus).
Environment (v2.57): FOLD_STAGED (default 1) builds the closure equations in batches of FOLD_BATCH column pairs;
PRESOLVE_RULE=markowitz picks substitutions by fill; BRANCH=1 uses fold_ideal.branch_presolve (case splits on factoring
equations, each leaf in a forked child, leaf ideals saturated by the leaf's inequations); CLOSURE_IN_SINGULAR=1 builds
the closure equations inside Singular.
Usage: python3 analysis/v2_55_pipeline.py QUIVER BLOCK OUT_PREFIX [--seed S]   (BLOCK: rows of V, e.g. 100,000,001; m for −1, e.g. 100m,0000,001m,0000)"""
import os, sys, re, random, subprocess, argparse, collections, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import numpy as np, sympy as sp
from sympy.polys.matrices import DomainMatrix
import fold_ideal as fi
import singular_nohilb; singular_nohilb.install()   # minAssGTZ without the overflowing Hilbert series
ap = argparse.ArgumentParser(); ap.add_argument('quiver'); ap.add_argument('block'); ap.add_argument('out')
ap.add_argument('--seed', type=int, default=55)
A = ap.parse_args(); rng = random.Random(A.seed); t0 = time.time()
V = [[-1 if ch == 'm' else int(ch) for ch in r] for r in A.block.split(',')]   # 'm' = −1 (v2.57)
S = fi.stratum(A.quiver, V)
R, P, ctx = S['R'], S['P'], S['ctx']
d, basis, arrows = ctx['d'], ctx['basis'], ctx['arrows']
M = np.asarray(ctx['M']); DL = np.asarray(ctx['DL'])
idx = {b: i for i, b in enumerate(basis)}
CLOSURE_IN_SINGULAR = os.environ.get('CLOSURE_IN_SINGULAR') == '1'
NZ = []                                                         # leaf inequations (BRANCH=1), saturated out below
LEAF = 0
import json
CERT = os.environ.get('CERT')                                   # v2.57: JSON-lines certificate file (appended)
def cert(rec):
    if not CERT: return
    rec = dict(quiver=A.quiver, block=A.block, leaf=LEAF, **rec)
    with open(CERT, 'a') as fh: fh.write(json.dumps(rec) + "\n")
def _leaf_parse(s_):
    gl_ = {str(g): g for g in R.gens}
    return R(sp.sympify(s_.replace('^', '**'), locals={k: sp.Symbol(k) for k in gl_}))
if os.environ.get('LEAF_LOAD'):
    # v2.66: one leaf of a case split written by LEAF_DUMP, run in a fresh process (forked leaf children of a large
    # parent were killed for memory); LEAVES must name exactly one leaf
    _D = json.load(open(os.environ['LEAF_LOAD'])); LEAF = int(os.environ['LEAVES'])
    _rem, _subs, _nz = _D['leaves'][LEAF - 1]
    rem = [_leaf_parse(e) for e in _rem]; subs = [(R.gens[[str(g) for g in R.gens].index(x)], _leaf_parse(v)) for x, v in _subs]
    NZ = [_leaf_parse(e) for e in _nz]; A.out += f'_L{LEAF}'
    fi.PIVOTS[:] = [tuple(p) for p in _D['pivots']]; fi.BRANCH_CONSTS |= {tuple(c) for c in _D['consts']}; fi.BRANCH_TAGS |= {tuple(c) for c in _D['tags']}
    print(f"leaf {LEAF}/{len(_D['leaves'])} ({len(NZ)} inequations, loaded):", flush=True)
elif os.environ.get('BRANCH') == '1':
    # v2.57: case-split presolve (fold_ideal.branch_presolve). The stratum is the union of the leaf varieties; each leaf
    # runs the rest of this script in a forked child, so every irreducible component of the stratum is among the
    # (saturated) leaf primes checked below. Leaf primes need not be components (one may lie in another leaf's).
    leaves = fi.branch_presolve(S)
    print(f"{A.quiver} V = {A.block}: {len(leaves)} leaves (case-split presolve); t = {time.time() - t0:.0f}s", flush=True)
    if os.environ.get('LEAF_DUMP'):
        json.dump(dict(leaves=[[[str(e) for e in r_], [[str(x), str(v)] for x, v in s_], [str(f) for f in n_]] for r_, s_, n_ in leaves],
                       pivots=[list(p) for p in fi.PIVOTS], consts=sorted(fi.BRANCH_CONSTS), tags=sorted(fi.BRANCH_TAGS)),
                  open(os.environ['LEAF_DUMP'], 'w'))
        print("leaves written to", os.environ['LEAF_DUMP'], flush=True); sys.exit(0)
    failed = []
    fi._FCACHE.clear(); import gc; gc.collect()       # v2.66: smaller parent before forking (children were OOM-killed)
    ONLY_LEAVES = {int(x) for x in os.environ.get('LEAVES', '').split(',') if x.strip()}   # v2.66: rerun these leaves only
    for k_, (rem_, subs_, nz_) in enumerate(leaves):
        if ONLY_LEAVES and k_ + 1 not in ONLY_LEAVES: continue
        sys.stdout.flush(); pid = os.fork()
        if pid == 0:
            rem, subs, NZ = rem_, subs_, nz_; A.out += f'_L{k_ + 1}'; LEAF = k_ + 1
            print(f"leaf {k_ + 1}/{len(leaves)} ({len(nz_)} inequations):", flush=True); break
        _, st_ = os.waitpid(pid, 0)
        if st_ != 0: failed.append((k_ + 1, st_))
    else:
        print(f"{A.quiver} V = {A.block}: all leaves done; " + (f"FAILED leaves {failed}" if failed else "every leaf exited cleanly"), flush=True)
        sys.exit(1 if failed else 0)
elif CLOSURE_IN_SINGULAR: rem, subs = fi.stage1_presolve(S)    # closure equations are generated by Singular below
elif S['eqs'] is None: rem, subs = fi.staged_presolve(S)        # v2.57: staged build + indexed presolve (default)
else: rem, subs = (fi.presolve if os.environ.get('PRESOLVE') == 'v255' else fi.presolve_fast)(S['eqs'], R)
elim = {str(x) for x, _ in subs}
vs = [n for n in S['names'] if n not in elim]                    # coordinates of the presolved stratum
if not vs:                       # every coordinate eliminated: the stratum is one point, or empty if a constant ≠ 0 remains
    empty = any(e != 0 for e in rem)
    print(f"{A.quiver} V = {A.block}: presolved to {'no point (inconsistent)' if empty else 'a single point'}; "
          f"{0 if empty else 1} components", flush=True)
    if not empty: print("  component 1: dim 0 (points)", flush=True)
    cert(dict(kind='stratum', vs=[], subs=[[str(x), str(v)] for x, v in subs], rem=[str(e) for e in rem],
              nz=[str(f) for f in NZ], pivots=fi.PIVOTS, branch_consts=sorted(fi.BRANCH_CONSTS), branch_tags=sorted(fi.BRANCH_TAGS), ncomp=0 if empty else 1))
    if not empty: cert(dict(kind='component', comp=1, dim=0, gens=[], indep=[], W=None, cuts=[], result='point'))
    sys.exit(0)
Pr = [[e for e in row] for row in P]
for x, val in subs: Pr = [[e.compose(x, val) for e in row] for row in Pr]
# class defects (ring arithmetic)
cls = collections.defaultdict(list)
for a_, st in enumerate(arrows): cls[st].append(a_)
cl = list(cls.values())
def DLclass(c):
    D_ = np.zeros((d, d, d), dtype=np.int64)
    for l, b in enumerate(basis):
        if b[0] == 'v': continue
        for k, x in enumerate(b[1]):
            if x not in c: continue
            pre = idx[('p', b[1][:k])] if k > 0 else idx[('v', arrows[x][0])]
            suf = idx[('p', b[1][k + 1:])] if k + 1 < len(b[1]) else idx[('v', arrows[x][1])]
            D_[l, pre, suf] += 1
    return D_
DLc = [DLclass(c) for c in cl]
def defect(Dl):
    out = []
    for l in range(d):
        PD = [[sum((Pr[a][x] * int(Dl[l][x, y]) for x in range(d) if Dl[l][x, y]), R(0)) for y in range(d)] for a in range(d)]
        for a in range(d):
            for b in range(d):
                lhs = sum((Pr[m][l] * int(Dl[m][a, b]) for m in range(d) if Dl[m][a, b]), R(0))
                out.append(lhs - sum((PD[a][y] * Pr[b][y] for y in range(d) if PD[a][y]), R(0)))
    return out
defects = [defect(D_) for D_ in DLc]
nz = [k for k in range(len(defects[0])) if any(defects[c][k] for c in range(len(cl)))]
ss = lambda f: str(f).replace('**', '^')
def closure_singular():
    """Singular code appending the closure equations of Pr (image and kernel closed under products) to I."""
    if not CLOSURE_IN_SINGULAR: return ""
    nzM = [(i + 1, j + 1, k + 1, int(M[i, j, k])) for i in range(d) for j in range(d) for k in range(d) if M[i, j, k]]
    lines = [f"matrix PP[{d}][{d}] = " + ",".join(ss(Pr[r][c]) if Pr[r][c] else "0" for r in range(d) for c in range(d)) + ";",
             f"matrix QQm[{d}][{d}] = unitmat({d}) - PP;",
             f"intmat MI[{len(nzM)}][4] = " + ",".join(f"{a},{b},{k},{c}" for a, b, k, c in nzM) + ";",
             "ideal CL; int ca; int cb; int ct;",
             f"for (ca = 1; ca <= {d}; ca++) {{ for (cb = 1; cb <= {d}; cb++) {{",
             f"  matrix u[{d}][1]; matrix w[{d}][1];",
             f"  for (ct = 1; ct <= {len(nzM)}; ct++) {{",
             "    u[MI[ct,3],1] = u[MI[ct,3],1] + MI[ct,4] * PP[MI[ct,1],ca] * PP[MI[ct,2],cb];",
             "    w[MI[ct,3],1] = w[MI[ct,3],1] + MI[ct,4] * QQm[MI[ct,1],ca] * QQm[MI[ct,2],cb]; }",
             "  CL = CL + ideal(PP * u - u) + ideal(PP * w); kill u; kill w; } }",
             "I = simplify(I + CL, 2 + 4 + 8);"]
    return "\n".join(lines)
src = f"""LIB "primdec.lib";
ring r = 0,({','.join(vs)}),dp;
ideal I = {','.join(ss(e) for e in rem) if rem else '0'};
{closure_singular()}
{''.join(f"I = sat(I, {ss(f)});" + chr(10) for f in NZ)}ideal SI = std(I); if (size(SI) > 0) {{ if (deg(SI[1]) == 0) {{ print("END"); quit; }} }} kill SI;
list L = minAssGTZ(I);
ideal DF = {','.join(ss(defects[c][k]) for c in range(len(cl)) for k in nz) if nz else '0'};
int i;
for (i = 1; i <= size(L); i++) {{
  ideal Gi = std(L[i]);
  print("PRIME " + string(i) + " DIM " + string(dim(Gi)));
  print("GENS " + string(Gi));
  intvec iv = indepSet(Gi); print("INDEP " + string(iv));
  ideal NFd = reduce(DF, Gi); print("NFD " + string(NFd));
  kill iv; kill NFd; kill Gi;
}}
print("END");
quit;
"""
open(A.out + '_base.sing', 'w').write(src)
res = subprocess.run(["Singular", "-q", "--no-rc", A.out + '_base.sing'], capture_output=True, text=True, timeout=20000)
out = res.stdout
assert "END" in out and "? error" not in out + res.stderr and "overflow" not in out + res.stderr, \
    "Singular decomposition failed: " + (out + res.stderr)[-800:]
loc = {v: sp.Symbol(v) for v in vs}; Sv = [sp.Symbol(v) for v in vs]
parse = lambda line: [sp.sympify(g.replace('^', '**'), locals=loc) for g in line.split(',') if g.strip()]
comps = []
for blk in out.split("PRIME ")[1:]:
    m_ = re.match(r"(\d+) DIM (-?\d+)", blk)
    comps.append((int(m_.group(1)), int(m_.group(2)), parse(re.search(r"GENS (.*)", blk).group(1)),
                  [int(x) for x in re.search(r"INDEP (.*)", blk).group(1).split(',')], parse(re.search(r"NFD (.*)", blk).group(1))))
# Der, Der_Δ, class torus derivations
QQ = sp.QQ
def nullsp(rows, n):
    if not rows: return [sp.Matrix([1 if i == j else 0 for i in range(n)]) for j in range(n)]
    Nm = DomainMatrix.from_list_sympy(len(rows), n, rows).convert_to(QQ).nullspace().to_Matrix()
    return [Nm.row(k).T for k in range(Nm.rows)]
rows = []
for i in range(d):
    for j in range(d):
        for m in range(d):
            row = [0] * (d * d)
            for k in range(d):
                if M[i, j, k]: row[m * d + k] += int(M[i, j, k])
            for r_ in range(d):
                if M[r_, j, m]: row[r_ * d + i] -= int(M[r_, j, m])
                if M[i, r_, m]: row[r_ * d + j] -= int(M[i, r_, m])
            if any(row): rows.append(row)
Der = [sp.Matrix(d, d, list(v)) for v in nullsp(rows, d * d)]
DLs = [sp.Matrix(DL[l].tolist()) for l in range(d)]
cond = [sum((list(sum((D[m, l] * DLs[m] for m in range(d)), sp.zeros(d, d)) - D * DLs[l] - DLs[l] * D.T) for l in range(d)), [])
        for D in Der]
cx = nullsp([list(r) for r in sp.Matrix(cond).T.tolist()], len(Der))
DerD = [sum((c[k] * Der[k] for k in range(len(Der))), sp.zeros(d, d)) for c in cx]
Ec = [sp.diag(*[0 if b[0] == 'v' else sum(1 for x in b[1] if x in c) for b in basis]) for c in cl]
PrE = None   # v2.57: P stays in ring arithmetic (sympy expressions of the presolved P ran out of memory on d = 18)
from sympy.polys.matrices import DomainMatrix as _DM
KD = R.to_domain()
gpos = {str(g): i for i, g in enumerate(R.gens)}
def ringify(D):                                        # rational d×d matrix → ring constants
    return [[R(sp.QQ(int(sp.fraction(D[a, b])[0]), int(sp.fraction(D[a, b])[1]))) if D[a, b] != 0 else R(0)
             for b in range(d)] for a in range(d)]
def commutator(D):                                     # [D, P] with P = Pr in the ring
    Dr = ringify(D)
    return [[sum((Dr[a][k] * Pr[k][b] for k in range(d) if Dr[a][k]), R(0)) -
             sum((Pr[a][k] * Dr[k][b] for k in range(d) if Dr[k][b]), R(0)) for b in range(d)] for a in range(d)]
def evaluate(f, pt):                                   # ring element at a point {Symbol: Rational} on the coordinates vs
    vals = [0] * len(R.gens)
    for sym, v in pt.items(): vals[gpos[str(sym)]] = sp.QQ(int(sp.fraction(v)[0]), int(sp.fraction(v)[1]))
    out = sp.QQ(0)                                     # direct sum over terms (f(*vals) builds a ring per variable)
    for mon, c in f.items():
        t = c
        for i, e in enumerate(mon):
            if e: t *= vals[i] ** e
        out += t
    return out
print(f"{A.quiver} V = {A.block}: presolved to {len(rem)} equations in {len(vs)} coordinates; {len(comps)} components; "
      f"dim Der = {len(Der)}, dim Der_Δ = {len(DerD)}, classes {len(cl)}; t = {time.time() - t0:.0f}s", flush=True)
def singular(lines, fn, tmo=int(os.environ.get('SING_TMO', '1800'))):   # v2.66: SING_TMO
    open(fn, 'w').write("\n".join(lines))
    try: res = subprocess.run(["Singular", "-q", "--no-rc", fn], capture_output=True, text=True, timeout=tmo)
    except subprocess.TimeoutExpired: return ""
    return "" if ("? error" in res.stdout + res.stderr or "overflow" in res.stdout + res.stderr) else res.stdout
def points(gl, iv, k):
    pts = []
    for t in range(12 * k):
        vals = {Sv[j]: sp.Integer(rng.randint(-5, 5)) for j in range(len(vs)) if iv[j] == 1}
        g2 = [g for g in (sp.expand(g.subs(vals)) for g in gl) if g != 0]
        unk = sorted(set().union(*[g.free_symbols for g in g2]), key=str) if g2 else []
        try: sols = sp.solve(g2, unk, dict=True) if g2 else [{}]
        except Exception: sols = []
        for so in sols:
            if len(so) < len(unk) or not all(x.is_rational for x in so.values()): continue
            pt = dict(vals); pt.update(so)
            for v in Sv: pt.setdefault(v, sp.Integer(rng.randint(-5, 5)))
            pts.append(pt); break
        if len(pts) >= k: break
    return pts
def points_numeric(gl, iv, k):                     # v2.66: algebraic points (complex floats) when no rational one is found;
    pts = []                                          # they only choose which minor to add, the cut itself is exact
    for t in range(12 * k):
        vals = {Sv[j]: sp.Integer(rng.randint(-5, 5)) for j in range(len(vs)) if iv[j] == 1}
        g2 = [g for g in (sp.expand(g.subs(vals)) for g in gl) if g != 0]
        unk = sorted(set().union(*[g.free_symbols for g in g2]), key=str) if g2 else []
        try: sols = sp.solve(g2, unk, dict=True) if g2 else [{}]
        except Exception: sols = []
        for so in sols:
            if len(so) < len(unk): continue
            try: pt = {x: complex(sp.N(v, 30)) for x, v in so.items()}
            except Exception: continue
            pt.update({x: complex(int(v)) for x, v in vals.items()})
            for v in Sv: pt.setdefault(v, complex(rng.randint(-5, 5)))
            pts.append(pt); break
        if len(pts) >= k: break
    return pts
def evaluate_num(f, pt):
    vals = [0j] * len(R.gens)
    for sym, v in pt.items(): vals[gpos[str(sym)]] = v
    out = 0j
    for mon, c in f.items():
        t = complex(float(c))
        for i, e in enumerate(mon):
            if e: t *= vals[i] ** e
        out += t
    return out
def numeric_pivots(Mn, r):                            # greedy column then row selection of an (r+1)-minor, numerically
    import numpy as _np
    tol = 1e-7 * max(1.0, _np.abs(Mn).max())
    rk = lambda A: _np.linalg.matrix_rank(A, tol=tol) if A.size else 0
    if rk(Mn) <= r: return None
    pc = []
    for j in range(Mn.shape[1]):
        if rk(Mn[:, pc + [j]]) > len(pc): pc.append(j)
        if len(pc) == r + 1: break
    pr = []
    for i in range(Mn.shape[0]):
        if rk(Mn[pr + [i]][:, pc]) > len(pr): pr.append(i)
        if len(pr) == r + 1: break
    return (pr, pc) if len(pr) == len(pc) == r + 1 else None
cert(dict(kind='stratum', vs=vs, subs=[[str(x), str(v)] for x, v in subs], rem=[str(e) for e in rem],
          nz=[str(f) for f in NZ], pivots=fi.PIVOTS, branch_consts=sorted(fi.BRANCH_CONSTS), branch_tags=sorted(fi.BRANCH_TAGS), ncomp=len(comps), DerD=[[str(x) for x in D] for D in DerD],
          Ec=[[int(x) for x in E.diagonal()] for E in Ec], classes=cl))
ONLY = {int(x) for x in os.environ.get('COMPS', '').split(',') if x.strip()}   # v2.66: run only these components
for i_, dim_, gens, indep, nf in comps:
    if ONLY and i_ not in ONLY: continue
    n_ = len(nz); wv = sp.symbols(f'w0:{len(cl)}')
    lin = []
    for k in range(n_):
        e = sp.expand(sum(wv[c] * nf[c * n_ + k] for c in range(len(cl))))
        if e != 0:
            for co in (sp.Poly(e, *Sv).coeffs() if Sv else [e]): lin.append([sp.expand(co).coeff(w) for w in wv])
    W = nullsp(lin, len(cl))
    h = DerD + [sum((w[c] * Ec[c] for c in range(len(cl))), sp.zeros(d, d)) for w in W]
    m = len(h)
    Mrows_all = [[e for row in commutator(D) for e in row] for D in h]
    cols = [j for j in range(d * d) if any(Mrows_all[i][j] for i in range(m))]
    Mrows = [[Mrows_all[i][j] for j in cols] for i in range(m)]
    # v2.66: prefer sparse minors: pivot columns and rows are chosen greedily in order of increasing term count
    nterm = lambda e: len(e.terms()) if e else 0
    corder = sorted(range(len(cols)), key=lambda j: (sum(nterm(Mrows[i][j]) for i in range(m)), j))
    rorder = sorted(range(m), key=lambda i: (sum(nterm(e) for e in Mrows[i]), i))
    if dim_ <= 0:
        print(f"  component {i_}: dim {dim_} (points)", flush=True)
        cert(dict(kind='component', comp=i_, dim=dim_, gens=[str(g) for g in gens], indep=indep,
                  W=[[str(x) for x in w] for w in W], cuts=[], result='point')); continue
    cuts = []
    res_r = []
    for r in range(dim_ - 1, -1, -1):
        F = []; FP = []; status = None
        for it in range(80):
            # v2.66: eliminate the variables that occur linearly first (presolve.lib elimpart: V(J) is the graph of the
            # substitution SB over V(J1), J1 in the remaining variables), then check the dimension before decomposing
            out = singular([f'LIB "primdec.lib"; LIB "presolve.lib";', f"ring R = 0,({','.join(vs)}),dp;",
                            f"ideal J = {','.join([str(g).replace('**', '^') for g in gens] + [str(f).replace('**', '^') for f in F]) or '0'};",
                            "list EP = elimpart(J); ideal J1 = EP[1]; ideal SB = EP[3]; int ne = size(EP[2]);",
                            f"ideal Js = std(J1); if (dim(Js) - ne <= {r}) {{ print(\"END\"); quit; }}",
                            "list L = minAssGTZ(Js); int i;",
                            f"for (i = 1; i <= size(L); i++) {{ ideal Gi = std(L[i] + SB); if (dim(Gi) > {r}) {{ print(\"COMP \" + string(dim(Gi))); "
                            "print(\"GENS \" + string(Gi)); intvec iv = indepSet(Gi); print(\"INDEP \" + string(iv)); kill iv; } kill Gi; }",
                            'print("END");', "quit;"], f"{A.out}_c{i_}.sing")
            if "END" not in out: status = 'singular timeout'; break
            above = [(int(b.split()[0]), parse(b.split("GENS ")[1].split("\n")[0]), [int(x) for x in b.split("INDEP ")[1].split("\n")[0].split(',')])
                     for b in out.split("COMP ")[1:]]
            if not above: status = 'ok'; break
            added = False
            for dm, gl, iv in above:
                pts = points(gl, iv, 4)
                if not pts:
                    import numpy as _np
                    for pt in points_numeric(gl, iv, 4):
                        Mn = _np.array([[evaluate_num(x, pt) for x in row] for row in Mrows])
                        sel = numeric_pivots(Mn[_np.ix_(rorder, corder)], r)
                        if sel is None: continue
                        pr, pc = [rorder[i] for i in sel[0]], [corder[j] for j in sel[1]]
                        F.append(_DM([[Mrows[i][j] for j in pc] for i in pr], (len(pr), len(pc)), KD).det())
                        FP.append(dict(rows=[int(i) for i in pr], cols=[int(cols[j]) for j in pc], det=str(F[-1]), numeric=True))
                        added = True; break
                    if not added: status = f'no rational point on a component of dim {dm}'; break
                    continue
                for pt in pts:
                    Mv = sp.Matrix([[sp.Rational(int(sp.QQ.numer(evaluate(x, pt))), int(sp.QQ.denom(evaluate(x, pt))))
                                     for x in row] for row in Mrows])
                    if Mv.rank() <= r: continue
                    _, pc = Mv[:, corder].rref(); pc = [corder[j] for j in list(pc)[:r + 1]]
                    _, pr = Mv[rorder, pc].T.rref(); pr = [rorder[i] for i in list(pr)[:r + 1]]
                    F.append(_DM([[Mrows[i][j] for j in pc] for i in pr], (len(pr), len(pc)), KD).det())
                    FP.append(dict(rows=[int(i) for i in pr], cols=[int(cols[j]) for j in pc], det=str(F[-1])))
                    added = True; break
                if not added and not status: status = f'rank ≤ {r} at every sampled point of a component of dim {dm}'
                if status: break
            if status: break
        res_r.append((r, 'ok' if status == 'ok' else status)); cuts.append(dict(r=r, status=status, minors=FP))
        if status != 'ok': break                    # the first failure is what matters
    ok = all(s_ == 'ok' for _, s_ in res_r)
    cert(dict(kind='component', comp=i_, dim=dim_, gens=[str(g) for g in gens], indep=indep,
              W=[[str(x) for x in w] for w in W], cuts=cuts, result='pass' if ok else 'fail'))
    print(f"  component {i_}: dim {dim_}; W_Z = {[list(w) for w in W]}; " +
          ("finitely many orbits: True" if ok else f"FAILS at r = {res_r[-1][0]}: {res_r[-1][1]}") + f"; t = {time.time() - t0:.0f}s", flush=True)
