"""v2.62 (Program 3): the relative orbit criterion for Q₀ → Q₁ over the 366 characteristic-0 census classes of Q₀.

For each class representative P̄ (logs/v2_62/q0_classes.jsonl, from v2_62_dedupe.py) and each of
  kind 'J' (P(J) ⊆ J: P = [[P̄, 0], [X, P_J]]) and kind 'Q' (P(kQ₀) ⊆ kQ₀: P = [[P̄, Y], [0, P_J]]), λ = 0 and 1,
the fibre F = {folds P of kQ₁ of that shape} is computed exactly (all other entries unknown; presolve; Singular
decomposition). The census classes of Q₁ lying over the class of P̄ (in these kinds) are the orbits of the stabiliser
Stab(P̄) ⊆ H_{Q₁} on F. On each positive-dimensional component C of F the v2.54 criterion is applied with
  𝔰_C = { h ∈ Der_Δ(Q₁) ⊕ E_{W_C} : [h|_{kQ₀}, P̄] = 0 },
W_C the class weightings of Q₁ compatible on all of C: dim {P ∈ C : rank [𝔰_C, P] ≤ r} ≤ r for every r < dim C, by
point-guided minors (as v2_55_pipeline.py). A pass means Stab(P̄)° has finitely many orbits on C (characteristic 0).
Output per (class, kind, λ): empty / points / components with dims and pass/fail; a tally.
Usage: python3 analysis/v2_62_relative.py OUT_PREFIX [--start I] [--stop J] [--timeout S]"""
import os, sys, json, subprocess, collections, argparse, random, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE)); sys.path.insert(0, HERE)
import numpy as np, sympy as sp
from sympy import QQ
from sympy.polys.rings import ring
from sympy.polys.matrices import DomainMatrix
import fold_ideal as fi, pcensus as pc, singular_nohilb
singular_nohilb.install()
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--start', type=int, default=0); ap.add_argument('--stop', type=int, default=10**9)
ap.add_argument('--timeout', type=int, default=1800); ap.add_argument('--classes', default=None)
A = ap.parse_args(); rng = random.Random(62)
root = os.path.dirname(HERE)
classes = [json.loads(l) for l in open(A.classes or os.path.join(root, 'logs', 'v2_62', 'q0_classes.jsonl'))]
c1 = pc.setup('pairs_in_short_tail', 3); c0 = pc.setup('pairs_in_short', 3)
b1, b0 = c1['basis'], c0['basis']; d1, d0 = len(b1), len(b0)
idx1 = {b: i for i, b in enumerate(b1)}; e4 = idx1[('v', 3)]
M1 = np.asarray(c1['M'], dtype=np.int64); DL1 = np.asarray(c1['DL'], dtype=np.int64)
nzM1 = [(i, j, k, int(M1[i, j, k])) for i in range(d1) for j in range(d1) for k in range(d1) if M1[i, j, k]]
Qx = [idx1[b] for b in b0]; Jx = [i for i in range(d1) if i not in Qx]
arrows = c1['arrows']
# Der_Δ(Q₁) and class scalings of Q₁ (as in v2_55_pipeline.py)
def nullsp(rows, n):
    if not rows: return [sp.Matrix([1 if i == j else 0 for i in range(n)]) for j in range(n)]
    Nm = DomainMatrix.from_list_sympy(len(rows), n, rows).convert_to(sp.QQ).nullspace().to_Matrix()
    return [Nm.row(k).T for k in range(Nm.rows)]
d = d1; M = M1; DL = DL1
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
cond = [sum((list(sum((D[m, l] * DLs[m] for m in range(d)), sp.zeros(d, d)) - D * DLs[l] - DLs[l] * D.T) for l in range(d)), []) for D in Der]
cx = nullsp([list(r) for r in sp.Matrix(cond).T.tolist()], len(Der))
DerD = [sum((c[k] * Der[k] for k in range(len(Der))), sp.zeros(d, d)) for c in cx]
cls = collections.defaultdict(list)
for a_, st in enumerate(arrows): cls[st].append(a_)
cl = list(cls.values())
Ec = [sp.diag(*[0 if b[0] == 'v' else sum(1 for x in b[1] if x in c) for b in b1]) for c in cl]
def DLclass(c):
    D_ = np.zeros((d, d, d), dtype=np.int64)
    for l, b in enumerate(b1):
        if b[0] == 'v': continue
        for k, x in enumerate(b[1]):
            if x not in c: continue
            pre = idx1[('p', b[1][:k])] if k > 0 else idx1[('v', arrows[x][0])]
            suf = idx1[('p', b[1][k + 1:])] if k + 1 < len(b[1]) else idx1[('v', arrows[x][1])]
            D_[l, pre, suf] += 1
    return D_
DLc = [DLclass(c) for c in cl]
ss = lambda e: str(e).replace('**', '^')
def singular(lines, fn, tmo=None):
    open(fn, 'w').write("\n".join(lines))
    try: rr = subprocess.run(["Singular", "-q", "--no-rc", fn], capture_output=True, text=True, timeout=tmo or A.timeout)
    except subprocess.TimeoutExpired: return None
    out = rr.stdout + rr.stderr
    return None if "? " in out else out
def parse_comps(out):
    res = []
    for blk in out.split("COMP ")[1:]:
        dim_ = int(blk.split()[0]); gens = blk.split("GENS ")[1].split("\n")[0].split(',')
        iv = [int(x) for x in blk.split("INDEP ")[1].split("\n")[0].split(',')]
        res.append((dim_, gens, iv))
    return res
def decompose(vs, gens, above=-1):
    out = singular([f'LIB "primdec.lib"; ring R = 0,({",".join(vs)}),dp;', f"ideal J = {','.join(g.replace('**', '^') for g in gens) or '0'};",
                    "list L = minAssGTZ(J); int i;",
                    f'for (i = 1; i <= size(L); i++) {{ ideal G = std(L[i]); if (dim(G) > {above}) {{ print("COMP " + string(dim(G))); print("GENS " + string(G)); intvec iv = indepSet(G); print("INDEP " + string(iv)); kill iv; }} kill G; }}',
                    'print("END"); quit;'], A.out + '_dec.sing')
    if out is None or "END" not in out: return None
    return parse_comps(out)
def points(Sv, gl, iv, k):
    pts = []
    for t in range(40 * k):
        vals = {Sv[j]: sp.Integer(rng.randint(-5, 5)) for j in range(len(Sv)) if iv[j] == 1}
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
tally = collections.Counter()
for ci in range(A.start, min(A.stop, len(classes))):
    C_ = classes[ci]
    Pb = sp.Matrix([[sp.Rational(x) for x in row] for row in C_['P']])
    for kind in ('J', 'Q'):
        for lam in (0, 1):
            cells = [(r, c) for c in Jx for r in Jx] + ([(r, c) for c in Jx for r in Qx] if kind == 'Q' else [(r, c) for c in Qx for r in Jx])
            cells.remove((e4, e4))
            R, *xs = ring(','.join(f"x{r}_{c}" for r, c in cells), QQ)
            Xv = dict(zip(cells, xs))
            P = [[R(0)] * d1 for _ in range(d1)]
            for i_, qi in enumerate(Qx):
                for j_, qj in enumerate(Qx): P[qi][qj] = R(QQ(int(sp.fraction(Pb[i_, j_])[0]), int(sp.fraction(Pb[i_, j_])[1])))
            for rc, x in Xv.items(): P[rc[0]][rc[1]] = x
            P[e4][e4] = R(lam)
            T0 = time.time()                    # staged: P² = P and Δ first, substitute, then the closure equations
            rem1, subs = fi.presolve_fast(fi._equations(P, R, d1, nzM1, DL1, ('idem', 'delta')), R)
            P1 = P
            for x, v in subs: P1 = [[e.compose(x, v) for e in row] for row in P1]
            rem, subs2 = fi.presolve_fast(set(rem1) | fi._equations(P1, R, d1, nzM1, DL1, ('closure',)), R)
            subs = subs + subs2; T2 = time.time()
            if os.environ.get('RTIME'): print(f"   timing: build + presolve {T2 - T0:.1f}s", flush=True)
            tag = f"class {ci} ({C_['block']}) kind {kind} λ = {lam}"
            if any(e.is_ground for e in rem): print(f"{tag}: empty", flush=True); tally[(kind, lam, 'empty')] += 1; continue
            for x, v in subs: P = [[e.compose(x, v) for e in row] for row in P]
            elim = {str(x) for x, _ in subs}; vs = [str(x) for x in xs if str(x) not in elim]
            if not vs: print(f"{tag}: one point", flush=True); tally[(kind, lam, 'one point')] += 1; continue
            comps = decompose(vs, [str(e) for e in rem]) if rem else [(len(vs), ['0'], [1] * len(vs))]
            if comps is None: print(f"{tag}: decomposition unresolved", flush=True); tally[(kind, lam, 'unresolved')] += 1; continue
            Sv = sp.symbols(vs); loc = {str(v): v for v in Sv}
            PE = sp.Matrix([[sp.sympify(str(e), locals=loc) for e in row] for row in P])
            verdicts = []
            for dim_, gens, iv in comps:
                if dim_ <= 0: verdicts.append(('0', 'points')); continue
                gl = [sp.sympify(g, locals=loc) for g in gens if g.strip() not in ('0', '')]
                pts0 = points(Sv, gl, iv, 3)
                if not pts0: verdicts.append((str(dim_), 'no rational point')); continue
                # W_C: weightings of Q₁ compatible at sampled points (linear in w); take the common kernel
                lin = []
                for pt in pts0:
                    Pp = np.array(PE.subs(pt).tolist(), dtype=object)
                    for D_ in DLc: pass
                    dfs = [sp.Matrix((np.einsum('ml,mab->lab', Pp, D_) - np.einsum('ax,lxy,by->lab', Pp, D_, Pp)).reshape(-1).tolist()) for D_ in DLc]
                    H_ = sp.Matrix.hstack(*dfs)
                    lin += [list(H_.row(i)) for i in range(H_.rows) if any(H_.row(i))]
                W = nullsp(lin, len(cl))
                hs = DerD + [sum((w[k] * Ec[k] for k in range(len(cl))), sp.zeros(d, d)) for w in W]
                # 𝔰 = {h : [h, P̄] = 0 on kQ₀}: h preserves kQ₀, so [h, P] restricted to rows/cols Qx
                rowsS = [list((Dh[np.ix_(Qx, Qx)] if False else Dh.extract(Qx, Qx)) * Pb - Pb * Dh.extract(Qx, Qx)) for Dh in hs]
                coef = nullsp([list(r) for r in sp.Matrix(rowsS).T.tolist()], len(hs))
                S_alg = [sum((c[k] * hs[k] for k in range(len(hs))), sp.zeros(d, d)) for c in coef]
                if not S_alg: verdicts.append((str(dim_), 'FAIL (trivial stabiliser)')); continue
                Mfun = sp.lambdify(Sv, sp.Matrix([list(Dh * PE - PE * Dh) for Dh in S_alg]), 'sympy')
                Msym = sp.Matrix([list(Dh * PE - PE * Dh) for Dh in S_alg])
                ok = True; why = ''
                for r in range(dim_ - 1, -1, -1):
                    F = []; status = None
                    for it in range(60):
                        above = decompose(vs, gens + F, r)
                        if above is None: status = 'unresolved'; break
                        if not above: status = 'ok'; break
                        added = False
                        for dm, gl2s, iv2 in above:
                            gl2 = [sp.sympify(g, locals=loc) for g in gl2s if g.strip() not in ('0', '')]
                            pts = points(Sv, gl2, iv2, 4)
                            if not pts: status = f'no rational point (dim {dm})'; break
                            for pt in pts:
                                Mv = sp.Matrix(Mfun(*[pt[v] for v in Sv]))
                                if Mv.rank() <= r: continue
                                _, pc_ = Mv.rref(); pc_ = list(pc_)[:r + 1]
                                _, pr_ = Mv[:, pc_].T.rref(); pr_ = list(pr_)[:r + 1]
                                F.append(str(sp.expand(Msym.extract(pr_, pc_).det())))
                                added = True; break
                            if not added and not status: status = f'rank ≤ {r} at every sampled point (dim {dm})'
                            if status: break
                        if status: break
                    if status != 'ok': ok = False; why = f'r = {r}: {status}'; break
                verdicts.append((str(dim_), 'pass' if ok else 'FAIL ' + why))
            print(f"{tag}: components {[v for v in verdicts]}", flush=True)
            for dm, v in verdicts: tally[(kind, lam, f'dim {dm}', v.split(' ')[0])] += 1
print("TALLY")
for k, v in sorted(tally.items(), key=str): print(" ", k, v)
