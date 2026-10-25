"""v2.54: finiteness of census orbits on each component of a vertex-block stratum of Q₀ = 1⇉2⇇3 + 1→3.

For each component Z (as in v2_53_orbits.py: the same presolve, primes, W_Z and 𝔥_Z = Der_Δ ⊕ E_{W_Z}) and each
r < dim Z, this computes an upper bound for dim Z_r, Z_r = {P ∈ Z : rank [𝔥_Z, P] ≤ r}:
    dim (Z ∩ V(f_1, …, f_N)),   f_k = an (r+1)-minor of M(P) on random rows and columns, reduced modulo the prime,
where M(P) is the m × (nonzero columns) matrix of the [h_i, P]. Each f_k vanishes on Z_r. N = 3(dim Z + 2).
(v2.54's first version used dense random combinations det(L·M·R); they were too slow for 5×5 minors.)
The criterion (logs/v2_54/PREDICTIONS.md): if the bound is ≤ r for every r < dim Z, then H_Z has finitely many
orbits on Z. When a bound is > r (or a Singular call exceeds 900 s), the r is retried with 8(dim Z + 2) minors before it is reported as failed.
Usage: python3 analysis/v2_54_orbit_strata.py BLOCK OUT_PREFIX [--seed S]"""
import os, sys, io, random, subprocess, contextlib, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser(); ap.add_argument('block'); ap.add_argument('out'); ap.add_argument('--seed', type=int, default=54)
A = ap.parse_args()
src = open(os.path.join(HERE, 'v2_53_orbits.py')).read().split("for i_, dim_, gens, indep, nf in comps:")[0]
sys.argv = ['v2_53_orbits.py', A.block, A.out + '_base']
G = {'__name__': 'orb', '__file__': os.path.join(HERE, 'v2_53_orbits.py')}
with contextlib.redirect_stdout(io.StringIO()):
    exec(src, G)
sp, comps, vs, cl, nz, DerD, Ec, Pr, d, nullsp = (G[k] for k in ('sp', 'comps', 'vs', 'cl', 'nz', 'DerD', 'Ec', 'Pr', 'd', 'nullsp'))
rng = random.Random(A.seed)
def sing_poly(e): return str(sp.expand(e)).replace('**', '^')
print(f"V = {A.block}: {len(comps)} components", flush=True)
summary = []
for i_, dim_, gens, indep, nf in comps:
    n_ = len(nz); wv = sp.symbols(f'w0:{len(cl)}')
    lin = []
    for k in range(n_):
        e = sp.expand(sum(wv[c] * nf[c * n_ + k] for c in range(len(cl))))
        if e != 0:
            for co in (sp.Poly(e, *[sp.Symbol(v) for v in vs]).coeffs() if vs else [e]): lin.append([sp.expand(co).coeff(w) for w in wv])
    W = nullsp(lin, len(cl))
    h = DerD + [sum((w[c] * Ec[c] for c in range(len(cl))), sp.zeros(d, d)) for w in W]
    m = len(h)
    rows = [list((D * Pr - Pr * D).applyfunc(sp.expand)) for D in h]
    cols = [j for j in range(d * d) if any(rows[i][j] != 0 for i in range(m))]
    if dim_ == 0: summary.append((i_, dim_, [])); print(f"  component {i_}: dim 0 (one point)", flush=True); continue
    res_r = []
    Sv = [sp.Symbol(v) for v in vs]
    def singular(lines, tmo=1200):
        fn = f"{A.out}_c{i_}.sing"; open(fn, 'w').write("\n".join(lines))
        try: return subprocess.run(["Singular", "-q", "--no-rc", fn], capture_output=True, text=True, timeout=tmo).stdout
        except subprocess.TimeoutExpired: return ""
    def components_above(F, r):
        """components of Z ∩ V(F) of dimension > r: list of (dim, generators, independent set)"""
        out = singular([f'LIB "primdec.lib";', f"ring R = 0,({','.join(vs)}),dp;",
                        f"ideal J = {','.join([sing_poly(g) for g in gens] + [sing_poly(f) for f in F]) or '0'};",
                        "list L = minAssGTZ(J); int i;",
                        f"for (i = 1; i <= size(L); i++) {{ ideal Gi = std(L[i]); if (dim(Gi) > {r}) {{ print(\"COMP \" + string(dim(Gi))); "
                        "print(\"GENS \" + string(Gi)); intvec iv = indepSet(Gi); print(\"INDEP \" + string(iv)); kill iv; } kill Gi; }",
                        'print("END");', "quit;"])
        if "END" not in out: return None
        res = []
        for blk in out.split("COMP ")[1:]:
            dm = int(blk.split()[0])
            gl = [sp.sympify(g.replace('^', '**'), locals={v: sp.Symbol(v) for v in vs})
                  for g in blk.split("GENS ")[1].split("\n")[0].split(',') if g.strip()]
            iv = [int(x) for x in blk.split("INDEP ")[1].split("\n")[0].split(',')]
            res.append((dm, gl, iv))
        return res
    def points(gl, iv, k):
        """up to k rational points of the variety of gl (random values on the independent set, solve the rest)"""
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
    Mrows = [[rows[i][j] for j in cols] for i in range(m)]
    for r in range(dim_ - 1, -1, -1):
        F = []; status = None
        for it in range(60):
            comps_above = components_above(F, r)
            if comps_above is None: status = 'singular timeout'; break
            if not comps_above: status = 'ok'; break
            added = False
            for dm, gl, iv in comps_above:
                pts = points(gl, iv, 4)
                if not pts: status = f'no rational point found on a component of dim {dm}'; break
                for pt in pts:
                    Mv = sp.Matrix([[x.subs(pt) for x in row] for row in Mrows])
                    if Mv.rank() <= r: continue
                    # a nonzero (r+1)-minor at pt: pivot columns of Mv, then pivot rows of those columns
                    _, pc = Mv.rref(); pc = list(pc)[:r + 1]
                    _, pr = Mv[:, pc].T.rref(); pr = list(pr)[:r + 1]
                    F.append(sp.expand(sp.Matrix([[Mrows[i][j] for j in pc] for i in pr]).det(method='berkowitz')))
                    added = True; break
                if not added and not status:
                    status = f'rank ≤ {r} at every sampled point of a component of dim {dm}'; break
            if status: break
        bound = r if status == 'ok' else None
        res_r.append((r, bound if bound is not None else status))
    ok = all(isinstance(b, int) and b <= r for r, b in res_r)
    summary.append((i_, dim_, res_r))
    print(f"  component {i_}: dim {dim_}; W_Z = {[list(w) for w in W]}; bounds dim Z_r (r: bound) = "
          f"{', '.join(f'{r}: {b}' for r, b in res_r)}; finitely many orbits: {ok}", flush=True)
