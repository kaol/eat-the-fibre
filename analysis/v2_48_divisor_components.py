"""v2.48: components of the rank-drop divisor D ⊂ V_A, and the behaviour test at a point of D.

(29) For each chart k: eliminate the c-variables from the incidence ideal (scripts written by v2_47_divisor.py),
     map to the ring of the z-variables, minAssGTZ; report the components' dimensions; dump the 6-dimensional primes
     and compare them across charts (mutual ideal containment via reduction modulo Gröbner bases).
(30) At P₁ = P₀ + N (the v2.44 line at t = 1): whether P₁ lies on the divisorial component, dim T_{P₁}D (exact rank
     of the Jacobian of the prime's generators), and dim(T_{P₁}D ∩ [Der(kQ), P₁]).
Usage: python3 analysis/v2_48_divisor_components.py REPS.pkl REP_INDEX LINES_LOG CHART_PREFIX CHARTS(e.g. 0,1,3)"""
import os, sys, io, re, ast, contextlib, subprocess
import sympy as sp
from sympy.polys.matrices import DomainMatrix
HERE = os.path.dirname(os.path.abspath(__file__))
rf, ri, logf, pref, charts = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4], [int(x) for x in sys.argv[5].split(',')]
QQ = sp.QQ
def rank_rows(rows): return DomainMatrix.from_list_sympy(len(rows), len(rows[0]), rows).convert_to(QQ).rank() if rows else 0
comps = {}
for k in charts:
    src = open(f"{pref}_chart{k}.sing").read()
    ring = re.search(r"ring r = 0,\((.*?)\),dp;", src).group(1).split(',')
    ideal = re.search(r"ideal I = (.*?);\n", src, re.S).group(1)
    cv = [v for v in ring if v.startswith('c')]; zv = [v for v in ring if v.startswith('z')]
    out = f"{pref}_comp{k}.txt"
    script = f"""LIB "primdec.lib";
ring r = 0,({','.join(ring)}),dp;
ideal I = {ideal};
ideal E = eliminate(I, {'*'.join(cv)});
ring rz = 0,({','.join(zv)}),dp;
ideal E2 = imap(r, E);
list L = minAssGTZ(E2);
int i; int j;
for (i = 1; i <= size(L); i++) {{
  ideal Gi = std(L[i]);
  print("comp " + string(i) + " dim " + string(dim(Gi)));
  if (dim(Gi) == 6) {{ write(":w {out}." + string(i), string(Gi)); }}
}}
quit;
"""
    open(f"{pref}_comp{k}.sing", 'w').write(script)
    res = subprocess.run(["Singular", "-q", "--no-rc", f"{pref}_comp{k}.sing"], capture_output=True, text=True, timeout=3000)
    print(f"chart {k}: " + res.stdout.strip().replace("\n", "; "), flush=True)
    zs = sp.symbols(zv)
    for f in sorted(os.listdir(os.path.dirname(out) or '.')):
        if os.path.basename(f).startswith(os.path.basename(out) + "."):
            gens = [sp.sympify(gstr.replace('^', '**'), locals={str(z): z for z in zs})
                    for gstr in open(os.path.join(os.path.dirname(out), f)).read().replace('\n', '').split(',') if gstr.strip()]
            comps[(k, f)] = gens
print(f"6-dimensional primes found: {len(comps)}", flush=True)
# (29) compare primes across charts by mutual containment of Gröbner bases (sympy over QQ)
keys = list(comps)
if len(keys) > 1:
    zs_all = sorted(set().union(*[set().union(*[g.free_symbols for g in comps[kk]]) for kk in keys]), key=str)
    GB = {kk: sp.groebner(comps[kk], *zs_all, order='grevlex') for kk in keys}
    for a in range(len(keys)):
        for b in range(a + 1, len(keys)):
            ab = all(GB[keys[b]].reduce(gg)[1] == 0 for gg in comps[keys[a]])
            ba = all(GB[keys[a]].reduce(gg)[1] == 0 for gg in comps[keys[b]])
            print(f"(29) prime {keys[a]} vs {keys[b]}: equal = {ab and ba}", flush=True)
# (30) the point P₁ and the behaviour test
sys.argv = ['family_curves.py', rf]
g = {'__name__': 'fc', '__file__': os.path.join(HERE, 'family_curves.py')}
with contextlib.redirect_stdout(io.StringIO()):
    exec(open(os.path.join(HERE, 'family_curves.py')).read().split("fams = []; nlift")[0], g)
d, M = g['d'], g['M']; I = sp.eye(d)
P0 = sp.Matrix(g['lifted'](ri).tolist())
Nd = [ast.literal_eval(m.group(1)) for m in re.finditer(r"N(?: \(nonzero entries\):| =) (\{.*\})", open(logf).read())][0]
N = sp.zeros(d, d)
for (i, j), v in Nd.items(): N[i, j] = sp.Rational(v)
BA = P0.columnspace(); BK = (I - P0).columnspace(); nA, nK = len(BA), len(BK)
B = sp.Matrix.hstack(*BA, *BK); Binv = B.inv(); RK = Binv[nA:, :]; BAm = sp.Matrix.hstack(*BA)
ZN = Binv[:nA, :] * N * sp.Matrix.hstack(*BK)
zsyms = [sp.Symbol(f'z{i}_{j}') for i in range(nA) for j in range(nK)]
pt = {sp.Symbol(f'z{i}_{j}'): ZN[i, j] for i in range(nA) for j in range(nK)}
P1 = P0 + N
# derivations over Q
Mi = [[[int(M[i, j, m]) for m in range(d)] for j in range(d)] for i in range(d)]
rows = []
for i in range(d):
    for j in range(d):
        for m in range(d):
            row = [0] * (d * d)
            for k_ in range(d):
                if Mi[i][j][k_]: row[m * d + k_] += Mi[i][j][k_]
            for r_ in range(d):
                if Mi[r_][j][m]: row[r_ * d + i] -= Mi[r_][j][m]
                if Mi[i][r_][m]: row[r_ * d + j] -= Mi[i][r_][m]
            if any(row): rows.append(row)
Dn = DomainMatrix.from_list_sympy(len(rows), d * d, rows).convert_to(QQ).nullspace().to_Matrix()
brs = [list(sp.Matrix(d, d, list(Dn.row(k_))) * P1 - P1 * sp.Matrix(d, d, list(Dn.row(k_)))) for k_ in range(Dn.rows)]
rB = rank_rows(brs)
for kk, gens in comps.items():
    on = all(sp.simplify(gg.subs(pt)) == 0 for gg in gens)
    if not on: print(f"(30) prime {kk}: P1 not on it"); continue
    J = [[sp.diff(gg, z).subs(pt) for z in zsyms] for gg in gens]
    rJ = rank_rows(J)
    # tangent vectors of D at P1 as d×d matrices N' = B_A Z' R_K
    Tz = DomainMatrix.from_list_sympy(len(J), len(zsyms), J).convert_to(QQ).nullspace().to_Matrix()
    Tm = [list(BAm * sp.Matrix(nA, nK, list(Tz.row(r_))) * RK) for r_ in range(Tz.rows)]
    inter = len(Tm) + rB - rank_rows(brs + Tm)
    print(f"(30) prime {kk}: P1 on it; dim T_P1 D = {len(zsyms) - rJ}; dim(T_P1 D ∩ [Der, P1]) = {inter}; "
          f"directions outside [Der, P1]: {len(Tm) - inter}", flush=True)
