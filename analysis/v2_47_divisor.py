"""v2.47: the rank-drop locus D = {P ∈ V_A : rank [𝔥′, P] ≤ dim 𝔥′ − 1} via the incidence variety
{(c, Z) : Z ∈ V_A, Σ cᵢ [Dᵢ, P(Z)] = 0}, one Singular run per affine chart c_k = 1.
𝔥′ = stabiliser of A = im P₀ in Der_Δ ⊕ Lie(T_W), W = admissible weights along the given line (as in v2_46_slice.py).
Usage: python3 analysis/v2_47_divisor.py REPS.pkl REP_INDEX LINES_LOG OUT_PREFIX"""
import os, sys, io, re, contextlib, subprocess
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
rf, ri, logf, outp = sys.argv[1:5]
# reuse v2_46_slice.py up to the rank printout: gives g, d, P0, Pt, t, DerD, Ec, W, h (𝔥′ basis) ...
src = open(os.path.join(HERE, 'v2_46_slice.py')).read().split("for tv in (1, 2, 3):")[0]
sys.argv = ['v2_46_slice.py', rf, ri, logf, '0']
G = {'__name__': 'slice', '__file__': os.path.join(HERE, 'v2_46_slice.py')}
exec(src, G)
d, P0, h, M, DL = G['d'], G['P0'], G['h'], G['M'], G['DL']
I = sp.eye(d)
BA = P0.columnspace(); BK = (I - P0).columnspace(); nA, nK = len(BA), len(BK)
B = sp.Matrix.hstack(*BA, *BK); Binv = B.inv(); RK = Binv[nA:, :]; BAm = sp.Matrix.hstack(*BA)
Z = sp.Matrix(nA, nK, lambda i, j: sp.Symbol(f'z{i}_{j}'))
P = (P0 + BAm * Z * RK).applyfunc(sp.expand); Qm = (I - P).applyfunc(sp.expand)
def mul(x, y): return [sp.expand(sum(int(M[i, j, m]) * x[i] * y[j] for i in range(d) for j in range(d) if M[i, j, m])) for m in range(d)]
eqs = set()
for a in range(d):
    for b in range(d):
        for e in (P * sp.Matrix(mul(list(Qm[:, a]), list(Qm[:, b])))):
            e = sp.expand(e)
            if e != 0: eqs.add(e)
DLs = [sp.Matrix(DL[m].tolist()) for m in range(d)]
for l in range(d):
    for e in (sum((P[m, l] * DLs[m] for m in range(d)), sp.zeros(d, d)) - P * DLs[l] * P.T):
        e = sp.expand(e)
        if e != 0: eqs.add(e)
nh = len(h)
cs = sp.symbols(f'c0:{nh}')
comm = (sum((cs[i] * h[i] for i in range(nh)), sp.zeros(d, d)))
inc = set(e for e in (comm * P - P * comm).applyfunc(sp.expand) if e != 0)
print(f"V_A: {len(eqs)} equations; dim 𝔥′ = {nh}; incidence equations {len(inc)}", flush=True)
zs = [str(v) for v in Z]
for k in [int(x) for x in os.environ.get('CHARTS', ','.join(map(str, range(nh)))).split(',')]:
    sub = {cs[k]: 1}
    eq_all = [str(e).replace('**', '^') for e in eqs] + [str(sp.expand(e.subs(sub))).replace('**', '^') for e in inc]
    eq_all = [e for e in eq_all if e != '0']
    cvars = [str(c) for i, c in enumerate(cs) if i != k]
    script = f"""ring r = 0,({','.join(cvars + zs)}),dp;
ideal I = {','.join(eq_all)};
ideal G = std(I);
print("chart {k}: dim " + string(dim(G)));
ideal E = eliminate(I, {'*'.join(cvars) if cvars else '1'});
print("chart {k}: dim of projection to V_A " + string(dim(std(E))));
quit;
"""
    open(f"{outp}_chart{k}.sing", 'w').write(script)
    res = subprocess.run(["Singular", "-q", "--no-rc", f"{outp}_chart{k}.sing"], capture_output=True, text=True, timeout=3000)
    print(res.stdout.strip() or res.stderr[-500:], flush=True)
