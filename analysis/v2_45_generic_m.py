"""v2.45: (19) do two lifted representatives have the same image? (20) m at generic points of the v2.44 line
P₀ + t·N (random t mod a 7-digit prime), with m = dim T_P V_A − dim([Der_A, P] ∩ T_P V_A), A = im P₀ (fixed along the
line). Usage: python3 analysis/v2_45_generic_m.py REPS.pkl REP_A REP_B LINES_LOG"""
import os, sys, io, re, ast, contextlib, random
import numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
rf, ra, rb, logf = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
PR = 9999991
sys.argv = ['family_curves.py', rf]
g = {'__name__': 'fc', '__file__': os.path.join(HERE, 'family_curves.py')}
with contextlib.redirect_stdout(io.StringIO()):
    exec(open(os.path.join(HERE, 'family_curves.py')).read().split("fams = []; nlift")[0], g)
sys.argv = ['fibre_rigidity.py']
fr = {'__name__': 'fr', '__file__': os.path.join(HERE, 'fibre_rigidity.py')}
exec(open(os.path.join(HERE, 'fibre_rigidity.py')).read().split("\ndef main():")[0], fr)
d = g['d']; M = np.array(g['M'], np.int64); DL = np.array(g['DL'], np.int64)
CP = fr['setup'](g['r']['name'], PR); rank, nullspace = fr['rank'], fr['nullspace']
PA = sp.Matrix(g['lifted'](ra).tolist()); PB = sp.Matrix(g['lifted'](rb).tolist())
ca, cb = PA.columnspace(), PB.columnspace()
same = sp.Matrix.hstack(*ca).rank() == sp.Matrix.hstack(*cb).rank() == sp.Matrix.hstack(*ca, *cb).rank()
print(f"(19) im P(rep {ra}) = im P(rep {rb}) over Q: {same}  (dims {len(ca)}, {len(cb)}, joint rank "
      f"{sp.Matrix.hstack(*ca, *cb).rank()})", flush=True)
I = np.eye(d, dtype=np.int64)
def m_at(P):
    Q = (I - P) % PR
    J = fr['jacobian'](P, dict(d=d, M=M, DL=DL), PR)
    Wv = np.array([(P[:, [a]] @ Q[[b], :] % PR).ravel() for a in range(d) for b in range(d)])
    Wv = Wv[[k for k in range(len(Wv)) if Wv[k].any()]]
    R, piv = fr['rref'](Wv, PR); Wb = R[:len(piv)]
    Cc = nullspace((J @ Wb.T) % PR, PR); TV = (Cc @ Wb) % PR
    cond = np.array([((Q @ D @ P) % PR).ravel() for D in CP['Der']]).T
    cA = nullspace(cond % PR, PR)
    br = np.array([((np.einsum('k,kij->ij', c, CP['Der']) % PR) @ P - P @ (np.einsum('k,kij->ij', c, CP['Der']) % PR)) % PR
                   for c in cA]).reshape(len(cA), -1) if len(cA) else np.zeros((0, d * d), np.int64)
    inter = rank(TV, PR) + rank(br, PR) - rank(np.concatenate([TV, br]), PR) if len(br) else 0
    return len(TV), len(cA), len(TV) - inter
Nd = [ast.literal_eval(m.group(1)) for m in re.finditer(r"N \(nonzero entries\): (\{.*\})", open(logf).read())][0]
N = np.zeros((d, d), np.int64)
for (i, j), v in Nd.items(): N[i, j] = int(sp.Rational(v))
P0 = np.array(g['lifted'](ra), dtype=np.int64) % PR
print("(20) at t = 0: (dim T_V, dim Der_A, m) =", m_at(P0), flush=True)
random.seed(1)
for _ in range(4):
    tv = random.randrange(2, PR - 2)
    print(f"     at random t: (dim T_V, dim Der_A, m) = {m_at((P0 + tv * N) % PR)}", flush=True)
