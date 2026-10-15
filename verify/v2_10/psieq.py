"""First-order chaos-last one-a corrections psi(h) = C(h a q0):  M psi(h) - psi(M h) = Phi((1-p)(q0 h)) psi(1).
Compare the solution space with the engine's tangent projected to one-a chaos-last coordinates."""
import sys
exec(open('psi.py').read().split("words = list")[0])
import fold as FD
from gen2 import GenSym, equations
import gen as G0
Lh, delta = int(sys.argv[1]), 1
G = Gen(S, {})            # for the power expansion table T (q_k = sum_j T[k][j] q0^j)
var = {}; names = []
def V(h, b):
    if (h, b) not in var: var[(h, b)] = len(names); names.append((h, b))
    return var[(h, b)]
Hl = [h for h in Hs if S.lev[h] <= Lh]
def psi_lin(vec):   # psi applied to an H-vector, as linear forms {(var): coef} per output basis
    out = {}
    for h, c in vec.items():
        if S.lev[h] > Lh: raise FD.OutOfWindow
        for b in Hs:
            if S.lev[b] <= S.lev[h] + delta:
                out.setdefault(b, {}); j = V(h, b); out[b][j] = (out[b].get(j, 0) + c) % p
    return out
def M_lin(lin):     # apply M to a vector of linear forms
    out = {}
    for b, form in lin.items():
        for k, v in Mop({b: 1}).items():
            out.setdefault(k, {})
            for j, c in form.items(): out[k][j] = (out[k].get(j, 0) + v * c) % p
    return out
rows = []
for h in Hl:
    try:
        lhs = M_lin(psi_lin({h: 1}))
        Mh = Mop({h: 1})
        rhs1 = psi_lin(Mh)
        # (1-p)(q0 h) in Q-basis -> powers of q0 -> Phi = sum beta_n M^n applied to psi(1)
        qpart = {k: v for k, v in S.mul(q0, h).items() if S.side[k] == 'Q'}
        beta = {}
        for qk, c in qpart.items():
            for j, t in G.T[qk].items(): beta[j] = (beta.get(j, 0) + c * t) % p
        cur = psi_lin({one: 1}); phi = {}
        for n in range(1, max(beta) + 1 if beta else 1):
            cur = M_lin(cur)
            if beta.get(n):
                for b, form in cur.items():
                    phi.setdefault(b, {})
                    for j, c in form.items(): phi[b][j] = (phi[b].get(j, 0) + beta[n] * c) % p
        keys = set(lhs) | set(rhs1) | set(phi)
        for b in keys:
            r = dict(lhs.get(b, {}))
            for j, c in rhs1.get(b, {}).items(): r[j] = (r.get(j, 0) - c) % p
            for j, c in phi.get(b, {}).items(): r[j] = (r.get(j, 0) - c) % p
            rows.append({j: c for j, c in r.items() if c})
    except FD.OutOfWindow: pass
K, _ = FD.kernel(rows, len(names))
low = [j for j, (h, b) in enumerate(names) if S.lev[h] <= 1]
E = FD.Echelon(); rk = sum(E.add({j: v[j] for j in low if j in v}) for v in K)
print(f"ψ-equation (inputs ≤ level {Lh}): solution space projected to inputs ≤ level 1: dim {rk}")
# engine comparison: tangent with chaos-last support, projected to one-a words h a q0 with lev(h) ≤ 1
words = list(G0.gen_words(S, 3, 3, q0))
Gs = GenSym(S, delta, lambda m: m[-1] == q0, 12); eqs, sk = equations(Gs, words)
Kt, _ = FD.kernel(FD.linear_part(eqs), len(Gs.names))
lowt = [j for j, (m, b) in enumerate(Gs.names) if len(m) == 2 and S.lev[m[0]] <= 1]
E2 = FD.Echelon(); rk2 = sum(E2.add({j: v[j] for j in lowt if j in v}) for v in Kt)
print(f"engine tangent (all a-degrees, chaos-last) projected to C(h a q0), lev(h) ≤ 1: dim {rk2}   (skips {sk})")
