"""One-a chaos-last folds C(h a q0) = psi(h), all other C = 0.  M = p∘(q0·) on H0;  test exactness of psi families."""
import sys, sympy as sp, random
from gen import Gen, gen_words, check
from fold import MOD, md
from tangent import spec
S = spec('dihedral', 50); p = MOD; u = sp.Symbol('u')
q0 = sorted([i for i in range(len(S.B)) if S.side[i] == 'Q'], key=lambda i: S.lev[i])[0]
Hs = sorted([i for i in range(len(S.B)) if S.side[i] == 'H'], key=lambda i: S.lev[i]); one = Hs[0]
rt2 = pow(2, (p + 1) // 4, p); roots = {'r+': (1 + rt2) % p, 'r-': (1 - rt2) % p}
def ev(i, r): return sum(md(c) * pow(r, m[0], p) for m, c in S.B[i][2].terms()) % p
LEV = 20
def Mop(vec):   # vec: dict H-index -> coef ; returns p(q0 * vec)
    out = {}
    for h, c in vec.items():
        for k, v in S.mul(q0, h).items():
            if S.side[k] == 'H': out[k] = (out.get(k, 0) + c * v) % p
    return {k: v for k, v in out.items() if v}
def add(*terms):
    out = {}
    for c, vec in terms:
        for k, v in vec.items(): out[k] = (out.get(k, 0) + c * v) % p
    return {k: v for k, v in out.items() if v}
def evv(vec, r): return sum(c * ev(k, r) for k, c in vec.items()) % p
def Fpoly(coeffs, vec):  # sum_n coeffs[n] M^n vec
    out = {}; cur = dict(vec)
    for n, c in enumerate(coeffs):
        if c: out = add((1, out), (c, cur))
        cur = Mop(cur)
    return out
def C_of(psi):
    return {(h, q0): d for h in Hs if S.lev[h] <= LEV for d in [psi({h: 1})] if d}
words = list(gen_words(S, 3, int(sys.argv[1]) if len(sys.argv) > 1 else 3, q0))
def test(nm, psi):
    b, t = check(Gen(S, C_of(psi)), words, 2); print(f"{nm:58s}: bad {b}/{t}", flush=True); return b
r = roots['r+']
dev = lambda vec, r=r: add((1, vec), (-evv(vec, r), {one: 1}))
test("root deviation  x − ev_r(x)", dev)
test("M ∘ (id − ev_r)   [= M, since M kills 1̄]", lambda v: Mop(dev(v)))
test("(1 + 2M + 3M²) ∘ (id − ev_r)", lambda v: Fpoly([1, 2, 3], dev(v)))
test("(id − ev_r) ∘ (1 + 2M)", lambda v: dev(Fpoly([1, 2], v)))
test("M alone", lambda v: Mop(v))
test("2·(id − ev_r+) + 3·M", lambda v: add((2, dev(v)), (3, Mop(v))))
test("constant: ψ(x) = ev_r(x)·1̄  (e₀-type)", lambda v: {one: evv(v, r)} if evv(v, r) else {})
test("(id − ev_r) with r = r−", lambda v: dev(v, roots['r-']))
