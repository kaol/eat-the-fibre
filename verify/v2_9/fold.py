"""Exact engine for Δ-compatible folds on k<a,u> up to a-degree 2, over a fixed splitting k[u] = H0 ⊕ Q0.
Adapted basis of k[u]: H-basis and Q-basis elements with a 'level'; p = coordinate projection onto H.
P(x a y)      = p(x) a p(y) + c(x,y)
P(x a y a z)  = p(x) a p(y) a p(z) + p(x) a c(y,z) + c(x,y) a p(z) + c2(x,y,z)
Unknowns c, c2 (values in k[u], output level <= input level + delta). Arithmetic mod a prime.
Equations: P∘P = P, P(PsPt) = PsPt, P((1-P)s(1-P)t) = 0 for monomial s,t (a-degree sum <= 2).
Coefficients are polynomials in the unknowns truncated at degree 2 (tangent + second order)."""
import sympy as sp, itertools
MOD = 2305843009213693951
def md(r):
    r = sp.Rational(r); return (int(r.p) * pow(int(r.q), MOD - 2, MOD)) % MOD
u = sp.symbols('u')

class Splitting:
    def __init__(self, Hb, Qb, Lmax):
        # Hb, Qb: lists of (poly, level); covering k[u] up to level Lmax
        self.B = [('H', l, sp.Poly(f, u)) for f, l in Hb if l <= Lmax] + [('Q', l, sp.Poly(f, u)) for f, l in Qb if l <= Lmax]
        self.Lmax = Lmax
        D = max(b[2].degree() for b in self.B)
        M = sp.Matrix([[b[2].coeff_monomial(u**m) for b in self.B] for m in range(D + 1)])
        assert M.shape[0] == M.shape[1] and M.det() != 0, (M.shape)
        self.Minv = M.inv(); self.D = D
        self.side = [b[0] for b in self.B]; self.lev = [b[1] for b in self.B]
        self.idx = {}
        self._mul = {}
    def decompose(self, f):
        f = sp.Poly(sp.expand(f), u)
        assert f.degree() <= self.D, f
        v = sp.Matrix([f.coeff_monomial(u**m) for m in range(self.D + 1)])
        c = self.Minv * v
        return {i: md(c[i]) for i in range(len(self.B)) if c[i] != 0}
    def mul(self, i, j):
        if (i, j) not in self._mul:
            if self.lev[i] + self.lev[j] > self.Lmax: raise OutOfWindow
            r = self.decompose(self.B[i][2].as_expr() * self.B[j][2].as_expr())
            for k in r: assert self.lev[k] <= self.lev[i] + self.lev[j]
            self._mul[(i, j)] = self._mul[(j, i)] = r
        return self._mul[(i, j)]

class OutOfWindow(Exception): pass

# ---- coefficient polynomials (dict monomial-tuple -> int), truncated at degree 2
def padd(a, b, s=1):
    r = dict(a)
    for k, v in b.items():
        r[k] = (r.get(k, 0) + s * v) % MOD
        if r[k] == 0: del r[k]
    return r
def pmul(a, b):
    r = {}
    for k1, v1 in a.items():
        for k2, v2 in b.items():
            k = tuple(sorted(k1 + k2))
            if len(k) > 2: continue
            r[k] = (r.get(k, 0) + v1 * v2) % MOD
    return {k: v for k, v in r.items() if v}
def pscale(a, s): return {k: v * s % MOD for k, v in a.items() if v * s % MOD}
ONE = {(): 1}

class Fold:
    def __init__(self, S, Lw, delta, fixed_c=None):
        self.S, self.Lw, self.delta = S, Lw, delta
        self.var = {}; self.names = []
        self.fixed_c = fixed_c  # optional: dict word1 -> {out: coefpoly}
    def v(self, key):
        if key not in self.var: self.var[key] = len(self.names); self.names.append(key)
        return {(self.var[key],): 1}
    def wlev(self, w): return sum(self.S.lev[i] for i in w)
    def outs(self, w):
        L = min(self.wlev(w) + self.delta, self.S.Lmax)
        return [o for o in range(len(self.S.B)) if self.S.lev[o] <= L]
    def c(self, w):
        if self.wlev(w) > self.Lw: raise OutOfWindow
        if self.fixed_c is not None: return self.fixed_c(w)
        return {o: self.v(('c', w, o)) for o in self.outs(w)}
    def c2(self, w):
        if self.wlev(w) > self.Lw: raise OutOfWindow
        return {o: self.v(('c2', w, o)) for o in self.outs(w)}
    def isH(self, i): return self.S.side[i] == 'H'
    def P(self, E):
        R = {}
        def add(word, coef):
            R[word] = padd(R.get(word, {}), coef)
            if not R[word]: del R[word]
        for w, cf in E.items():
            if len(w) == 1:
                if self.isH(w[0]): add(w, cf)
            elif len(w) == 2:
                x, y = w
                if self.isH(x) and self.isH(y): add(w, cf)
                for o, val in self.c(w).items(): add((o,), pmul(cf, val))
            else:
                x, y, z = w
                if self.isH(x) and self.isH(y) and self.isH(z): add(w, cf)
                if self.isH(x):
                    for o, val in self.c((y, z)).items(): add((x, o), pmul(cf, val))
                if self.isH(z):
                    for o, val in self.c((x, y)).items(): add((o, z), pmul(cf, val))
                for o, val in self.c2(w).items(): add((o,), pmul(cf, val))
        return R
    def mul(self, A, B):
        R = {}
        for w1, c1 in A.items():
            for w2, c2 in B.items():
                if len(w1) + len(w2) - 2 > 2: continue
                cf = pmul(c1, c2)
                if not cf: continue
                for k, m in self.S.mul(w1[-1], w2[0]).items():
                    w = w1[:-1] + (k,) + w2[1:]
                    R[w] = padd(R.get(w, {}), pscale(cf, m))
        return {w: c for w, c in R.items() if c}
    def sub(self, A, B):
        R = dict(A)
        for w, c in B.items():
            R[w] = padd(R.get(w, {}), c, -1)
            if not R[w]: del R[w]
        return R

def monomials(S, Lw, adeg):
    n = len(S.B)
    for w in itertools.product(range(n), repeat=adeg + 1):
        if sum(S.lev[i] for i in w) <= Lw: yield w

SKIPS = [0]
def equations(F, Leq):
    S = F.S; eqs = []; SKIPS[0] = 0
    def emit(E):
        for w, cf in E.items():
            if cf: eqs.append(cf)
    mons = {d: list(monomials(S, Leq, d)) for d in range(3)}
    for d in range(3):
        for s in mons[d]:
            try: Ps = F.P({s: ONE}); emit(F.sub(F.P(Ps), Ps))
            except OutOfWindow: SKIPS[0] += 1
    for d1 in range(3):
        for d2 in range(3 - d1):
            for s in mons[d1]:
                for t in mons[d2]:
                    if F.wlev(s) + F.wlev(t) > Leq: continue
                    try:
                        Ps, Pt = F.P({s: ONE}), F.P({t: ONE})
                        PP = F.mul(Ps, Pt); emit(F.sub(F.P(PP), PP))
                        Qs, Qt = F.sub({s: ONE}, Ps), F.sub({t: ONE}, Pt)
                        emit(F.P(F.mul(Qs, Qt)))
                    except OutOfWindow: SKIPS[0] += 1
    return eqs

# ---- sparse linear algebra mod p
class Echelon:
    def __init__(self): self.piv = {}  # col -> row(dict), row normalized with row[col]=1
    def reduce(self, r):
        r = dict(r)
        while True:
            cols = [c for c in r if c in self.piv]
            if not cols: return r
            c = cols[0]; f = r[c]
            for k, v in self.piv[c].items():
                r[k] = (r.get(k, 0) - f * v) % MOD
                if r[k] == 0: del r[k]
    def add(self, r):
        r = self.reduce(r)
        if not r: return False
        c = min(r); inv = pow(r[c], MOD - 2, MOD)
        r = {k: v * inv % MOD for k, v in r.items()}
        self.piv[c] = r; return True
    def full_reduce(self):
        # back-substitute to RREF
        for c in sorted(self.piv, reverse=True):
            r = self.piv[c]
            for c2 in list(self.piv):
                if c2 != c and c in self.piv[c2]:
                    f = self.piv[c2][c]
                    nr = dict(self.piv[c2])
                    for k, v in r.items():
                        nr[k] = (nr.get(k, 0) - f * v) % MOD
                        if nr[k] == 0: del nr[k]
                    self.piv[c2] = nr

def linear_part(eqs):
    return [{k[0]: v for k, v in e.items() if len(k) == 1} for e in eqs]

def kernel(rows, nvars):
    E = Echelon()
    for r in rows:
        if r: E.add(r)
    E.full_reduce()
    free = [j for j in range(nvars) if j not in E.piv]
    K = []
    for f in free:
        vec = {f: 1}
        for c, r in E.piv.items():
            if f in r: vec[c] = (-r[f]) % MOD
        K.append(vec)
    return K, E

def rank(vecs):
    E = Echelon(); return sum(E.add(v) for v in vecs if v)
