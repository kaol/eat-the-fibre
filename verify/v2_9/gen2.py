"""gen.py with symbolic C: coefficients are polynomials in unknowns (truncated at degree 2), variables created lazily."""
import itertools, functools, sympy as sp
import fold as FD
from fold import MOD, padd, pmul, pscale
import gen as G0
class GenSym(G0.Gen):
    def __init__(self, S, delta, support, maxlev_vars):
        super().__init__(S, {}); self.delta = delta; self.support = support; self.maxlev_vars = maxlev_vars
        self.var = {}; self.names = []
    def Cval(self, m):
        if not self.support(m): return {}
        lev = sum(self.S.lev[x] for x in m)
        if lev > self.maxlev_vars: raise FD.OutOfWindow
        out = {}
        for b in range(len(self.S.B)):
            if self.S.side[b] == 'H' and self.S.lev[b] <= lev + self.delta:
                key = (m, b)
                if key not in self.var: self.var[key] = len(self.names); self.names.append(key)
                out[b] = {(self.var[key],): 1}
        return out
    # override arithmetic: coefficients are polys
    def mul(self, A, B):
        out = {}
        for w1, c1 in A.items():
            for w2, c2 in B.items():
                cc = pmul(c1, c2)
                if not cc: continue
                for w, m in self.mulw(w1, w2).items(): out[w] = padd(out.get(w, {}), pscale(cc, m))
        return {w: c for w, c in out.items() if c}
    def kappa(self, f):
        i = [t for t, x in enumerate(f) if self.isQ(x)][0]; n = len(f) - 1; out = {}
        for j in [-1] + list(range(i)):
            for k in [n] + list(range(i, n)):
                m = f[j + 1:k + 1]
                if len(m) < 2: continue
                for b, v in self.Cval(m).items():
                    w = f[:j + 1] + (b,) + f[k + 1:]; out[w] = padd(out.get(w, {}), v)
        return {w: c for w, c in out.items() if c}
    def g(self, f):
        out = {f: {(): 1}}
        for w, c in self.kappa(f).items(): out[w] = padd(out.get(w, {}), c, -1)
        return {w: c for w, c in out.items() if c}
    def P(self, w):
        if w in self._P: return self._P[w]
        qs = [t for t, x in enumerate(w) if self.isQ(x)]
        if len(w) == 1: res = {} if qs else {w: {(): 1}}
        elif not qs: res = {w: {(): 1}}
        else:
            res = {}
            for combo in itertools.product(*[list(self.T[w[t]].items()) for t in qs]):
                coef = 1; powers = {}
                for t, (j, c) in zip(qs, combo): coef = coef * c % MOD; powers[t] = j
                facs = self.factors(w, powers)
                pf = functools.reduce(self.mul, [{f: {(): 1}} for f in facs])
                pg = functools.reduce(self.mul, [self.g(f) for f in facs])
                diff = dict(pf)
                for ww, c in pg.items(): diff[ww] = padd(diff.get(ww, {}), c, -1)
                for ww, c in diff.items():
                    if c:
                        for k, v in self.P(ww).items(): res[k] = padd(res.get(k, {}), pscale(pmul(c, v), coef))
            res = {k: v for k, v in res.items() if v}
        self._P[w] = res; return res
    def Pel(self, E):
        out = {}
        for w, c in E.items():
            for k, v in self.P(w).items(): out[k] = padd(out.get(k, {}), pmul(c, v))
        return {k: v for k, v in out.items() if v}
def equations(Gs, words, depth=2):
    eqs = []; skipped = 0
    for combo in itertools.product(words, repeat=depth):
        try:
            r = Gs.Pel(functools.reduce(Gs.mul, [Gs.g(f) for f in combo]))
            eqs += [c for c in r.values() if c]
        except FD.OutOfWindow: skipped += 1
    return eqs, skipped
