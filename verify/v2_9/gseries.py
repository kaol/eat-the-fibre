"""Order-by-order integration in C-coordinates, ALL a-degrees (window = generator words / variable levels)."""
import itertools, functools, random, sys, time
import fold as FD
from fold import MOD
import gen as G0
K_ORDER = [1]
DCORR = [3]
LAST = [None]
def sadd(a, b, s=1):
    r = dict(a)
    for k, v in b.items():
        r[k] = (r.get(k, 0) + s * v) % MOD
        if r[k] == 0: del r[k]
    return r
def smul(a, b):
    K = K_ORDER[0]; r = {}
    for (o1, v1), x1 in a.items():
        for (o2, v2), x2 in b.items():
            o = o1 + o2
            if o > K or (v1 is not None and v2 is not None): continue
            v = v1 if v1 is not None else v2
            r[(o, v)] = (r.get((o, v), 0) + x1 * x2) % MOD
    return {k: v for k, v in r.items() if v}
def sscale(a, s): return {k: v * s % MOD for k, v in a.items() if v * s % MOD}
ONE = {(0, None): 1}
class GSer(G0.Gen):
    def __init__(self, S, delta, Lvars, known, fixed_first=None, Llow=99):
        super().__init__(S, {}); self.delta, self.Lvars, self.known = delta, Lvars, known
        self.var = {}; self.names = []; self.fixed_first = fixed_first; self.Llow = Llow
    def Cval(self, m):
        lev = sum(self.S.lev[x] for x in m)
        if lev > self.Lvars: raise FD.OutOfWindow
        K = K_ORDER[0]; out = {}
        for b in range(len(self.S.B)):
            if self.S.side[b] == 'H' and self.S.lev[b] <= lev + (self.delta if K_ORDER[0] == 1 else DCORR[0]):
                key = (m, b); val = {(o, None): x for o, x in self.known.get(key, {}).items() if o < K and x}
                if not (K == 1 and self.fixed_first is not None and lev <= self.Llow):
                    if key not in self.var: self.var[key] = len(self.names); self.names.append(key)
                    val[(K, self.var[key])] = 1
                else:
                    x = self.known.get(key, {}).get(1, 0)
                    if x: val[(1, None)] = x
                if val: out[b] = val
        return out
    def mul(self, A, B):
        out = {}
        for w1, c1 in A.items():
            for w2, c2 in B.items():
                cc = smul(c1, c2)
                if not cc: continue
                for w, m in self.mulw(w1, w2).items(): out[w] = sadd(out.get(w, {}), sscale(cc, m))
        return {w: c for w, c in out.items() if c}
    def kappa(self, f):
        i = [t for t, x in enumerate(f) if self.isQ(x)][0]; n = len(f) - 1; out = {}
        for j in [-1] + list(range(i)):
            for k in [n] + list(range(i, n)):
                m = f[j + 1:k + 1]
                if len(m) < 2: continue
                for b, v in self.Cval(m).items():
                    w = f[:j + 1] + (b,) + f[k + 1:]; out[w] = sadd(out.get(w, {}), v)
        return {w: c for w, c in out.items() if c}
    def g(self, f):
        out = {f: dict(ONE)}
        for w, c in self.kappa(f).items(): out[w] = sadd(out.get(w, {}), c, -1)
        return {w: c for w, c in out.items() if c}
    def P(self, w):
        if w in self._P: return self._P[w]
        qs = [t for t, x in enumerate(w) if self.isQ(x)]
        if len(w) == 1: res = {} if qs else {w: dict(ONE)}
        elif not qs: res = {w: dict(ONE)}
        else:
            res = {}
            for combo in itertools.product(*[list(self.T[w[t]].items()) for t in qs]):
                coef = 1; powers = {}
                for t, (j, c) in zip(qs, combo): coef = coef * c % MOD; powers[t] = j
                facs = self.factors(w, powers)
                pf = functools.reduce(self.mul, [{f: dict(ONE)} for f in facs])
                pg = functools.reduce(self.mul, [self.g(f) for f in facs])
                diff = dict(pf)
                for ww, c in pg.items(): diff[ww] = sadd(diff.get(ww, {}), c, -1)
                for ww, c in diff.items():
                    if c:
                        for k, v in self.P(ww).items(): res[k] = sadd(res.get(k, {}), sscale(smul(c, v), coef))
            res = {k: v for k, v in res.items() if v}
        self._P[w] = res; return res
    def Pel(self, E):
        out = {}
        for w, c in E.items():
            for k, v in self.P(w).items(): out[k] = sadd(out.get(k, {}), smul(c, v))
        return {k: v for k, v in out.items() if v}
CONST = 10**12
def integrate(S, words, delta, Lvars, tau, N, depth=2, Llow=99):
    """tau: dict (m, b) -> value: first-order C. Returns (order reached, status)."""
    known = {k: {1: v % MOD} for k, v in tau.items()}
    for K in range(1, N + 1):
        K_ORDER[0] = K
        Gs = GSer(S, delta, Lvars, known, fixed_first=(K == 1) or None, Llow=Llow)
        E = FD.Echelon(); skipped = 0; bad = False
        for combo in itertools.product(words, repeat=depth):
            try: r = Gs.Pel(functools.reduce(Gs.mul, [Gs.g(f) for f in combo]))
            except FD.OutOfWindow: skipped += 1; continue
            for c in r.values():
                row = {}
                for (o, v), x in c.items():
                    if o < K and x: return K, f"inconsistent lower order"
                    if o == K: row[CONST if v is None else v] = x
                rr = E.reduce(row)
                if rr and min(rr) == CONST: bad = True; break
                if rr: E.add(rr)
            if bad: break
        if skipped: return K, f"window too small ({skipped} skipped)"
        if bad: return K, "obstructed"
        E.full_reduce()
        for c, r in E.piv.items():
            if c != CONST and CONST in r: known.setdefault(Gs.names[c], {})[K] = (-r[CONST]) % MOD
    LAST[0] = known
    return N, "ok"
