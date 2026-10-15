"""Order-by-order (in ε) integration of a first-order fold deformation, a-degree <= 2.
Coefficients: dict (order, var|None) -> int mod p, truncated at order K; the unknown x_K appears only at order K."""
import fold as FD
from fold import MOD, Echelon, md
from tangent import spec
K_ORDER = [1]
LAST = [None]
def padd(a, b, s=1):
    r = dict(a)
    for k, v in b.items():
        r[k] = (r.get(k, 0) + s * v) % MOD
        if r[k] == 0: del r[k]
    return r
def pmul(a, b):
    K = K_ORDER[0]; r = {}
    for (o1, v1), x1 in a.items():
        for (o2, v2), x2 in b.items():
            o = o1 + o2
            if o > K or (v1 is not None and v2 is not None): continue
            v = v1 if v1 is not None else v2
            r[(o, v)] = (r.get((o, v), 0) + x1 * x2) % MOD
    return {k: v for k, v in r.items() if v}
def pscale(a, s): return {k: v * s % MOD for k, v in a.items() if v * s % MOD}
FD.padd, FD.pmul, FD.pscale = padd, pmul, pscale
FD.ONE = {(0, None): 1}
class SFold(FD.Fold):
    def __init__(self, S, Lw, D, known):
        super().__init__(S, Lw, D); self.known = known
    def _val(self, key):
        K = K_ORDER[0]; val = {(o, None): x for o, x in self.known.get(key, {}).items() if o < K and x}
        if K == 1 and key[0] == 'c':      # order 1: the a-degree-1 part is the given direction
            return val if val else {(1, None): 0} if False else {(o, None): x for o, x in self.known.get(key, {}).items() if x}
        if key not in self.var: self.var[key] = len(self.names); self.names.append(key)
        val[(K, self.var[key])] = 1
        return val
    def c(self, w):
        if self.wlev(w) > self.Lw: raise FD.OutOfWindow
        return {o: self._val(('c', w, o)) for o in self.outs(w)}
    def c2(self, w):
        if self.wlev(w) > self.Lw: raise FD.OutOfWindow
        return {o: self._val(('c2', w, o)) for o in self.outs(w)}
CONST = 10**12
def integrate(S, Lw, D, tau, N, Leq=None):
    Leq = Lw if Leq is None else Leq
    """tau: dict key -> value (order-1 term, keys ('c',w,o)). Returns highest order reached (<= N) and status."""
    known = {k: {1: md(v)} for k, v in tau.items()}
    for K in range(1, N + 1):
        K_ORDER[0] = K
        F = SFold(S, Lw, D, known); eqs = FD.equations(F, Leq)
        if FD.SKIPS[0]: return K, f'window too small ({FD.SKIPS[0]} equations skipped)'
        E = Echelon(); bad = False; lower = 0
        for e in eqs:
            row = {}
            for (o, v), x in e.items():
                if o < K: lower += 1
                elif o == K: row[CONST if v is None else v] = x
            r = E.reduce(row)
            if r and min(r) == CONST: bad = True; break
            if r: E.add(r)
        if lower: return K, f"inconsistent lower-order terms ({lower})"
        if bad: return K, "obstructed" if K > 1 else "not a tangent vector"
        E.full_reduce()
        for c, r in E.piv.items():
            if c != CONST and CONST in r:
                known.setdefault(F.names[c], {})[K] = (-r[CONST]) % MOD
    LAST[0] = known
    return N, "ok"
