"""Degrees of laws with a given c.e. support U (exact rational arithmetic on toy c.e. trees).
(a) Stage-weight law: every c.e. tree without dead ends is the positive set of a COMPUTABLE μ ∈ S.
    weight(v) = 2^{-s(v)} (s = enumeration stage), c_w = weight(w1)/(weight(w0)+weight(w1)).
    Lower bounds at stage s: ratio_s − 2^{m−s} (m = stage of the first child), clipped at 0, multiplied along w.
(b) Coding: with designated branching nodes d_0, d_1, … (c.e., infinitely many), put y_i at d_{2i} (scaled into [¼,¾])
    and the original c_{d_j} at d_{2j+1}; μ' ≡ μ ⊕ y, and the bounds of S survive."""
import random
from fractions import Fraction as Fr
from itertools import product
random.seed(7)
DEPTH = 9
def random_ce_tree():
    """a tree without dead ends up to DEPTH, with enumeration stages; some branches enumerated very late"""
    stage = {'': 0}; frontier = ['']
    for n in range(DEPTH):
        nf = []
        for w in frontier:
            kids = random.choice([['0'], ['1'], ['0', '1'], ['0', '1']])
            for b in kids:
                stage[w + b] = max(stage[w], 0) + random.choice([1, 2, 3, 5, 8, 13, 40])
                nf.append(w + b)
        frontier = nf
    return stage
def weight(stage, v, s=None):
    if v not in stage: return Fr(0)
    if s is not None and stage[v] > s: return Fr(0)
    return Fr(1, 2 ** stage[v])
def ratio(stage, w, b, s=None):
    x, y = weight(stage, w + b, s), weight(stage, w + ('1' if b == '0' else '0'), s)
    return x / (x + y) if x + y else Fr(0)
def mu(stage, w):
    p = Fr(1)
    for i in range(len(w)): p *= ratio(stage, w[:i], w[i])
    return p
def lower(stage, w, s):
    p = Fr(1)
    for i in range(len(w)):
        v = w[:i]; m = min((stage[v + b] for b in '01' if v + b in stage), default=None)
        if m is None or m > s: return Fr(0)
        p *= max(Fr(0), ratio(stage, v, w[i], s) - Fr(2 ** m, 2 ** s))
    return p
ok_pos = ok_bound = ok_low = ok_conv = True
for trial in range(6):
    st = random_ce_tree()
    words = [''.join(t) for n in range(DEPTH + 1) for t in product('01', repeat=n)]
    for w in words:
        m = mu(st, w)
        ok_pos &= ((m > 0) == (w in st))                                     # positive set = U exactly
        if w in st and w:
            ok_bound &= ratio(st, w[:-1], w[-1]) >= Fr(1, 2 ** (st[w] + 1))  # S-bound with ε = ½, L = s + 1
    for s in (10, 30, 60, 200):
        for w in words[:200]:
            lo = lower(st, w, s); ok_low &= lo <= mu(st, w)
    tot = sum(lower(st, ''.join(t), 1000) for t in product('01', repeat=DEPTH))
    ok_conv &= (1 - tot) < Fr(1, 2 ** 500)                                   # lower bounds exhaust the mass
print("(a) positive set is exactly U:", ok_pos)
print("(a) S-lower bound c ≥ 2^{-(s+1)} holds:", ok_bound)
print("(a) stagewise lower approximations never exceed μ:", ok_low)
print("(a) lower approximations at stage 1000 capture all but < 2^-500 of the mass at depth 9:", ok_conv)
# (b) coding
st = random_ce_tree()
branching = sorted([w for w in st if len(w) < DEPTH and w + '0' in st and w + '1' in st], key=lambda w: (max(st[w + '0'], st[w + '1']), w))
dn = branching                                                                   # designated nodes, in discovery order
y = [Fr(random.randint(0, 999), 999) for _ in range(len(dn) // 2 + 1)]
c = {w: ratio(st, w, '1') for w in st if len(w) < DEPTH}
c2 = dict(c)
for i, w in enumerate(dn):
    if i % 2 == 0: c2[w] = Fr(1, 4) + y[i // 2] / 2                             # y_i ∈ [0,1] ↦ [¼,¾]
    else: c2[w] = c[dn[(i - 1) // 2]]                                            # original of d_j stored at d_{2j+1}
dec_y = [ (c2[dn[2 * i]] - Fr(1, 4)) * 2 for i in range((len(dn) + 1) // 2)]
dec_c = {dn[j]: c2[dn[2 * j + 1]] for j in range(len(dn)) if 2 * j + 1 < len(dn)}
ok_dec = dec_y == y[:len(dec_y)] and all(dec_c[w] == c[w] for w in dec_c)
ok_rest = all(c2[w] == c[w] for w in c if w not in dn)
ok_int = all(0 < c2[w] < 1 for w in dn)
print(f"(b) {len(dn)} designated nodes: y and the displaced conditionals decode exactly: {ok_dec}; other nodes untouched: {ok_rest}; "
      f"designated conditionals stay in (0,1): {ok_int}")
