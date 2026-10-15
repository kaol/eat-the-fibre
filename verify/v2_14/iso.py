"""k⟨a,u⟩ with a-deletion Δ  ≅  path algebra of the one-vertex 2-loop quiver k⟨α,β⟩ with arrow deletion,
via φ(a) = α, φ(u) = α − β.  Exact check on all words of length ≤ N: φ is multiplicative by construction, bijective
(unitriangular in the letter-count filtration), and Δφ = (φ⊗φ)Δ."""
import itertools
from collections import defaultdict
N = 6
def phi(word):                                   # word over 'au' -> dict over words in 'xy' (x=α, y=β)
    out = {'': 1}
    for ch in word:
        nxt = defaultdict(int)
        for w, c in out.items():
            if ch == 'a': nxt[w + 'x'] += c
            else: nxt[w + 'x'] += c; nxt[w + 'y'] -= c
        out = {w: c for w, c in nxt.items() if c}
    return out
def delta_a(word):                               # delete one 'a'
    out = defaultdict(int)
    for i, ch in enumerate(word):
        if ch == 'a': out[(word[:i], word[i + 1:])] += 1
    return out
def delta_xy(word):                              # delete any letter
    out = defaultdict(int)
    for i in range(len(word)): out[(word[:i], word[i + 1:])] += 1
    return out
ok = True
for n in range(N + 1):
    for t in itertools.product('au', repeat=n):
        w = ''.join(t)
        lhs = defaultdict(int)
        for v, c in phi(w).items():
            for k, e in delta_xy(v).items(): lhs[k] += c * e
        rhs = defaultdict(int)
        for (x, y), e in delta_a(w).items():
            for vx, cx in phi(x).items():
                for vy, cy in phi(y).items(): rhs[(vx, vy)] += e * cx * cy
        ok &= {k: v for k, v in lhs.items() if v} == {k: v for k, v in rhs.items() if v}
# bijectivity: the image of the 2^n words of length n spans the 2^n words in x,y (determinant ±1)
from sympy import Matrix
bij = True
for n in range(1, 6):
    words = [''.join(t) for t in itertools.product('au', repeat=n)]; targets = [''.join(t) for t in itertools.product('xy', repeat=n)]
    M = Matrix([[phi(w).get(v, 0) for v in targets] for w in words]); bij &= abs(M.det()) == 1
print(f"Δφ = (φ⊗φ)Δ on all words of length ≤ {N}: {ok}; φ bijective in each length ≤ 5 (det ±1): {bij}")
