# The full-support box B(2^{<ω}, ¼) is computably homeomorphic to the Hilbert-cube box [¼,¾]^{2^{<ω}}:
#   x ↦ μ_x(w) = ∏_{v ⊏ w} (x_v or 1−x_v),   μ ↦ c_w = μ(w1)/μ(w)  (μ(w) ≥ 4^{-|w|}, so division has a computable modulus).
# Hence a computes μ iff a computes x, and Name(μ) = cone(μ) inherits Miller's non-total continuous degrees.
import random
from fractions import Fraction as Fr
from itertools import product
random.seed(12)
N = 7
words = [''.join(t) for n in range(N + 1) for t in product('01', repeat=n)]
ok_round = ok_mod = True
for trial in range(30):
    x = {w: Fr(random.randint(1, 3 * 64), 4 * 64) for w in words if len(w) < N}
    x = {w: min(max(v, Fr(1, 4)), Fr(3, 4)) for w, v in x.items()}
    mu = {'': Fr(1)}
    for w in words[1:]:
        mu[w] = mu[w[:-1]] * (x[w[:-1]] if w[-1] == '1' else 1 - x[w[:-1]])
    ok_round &= all(mu[w + '1'] / mu[w] == x[w] for w in x)
    ok_round &= all(mu[w] == mu[w + '0'] + mu[w + '1'] for w in x)           # additivity
    # modulus: perturb every μ-value by at most η; the recovered c_w moves by at most 2η·4^{|w|}/(1 − η·4^{|w|})
    eta = Fr(1, 10 ** 9)
    pert = {w: mu[w] + Fr(random.randint(-1000, 1000), 1000) * eta for w in mu}
    for w in x:
        b = 4 ** len(w)
        ok_mod &= abs(pert[w + '1'] / pert[w] - x[w]) <= 2 * eta * b / (1 - eta * b)
print("round trip x -> mu -> x exact, mu additive:", ok_round)
print("inverse map has the computable modulus 2η·4^|w|/(1−η·4^|w|):", ok_mod)
