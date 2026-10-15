"""Point count N_m(q) of the Kronecker fold variety predicted by the hand classification (14m − 8 orbit types)."""
import sys
from sympy import symbols, expand, Poly
def gb(n, k, q):
    if k < 0 or k > n: return 0
    num = 1; den = 1
    for i in range(k): num *= (q ** (n - i) - 1); den *= (q ** (i + 1) - 1)
    return num // den if isinstance(q, int) else expand(num / den)
def N(m, q):
    caseI = sum(gb(m - 1, m - k, q) * q ** (k * (m - k)) * (2 * q ** (m - k) - 1) for k in range(1, m + 1))
    IIaA0 = sum(gb(m - 1, j, q) * q ** (j * (m - j)) * (2 * q ** j - 1) for j in range(0, m))
    Ske2 = sum(gb(m - 1, j, q) * q ** (j * (m - j)) * q ** (m - 1 - j) for j in range(0, m))
    Sk1 = sum(gb(m - 1, j, q) * q ** (j * (m - j)) for j in range(0, m))
    SA0 = sum(gb(m - 1, j, q) * (q ** (m - 1 - j) - 1) * q ** ((j + 1) * (m - 1 - j)) for j in range(0, m - 1))
    return caseI + IIaA0 + 2 * (Ske2 + Sk1 + SA0)
if __name__ == '__main__':
    for m in range(1, 6):
        print(f"m={m}: N(3)={N(m, 3)}, N(5)={N(m, 5)}, N(7)={N(m, 7)}")
