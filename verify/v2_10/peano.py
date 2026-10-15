import numpy as np
def hilbert_d2xy(order, d):
    n = 2 ** order; x = y = 0; t = d; s = 1
    while s < n:
        rx = 1 & (t // 2); ry = 1 & (t ^ rx)
        if ry == 0:
            if rx == 1: x, y = s - 1 - x, s - 1 - y
            x, y = y, x
        x += s * rx; y += s * ry; t //= 4; s *= 2
    return x, y
ORDER = 11; N = 4 ** ORDER; n = 2 ** ORDER
X = np.empty(N); Y = np.empty(N)
for d in range(N):
    x, y = hilbert_d2xy(ORDER, d); X[d] = (x + 0.5) / n; Y[d] = (y + 0.5) / n
def H(t):   # continuous Hilbert curve [0,1] -> [0,1]^2 (piecewise-constant approximation at this order)
    i = np.minimum((t * N).astype(np.int64), N - 1); return X[i], Y[i]
ts = (np.arange(N) + 0.5) / N
c1, r = H(ts); c2, r2 = H(r); c3, _ = H(r2)
for k, coords in ((2, (c1, c2)), (3, (c1, c2, c3))):
    res = 8 if k == 2 else 4
    cells = set(tuple(np.minimum((np.array(cc) * res).astype(int), res - 1)) for cc in zip(*coords))
    print(f"projection to first {k} coordinates of x ↦ (H₁x, H₁H₂x, H₁H₂H₂x, …): covers {len(cells)}/{res**k} cells of a {res}^{k} grid")
