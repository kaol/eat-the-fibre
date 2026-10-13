"""solve_affine ≡ brute force: the solution set of E κ = r over F_p, for random small systems (hypothesis)."""
import itertools, os, sys
import numpy as np, torch
from hypothesis import given, settings, strategies as st
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import fibre_torch as ft


def solutions(cons, k0, N, free, p):
    if not cons: return set()
    fs = [f for f in range(len(free)) if free[f]]
    return {tuple(int(x) for x in (k0 + sum((c * N[f] for c, f in zip(cs, fs)), np.zeros_like(k0))) % p)
            for cs in itertools.product(range(p), repeat=len(fs))}


@settings(max_examples=300, deadline=None)
@given(st.sampled_from([2, 3, 5, 7]), st.integers(1, 4), st.integers(0, 7), st.integers(0, 2**32 - 1),
       st.floats(0, 1))
def test_solve_affine(p, k, R, seed, density):
    g = np.random.default_rng(seed)
    n = 3
    E = g.integers(0, p, (n, R, k)) * (g.random((n, R, k)) < density)
    if R and g.random() < 0.5: E[:, -1] = (2 * E[:, 0]) % p                 # dependent rows
    x = g.integers(0, p, (n, k))
    r = np.einsum('nrk,nk->nr', E, x) % p
    if g.random() < 0.3: r = (r + (g.random((n, R)) < 0.3)) % p             # sometimes inconsistent
    F = ft.Fp(p, torch.device('cpu'))
    Ab = torch.as_tensor(np.concatenate([E, r[:, :, None]], 2), dtype=torch.float32)
    cons, k0, N, free = (t.numpy() for t in ft.solve_affine(F, Ab))
    for i in range(n):
        want = {kap for kap in itertools.product(range(p), repeat=k)
                if ((E[i] @ np.array(kap) - r[i]) % p == 0).all()}
        got = solutions(cons[i], k0[i].astype(np.int64), N[i].astype(np.int64), free[i], p)
        assert got == want
