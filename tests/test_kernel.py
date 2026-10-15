"""Kernel stage peers: kernel_torch.kernel_stage ≡ (pcensus.kernel_folds, pcensus.torus_reps) as LISTS (same folds, same
order, same reps and l-lists), since fibre indices in checkpoints and representative-dependent fields depend on the order.
Plus the exact Q-rank (rank_Q ≡ sympy) and the kernel cache (hit, and digest mismatch → recompute)."""
import os, sys
import numpy as np, torch
import pytest
from hypothesis import given, settings, strategies as st
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pcensus as pc

CASES = [('K2', 3), ('K2', 5), ('pair_after', 3), ('pair_srcout', 5), ('pair_tgtfork', 3), ('K3_after', 3), ('two_pairs', 3)]


def same(A, B):
    (Ka, ra), (Kb, rb) = A, B
    return (np.asarray(Ka).shape == np.asarray(Kb).shape and (np.asarray(Ka) == np.asarray(Kb)).all() and len(ra) == len(rb)
            and all((x[0] == y[0]).all() and x[1] == y[1] for x, y in zip(ra, rb)))


@pytest.mark.parametrize('q,p', CASES)
def test_kernel_torch_equals_ref(q, p):
    ctx = pc.setup(q, p); phis = pc.grading(ctx)
    assert same(pc.kernel_stage(ctx, phis, 'torch'), pc.kernel_stage(ctx, phis, 'ref'))


def test_kernel_cache(tmp_path):
    ctx = pc.setup('pair_after', 3); phis = pc.grading(ctx); f = str(tmp_path / 'k.npz')
    a = pc.kernel_stage(ctx, phis, 'torch', cache=f)
    assert os.path.exists(f)
    assert same(a, pc.kernel_stage(ctx, phis, 'ref', cache=f))              # served from the cache
    other = pc.setup('K2', 3)                                                # same file, other seed: digest mismatch
    assert same(pc.kernel_stage(other, pc.grading(other), 'torch', cache=f), pc.kernel_stage(other, pc.grading(other), 'ref'))


@settings(max_examples=300, deadline=None)
@given(st.integers(0, 5), st.integers(1, 7), st.integers(0, 2**32 - 1))
def test_rank_Q(R, n, seed):
    import sympy as sp
    g = np.random.default_rng(seed)
    A = g.integers(-3, 4, (R, n)) * (g.random((R, n)) < 0.6)
    if R > 2: A[-1] = 2 * A[0] - 3 * A[1]                                    # a dependent row
    want = sp.Matrix(A.tolist()).rank() if R else 0
    assert pc.rank_Q(A.tolist()) == want
