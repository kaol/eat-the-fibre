"""q123s.py with a cheaper stage check. In a batch only the new column varies, so products of two old columns are computed
once (unbatched); only pairs involving the new column are batched. Same masked conditions as q123v.stage_ok."""
import sys, itertools, time, collections
import numpy as np
src = open("q123s.py").read()
head, tail = src.split("\nt0 = time.time(); done")
exec(head)
def stage_ok_new(Pb_all, upto, new):
    idx = np.arange(Pb_all.shape[0]); P0 = Pb_all[0]
    def shrink(ok):
        nonlocal idx, Pb
        idx = idx[ok]; Pb = Pb[ok]
    Pb = Pb_all
    unknown = lambda v: v[..., upto + 1:].any(axis=-1)
    cols = range(upto + 1)
    for i in cols:                                   # idempotency P v = v (masked if v needs unknown columns)
        v = Pb[:, :, i] if i == new else P0[:, i]
        if v.ndim == 1:
            if unknown(v): continue
            shrink(~((np.einsum('nij,j->ni', Pb, v) - v) % p).any(axis=1))
        else:
            shrink(~((((np.einsum('nij,nj->ni', Pb, v) - v) % p).any(axis=1)) & ~unknown(v)))
        if not len(idx): return idx
    for i in cols:
        for j in cols:
            for comp in (False, True):               # comp: complement Q = (1−P)A is a subalgebra: P q = 0
                xi = (I[:, i] - (Pb[:, :, i] if i == new else P0[:, i])) % p if comp else (Pb[:, :, i] if i == new else P0[:, i])
                yj = (I[:, j] - (Pb[:, :, j] if j == new else P0[:, j])) % p if comp else (Pb[:, :, j] if j == new else P0[:, j])
                if xi.ndim == 1 and yj.ndim == 1: h = mul(xi, yj)
                elif xi.ndim == 1: h = np.einsum('i,nj,ijk->nk', xi, yj, M) % p
                elif yj.ndim == 1: h = np.einsum('ni,j,ijk->nk', xi, yj, M) % p
                else: h = np.einsum('ni,nj,ijk->nk', xi, yj, M) % p
                if h.ndim == 1:
                    if unknown(h): continue
                    r = np.einsum('nij,j->ni', Pb, h) % p
                else:
                    r = np.einsum('nij,nj->ni', Pb, h) % p
                bad = (r if comp else (r - h) % p).any(axis=1)
                if h.ndim == 2: bad &= ~unknown(h)
                shrink(~bad)
                if not len(idx): return idx
    return idx
def extend(P, col, L, upto, extra=None):
    Pb = np.repeat(P[None], len(allk), axis=0)
    Pb[:, :, col] = (L + allk) % p
    if extra is not None: extra(Pb)
    # with an extra hook (α2 = α1 − δ) two columns vary: fall back to the full check
    return Pb[stage_ok(Pb, list(range(upto + 1)))] if extra is not None else Pb[stage_ok_new(Pb, upto, col)]
exec("t0 = time.time(); done" + tail)
