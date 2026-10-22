"""v2.42: vertex blocks εPι (symmetric residues) of the merged orbits in given signature groups, and which blocks
grow from one prime to the next. The vertex block is invariant under G, class weightings and the shear
c ↦ c + s·δ_aδ_b (all fix vertex idempotents mod rad and preserve rad).
Usage: python3 analysis/v2_42_vblock.py REPS_q1.pkl REPS_q2.pkl SIGNATURE_SUBSTRING..."""
import os, sys, io, contextlib, collections
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
files = sys.argv[1:3]; targets = sys.argv[3:]
out = {}
for f in files:
    sys.argv = ['family_curves.py', f]
    g = {'__name__': 'fc', '__file__': os.path.join(HERE, 'family_curves.py')}
    with contextlib.redirect_stdout(io.StringIO()):
        exec(open(os.path.join(HERE, 'family_curves.py')).read().split("fams = []; nlift")[0], g)
    r, sig, p = g['r'], g['sig_of_loc'], g['p']
    nv = sum(1 for b in g['basis'] if b[0] == 'v')
    vb = {}
    for i, l in enumerate(r['loc']):
        s = sig[l]
        if not any(tg in str(s) for tg in targets): continue
        P = r['greps'][i].astype(np.int64)
        B = tuple(int(x) for x in ((P[:nv, :nv] % p + p // 2) % p - p // 2).flatten())
        if l in vb: assert vb[l] == B, "vertex block not constant on an orbit"
        vb[l] = B
    out[p] = collections.Counter(vb.values())
    print(f"q={p}: {len(vb)} orbits, {len(out[p])} vertex blocks")
q1, q2 = sorted(out)
for B in sorted(set(out[q1]) | set(out[q2])):
    if out[q1][B] != out[q2][B]: print("  grows:", B, out[q1][B], "->", out[q2][B])
