"""v2.42: v2.34's infinitesimal invariants (b, e) on the representatives of given signature groups (merged labels).
b = dim T_P Folds − dim([Der, P] ∩ T): tangent directions not produced by derivations (splitting deformations);
e = dim([Der, P] ∩ T) − dim(census symmetry directions): derivation directions outside the census group.
Usage: python3 analysis/v2_42_be.py REPS.pkl SIGNATURE_SUBSTRING..."""
import os, sys, io, contextlib, collections, pickle
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
rf = sys.argv[1]; targets = sys.argv[2:]
sys.argv = ['family_curves.py', rf]
g = {'__name__': 'fc', '__file__': os.path.join(HERE, 'family_curves.py')}
with contextlib.redirect_stdout(io.StringIO()):
    exec(open(os.path.join(HERE, 'family_curves.py')).read().split("fams = []; nlift")[0], g)
sys.argv = ['fibre_rigidity.py']
fr = {'__name__': 'fr', '__file__': os.path.join(HERE, 'fibre_rigidity.py')}
exec(open(os.path.join(HERE, 'fibre_rigidity.py')).read().split("\ndef main():")[0], fr)
r = g['r']; p = r['p']; sig = g['sig_of_loc']
C = fr['setup'](r['name'], p)
print(f"q={p}: dim Der = {len(C['Der'])}, dim Der_Δ = {len(C['DerD'])}")
for tg in targets:
    locs = {l for l, s in sig.items() if tg in str(s)}
    h = collections.Counter(); hl = collections.defaultdict(set)
    for i, l in enumerate(r['loc']):
        if l in locs:
            a = fr['at'](r['greps'][i], C, p); h[(a['b'], a['e'])] += 1; hl[(a['b'], a['e'])].add(l)
    print(f"group …{tg[-40:]}: {len(locs)} orbits; (b, e) over G-reps: {dict(h)}; orbits per (b, e): "
          f"{ {k: len(v) for k, v in hl.items()} }", flush=True)
