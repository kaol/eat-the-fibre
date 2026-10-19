"""Audit (v2.28): for every quiver in the census table, dim Lie(Aut_Δ) (all Δ-derivations) against dim Lie G (the
census group's substitutions). Where they differ, census orbit counts are counts for a group smaller than Aut_Δ°.
Usage: python3 analysis/aut_audit.py [max basis size, default 14]"""
import os as _os, sys as _sys, runpy
ARGS = list(_sys.argv)
AN = _os.path.dirname(_os.path.abspath(__file__))
T3 = runpy.run_path(_os.path.join(AN, 'theorem3.py'), run_name='t3')
dmax = int(ARGS[1]) if len(ARGS) > 1 else 14
for name, (nv, arrows) in sorted(T3['QUIVERS'].items()):
    basis, M, DL = T3['build'](nv, arrows)
    if len(basis) > dmax: print(f"{name:22s} d = {len(basis):3d}: skipped"); continue
    n, g, k = T3['lie_aut_delta'](name)
    print(f"{name:22s} d = {len(basis):3d}: dim Lie(Aut_Δ) = {n}, dim Lie G = {g}, vertices fixed: {k}" + ("   <-- larger" if n != g else ""), flush=True)
