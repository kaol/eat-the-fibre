"""symmetry_plan (v2.28: streaming over H, orbits from generators) against symmetry_plan_ref: same sources, same
source for every transported rep, same N; and every transport element h is valid: h·PK_source·h⁻¹ = PK_r on K,
Rh = Binv h⁻¹ B. (The h itself may differ from the reference's choice; the fibre it produces is the same set.)
Usage: python3 tests/test_symmetry_plan.py [QUIVER:p ...]"""
import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); os.chdir(ROOT); sys.path.insert(0, ROOT)
ARGS = sys.argv[1:] or ['pairs_in:3', 'two_pairs_short:3', 'pairs_in_tail:3', 'two_pairs:3', 'K3_after:3', 'triple_pair_in_short:3']
sys.argv = ['x']; __file__ = os.path.join(ROOT, 'pcensus.py')
exec(open(__file__).read().split("if __name__ == '__main__':")[0])
import numpy as np
sys.path.insert(0, ROOT)
ok_all = True
for qp in ARGS:
    name, p = qp.split(':'); p = int(p)
    ctx = setup(name, p); C.update(ctx); phis = grading(ctx)
    Kf, reps = kernel_stage(ctx, phis, 'auto', 'cpu', None)
    s1, t1, N1 = symmetry_plan(ctx, reps, phis)
    s2, t2, N2 = symmetry_plan_ref(ctx, reps, phis)
    Kinv = ctx['Binv'][:ctx['kdim']]; Kmat = ctx['Kmat'] % p
    def valid(r, s_, h, Rh):
        hK = Kinv @ h @ Kmat % p; hi = ctx['B'] @ Rh @ ctx['Binv'] % p                 # Rh = Binv h⁻¹ B
        hKi = Kinv @ hi @ Kmat % p
        return ((h @ hi) % p == np.eye(ctx['d'], dtype=np.int64)).all() and \
               ((hK @ reps[s_][0] @ hKi - reps[r][0]) % p == 0).all()
    ok = s1 == s2 and t1.keys() == t2.keys() and all(t1[r][0] == t2[r][0] and valid(r, *t1[r]) for r in t1)
    ok &= N1[0].shape == N2[0].shape and (N1[0] == N2[0]).all() and (N1[1] == N2[1]).all()
    ok_all &= ok
    print(f"{name} F_{p}: {len(reps)} reps, {len(s1)} sources, {len(t1)} transported, |N| = {len(N1[0])}: identical = {ok}", flush=True)
print("ALL IDENTICAL" if ok_all else "MISMATCH")
