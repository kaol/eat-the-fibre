"""Toy check of the verification-with-budget construction for undecidable supports (exact rational arithmetic).
Node n: child 0 always in U; child 1 in U iff 'machine n halts', at time T_n (toy table; None = never).
Target conditional probability of child 1 at node n: c_n in [delta, 1 - delta]."""
from fractions import Fraction as Fr
import random
rnd = random.Random(1)
T = {n: (None if n % 3 == 1 else rnd.randint(1, 40) * (n + 1)) for n in range(12)}       # halting times
delta = Fr(1, 20)
c = {n: Fr(rnd.randint(1, 18), 20) for n in T}                                           # targets in [1/20, 18/20]
def run_node(n, schedule, cap_steps=4000):
    """exact P(emit 1), P(emit 0), P(still running after cap) for the budgeted verification at node n."""
    alive = Fr(1); p1 = Fr(0); p0 = Fr(0)
    for t in range(cap_steps):
        if T[n] is not None and t == T[n]:                     # verification succeeds now
            rho = alive
            ratio = c[n] / rho if c[n] <= rho else Fr(1)
            p1 += alive * ratio; p0 += alive * (1 - ratio); alive = Fr(0); break
        q = schedule(n, t)
        p0 += alive * q; alive *= (1 - q)                     # give up -> emit the known child 0
    return p1, p0, alive
# (i) encoded deadline T*_n (from the real parameter): tiny q before the deadline, 1/2 after
Tstar = {n: (T[n] if T[n] is not None else 0) for n in T}
enc = lambda n, t: (delta / 2) * Fr(1, 2 ** (t + 1)) if t < Tstar[n] else Fr(1, 2)
ok_enc = True
for n in T:
    p1, p0, alive = run_node(n, enc)
    want = c[n] if T[n] is not None else Fr(0)
    ok_enc &= (p1 == want) and (alive < Fr(1, 2 ** 100))
print("encoded deadlines: P(emit 1) = c_n exactly on halting nodes, 0 elsewhere, and the run stops (mass left < 2^-100):", ok_enc)
# (ii) any FIXED computable schedule (no access to T*): productivity forces prod(1 - q_t) -> 0, which caps P(emit 1) at rho_{T_n}
fixed = lambda n, t: Fr(1, t + 2)
caps = {n: run_node(n, fixed)[0] for n in T if T[n] is not None}
fails = [n for n in caps if caps[n] < c[n]]
print(f"fixed schedule q_t = 1/(t+2): achieved P(emit 1) < target on {len(fails)}/{len(caps)} halting nodes (cap = 1/(T_n+1))")
