"""Parameter complexity: deadline functions, exact arithmetic on toy halting tables.
Node n: child 0 always in U; child 1 in U iff machine n halts at time T_n.  Branch mass c_n = mu(w1)/mu(w)."""
from fractions import Fraction as Fr
import random
rnd = random.Random(4)
T = {n: (None if n % 4 == 1 else rnd.randint(1, 30) * (n + 2)) for n in range(16)}
# ---------- necessity: extract d from a NAMING PROGRAM (fair-coin halting program: budget s ~ Geometric(1/2), emit fair bit if halted by s)
def program_histories(n, maxlen):
    """histories at node n as (prob, emitted symbol, enumeration stage at which the history is found)."""
    out = []
    for s in range(maxlen):                     # budget s with prob 2^{-(s+1)}
        pb = Fr(1, 2 ** (s + 1))
        if T[n] is not None and s >= T[n]:
            out.append((pb / 2, 1, s)); out.append((pb / 2, 0, s))
        else:
            out.append((pb, 0, s))
    return out
def mu1(n):  # exact branch mass of child 1 for the fair-coin program
    return Fr(0) if T[n] is None else Fr(1, 2 ** (T[n] + 1))
def d_from_program(n, k, maxlen=400):
    """explore histories in order until unexplored mass < 2^{-k-1}; return the stage reached."""
    explored = Fr(0)
    for (pr, sym, stage) in program_histories(n, maxlen):
        explored += pr
        if 1 - explored < Fr(1, 2 ** (k + 1)): return stage
    return maxlen
ok_nec = True
for n in T:
    for k in range(0, 60):
        if mu1(n) >= Fr(1, 2 ** k):                        # branch heavy at scale k ...
            ok_nec &= T[n] is not None and T[n] <= d_from_program(n, k)   # ... => certified by stage d(n,k)
print("necessity: the deadline extracted from the naming program certifies every heavy branch:", ok_nec)
# ---------- sufficiency: piecewise-halving survival built from a deadline function d
def run_node(n, c, d, cap=6000):
    """give-up schedule: survival rho(t) = 1 - delta/2 up to D_1, then 2^{-(k-1)} on (D_k, D_{k+1}]  (D_k = d(n,k))."""
    delta = Fr(1, 20)
    D = [d(n, k) for k in range(1, 40)]
    def rho_target(t):
        if t <= D[0]: return 1 - delta / 2
        for k in range(1, len(D)):
            if t <= D[k]: return Fr(1, 2 ** (k - 1))
        return Fr(1, 2 ** (len(D)))
    alive = Fr(1); p1 = Fr(0)
    for t in range(cap):
        if T[n] is not None and t == T[n]:
            ratio = c / alive
            if ratio > 1: return None                    # cannot hit the target
            return alive * ratio
        nxt = rho_target(t + 1); q = 1 - nxt / alive if alive else Fr(0)
        alive = nxt
    return p1
targets = {n: (Fr(rnd.randint(1, 18), 20) if T[n] is not None else Fr(0)) for n in T}   # halting-type: c bounded below
dgood = lambda n, k: (T[n] if T[n] is not None else 0) + k                                # dominates verification time
dbad = lambda n, k: 5 * k                                                                  # computable, does not dominate
good = all(run_node(n, targets[n], dgood) == targets[n] for n in T if T[n] is not None)
bad = [n for n in T if T[n] is not None and run_node(n, targets[n], dbad) != targets[n]]
print("sufficiency: with a dominating deadline function every target branch mass is hit exactly:", good)
print(f"           with a non-dominating computable d, targets missed on {len(bad)} halting nodes")
# fair-coin measure: masses 2^{-(T+1)} are hit with the COMPUTABLE deadline d(n,k) = k (no oracle needed)
dcomp = lambda n, k: k
fc = all(run_node(n, mu1(n), dcomp) == mu1(n) for n in T if T[n] is not None)
print("the fair-coin halting measure is hit exactly with the computable deadline d(n,k) = k:", fc)
# ---------- the deadline is computable from mu itself: compare branch mass to 3·2^{-k-2}, search only when positive
def approx(x, k):   # an approximation of x within 2^{-k-2} (here: truncation of the exact value, as an oracle would give)
    s = Fr(1, 2 ** (k + 2)); return Fr(int(x / s)) * s
def d_from_mu(c_of, n, k):
    if approx(c_of(n), k) >= Fr(3, 2 ** (k + 2)):
        assert T[n] is not None                  # positive mass => in U => the enumeration search halts
        return T[n]                              # stage at which the branch is enumerated
    return 0
for nm, c_of in (("halting-type targets", lambda n: targets[n]), ("fair-coin masses", mu1)):
    ok = all((T[n] is not None and T[n] <= d_from_mu(c_of, n, k)) for n in T for k in range(0, 60) if c_of(n) >= Fr(1, 2 ** k))
    hit = all(run_node(n, c_of(n), lambda m, k: d_from_mu(c_of, m, k)) == c_of(n) for n in T if T[n] is not None)
    print(f"deadline computed from mu alone ({nm}): certifies every heavy branch: {ok}; reproduces every branch mass: {hit}")
