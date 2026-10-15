"""Key lemma for language-level Borelness, checked on random finite-state coin programs.
State q flips a coin (bias p or 1/2); heads/tails each lead to a state and emit 0, 1 or nothing."""
import random, itertools
from fractions import Fraction as Fr
import sympy as sp
def rand_prog(rnd, s):
    prog = []
    for q in range(s):
        coin = rnd.choice(['p', 'f'])
        prog.append((coin, [(rnd.randrange(s), rnd.choice([0, 1, None])) for _ in range(2)]))
    return prog
def productive(prog):
    s = len(prog)
    # reachable states from 0
    reach = {0}; st = [0]
    while st:
        q = st.pop()
        for (q2, e) in prog[q][1]:
            if q2 not in reach: reach.add(q2); st.append(q2)
    # from every reachable state, some emission is reachable with prob 1: no closed set of states with only silent moves
    # simple sufficient check: every reachable state can reach an emitting transition, and no silent closed class
    silent = {q: [q2 for (q2, e) in prog[q][1] if e is None] for q in range(s)}
    for q in reach:
        if all(e is None for (_, e) in prog[q][1]):
            # all moves silent: need that every silent path eventually hits a state with an emission
            seen = set(); st = [q]; closed = True
            while st:
                x = st.pop()
                if x in seen: continue
                seen.add(x)
                if any(e is not None for (_, e) in prog[x][1]): closed = False; continue
                st += silent[x]
            if closed: return False
    return True
def cyl_prob(prog, w, pval):
    """P(output starts with w) by solving the linear system over (state, matched) for silent moves."""
    s = len(prog); n = len(w)
    X = sp.symbols(f'x0:{s * n}') if n else []
    idx = lambda q, k: q * n + k
    eqs = []
    def val(q, k): return 1 if k == n else X[idx(q, k)]
    for q in range(s):
        coin, moves = prog[q]
        pr = [pval, 1 - pval] if coin == 'p' else [Fr(1, 2), Fr(1, 2)]
        for k in range(n):
            rhs = 0
            for (q2, e), c in zip(moves, pr):
                if e is None: rhs += c * val(q2, k)
                elif e == w[k]: rhs += c * val(q2, k + 1)
            eqs.append(sp.Eq(X[idx(q, k)], rhs))
    if n == 0: return Fr(1)
    sol = sp.solve(eqs, X, dict=True)[0]
    return sol[X[idx(0, 0)]]
def shortest_history(prog, w):
    from collections import deque
    n = len(w); start = (0, 0); dist = {start: 0}; dq = deque([start])
    while dq:
        q, k = dq.popleft()
        if k == n: return dist[(q, k)]
        for (q2, e) in prog[q][1]:
            k2 = k if e is None else (k + 1 if e == w[k] else None)
            if k2 is None: continue
            if (q2, k2) not in dist: dist[(q2, k2)] = dist[(q, k)] + 1; dq.append((q2, k2))
    return None
rnd = random.Random(3); tested = 0; ok_pattern = True; ok_bound = True; ok_ce = True
pvals = [Fr(1, 3), Fr(1, 7), Fr(5, 6)]
while tested < 12:
    prog = rand_prog(rnd, rnd.randint(2, 4))
    if not productive(prog): continue
    tested += 1
    for L in range(1, 5):
        for w in itertools.product([0, 1], repeat=L):
            probs = [sp.nsimplify(cyl_prob(prog, list(w), pv)) for pv in pvals]
            pos = [pr > 0 for pr in probs]
            ok_pattern &= len(set(pos)) == 1
            Lw = shortest_history(prog, list(w))
            ok_ce &= (Lw is not None) == pos[0]
            if Lw is not None:
                for pv, pr in zip(pvals, probs):
                    eps = min(pv, 1 - pv, Fr(1, 2)); ok_bound &= pr >= sp.Rational(eps.numerator, eps.denominator)**Lw
print(f"{tested} productive programs, all words of length ≤ 4, p ∈ {{1/3, 1/7, 5/6}}:")
print("  positive set independent of p:", ok_pattern)
print("  positive ⇔ some finite history emits w (c.e. by search):", ok_ce)
print("  μ_p(w) ≥ ε^{L(w)}, ε = min(p, 1−p, 1/2), L = shortest emitting history:", ok_bound)
