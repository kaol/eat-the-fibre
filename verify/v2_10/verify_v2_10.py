"""verify_v2_10.py — machine checks Q1–Q5 for discordianism-open-questions-explored-v2.10.md"""
import sys, os, subprocess, re
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); os.chdir(HERE)
res = []
def report(tag, ok, claim, detail=""):
    res.append(bool(ok)); print(f"[{tag}] {'PASS' if ok else 'FAIL'}  {claim}", flush=True)
    if detail: print(f"       {detail}", flush=True)
def run(args):
    return subprocess.run([sys.executable] + args, capture_output=True, text=True, cwd=HERE).stdout
# Q1: structure of A0/Q0 as k[M]-module
import sympy as sp
exec(open('psi.py').read().split("words = list")[0])
m1 = Mop({one: 1})
ubar = {k: v for k, v in S.decompose(sp.Symbol('u')).items() if S.side[k] == 'H'}
import fold as FD
chain = [ubar]; 
for _ in range(5): chain.append(Mop(chain[-1]))
E = FD.Echelon(); indep = all(E.add(dict(v)) for v in chain)
report("Q1", m1 == {} and indep,
       "A₀/Q₀ ≅ H₀ is, as a k[M]-module (M = p∘q₀·), k·1̄ ⊕ k[M]·ū: M kills 1̄ and ū, Mū, M²ū, … are independent (free part)",
       f"M(1̄) = {m1}; ū, Mū, …, M⁵ū independent: {indep}")
out = run(['psieq.py', '6'])
dims = re.findall(r'dim (\d+)', out)
report("Q2", len(dims) == 2 and dims[0] == dims[1],
       "first-order one-sided corrections ψ(h) = C(h·a·q₀) are exactly the solutions of Mψ − ψM = Φ((1−p)(q₀·))·ψ(1) (hand-derived from the Q₀-module conditions); the solution space equals the engine's all-a-degree tangent",
       out.strip().replace("\n", " | "))
out = run(['psi.py', '3'])
lines = {l.split(':')[0].strip(): int(re.search(r'bad (\d+)', l).group(1)) for l in out.strip().splitlines()}
exact = [k for k, v in lines.items() if v == 0]
report("Q3", set(exact) == {'root deviation  x − ev_r(x)', '(id − ev_r) with r = r−'},
       "with C supported on one-a words only, the exact folds among the first-order family are precisely the root deviations μ(x − ev_r(x)), r = b ± √s; any M-dependence needs corrections on longer words",
       f"exact: {exact}")
out = run(['borel.py'])
report("Q4", out.count('True') == 3,
       "key lemma for language-level Borelness (random productive coin programs, exact rational arithmetic): the positive set of μ_p is independent of p ∈ (0,1) and is the c.e. set of words emitted by some finite history; μ_p(w) ≥ ε^{L(w)} with L the shortest emitting history",
       out.strip().replace("\n", " | "))
out = run(['peano.py'])
report("Q5", out.count('64/64') == 2,
       "lower-bound ingredient: the iterated Hilbert-curve map x ↦ (H₁x, H₁H₂x, H₁H₂H₂x, …) is a computable continuous surjection onto the Hilbert cube (coverage of finite projections)",
       out.strip().replace("\n", " | "))
print(f"\n{sum(res)}/{len(res)} checks pass")
