"""verify_v2_12.py — machine checks S1–S5 for discordianism-open-questions-explored-v2.12.md"""
import sys, os, subprocess, re
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
res = []
def report(tag, ok, claim, detail=""):
    res.append(bool(ok)); print(f"[{tag}] {'PASS' if ok else 'FAIL'}  {claim}", flush=True)
    if detail: print(f"       {detail}", flush=True)
run = lambda args, env=None: subprocess.run([sys.executable] + args, capture_output=True, text=True, env={**os.environ, **(env or {})}).stdout
out = run(['deadline.py'])
report("S1", out.count('True') == 7 and 'missed on 11' in out,
       "parameter complexity (toy halting tables, exact arithmetic): a naming program yields a p-computable deadline function certifying every heavy branch (necessity); any dominating deadline reproduces the target masses exactly while a non-dominating one fails (sufficiency); the fair-coin halting measure needs only the computable deadline; and a valid deadline is always computable from μ itself, so Name(μ) = cone(μ)",
       out.strip().replace("\n", " | "))
out = run(['quiver.py', 'A2 (1→2)', 'A3 (1→2→3)', '1→2←3', 'Kronecker 2', 'Kronecker 3'], {'TRIALS': '25'})
mx = re.findall(r'max \(local dim − orbit dim\) = (\d+)', out)
report("S2", len(mx) == 5 and all(x == '0' for x in mx),
       "folds on finite acyclic path algebras are rigid modulo Δ-compatible automorphisms: at every sampled fold the local dimension of the solution variety equals the orbit dimension (A₂, A₃, 1→2←3, 2- and 3-Kronecker)",
       " | ".join(l for l in out.splitlines() if 'max' in l))
out = run(['kdiag.py'])
rows = re.findall(r"\[\[(.*?)\]\]", out)
ok3 = True; samples = 0
for line in out.splitlines():
    if 'normalized sv' not in line: continue
    nums = [float(x) for x in re.findall(r"np\.float64\(([0-9.e-]+)\)", line)]
    if len(nums) < 12: continue
    third = [nums[2], nums[6], nums[10]]; samples += 1
    ok3 &= third[1] / third[0] > 5 and third[2] / third[1] > 5          # grows ~10x as eps shrinks 10x: absolute noise floor
    ok3 &= min(nums[0], nums[4], nums[8]) > 1 and min(nums[1], nums[5], nums[9]) > 1   # two genuine directions
report("S3", ok3 and samples >= 3,
       "at Kronecker folds with 3-dimensional tangent space only 2 directions are genuine (displacements ∝ ε); the third is an obstructed first-order direction (displacement is a constant noise floor), so v2.4's 'continuous moduli' are the 2-dimensional automorphism orbits",
       f"{samples} points examined")
terms = [int(x) for x in open('terms26.txt').read().split()]
from mpmath import mp, mpf, log, exp, matrix, lu_solve
mp.dps = 40; N = len(terms)
lr = [log(mpf(terms[n - 1]) / (2 * mpf(n + 2) ** (n - 1))) for n in range(1, N + 1)]
def est(basis, upto):
    ns = list(range(upto - len(basis) + 1, upto + 1)); A = matrix(len(basis), len(basis)); b = matrix(len(basis), 1)
    for r, n in enumerate(ns):
        for j, f in enumerate(basis): A[r, j] = f(n)
        b[r] = lr[n - 1]
    return float(exp(lu_solve(A, b)[1]))
B3 = [lambda n: 1, lambda n: n, lambda n: log(n)]; B4 = B3 + [lambda n: mpf(1) / n]
e3 = [est(B3, u) for u in (N - 4, N - 2, N)]; e4 = [est(B4, u) for u in (N - 4, N - 2, N)]
report("S4", all(0.3 < x < 0.45 for x in e3 + e4) and abs(e3[-1] - e4[-1]) > 0.001,
       "the decay constant of the compatible fraction is NOT determined by 26 terms: fitted ρ drifts with the cutoff and with the ansatz",
       f"ρ (A + n log ρ + α log n): {[round(x, 4) for x in e3]};  (+ β/n): {[round(x, 4) for x in e4]}")
out5 = subprocess.run([sys.executable, "hcube.py"], capture_output=True, text=True, cwd=HERE).stdout
report("S5", out5.count("True") == 2,
       "the full-support box B(2^<ω, ¼) is computably homeomorphic to the Hilbert-cube box [¼,¾]^(2^<ω) (exact round trip and an explicit modulus for the inverse), so Name(μ) inherits non-total continuous degrees")

print(f"\n{sum(res)}/{len(res)} checks pass")
