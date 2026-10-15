"""verify_v2_11.py — machine checks R1–R5 for discordianism-open-questions-explored-v2.11.md"""
import sys, os, subprocess, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); os.chdir(HERE)
res = []
def report(tag, ok, claim, detail=""):
    res.append(bool(ok)); print(f"[{tag}] {'PASS' if ok else 'FAIL'}  {claim}", flush=True)
    if detail: print(f"       {detail}", flush=True)
from count import a
# R1: formula vs brute force
brute = []
try:
    subprocess.run(["gcc", "-O2", "-o", "cnt", "cnt.c"], check=True)
    for n in range(1, 8):
        out = subprocess.run(["./cnt", str(n)], capture_output=True, text=True).stdout
        brute.append(int(out.split("compatible=")[1]))
except Exception as ex: brute = [str(ex)]
cached = [2, 6, 23, 100, 563, 3734, 28777, 246164, 2395445, 25665390]
form = [a(n) for n in range(1, 11)]
report("R1", form == cached and brute == cached[:7],
       "polynomial-time formula for arrow-compatible splittings (leaf-type rule + K_{1,x,y} spanning-tree count + inclusion–exclusion + matched P1-run/HL-block DP) agrees with exhaustive enumeration",
       f"formula n=1..10: {form}; brute force (C, recomputed) n ≤ 7: {brute}; n = 8..10 from the earlier exhaustive runs")
t0 = time.time(); more = [a(n) for n in range(11, 16)]
terms = [int(x) for x in open('terms25.txt').read().split()]
report("R2", more == terms[10:15], "the formula extends the sequence (terms to n = 25 in terms25.txt; n = 11..15 recomputed here)",
       f"a(11..15) = {more}  [{time.time()-t0:.0f}s]")
import sympy as sp
from math import factorial
A = {n + 1: v for n, v in enumerate(terms)}; N = len(terms); found = []
for lab, seq in (("a(n)", lambda n: sp.Integer(A[n])), ("a(n)/n!", lambda n: sp.Rational(A[n], factorial(n)))):
    for order in range(1, 5):
        for deg in range(0, 4):
            unk = (order + 1) * (deg + 1)
            if N - order < unk + 3: continue
            M = sp.Matrix([[sp.Integer(n) ** d * seq(n - k) for k in range(order + 1) for d in range(deg + 1)] for n in range(1 + order, N + 1)])
            if M.nullspace(): found.append((lab, order, deg))
report("R3", not found, "no linear recurrence with polynomial coefficients (order ≤ 4, degree ≤ 3, ≥ 3 surplus equations) fits a(n) or a(n)/n! on 25 terms",
       f"found: {found or 'none'}")
out = subprocess.run([sys.executable, "budget.py"], capture_output=True, text=True).stdout
report("R4", "True" in out and "8/8" in out,
       "undecidable supports: verification with give-up coins whose deadline is ENCODED in a real parameter reproduces any target conditional probability exactly and stays productive; any fixed computable schedule is capped by ρ_T → 0 (toy halting table, exact arithmetic)",
       out.strip().replace("\n", " | "))
out = subprocess.run([sys.executable, "borel.py"], capture_output=True, text=True).stdout
report("R5", out.count("True") == 3, "upper-bound lemma re-checked (positive set p-independent and c.e.; μ_p(w) ≥ ε^{L(w)})", out.strip().replace("\n", " | "))
print(f"\n{sum(res)}/{len(res)} checks pass")
