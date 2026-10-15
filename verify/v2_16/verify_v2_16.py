"""verify_v2_16.py — machine checks W1–W6 for discordianism-open-questions-explored-v2.16.md
Install:  pip install -r requirements.txt      SLOW=1 adds the 5-vertex DAG sample and the n = 7 star."""
import sys, os, subprocess, re
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
SLOW = os.environ.get('SLOW') == '1'
missing = []
for mod, pkg in (('sympy', 'sympy'), ('numpy', 'numpy'), ('scipy', 'scipy'), ('sage.all__sagemath_singular', 'passagemath-singular')):
    try: __import__(mod)
    except Exception as e: missing.append(f"{pkg} ({type(e).__name__}: {e})")
if missing:
    print("Missing dependencies:\n  " + "\n  ".join(missing) + "\nInstall with:  pip install -r requirements.txt"); sys.exit(2)
res = []
def report(tag, ok, claim, detail=""):
    res.append(bool(ok)); print(f"[{tag}] {'PASS' if ok else 'FAIL'}  {claim}", flush=True)
    if detail: print(f"       {detail}", flush=True)
def run(args, env=None):
    r = subprocess.run([sys.executable] + args, capture_output=True, text=True, env={**os.environ, **(env or {})})
    if r.returncode != 0:
        print(f"  !! {' '.join(args)} exited with code {r.returncode}; stderr ends with:\n     " + "\n     ".join(r.stderr.strip().splitlines()[-6:]), flush=True)
    return r.stdout
allT = lambda o: o.count('True') >= 4 and 'False' not in o
o1 = run(['proofcheck.py', 'A3;A4;diamond;D4;A4 →→←;swapfail']); o2 = run(['runcase.py'])
report("W1", allT(o1) and 'P1..P4: (True, True, True, True)' in o2 and 'nullities by length: [0, 0, 0]' in o2,
       "proof steps P1–P4 (forest data; Hodge expansion path-supported; diagonal endpoints Hodge; run decomposition exact with path-supported, T-killed factors) on named quivers and on a hand-built run of two consecutive Podge vertices",
       o1.strip() + " | " + " ".join(l for l in o2.splitlines() if 'proof check' in l or 'nullities' in l))
o = run(['proofcheck.py', 'dags', '4'])
report("W2", allT(o), "proof steps P1–P4 on all nonempty DAGs on 4 vertices", o.strip())
if SLOW:
    o = run(['sample5.py', '25', '7'])
    report("W3", allT(o.strip().splitlines()[-1] if o.strip() else ''), "proof steps P1–P4 on 25 random 5-vertex DAGs", o.strip().splitlines()[-1] if o.strip() else '')
o = run(['randq.py', '8', '4', '3'])
report("W4", 'holds on all: True' in o, "independent cross-check: linear uniqueness certificate on random 4-vertex quivers", o.strip().splitlines()[-1])
N = 7 if SLOW else 6
code = ("from forestcount import count\n"
        f"print([count(n, [(0, i) for i in range(1, n)]) for n in range(1, {N} + 1)])\n"
        f"print([count(n, [(i, i + 1) for i in range(n - 1)]) for n in range(1, 7)])")
o = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True).stdout.strip().splitlines()
star = eval(o[0]); lin = eval(o[1])
formula = [(n + 1) ** (n - 1) + 2 * (n + 1) ** (n - 2) + 2 ** (n - 1) - 1 if n >= 2 else 2 for n in range(1, N + 1)]
terms = [int(x) for x in open('terms26.txt').read().split()]
report("W5", star == formula and lin == terms[:len(lin)],
       "fold counts: star = (n+1)^{n−1} + 2(n+1)^{n−2} + 2^{n−1} − 1 (n ≤ %d), and the linear Aₙ count equals v2.11's a(n)" % N,
       f"star {star} | formula {formula} | linear {lin}")
o = run(['forestcount.py', 'zigzag', '6']) + subprocess.run([sys.executable, '-c', "from arrowrule import count\nprint([count(n, [(i, i + 1) if i % 2 == 0 else (i + 1, i) for i in range(n - 1)]) for n in range(1, 7)])"], capture_output=True, text=True).stdout
vals = re.findall(r'\[[\d, ]+\]', o)
report("W6", len(vals) == 2 and vals[0] == vals[1], "the explicit forest description S_v = subtree(v) / subtree(v)∖{v} reproduces the matrix-based arrow-rule counts (zigzag Aₙ, n ≤ 6)", " vs ".join(vals))
print(f"\n{sum(res)}/{len(res)} checks pass")
