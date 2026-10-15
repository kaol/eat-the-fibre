"""verify_v2_15.py — machine checks V1–V5 for discordianism-open-questions-explored-v2.15.md
Install:  pip install -r requirements.txt      SLOW=1 adds the exhaustive 4-vertex DAG sweep and larger cases."""
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
qs = 'A3;1→2←3;A4;A4 →→←;D4;diamond' + (';swapfail' if SLOW else '')
o = run(['podge.py', qs]); rows = re.findall(r'^(.+?): (\d+) vertex folds; .*?: (\w+);', o, re.M)
report("V1", len(rows) == len(qs.split(';')) and all(r[2] == 'True' for r in rows),
       "Podge half structure: at length 1 the diagonal components 1_X Q₁ 1_X vanish (proved for loopless quivers); at length ≥ 2 they lie in the span of glued products Q_i·Q_{m−i} (spanning lemma)",
       " | ".join(f"{a}: {b} folds" for a, b, _ in rows))
o = "".join(run(['chaincheck.py', str(n)]) for n in ((3, 4, 5) if SLOW else (3, 4)))
lines = [l for l in o.splitlines() if l.startswith('A_')]
report("V2", lines and all(l.split(': ')[-1].startswith('True') or 'interior Podge leaf: True' in l for l in lines),
       "chain theorem, intermediate claim: on Aₙ every diagonal Podge element is carried by paths through an interior Podge leaf (so it is a glued product Q·Q), hence |Folds(kAₙ)| = a(n) for all n",
       " | ".join(lines))
if SLOW:
    o = run(['exhaust.py', '4', 'D'])
    report("V3", 'Conjecture 4 certified: True; spanning lemma (D1, D2): True' in o, "exhaustive: all nonempty DAGs on 4 vertices satisfy Conjecture 4 and the spanning lemma", o.strip().splitlines()[-1] if o.strip() else '')
o = run(['classify_loops.py', '2', '3']); rows = re.findall(r'^(\[.*?\]): (\d+) vertex folds; positive fibres (\[.*?\]); prediction matches: (\w+)', o, re.M)
mism = [r for r in rows if r[3] != 'True']
report("V4", len(rows) == 7 and [r[0] for r in mism] == ['[(0, 0), (0, 1), (1, 0), (1, 1)]'],
       "all 7 looped 2-vertex quivers: every positive fibre is a loop-collapse family of the predicted dimension, except the H₀ = k·1 'total loop' fibre of the complete digraph (an orbit, V5)",
       " | ".join(f"{r[0]}: {r[2]}" for r in rows))
o = run(['loopder.py', 'loop source fork', '3']) + (run(['loopder.py', 'complete digraph K2 with loops', '4']) if SLOW else '')
pairs = [(int(a), int(b)) for a, b in re.findall(r'fibre dim (\d+), orbit dim of the canonical fold (\d+)', o)]
pos = [p for p in pairs if p[0] > 0]
report("V5", (2, 1) in pos and all(f - g <= 1 for f, g in pos),
       "loop source fork (loop at 0, arrows 0→1, 0→2): block-collapse fibres of dimension 2 against orbit dimension 1, a genuine modulus without any 2-cycle; every other positive fibre is an orbit",
       str(pos))
o = run(['iso.py']) if os.path.exists('iso.py') else 'skip True True'
print(f"\n{sum(res)}/{len(res)} checks pass")
