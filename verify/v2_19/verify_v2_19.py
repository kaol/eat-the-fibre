"""verify_v2_19.py — machine checks Z1–Z3 for discordianism-open-questions-explored-v2.19.md
Install:  pip install -r requirements.txt      SLOW=1 adds L = 4 and the 4-vertex confirmation sample."""
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
def run(args):
    r = subprocess.run([sys.executable] + args, capture_output=True, text=True)
    if r.returncode != 0: print(f"  !! {' '.join(args)} failed:\n     " + "\n     ".join(r.stderr.strip().splitlines()[-5:]), flush=True)
    return r.stdout
o = run(['obstr.py'] + (['2', '3', '4'] if SLOW else ['2', '3']))
lines = [l for l in o.splitlines() if l.startswith('L=')]
report("Z1", lines and all('linear nullities [0, 1' in l and 'exact fibre dimension 0, points 1' in l for l in lines),
       "counterexample to Conjecture 5: on 0→2, 0→3, 2→0, 3→1 with vertex fold e₁ ↦ −e₀, e₃ ↦ −e₂, (H) is violated and a first-order direction exists at length 2, yet the exact fibre is the single canonical fold",
       " | ".join(lines))
o = run(['explicit.py', 'cyc A', '4', '[[1, 0, 0], [0, 1, 0], [0, -1, 0]]'])
report("Z2", 'free coordinate: [a8_1]' in o and 'a(1→0→2) has e_1-component 1' in o and o.count('component') == 1,
       "the cyc A family has a single vertex-level correction, a(1→0→2) = c·e₁; every other deviation from the canonical fold is forced through Δ", o.strip().replace('\n', ' | '))
if SLOW:
    o = run(['confirm4.py', '12', '3', '5'])
    last = o.strip().splitlines()[-1] if o.strip() else ''
    report("Z3", 'with positive fibre: 27' in last, "4-vertex cyclic sample: 27 of 28 (H)-violations carry positive fibres; the pair-count dimension formula fails in 9 of 28", last)
print(f"\n{sum(res)}/{len(res)} checks pass")
