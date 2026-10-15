"""verify_v2_17.py — machine checks X1–X6 for discordianism-open-questions-explored-v2.17.md
Install:  pip install -r requirements.txt      SLOW=1 adds the looped and 4-vertex cyclic comparisons and longer truncations."""
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
o = run(['cyclic.py', '3', '4']); sus = re.findall(r'SUSPICIOUS \((\[.*?\]), (\[\[.*?\]\]), ', o)
o2 = run(['hcheck.py'])
report("X1", len(sus) == 4 and 'violated by 4 vertex folds; predictor flags 4; same set: True' in o2,
       "all 10 loopless cyclic 3-vertex quivers: the linear certificate fails for exactly 4 vertex folds, and these are exactly the (H)-violations (= the closed-walk predictor)",
       f"{len(sus)} suspicious; " + o2.strip())
Ls = (3, 4, 5) if SLOW else (3, 4)
dims = {}
for q in ('cyc A', 'cyc B', 'cyc C'):
    for L in (Ls if q != 'cyc C' else Ls[:2]):
        o = run(['cfibre.py', q, str(L)]); dims[(q, L)] = sorted(int(x) for x in re.findall(r'fibre dimension (\d+)', o) if int(x) > 0)
report("X2", all(v == ([1] if k[0] != 'cyc C' else [2, 2]) for k, v in dims.items()),
       "the (H)-violations carry genuine fold families, stable in the truncation: dimension 1 for cyc A (0⇄1, 0→2) and cyc B (0⇄1, 2→0), 2 for cyc C (0⇄1, 0⇄2)",
       str(dims))
o = run(['loopder.py', 'cyc A', '4'])
report("X3", 'derivations (on generators) = 0' in o and 'fibre dim 1, orbit dim of the canonical fold 0' in o,
       "cyc A has no Δ-compatible derivations: its family is a genuine 1-parameter modulus (fibre 1, orbit 0)", " | ".join(l.strip() for l in o.splitlines() if 'dim' in l and ('derivations' in l or 'fibre dim 1' in l)))
code = ("import itertools\nfrom dags import dags\nfrom hcheck import violates_H\nfrom extend import vertex_folds\nbad = 0; n = 0\n"
        "for arrows in [a for a in dags(4) if a]:\n    A = set(arrows)\n    for P0 in vertex_folds(4):\n"
        "        if not all(not [a for a in range(4) if P0[a, s]] or not [b for b in range(4) if P0[b, t]] or all((a, b) in A for a in range(4) if P0[a, s] for b in range(4) if P0[b, t]) for s, t in arrows): continue\n"
        "        n += 1; bad += violates_H(4, arrows, P0, 4)\nprint(n, bad)")
o = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True).stdout.split()
report("X4", len(o) == 2 and o[1] == '0', "acyclic quivers satisfy (H) automatically: no violation among the arrow-rule vertex folds of all 4-vertex DAGs", f"{o[0] if o else '?'} vertex folds checked")
if SLOW:
    o = run(['hloops.py', '3'])
    bad = [l for l in o.splitlines() if 'positive fibre' in l and ': False' in l]
    report("X5", len(bad) == 1 and bad[0].startswith('[(0, 0), (0, 1), (1, 0), (1, 1)]') and '[[0, 0], [1, 1]]' in bad[0] and '[[1, 1], [0, 0]]' in bad[0],
           "looped 2-vertex quivers: every (H)-violation has a positive fibre, and the only other positive fibres are the two H₀ = k·1 'total loop' translation orbits of the complete digraph (looped Hodge block)",
           " | ".join(bad))
    o = run(['cyc4.py', '12', '4', '5'])
    report("X6", 'every positive nullity is predicted: True' in o, "12 random loopless cyclic 4-vertex quivers: every positive linear nullity is an (H)-type violation", o.strip().splitlines()[-1])
print(f"\n{sum(res)}/{len(res)} checks pass")
