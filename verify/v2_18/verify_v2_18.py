"""verify_v2_18.py — machine checks Y1–Y3 for discordianism-open-questions-explored-v2.18.md
Install:  pip install -r requirements.txt  (numpy, sympy).   SLOW=1 adds K_2 over F_7 and K_3 over F_3."""
import sys, os, subprocess, re
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
SLOW = os.environ.get('SLOW') == '1'
for mod in ('numpy', 'sympy'):
    try: __import__(mod)
    except Exception as e: print(f"Missing dependency {mod}: {e}\nInstall with:  pip install -r requirements.txt"); sys.exit(2)
from kronformula import N
res = []
def report(tag, ok, claim, detail=""):
    res.append(bool(ok)); print(f"[{tag}] {'PASS' if ok else 'FAIL'}  {claim}", flush=True)
    if detail: print(f"       {detail}", flush=True)
cases = [(1, 3), (2, 3), (2, 5)] + ([(2, 7), (3, 3)] if SLOW else [])
rows = []
for m, p in cases:
    r = subprocess.run([sys.executable, 'kronff.py', str(m), str(p)], capture_output=True, text=True)
    if r.returncode != 0: print("  !! kronff.py failed:", r.stderr.strip().splitlines()[-3:])
    f = re.search(r'(\d+) folds, \|G\| = (\d+), orbits = (\d+)', r.stdout)
    rows.append((m, p, int(f.group(1)), int(f.group(3))) if f else (m, p, None, None))
report("Y1", all(o == 14 * m - 8 for m, p, n, o in rows),
       "exhaustive census over F_p: the Δ-compatible folds of the Kronecker algebra K_m form exactly 14m − 8 orbits of Aut_Δ = {g ∈ GL_m : σ∘g = σ}",
       " | ".join(f"K_{m}/F_{p}: {o} orbits" for m, p, n, o in rows))
report("Y2", all(n == N(m, p) for m, p, n, o in rows),
       "the census fold counts equal the point count N_m(q) predicted by the hand classification",
       " | ".join(f"K_{m}/F_{p}: {n} = {N(m, p)}" for m, p, n, o in rows))
from sympy import symbols, Poly, cancel, expand
q = symbols('q')
degs = [Poly(expand(cancel(N(m, q))), q).degree() for m in range(2, 7)]
report("Y3", degs == [m * m // 2 for m in range(2, 7)],
       "dim of the fold variety (degree of N_m) = ⌊m²/2⌋ < dim Aut_Δ = m² − m for m ≥ 3", f"degrees m=2..6: {degs}")
print(f"\n{sum(res)}/{len(res)} checks pass")
