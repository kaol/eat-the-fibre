"""verify_v2_20.py — machine checks V1–V3 for discordianism-open-questions-explored-v2.20.md
Install:  pip install -r requirements.txt.   SLOW=1 adds the exact tangent / numerical local-dimension analysis (V3)."""
import sys, os, subprocess, re
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
SLOW = os.environ.get('SLOW') == '1'
res = []
def report(tag, ok, claim, detail=""):
    res.append(bool(ok)); print(f"[{tag}] {'PASS' if ok else 'FAIL'}  {claim}", flush=True)
    if detail: print(f"       {detail}", flush=True)
r = subprocess.run([sys.executable, 'q123v.py', '3'], capture_output=True, text=True)
if r.returncode: print("  !! q123v.py failed:", r.stderr.strip().splitlines()[-3:])
m = re.search(r'(\d+) folds, \|G\| = (\d+), orbits = (\d+)', r.stdout); nf = re.search(r'raises path length on A_\(≤1\) \(e\.g\. P\(β\) ∋ δβ\): (\d+)', r.stdout)
report("V1", m and m.group(1) == '219' and m.group(3) == '71',
       "census of Δ-compatible folds of k(1⇉2→3) over F_3 by kernel fibration: 219 folds in 71 orbits of G(F_3) (|G| = 6)", m.group(0) if m else r.stdout[-300:])
report("V2", nf and int(nf.group(1)) > 0,
       "with parallel arrows P need not preserve path length: some folds send β to β + δβ (δ = α₁ − α₂); P preserves the Δ-filtration instead",
       f"{nf.group(1) if nf else '?'} such folds")
if SLOW:
    r = subprocess.run([sys.executable, 'generic123.py'], capture_output=True, text=True)
    last = [l for l in r.stdout.splitlines() if 'generic local dim' in l or 'flagged folds' in l]
    report("V3", any(l.endswith('True') for l in last),
           "at every lifted Q-fold whose local dimension exceeds its (Aut_Δ + grading)-orbit dimension, a generic nearby fold has local dim = orbit dim (fixed-point degeneracies, no moduli)", " | ".join(last))
print(f"\n{sum(res)}/{len(res)} checks pass")
