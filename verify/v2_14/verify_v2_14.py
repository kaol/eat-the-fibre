"""verify_v2_14.py — machine checks U1–U7 for discordianism-open-questions-explored-v2.14.md
Install:  pip install -r requirements.txt      SLOW=1 adds larger cases."""
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
# U1 two-stage structure of the uniqueness system
qs = ['A4', 'diamond', 'A4 →→←'] + (['D4', 'A4 alt'] if SLOW else [])
outs = [run(['perlen.py', q, 'ZQ,HL,HR,HH']) + run(['perlen.py', q, 'ZH,QL,QR,QQ,HL,HR']) for q in qs]
nul = re.findall(r'max nullity per length: (\{.*\})', "".join(outs))
report("U1", len(nul) == 2 * len(qs) and all(set(eval(x).values()) <= {0} for x in nul),
       "Hodge half (proved): with a = 0 on the Podge part, HL, HR, HH alone force a = 0; Podge half: with a = 0 on the Hodge part, Q₀-multiples and Q·Q force a = 0 (" + ", ".join(qs) + ")",
       " | ".join(nul))
# U2 random acyclic quivers
o = run(['randq.py', '10', '4', '1']) + (run(['randq.py', '12', '5', '2']) if SLOW else '')
fin = re.findall(r'random quivers on (\d) vertices: Conjecture 4 holds on all: (\w+)', o)
report("U2", fin and all(x[1] == 'True' for x in fin), "Conjecture 4 (unique extension) on random acyclic quivers without parallel arrows (linear certificate, all lengths)",
       " | ".join(f"{n} vertices: {v}" for n, v in fin) + f" | {len(re.findall('max nullity 0', o))} quivers certified")
# U3 isomorphism
o = run(['iso.py'])
report("U3", o.count('True') == 2, "k⟨a,u⟩ with a-deletion ≅ path algebra of the 2-loop quiver (a ↦ α, u ↦ α − β) as ε-bialgebras", o.strip())
# U4 cycles rigid
o = "".join(run(['cfibre.py', '2-cycle', str(L)]) for L in ((2, 3, 4, 5) if SLOW else (2, 3, 4))) + run(['cfibre.py', '3-cycle', '3'])
lines = [l for l in o.splitlines() if 'fibre dimension' in l]
report("U4", lines and all('fibre dimension 0, points = 1' in l for l in lines),
       "cycles without loops are rigid: every fibre over an arrow-rule vertex fold is a single point (2-cycle up to L = 4/5, 3-cycle L = 3)", f"{len(lines)} fibres, all single points")
# U5 looped fibres = collapse-family dimension
exp = {'loop+arrow': 1, 'arrow+loop at target': 1, 'two loops+arrow': 2, '2-cycle+loop': 2}
det = []
okk = True
for q, e in exp.items():
    for L in ((3, 4) if SLOW else (3,)):
        o = run(['cfibre.py', q, str(L)])
        dims = [int(x) for x in re.findall(r'fibre dimension (\d+)', o)]
        pos = [x for x in dims if x > 0]
        okk &= bool(dims) and pos and all(x == e for x in pos)
        det.append(f"{q} L={L}: positive fibres {pos}")
report("U5", okk, "loops create positive-dimensional fibres, of exactly the dimension of the loop-collapse family (parameters = arrows among {s, s′} other than u)", " | ".join(det))
# U6 collapse theorem
o = run(['collapse.py', '4' if SLOW else '3'])
report("U6", o.strip().endswith('True'), "loop-collapse folds P = χ(·)e_s with σ = (1, −1) on {s, s′} and monic linear weights are Δ-compatible folds (symbolic parameters); three-vertex support, equal signs, non-monic weights and a translated loop at s all fail",
       " | ".join(l for l in o.splitlines() if ':' in l))
# U7 orbits: Δ-compatible derivations (loop translations) vs fibre dimensions
o = "".join(run(['loopder.py', q, '4']) for q in ('loop+arrow', 'two loops+arrow', '2-cycle+loop'))
pairs = re.findall(r'fibre dim (\d+), orbit dim of the canonical fold (\d+)', o)
pos = [(int(a), int(b)) for a, b in pairs if int(a) > 0]
report("U7", pos == [(1, 1), (2, 2), (2, 2), (2, 1)],
       "loop translations u ↦ u − c·e_s are Δ-compatible; the loop+arrow and two-loops+arrow families are single orbits, the 2-cycle+loop family has a genuine 1-dimensional modulus (fibre 2, orbit 1)",
       " | ".join(l.strip() for l in o.splitlines() if 'derivations' in l or 'fibre dim' in l and not l.strip().endswith('orbit dim of the canonical fold 0')))
print(f"\n{sum(res)}/{len(res)} checks pass")
