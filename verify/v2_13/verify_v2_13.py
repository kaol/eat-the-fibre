"""verify_v2_13.py — machine checks T1–T10 for discordianism-open-questions-explored-v2.13.md
Requires:  pip install -r requirements.txt   (sympy, numpy, scipy, mpmath, passagemath-singular)
SLOW=1 also runs the slow checks (3-Kronecker orbit decomposition, A4/diamond fibres, numeric cross-check)."""
import sys, os, subprocess, re
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
SLOW = os.environ.get('SLOW') == '1'
res = []
missing = []
for mod, pkg in (('sympy', 'sympy'), ('numpy', 'numpy'), ('scipy', 'scipy'), ('mpmath', 'mpmath'),
                 ('sage.all__sagemath_singular', 'passagemath-singular')):
    try: __import__(mod)
    except Exception as e: missing.append(f"{pkg} ({type(e).__name__}: {e})")
if missing:
    print("Missing dependencies:\n  " + "\n  ".join(missing) + "\nInstall with:  pip install -r requirements.txt")
    sys.exit(2)
def report(tag, ok, claim, detail=""):
    res.append(bool(ok)); print(f"[{tag}] {'PASS' if ok else 'FAIL'}  {claim}", flush=True)
    if detail: print(f"       {detail}", flush=True)
def run(args, env=None):
    r = subprocess.run([sys.executable] + args, capture_output=True, text=True, env={**os.environ, **(env or {})})
    if r.returncode != 0:                       # never hide a crashed sub-script
        print(f"  !! {' '.join(args)} exited with code {r.returncode}; its stderr ends with:", flush=True)
        print("     " + "\n     ".join(r.stderr.strip().splitlines()[-6:]), flush=True)
    return r.stdout
out = run(['degrees.py']); L = out.splitlines()
report("T1", len(L) >= 4 and all(l.endswith('True') for l in L[:4]),
       "every c.e. tree without dead ends carries a COMPUTABLE law in S (stage weights): positive set = U, S-bound, lower approximations sound and exhaustive",
       " | ".join(L[:4]))
report("T2", len(L) >= 5 and L[4].count('True') == 3,
       "coding at designated branching nodes: y and every displaced conditional decode exactly (μ′ ≡ μ ⊕ y)", L[4] if len(L) > 4 else out)
out = run(['foldcount.py', 'A2', 'A3', '1→2←3'])
f = re.findall(r'folds = (\d+) \(all satisfy the fold equations exactly: (\w+)\); positive-dim Q pieces admit no complement: (\w+)', out)
report("T3", [x[0] for x in f] == ['6', '23', '27'] and all(x[1] == 'True' and x[2] == 'True' for x in f),
       "Schubert-cell classification: A₂, A₃, 1→2←3 have exactly 6, 23, 27 folds; every reconstructed P satisfies the fold equations; positive-dimensional Q-pieces have no complement",
       " | ".join(l for l in out.splitlines()))
qs = ['A2', 'K2', 'A3', '1→2←3'] + (['K3'] if SLOW else [])
out = run(['orbits.py'] + qs)
fin = re.findall(r'^(\S+): nonempty chart pairs = (\d+), chart dims = (\[.*?\]); fold variety is a finite union of G-orbits: (\w+)', out, re.M)
report("T4", len(fin) == len(qs) and all(x[3] == 'True' for x in fin),
       "open-chart pairs + minimal primes + orbit-rank recursion: the fold variety is a finite union of Aut_Δ-orbits for " + ", ".join(qs) + ("" if SLOW else " (SLOW=1 adds the 3-Kronecker quiver)"),
       " | ".join(" ".join(x) for x in fin))
fast = ['A1', 'A2', 'A3', '1→2←3', '1←2→3', 'D4', '1→2→3, 1→3'] + (['A4', 'A4 alt', 'diamond'] if SLOW else [])
out_e = run(['extend.py'] + fast); out_a = run(['arrowrule.py'] + fast)
ex = dict(re.findall(r'^(.+?): dim A = \d+; folds = (\d+)', out_e, re.M))
can = re.findall(r"^(.+?): .*?\{1: (\d+), 'the point is a = 0 \(canonical tensor extension\)': (\d+)\}", out_e, re.M)
ar = dict(l.rsplit(' ', 1) for l in out_a.strip().splitlines())
report("T5", len(ex) == len(fast) and all(ex[q] == ar[q] for q in fast) and len(can) == len(fast) and all(a == b for _, a, b in can),
       "no parallel arrows: every vertex fold passing the arrow rule extends to exactly ONE fold, the canonical tensor extension P₀^{⊗•}; fold counts = arrow-rule counts",
       " | ".join(f"{q}: {ex.get(q)} = {ar.get(q)}" for q in fast))
out = run(['kerdelta.py'])
report("T6", out.strip().endswith('True'), "ker Δ = A₀ exactly when there are no parallel arrows (11 quivers, exact rank)", (out.strip().splitlines() or ['(no output)'])[-1])
out = run(['gdim.py'])
report("T7", out.strip().endswith('True'), "dim Lie Aut_Δ = Σ m_ij(m_ij − 1) on 11 quivers (exact nullspace)", " | ".join(out.strip().splitlines()[:-1]))
terms = [int(x) for x in open('terms26.txt').read().split()]
out = run(['arrowrule.py', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6'])
vals = [int(l.split()[-1]) for l in out.strip().splitlines()]
report("T8", vals == terms[:6], "the arrow-rule count on Aₙ reproduces v2.11's sequence a(1…6) = 2, 6, 23, 100, 563, 3734, so a(n) ≤ |Folds(kAₙ)| (canonical extension), with equality wherever T10 certifies uniqueness (n ≤ 6 with SLOW=1)", str(vals))
lq = 'A2;A3;1→2←3;D4;1→2→3, 1→3' + (';A4;A4 alt;diamond;A5;A6' if SLOW else '')
out = run(['lin.py', lq])
nl = re.findall(r'^(.+?): (\d+) arrow-rule vertex folds; max nullity .*?: (\d+)$', out, re.M)
report("T10", len(nl) == len(lq.split(';')) and all(x[2] == '0' for x in nl),
       "linear uniqueness certificate (conicity ⇒ the fibre equations at each length are linear): the vertex-product, H·H and Q·Q equations alone force a = 0, so every fold is the canonical extension"
       + ("" if SLOW else " (SLOW=1 adds A₄, 1→2←3→4, diamond, A₅, A₆)"),
       " | ".join(f"{x[0]}: {x[1]} vertex folds, nullity {x[2]}" for x in nl))
if SLOW:
    out = run(['crosscheck.py', 'A2', 'A3', '1→2←3'], {'TRIALS': '40'})
    report("T9", out.strip().endswith('True'), "every fold found by v2.12's numerical sampler is one of the exact folds", " | ".join(l for l in out.splitlines() if 'numerical' in l))
print(f"\n{sum(res)}/{len(res)} checks pass")
