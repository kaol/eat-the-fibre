"""verify_v2_22.py — machine checks X1–X10 for discordianism-open-questions-explored-v2.22.md (a parallel pair plus more arrows)
Install:  pip install -r requirements.txt.
Default: X1–X6 (a few minutes).  SLOW=1 adds F_7 for 1⇉2,1→3 (X7) and the 4-vertex quivers over F_3 (X8, X9).
SLOW=2 adds the 4-vertex quivers over F_5 (X10; about an hour or more)."""
import sys, os, subprocess, re, ast
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
SLOW = int(os.environ.get('SLOW', '0') or 0)
res = []
def report(tag, ok, claim, detail=""):
    res.append(bool(ok)); print(f"[{tag}] {'PASS' if ok else 'FAIL'}  {claim}", flush=True)
    if detail: print(f"       {detail}", flush=True)
def run(*args):
    r = subprocess.run([sys.executable, '-u', *args], capture_output=True, text=True)
    if r.returncode: print(f"  !! {' '.join(args)} failed:", *r.stderr.strip().splitlines()[-3:], sep="\n     ")
    return r.stdout
def census(q, p):
    out = run('gcensus.py', q, str(p)); line = next((l for l in out.splitlines() if l.startswith('{')), None)
    return ast.literal_eval(line) if line else None
def graded(q, p):
    m = re.search(r'G-orbits (\d+); \(G × grading\)-orbits (\d+)', run('graded.py', q, str(p)))
    return (int(m.group(1)), int(m.group(2))) if m else None
inv = lambda c: (c['N'], c['orbits'], c['sizes'], c['fib']) if c else None

k3, k5, a3 = census('K2', 3), census('K2', 5), census('pair_after', 3)
report("X1", k3 and k5 and a3 and (k3['N'], k3['orbits'], k5['N'], k5['orbits'], a3['N'], a3['orbits']) == (64, 20, 164, 20, 219, 71)
       and k3['sizes'] == {1: 6, 2: 2, 3: 6, 6: 6},
       "the general census (gcensus.py) reproduces v2.18 (K₂: 64 / 164 folds, 20 orbits) and v2.21 (1⇉2→3: 219 folds, 71 orbits); "
       "K₂ has orbit sizes 6×1, 2×(q−1), 6×q, 6×q(q−1), summing to v2.18's 6q² + 2q + 4",
       f"K2/F3 {k3 and k3['sizes']}  K2/F5 {k5 and k5['sizes']}")
b3, s3, t3 = census('pair_before', 3), census('pair_srcout', 3), census('pair_tgtin', 3)
report("X2", inv(b3) == inv(a3) and inv(t3) == inv(s3),
       "opposite quivers have identical censuses (kQ^op has the same subalgebras, and the flip τ commutes with P⊗P): "
       "1→2⇉3 ≡ 1⇉2→3 and 1⇉2,3→2 ≡ 1⇉2,1→3", f"{b3 and b3['N']} = {a3 and a3['N']}, {t3 and t3['N']} = {s3 and s3['N']}")
out = run('symalg.py', 'K2', 'pair_after', 'pair_srcout')
report("X3", out.count("conformal symmetry algebra) = 3, dim Lie Aut_Δ (c = 0) = 2") == 3,
       "for K₂, 1⇉2→3 and 1⇉2,1→3 the symmetry algebra of the fold equations (derivations with ΔD − (D⊗1+1⊗D)Δ = cΔ) is "
       "exactly Lie(G) ⊕ k·E: Aut_Δ° = G (dim 2) plus the path-length grading", out.strip().replace("\n", " | "))
s5 = census('pair_srcout', 5); g3, g5 = graded('pair_srcout', 3), graded('pair_srcout', 5)
cub = lambda q: 4 * q**3 + 25 * q**2 + 10 * q + 15
report("X4", s3 and s5 and (s3['N'], s5['N']) == (cub(3), cub(5)) and g3 == (105, 101) and g5 == (113, 101),
       "1⇉2, 1→3: N(q) = 4q³ + 25q² + 10q + 15 (378, 1190) — degree 3 > dim G = 2; G-orbits grow (105, 113) while "
       "(G × grading)-orbits stay 101: Conjecture 3′ fails, Conjecture 3″ holds here", f"G / G×grading orbits: F3 {g3}, F5 {g5}")
out = run('family_q.py')
report("X5", "all four families are genuine Aut_Δ-moduli removed by the grading: True" in out,
       "over Q(t): the four free (G × grading)-orbits are families P_t of folds whose G-invariant τ(P_t) (vertex part of P_t(γ)) "
       "equals t — a modulus for Aut_Δ that the grading rescales away", " | ".join(l for l in out.splitlines() if l.startswith('family')))
ga5 = graded('pair_after', 5); gk5 = graded('K2', 5)
report("X6", ga5 == (71, 71) and gk5 == (20, 20),
       "for K₂ and 1⇉2→3 the grading identifies no G-orbits (20 and 71 orbits either way)", f"{gk5}, {ga5}")
if SLOW >= 1:
    s7 = census('pair_srcout', 7); g7 = graded('pair_srcout', 7)
    report("X7", s7 and s7['N'] == cub(7) and g7 == (121, 101),
           "1⇉2, 1→3 over F_7: 2682 = 4·343 + 25·49 + 70 + 15 folds (predicted before the run), 121 = 4q + 93 G-orbits, 101 (G × grading)-orbits",
           f"{s7 and s7['N']} folds, orbits {g7}")
    out8 = {q: census(q, 3) for q in ('pair_after2', 'pair_mid')}
    report("X8", out8['pair_after2'] and (out8['pair_after2']['N'], out8['pair_after2']['orbits']) == (1057, 334)
           and out8['pair_mid'] and (out8['pair_mid']['N'], out8['pair_mid']['orbits']) == (696, 257),
           "4-vertex chains over F_3: 1⇉2→3→4 has 1057 folds in 334 G-orbits; 1→2⇉3→4 has 696 in 257",
           "; ".join(f"{q}: N={c['N']} orbits={c['orbits']} sizes={c['sizes']}" for q, c in out8.items() if c))
    gt, gb = graded('pair_tgtfork', 3), graded('pair_srcboth', 3)
    report("X9", gt == (428, 428) and gb == (511, 485),
           "fork test over F_3: arrows leaving the pair's target (1⇉2, 2→3, 2→4) give no families (428 = 428); "
           "a same-orientation arrow at both ends (1⇉2, 1→3, 4→2) gives 26 families (511 G-orbits, 485 with grading)", f"{gt}, {gb}")
if SLOW >= 2:
    X10 = {}
    for q in ('pair_after2', 'pair_mid'):
        c = census(q, 5); X10[q] = (c['N'], c['orbits']) if c else None
    gb5 = graded('pair_srcboth', 5)
    report("X10", X10 == {'pair_after2': (2749, 334), 'pair_mid': (1742, 257)} and gb5 == (563, 485),
           "over F_5: 1⇉2→3→4 has 2749 = 102q²+30q+49 folds in 334 orbits and 1→2⇉3→4 has 1742 = 61q²+35q+42 in 257 (constant, as "
           "predicted from F_3); the double fork has 563 = 26q+433 G-orbits but still 485 (G × grading)-orbits", f"{X10}, srcboth {gb5}")
print(f"\n{sum(res)}/{len(res)} checks pass")
