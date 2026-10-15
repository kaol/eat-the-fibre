"""verify_v2_23.py — machine checks Z1–Z8 for discordianism-open-questions-explored-v2.23.md
Install:  pip install -r requirements.txt.   JOBS=n sets worker processes (default: all cores).
Default: Z1–Z4 (minutes).  SLOW=1 adds Z5–Z8 (K₃ + arrow and 1⇉2⇉3 over F_3, exact tangent spaces and the numerical
local-dimension test; about 1–2 hours on 2 cores, mostly Z8, which is single-threaded).
The long censuses for larger primes are in run_batch.py, not here."""
import sys, os, subprocess, json
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE); os.makedirs('logs', exist_ok=True)
SLOW = int(os.environ.get('SLOW', '0') or 0); JOBS = os.environ.get('JOBS', str(os.cpu_count()))
res = []
def report(tag, ok, claim, detail=""):
    res.append(bool(ok)); print(f"[{tag}] {'PASS' if ok else 'FAIL'}  {claim}", flush=True)
    if detail: print(f"       {detail}", flush=True)
def run(*args):
    r = subprocess.run([sys.executable, '-u', *args], capture_output=True, text=True)
    if r.returncode: print(f"  !! {' '.join(args)} failed:", *r.stderr.strip().splitlines()[-3:], sep="\n     ")
    return r.stdout
def census(q, p, save=False):
    args = ['pcensus.py', q, str(p), '--jobs', JOBS] + (['--save-folds', f'logs/folds_{q}_{p}.pkl'] if save else [])
    line = next((l for l in run(*args).splitlines() if l.startswith('{')), None)
    return json.loads(line) if line else None
sig = lambda r: (r['N'], r['orbits_G'], r['orbits_GxGr'], r['orbits_GxT'], r['transversal'].get('transversal_modp', 0)) if r else None

k2, pa, ps = census('K2', 3), census('pair_after', 3, save=True), census('pair_srcout', 3, save=True)
report("Z1", sig(k2) == (64, 20, 20, 20, 0) and sig(pa) == (219, 71, 71, 71, 0) and sig(ps) == (378, 105, 101, 101, 4),
       "parallel census (pcensus.py) reproduces v2.18/v2.21/v2.22 over F_3; orbits under G, G×grading, G×T_Q; transversal "
       "(one-parameter-family) orbits 0, 0, 4 for K₂, 1⇉2→3, 1⇉2,1→3", f"{sig(k2)}, {sig(pa)}, {sig(ps)}  (N, G, G×Gr, G×T, transversal)")
k3 = census('K3', 3)
report("Z2", k3 and (k3['N'], k3['orbits_G'], k3['G']) == (1300, 34, 432),
       "K₃ over F_3: 1300 folds in 34 orbits of the order-432 column-sum-1 group (v2.18)", str(sig(k3)))
pa5, ps5 = census('pair_after', 5), census('pair_srcout', 5)
report("Z3", pa5 and ps5 and (pa5['orbits_GxT'], pa5['torus_preserving'], pa5['torus_all']) == (71, 16, 16)
       and (ps5['orbits_G'], ps5['orbits_GxT'], ps5['torus_preserving'], ps5['torus_all']) == (113, 101, 4, 16),
       "over F_5 the class-weight torus preserves all folds of 1⇉2→3 (16/16 weightings; still 71 orbits) but only the grading "
       "(4/16) for 1⇉2,1→3 (101 orbits)", f"{sig(pa5)} torus {pa5 and pa5['torus_preserving']}; {sig(ps5)} torus {ps5 and ps5['torus_preserving']}")
o1, o2, o3 = run('splitdelta.py', 'pair_after', '3'), run('splitdelta.py', 'pair_srcout', '3'), run('classmix.py')
t4 = run('tangent.py', 'pair_srcout', '3')
report("Z4", "class (0, 1): 219/219" in o1 and "class (1, 2): 219/219" in o1 and "class (0, 1): 346/378" in o2 and "[((2, 3, 3), 32)]" in t4 and "family folds all split compatible: True; split failure ⇔ mixing beyond ker Δ: True" in o3,
       "split compatibility: every fold of 1⇉2→3 respects each class deletion Δ_c; for 1⇉2,1→3, 32 do not — exactly the folds mixing the two classes beyond ker Δ; all family folds are split compatible; "
       "at the transversal family Q-folds of 1⇉2,1→3 the tangent dimension equals r_H = 3", " | ".join(l for l in (o1 + o2).splitlines() if 'class' in l) + " | " + t4.strip() + " | " + o3.strip().splitlines()[-1])
if SLOW >= 1:
    ka, ks = census('K3_after', 3, save=True), census('K3_srcout', 3)
    o5 = run('splitdelta.py', 'K3_after', '3')
    report("Z5", sig(ka) and sig(ka)[:3] == (5439, 120, 120) and sig(ka)[4] == 0 and sig(ks) and sig(ks)[:3] == (10410, 191, 183)
           and sig(ks)[4] == 8 and "class (0, 1): 5439/5439" in o5,
           "K₃ + one arrow over F_3: 1⇛2→3 has 5439 folds, 120 orbits, no families, and is split compatible; 1⇛2,1→3 has "
           "10410 folds, 191 G-orbits, 183 with grading and 8 transversal orbits (the fork pattern)", f"{sig(ka)}, {sig(ks)}")
    tp = census('two_pairs', 3, save=True)
    report("Z6", tp and sig(tp)[:3] == (16376, 814, 728) and sig(tp)[3] == 716 and sig(tp)[4] == 82 and tp['transversal'].get('lifted_transversal', 0) > 0,
           "1⇉2⇉3 over F_3: 16376 folds; 814 G-orbits, 728 with grading, 716 with the class torus; 82 transversal orbits (some exact over Q)",
           f"{sig(tp)}, torus {tp and (tp['torus_preserving'], tp['torus_all'])}, {tp and tp['transversal']}")
    o7 = run('splitdelta.py', 'two_pairs', '3')
    report("Z7", "class (0, 1): 16376/16376" in o7 and "class (1, 2): 16376/16376" in o7,
           "every fold of 1⇉2⇉3 respects each class deletion separately, so the class torus 𝔾_m² acts (Lemma)",
           " | ".join(l for l in o7.splitlines() if 'class' in l))
    o8 = run('torus.py', 'two_pairs', '3', '24')
    gen = [l for l in o8.splitlines() if 'generic nearby verdicts' in l]
    ok8 = "((4, 5, 6, 6, True), 576)" in o8 and gen and all(
        seg.split('local ')[1].split(',')[0] == seg.split('r_T ')[1].split(',')[0] for seg in gen[0].split("generic: ")[1:])
    report("Z8", ok8,
           "1⇉2⇉3: at 576 transversal Q-folds r_T = tangent = 6 (exact): a smooth 6-dimensional locus, so the 5-dimensional "
           "Aut_Δ ⋊ grading has infinitely many orbits (Conjecture 3″ fails); at generic points near every other class, "
           "local dim = r_T (numerical), so Aut_Δ ⋊ T_Q accounts for all local families found", o8.strip())
print(f"\n{sum(res)}/{len(res)} checks pass")
