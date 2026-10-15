"""refinish.py — recompute the orbit statistics of finished batch jobs from their checkpoints (no fibre is recomputed).

v2.24 adds orbits_GxTloc (G-orbits merged along local class-torus moves). For every complete checkpoint in checkpoints/,
this reruns pcensus.py on it (the engine finds all fibres done and goes straight to the orbit stage) and appends the new
JSON line to results_v224.jsonl. Predictions made before the rerun are printed next to the results.

  python3 refinish.py                                # all checkpoints in checkpoints/
  python3 refinish.py two_pairs                      # only checkpoints whose name starts with this
  python3 refinish.py --dir data/checkpoints --out X # the committed raw data of the v2.23 batch (see data/README.md)"""
import sys, os, glob, json, subprocess, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
args = [a for a in sys.argv[1:]]
def opt(flag, default):
    if flag in args:
        i = args.index(flag); v = args[i + 1]; del args[i:i + 2]; return v
    return default
CKDIR = opt('--dir', 'checkpoints'); OUT = opt('--out', 'results_v224.jsonl')
PRED = {('pair_tgtfork', 5): {'orbits_G': 456, 'orbits_GxTloc': 414},     # F_3: 428 G-orbits, 414 with the local torus
        ('pair_srcout', 11): {'orbits_GxTloc': 101}, ('pair_srcboth', 7): {'orbits_GxTloc': 485},
        ('pair_after2', 7): {'orbits_GxTloc': 334}, ('pair_mid', 7): {'orbits_GxTloc': 257},
        ('K3_after', 5): {'orbits_GxTloc': 120}, ('K3_srcout', 5): {'orbits_GxTloc': 183}}
done = set()
if os.path.exists(OUT):
    for l in open(OUT):
        try: r = json.loads(l); done.add((r['name'], r['p']))
        except Exception: pass
sel = args[0] if args else ''
for ck in sorted(glob.glob(os.path.join(CKDIR, '*.pkl'))):
    base = os.path.basename(ck)[:-4]; name, p = base.rsplit('_', 1); p = int(p)
    if not name.startswith(sel) or (name, p) in done: continue
    c = pickle.load(open(ck, 'rb'))
    if len(c['fibres']) < c['nreps']: print(f"{base}: incomplete ({len(c['fibres'])}/{c['nreps']}), skipped", flush=True); continue
    print(f"=== {name} over F_{p}: re-finishing from {ck}", flush=True)
    r = subprocess.run([sys.executable, '-u', 'pcensus.py', name, str(p), '--engine', 'pool', '--jobs', '1',
                        '--checkpoint', ck, '--out', OUT], capture_output=True, text=True)
    if r.returncode: print(r.stderr[-1500:], flush=True); continue
    line = next(l for l in r.stdout.splitlines() if l.startswith('{')); res = json.loads(line)
    keys = ('N', 'orbits_G', 'orbits_GxGr', 'orbits_GxT', 'orbits_GxTloc')
    pr = PRED.get((name, p), {})
    print("    " + " ".join(f"{k}={res.get(k)}" + ("" if k not in pr else (" ok" if res.get(k) == pr[k] else f" (predicted {pr[k]}) MISMATCH")) for k in keys)
          + f"; local torus sizes {res.get('local_torus_sizes')}; torus {res.get('torus_preserving')}/{res.get('torus_all')}", flush=True)
