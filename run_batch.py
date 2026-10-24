"""run_batch.py — the long censuses (v2.23; jobs 12–15 added in v2.28, job 16 in v2.42).

Each job runs pcensus.py with the torch fibre engine (fibre_torch.py: GPU via ROCm/CUDA when torch sees one, else CPU),
a per-job checkpoint in checkpoints/, and appends one JSON line to results.jsonl. Rerunning skips finished jobs and
resumes interrupted ones, so it is safe to stop at any time (Ctrl-C) and start again. Send back results.jsonl (and
batch.log). Checkpoints are engine-independent: a job started with one engine can be resumed with the other.

Env: ENGINE=torch|pool (default torch), DEVICE=cuda|cpu|... (torch engine; default: GPU if available),
     JOBS=n (pool engine: worker processes, default all cores), ETF_MEM / ETF_GROUP (see fibre_torch.py).

  python3 run_batch.py              # all jobs, cheapest first
  python3 run_batch.py 3 5          # only jobs 3 and 5 (numbers as in the listing printed at start)
  python3 run_batch.py --list       # just list the jobs and predictions
  python3 run_batch.py 16 --save-reps   # also write logs/reps_QUIVER_p.pkl (G-orbit reps + local-torus labels); a job
                                        # counts as done only when its reps file exists too, and a finished job is
                                        # redone from its checkpoint (fibres are not recomputed)

Predictions (made before these runs, from smaller primes) are printed next to each result: "ok" means it matches."""
import sys, os, json, subprocess, time
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
JOBS = os.environ.get('JOBS', str(os.cpu_count()))
ENGINE = os.environ.get('ENGINE', 'torch')
DEVICE = os.environ.get('DEVICE')
def engine_args():
    if ENGINE == 'pool': return ['--engine', 'pool', '--jobs', JOBS]
    return ['--engine', 'torch'] + (['--device', DEVICE] if DEVICE else [])
# (quiver, p, predictions {field: value} or {}, note[, extra env for pcensus.py])
BATCH = [
    ('pair_tgtfork', 5, {'orbits_G': 428, 'orbits_GxGr': 428}, "1⇉2, 2→3, 2→4: stalled here at q=5; fork rule predicts no families"),
    ('pair_srcboth', 7, {'N': 14928, 'orbits_G': 615, 'orbits_GxGr': 485}, "double fork: 26q³+115q²+45q+60, G-orbits 26q+433, 485"),
    ('pair_after2', 7, {'N': 5257, 'orbits_G': 334}, "1⇉2→3→4: 102q²+30q+49"),
    ('pair_mid', 7, {'N': 3276, 'orbits_G': 257}, "1→2⇉3→4: 61q²+35q+42"),
    ('pair_srcout', 11, {'N': 8474, 'orbits_G': 137, 'orbits_GxGr': 101}, "1⇉2, 1→3: 4q³+25q²+10q+15, 4q+93, 101"),
    ('K3_srcout', 5, {}, "1⇛2, 1→3: 8 families at q=3"),
    ('K3_after', 5, {}, "1⇛2→3: no families at q=3 (120 orbits)"),
    ('two_pairs_sep', 3, {}, "1⇉2→3⇉4: new"),
    ('pair_after3', 3, {}, "1⇉2→3→4→5: new"),
    ('two_pairs', 5, {}, "1⇉2⇉3: heaviest job (dim K = 6, ~197k kernel folds; tens of core-hours) — decides whether its families persist"),
    ('two_pairs', 7, {}, "1⇉2⇉3 over F_7 (v2.24): G×T orbits 716 = F_3 count if the extra F_5 orbits are defined over Q(i); "
                     "792 if they are one-parameter families (Conjecture 3*ˡ fails). Kernel stage: 62M subspaces, now batched on the GPU"),
    # v2.28 (predictions committed in docs/explorations/...-v2.28.md before these runs; jobs 1–11 are already done)
    ('triple_pair_in_short', 5, {'orbits_GxTloc': 752}, "P2: 1⇛2⇇3 + 1→3; local-torus count as over F_3 (752). OOM in the "
                     "container (~5.7e5 torus orbits, 5080 H-orbits after symmetry)"),
    ('pairs_in_short_tail', 3, {}, "P4 (part 1): 1⇉2⇇3 + 1→3 + 2→4, d = 18; too slow in the container"),
    ('pairs_in_short_tail', 5, {}, "P4 (part 2): orbits_GxTloc should equal job 13's. v2.51: per-fibre column-order portfolio "
                     "(the default order ran out of memory inside one fibre; no single order suits every fibre)",
     {'ETF_CORDERS': '4,6,8,14,15,9,10,12,11,13;4,6,8,14,15,9,12,16,13,17,10,11'}),
    ('pairs_in_short', 11, {'orbits_GxTloc': 366}, "1⇉2⇇3 + 1→3: 366 over F_3, F_5, F_7; one more prime"),
    # v2.42 (prediction committed in logs/v2_42/PREDICTIONS.md before this run): needs --save-reps for the shear merge
    ('two_pairs_short', 7, {'orbits_GxTloc': 2012}, "v2.42: 1⇉2⇉3 + 1→3 over F_7, linear in q (1668, 1840, 2012); "
                     "merged with c ↦ c + s·δ_aδ_b, 873 + 22q = 1027. OOM in the container (est. 1e7 folds)"),
]
SAVE_REPS = '--save-reps' in sys.argv          # (before v2.43 this flag was accepted but silently not passed on)
def reps_path(q, p): return os.path.join('logs', f'reps_{q}_{p}.pkl')
def complete(q, p, res): return (q, p) in res and (not SAVE_REPS or os.path.exists(reps_path(q, p)))
def done():
    out = {}
    if os.path.exists('results.jsonl'):
        for line in open('results.jsonl'):
            try: r = json.loads(line); out[(r['name'], r['p'])] = r
            except Exception: pass
    return out
def show(i, q, p, pred, note, r=None):
    s = f"[{i}] {q} over F_{p}: {note}"
    if r:
        got = {k: r.get(k) for k in ('N', 'orbits_G', 'orbits_GxGr', 'orbits_GxT', 'orbits_GxTloc') if k in r}
        chk = " ".join(f"{k}={got[k]}" + ("" if k not in pred else (" ok" if got[k] == pred[k] else f" (predicted {pred[k]}) MISMATCH")) for k in got)
        s += f"\n      {chk}; transversal={r.get('transversal')}; t={r.get('t')}s"
    else:
        s += "  — predicted " + (", ".join(f"{k}={v}" for k, v in pred.items()) or "nothing")
    print(s, flush=True)
if __name__ == '__main__':
    sel = [int(x) for x in sys.argv[1:] if x.isdigit()]
    res = done()
    for i, (q, p, pred, note, *_) in enumerate(BATCH, 1):
        show(i, q, p, pred, note, res.get((q, p)))
    if '--list' in sys.argv: sys.exit()
    os.makedirs('checkpoints', exist_ok=True)
    for i, (q, p, pred, note, *rest) in enumerate(BATCH, 1):
        env = dict(os.environ, **(rest[0] if rest else {}))
        if sel and i not in sel: continue
        if complete(q, p, done()): continue
        print(f"\n=== job {i}: {q} over F_{p} ({time.strftime('%H:%M:%S')})"
              + "".join(f" {k}={v}" for k, v in (rest[0] if rest else {}).items()) + " ===", flush=True)
        with open('batch.log', 'a') as log:
            pr = subprocess.Popen([sys.executable, '-u', 'pcensus.py', q, str(p), *engine_args(),
                                   '--checkpoint', f'checkpoints/{q}_{p}.pkl', '--out', 'results.jsonl',
                                   *(['--save-reps', reps_path(q, p)] if SAVE_REPS else [])],
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env)
            for line in pr.stdout:
                log.write(line); log.flush()
                if not line.startswith('{'): print("   ", line, end="", flush=True)
            if pr.wait(): print(f"    job {i} exited with code {pr.returncode} (see batch.log)", flush=True)
        show(i, q, p, pred, note, done().get((q, p)))
