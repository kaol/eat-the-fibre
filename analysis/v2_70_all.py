"""v2.70: all certificates for one quiver: per block tree of logs/v2_67/count/<QUIVER>/, run
  v2_70_gcert.py (G1, S, W, separation, point nodes), v2_68_uniform.py (C1–C4), v2_68_cover.py (C5),
J processes at a time, into logs/v2_70/<QUIVER>/; then tally the candidate primes.
Usage: python3 analysis/v2_70_all.py QUIVER CERT [J]"""
import sys, os, json, glob, subprocess, time, collections
Q, cert = sys.argv[1], sys.argv[2]; J = int(sys.argv[3]) if len(sys.argv) > 3 else 3
outd = f'logs/v2_70/{Q}'; os.makedirs(outd, exist_ok=True)
SCRIPTS = {'gc': 'analysis/v2_70_gcert.py', 'u': 'analysis/v2_68_uniform.py', 'cv': 'analysis/v2_68_cover.py'}
todo = []
for f in sorted(glob.glob(f'logs/v2_67/count/{Q}/*.json')):
    b = os.path.basename(f)[:-5]; blk = json.load(open(f))['block']
    for tag, scr in SCRIPTS.items():
        o = f'{outd}/{b}.{tag}'
        if not os.path.exists(o + '.json'): todo.append((scr, blk, f, o))
run = {}
while todo or run:
    while todo and len(run) < J:
        scr, blk, f, o = todo.pop(0)
        run[o] = subprocess.Popen([sys.executable, scr, cert, blk, f, o], stdout=open(o + '.txt', 'w'), stderr=subprocess.STDOUT)
    time.sleep(0.5)
    for o in list(run):
        if run[o].poll() is not None: del run[o]
cands = collections.Counter(); missing = []; problems = []; tally = collections.Counter()
for f in sorted(glob.glob(f'logs/v2_67/count/{Q}/*.json')):
    b = os.path.basename(f)[:-5]
    for tag in SCRIPTS:
        o = f'{outd}/{b}.{tag}.json'
        if not os.path.exists(o): missing.append(f'{b}.{tag}'); continue
        d = json.load(open(o))
        cs = [str(c) for c in d.get('candidates', [])]
        for c in cs: cands[(tag, c)] += 1
        if tag == 'gc': tally.update(d['tally'])
        if tag == 'u' and d.get('problems'): problems.append((b, 'u', d['problems'][:2]))
        if tag == 'cv' and (d.get('power_fail') and d.get('rab_fail')): problems.append((b, 'cv', d))
print(f"{Q}: gcert tally {dict(tally)}; candidates (script, prime): blocks {dict(sorted(cands.items()))}; "
      f"primes {sorted(set(c for _, c in cands))}; problems {len(problems)} {problems[:3]}; missing {len(missing)} {missing[:5]}")
