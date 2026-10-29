"""v2.69: run v2_69_stab.py on every block tree of logs/v2_67/count/<QUIVER>/ (J at a time) and tally.
Usage: python3 analysis/v2_69_stab_all.py QUIVER CERT [J]"""
import sys, os, json, glob, subprocess, time, collections
Q, cert = sys.argv[1], sys.argv[2]; J = int(sys.argv[3]) if len(sys.argv) > 3 else 2
outd = f'logs/v2_69/stab/{Q}'; os.makedirs(outd, exist_ok=True)
todo = []
for f in sorted(glob.glob(f'logs/v2_67/count/{Q}/*.json')):
    b = os.path.basename(f)[:-5]; d = json.load(open(f))
    if not any(n.get('result') == 'dense orbit' for n in d['nodes'].values()): continue
    if os.path.exists(f'{outd}/{b}.json'): continue
    todo.append((d['block'], f, f'{outd}/{b}'))
run = {}
while todo or run:
    while todo and len(run) < J:
        blk, f, o = todo.pop(0)
        run[o] = subprocess.Popen([sys.executable, 'analysis/v2_69_stab.py', cert, blk, f, o], stdout=open(o + '.txt', 'w'), stderr=subprocess.STDOUT)
    time.sleep(0.5)
    for o in list(run):
        if run[o].poll() is not None: del run[o]
tot = collections.Counter(); bad = []
for f in glob.glob(f'{outd}/*.json'):
    d = json.load(open(f)); tot.update(d['tally'])
    for n in d['nodes']:
        if n.get('stab_components', 1) != 1 or n.get('stab_dims') and max(n['stab_dims']) != n['expected'] or 'result' in n: bad.append((d['block'], n))
missing = [os.path.basename(f) for f in glob.glob(f'logs/v2_67/count/{Q}/*.json') if not os.path.exists(f'{outd}/' + os.path.basename(f))]
print(f"{Q}: {dict(tot)}; problems {len(bad)} {bad[:3]}; blocks without output {len(missing)}")
