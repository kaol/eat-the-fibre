"""v2.71: for one quiver, run v2_71_gcert.py on every block tree (separation per point, direct checks at every residual
prime ≥ 3), and v2_68_uniform.py --p 3 on the trees whose ℤ Lemma C claims (logs/v2_70/<Q>/<b>.u.json) list 3.
Outputs into logs/v2_71/<QUIVER>/. The ℤ claims C1–C5 themselves are those of v2.70 (unchanged).
Usage: python3 analysis/v2_71_all.py QUIVER CERT [J]"""
import sys, os, json, glob, subprocess, time
Q, cert = sys.argv[1], sys.argv[2]; J = int(sys.argv[3]) if len(sys.argv) > 3 else 3
outd = f'logs/v2_71/{Q}'; os.makedirs(outd, exist_ok=True)
todo = []
for f in sorted(glob.glob(f'logs/v2_67/count/{Q}/*.json')):
    b = os.path.basename(f)[:-5]; blk = json.load(open(f))['block']
    if not os.path.exists(f'{outd}/{b}.gc.json'):
        todo.append(([sys.executable, 'analysis/v2_71_gcert.py', cert, blk, f, f'{outd}/{b}.gc'], f'{outd}/{b}.gc'))
    u = json.load(open(f'logs/v2_70/{Q}/{b}.u.json'))
    if 3 in u.get('candidates', []) and not os.path.exists(f'{outd}/{b}.u3.json'):
        todo.append(([sys.executable, 'analysis/v2_68_uniform.py', cert, blk, f, f'{outd}/{b}.u3', '--p', '3'], f'{outd}/{b}.u3'))
run = {}
while todo or run:
    while todo and len(run) < J:
        cmd, o = todo.pop(0)
        run[o] = subprocess.Popen(cmd, stdout=open(o + '.txt', 'w'), stderr=subprocess.STDOUT)
    time.sleep(0.5)
    for o in list(run):
        if run[o].poll() is not None: del run[o]
print(f"{Q}: done")
