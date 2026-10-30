"""v2.72: as v2_71_all.py for a prime P: v2_71_gcert.py with no exempt prime (GCERT_EXEMPT empty: direct checks at
every residual prime, P included), and v2_68_uniform.py --p P on the trees whose ℤ claims list P. Outputs into
logs/v2_72/<QUIVER>/.
Usage: python3 analysis/v2_72_all.py QUIVER CERT P [J]"""
import sys, os, json, glob, subprocess, time
Q, cert, P = sys.argv[1], sys.argv[2], int(sys.argv[3]); J = int(sys.argv[4]) if len(sys.argv) > 4 else 3
os.environ['GCERT_EXEMPT'] = ''
outd = f'logs/v2_72/{Q}'; os.makedirs(outd, exist_ok=True)
todo = []
for f in sorted(glob.glob(f'logs/v2_67/count/{Q}/*.json')):
    b = os.path.basename(f)[:-5]; blk = json.load(open(f))['block']
    if not os.path.exists(f'{outd}/{b}.gc.json'):
        todo.append(([sys.executable, 'analysis/v2_71_gcert.py', cert, blk, f, f'{outd}/{b}.gc'], f'{outd}/{b}.gc'))
    u = json.load(open(f'logs/v2_70/{Q}/{b}.u.json'))
    if P in u.get('candidates', []) and not os.path.exists(f'{outd}/{b}.u{P}.json'):
        todo.append(([sys.executable, 'analysis/v2_68_uniform.py', cert, blk, f, f'{outd}/{b}.u{P}', '--p', str(P)], f'{outd}/{b}.u{P}'))
run = {}
while todo or run:
    while todo and len(run) < J:
        cmd, o = todo.pop(0)
        run[o] = subprocess.Popen(cmd, stdout=open(o + '.txt', 'w'), stderr=subprocess.STDOUT)
    time.sleep(0.5)
    for o in list(run):
        if run[o].poll() is not None: del run[o]
print(f"{Q}: done")
