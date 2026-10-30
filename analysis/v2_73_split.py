"""v2.73: ℤ-certificates for a case-split block (A–E of logs/v2_73/PREDICTIONS.md).
Roots: the H-stable leaf primes (the dense roots of the v2.67 leaf trees, by full-coordinate key); one representative
leaf chart per root. For every leaf: v2_73_leafcheck.py (A, B, C, point leaves). For every representative chart: the
root's subtree as a tree file, v2_68_uniform.py (C1–C4) and v2_71_gcert.py (node certificates, no exempt prime), both in
leaf mode.
Usage: python3 analysis/v2_73_split.py CERT BLOCK LEAFDIR OUTDIR [J]"""
import sys, os, json, glob, subprocess, time
cert, block, ld, outd = sys.argv[1:5]; J = int(sys.argv[5]) if len(sys.argv) > 5 else 2
os.makedirs(outd, exist_ok=True)
leaves = {}
for l in open(cert):
    if f'"block": "{block}"' not in l: continue
    r = json.loads(l)
    if r['kind'] == 'stratum' and r['block'] == block: leaves[r['leaf']] = r
roots = []; rep = {}
for f in sorted(glob.glob(f'{ld}/L*.json'), key=lambda x: int(os.path.basename(x)[1:-5]) if os.path.basename(x)[1:-5].isdigit() else 10**9):
    k = os.path.basename(f)[1:-5]
    if not k.isdigit(): continue
    T = json.load(open(f))
    for rk in T['roots']:
        nd = T['nodes'].get(rk) if rk is not None else None          # '' is a root too: the zero ideal of its chart
        if not nd or nd.get('result') != 'dense orbit': continue
        if nd['fullkey'] not in roots: roots.append(nd['fullkey']); rep[nd['fullkey']] = (int(k), rk)
json.dump(dict(roots=roots), open(f'{outd}/roots.json', 'w'))
print(f"{block}: {len(leaves)} leaves, {len(roots)} H-stable roots (representative leaves {sorted(set(v[0] for v in rep.values()))})", flush=True)
jobs = []
for k in sorted(leaves):
    o = f'{outd}/leaf{k}'
    if not os.path.exists(o + '.json'):
        jobs.append(([sys.executable, 'analysis/v2_73_leafcheck.py', cert, block, str(k), f'{outd}/roots.json', o], o, {}))
for key, (k, rk) in rep.items():
    T = json.load(open(f'{ld}/L{k}.json'))
    keep = {}; stack = [rk]
    while stack:
        x = stack.pop()
        if x in keep or x not in T['nodes']: continue
        keep[x] = T['nodes'][x]; stack += T['nodes'][x].get('children', [])
    tf = f'{outd}/tree_L{k}_{list(rep).index(key)}.json'
    json.dump(dict(nodes=keep, roots=[rk]), open(tf, 'w'))
    for tag, scr in (('u', 'analysis/v2_68_uniform.py'), ('gc', 'analysis/v2_71_gcert.py')):
        o = tf[:-5] + '.' + tag
        if not os.path.exists(o + '.json'):
            jobs.append(([sys.executable, scr, cert, block, tf, o], o, {'STAB_LEAF': str(k), 'GCERT_EXEMPT': ''}))
run = {}
while jobs or run:
    while jobs and len(run) < J:
        cmd, o, env = jobs.pop(0)
        run[o] = subprocess.Popen(cmd, stdout=open(o + '.txt', 'w'), stderr=subprocess.STDOUT, env=dict(os.environ, **env))
    time.sleep(0.5)
    for o in list(run):
        if run[o].poll() is not None: del run[o]
print("done", flush=True)
