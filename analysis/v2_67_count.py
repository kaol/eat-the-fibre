"""v2.67: geometric H-orbit count per vertex block (v2_67_orbits.py on all components of the block) against the census
local-torus classes per block (v2_53_block_census.py output at one prime).
Usage: python3 analysis/v2_67_count.py CERT BLOCKCENSUS OUTDIR [J]"""
import sys, os, json, lzma, re, subprocess, collections, time
cert, bc, out = sys.argv[1:4]; J = int(sys.argv[4]) if len(sys.argv) > 4 else 2
os.makedirs(out, exist_ok=True)
op = lzma.open if cert.endswith('.xz') else open
blocks = collections.OrderedDict()
for l in op(cert, 'rt'):
    r = json.loads(l)
    blocks.setdefault(r['block'], set()).add(r['leaf'])
census = {}
for l in open(bc).readlines()[1:-1]:
    v, cells = l.split('|')[0].strip(), l.split('|')[1].split()
    census[','.join(row.replace('-1', 'm') for row in v.split())] = int(cells[2])
todo = [b for b in blocks if blocks[b] == {0}]
skipped = [b for b in blocks if blocks[b] != {0}]
run = {}; res = {}
def result(b):
    f = os.path.join(out, b.replace(',', '_') + '.txt')
    m = re.search(r"RESULT (.*?): .*geometric H-orbits (\w+)", open(f).read()) if os.path.exists(f) else None
    return (m.group(1), m.group(2)) if m else None
for b in todo:
    if result(b): res[b] = result(b)
todo = [b for b in todo if b not in res]
while todo or run:
    while todo and len(run) < J:
        b = todo.pop(0); f = os.path.join(out, b.replace(',', '_'))
        run[b] = subprocess.Popen([sys.executable, 'analysis/v2_67_orbits.py', cert, b, '0', 'all', f], stdout=open(f + '.txt', 'w'),
                                  stderr=subprocess.STDOUT)
    time.sleep(1)
    for b in list(run):
        if run[b].poll() is not None:
            res[b] = result(b) or ('ERROR', 'None'); del run[b]
tot_o = tot_c = 0; bad = []
for b in blocks:
    if b in skipped: print(f"{b:24} case split: skipped"); continue
    st, n = res[b]; c = census.get(b)
    ok = st == 'finitely many orbits' and n != 'None' and c is not None and int(n) == c
    if not ok: bad.append(b)
    tot_o += int(n) if n != 'None' else 0; tot_c += c or 0
    print(f"{b:24} {st:22} orbits {n:>4}  census {c}  {'ok' if ok else 'MISMATCH'}")
print(f"TOTAL orbits {tot_o} census {tot_c}; blocks {len(blocks)}, skipped {len(skipped)}, mismatches {len(bad)} {bad}")
