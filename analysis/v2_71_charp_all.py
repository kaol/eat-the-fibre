"""v2.71: for every block of a quiver that fails at p (logs/v2_71/tally.json: C1–C4 directly at p, or node certificates),
build the char-p tree (v2_68_charp.py) and run Lemma L directly at p on it (v2_71_charp_cert.py); compare the
certified count with the char 0 count of the block (dense nodes + points of the v2.67 tree).
Usage: python3 analysis/v2_71_charp_all.py QUIVER CERT P"""
import sys, os, json, subprocess
Q, cert, P = sys.argv[1], sys.argv[2], int(sys.argv[3])
tal = json.load(open('logs/v2_71/tally.json'))[Q]
outd = f'logs/v2_71/charp/{Q}'; os.makedirs(outd, exist_ok=True)
rows = []
for b in sorted({x[0] for x in tal['fail3']} | set(tal.get('resolved3', []))):
    f = f'logs/v2_67/count/{Q}/{b}.json'; T0 = json.load(open(f)); blk = T0['block']
    c0 = sum(1 for n in T0['nodes'].values() if n.get('result') == 'dense orbit') + \
         sum(n.get('geometric', 1) for n in T0['nodes'].values() if n.get('result') == 'points')
    tr, ce = f'{outd}/{b}_p{P}', f'{outd}/{b}_p{P}_cert'
    if not os.path.exists(tr + '.json'):
        subprocess.run([sys.executable, 'analysis/v2_71_charp_tree.py', cert, blk, str(P), f, tr], stdout=open(tr + '.txt', 'w'), stderr=subprocess.STDOUT)
    if not os.path.exists(ce + '.json') and os.path.exists(tr + '.json'):
        subprocess.run([sys.executable, 'analysis/v2_71_charp_cert.py', cert, blk, str(P), tr + '.json', ce], stdout=open(ce + '.txt', 'w'), stderr=subprocess.STDOUT)
    r = json.load(open(ce + '.json'))['tally'] if os.path.exists(ce + '.json') else None
    st = json.load(open(tr + '.json'))['stats'] if os.path.exists(tr + '.json') else None
    rows.append(dict(block=b, char0=c0, tree=st, cert=r, ok=bool(r) and r['failed'] == 0 and r['certified'] == c0))
    print(f"  {b}: char 0 count {c0}; char-{P} tree {st}; certificates {r}; {'AGREES' if rows[-1]['ok'] else 'CHECK'}", flush=True)
json.dump(rows, open(f'{outd}/summary_p{P}.json', 'w'), indent=1)
print(f"{Q} at p = {P}: {sum(r['ok'] for r in rows)}/{len(rows)} blocks certified with the char 0 count")
