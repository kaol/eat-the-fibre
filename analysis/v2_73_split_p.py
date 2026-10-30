"""v2.73: a case-split block at a prime P where the char 0 root trees do not specialise: for each representative chart
(an H-stable root R, leaf k), the char-P tree of R (roots: minimal primes over 𝔽_P of R's leaf prime in the chart) and
Lemma L at its nodes (v2_71_charp_cert.py, Lemma N); then the certified nodes of all charts, deduplicated in full
coordinates (two certified nodes are one when each one's point lies on the other's closure), plus the fixed point leaves.
Usage: python3 analysis/v2_73_split_p.py CERT BLOCK OUTDIR P"""
import sys, os, json, glob, subprocess
cert, block, outd, P = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
comps = {}
for l in open(cert):
    if f'"block": "{block}"' not in l: continue
    r = json.loads(l)
    if r['kind'] == 'component' and r['block'] == block: comps[r['leaf']] = r
for tf in sorted(glob.glob(f'{outd}/tree_L*_*.json')):
    if tf.endswith(('.gc.json', '.u.json', '.full.json', '.u2.json')) or '_p' in os.path.basename(tf): continue
    base = tf[:-5]; k = int(os.path.basename(base).split('_')[1][1:])
    gens = [g for g in comps[k]['gens'] if g.strip() not in ('', '0')] or ['0']
    env = dict(os.environ, STAB_LEAF=str(k), CHARP_ROOTGENS=json.dumps(gens))
    tr, ce = f'{base}_p{P}', f'{base}_p{P}_cert'
    if not os.path.exists(tr + '.json'):
        subprocess.run([sys.executable, 'analysis/v2_71_charp_tree.py', cert, block, str(P), tf, tr], stdout=open(tr + '.txt', 'w'), stderr=subprocess.STDOUT, env=env)
    if not os.path.exists(ce + '.json'):
        subprocess.run([sys.executable, 'analysis/v2_71_charp_cert.py', cert, block, str(P), tr + '.json', ce], stdout=open(ce + '.txt', 'w'), stderr=subprocess.STDOUT, env=env)
    print(f"  {os.path.basename(base)}: {json.load(open(ce + '.json'))['tally'] if os.path.exists(ce + '.json') else 'no output'}", flush=True)
