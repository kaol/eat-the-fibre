"""v2.69: stabilisers on the leaf trees of a case-split block (logs/v2_67/leafunion/<DIR>/L*.json), H-stable roots
only, each orbit closure once (by its full-coordinate key).
Usage: python3 analysis/v2_69_stab_leaves.py CERT BLOCK LEAFDIR OUTDIR"""
import sys, os, json, glob, subprocess, collections
cert, block, ld, outd = sys.argv[1:5]; os.makedirs(outd, exist_ok=True)
seen = set(); tot = collections.Counter(); bad = []
for f in sorted(glob.glob(f'{ld}/L*.json'), key=lambda x: int(os.path.basename(x)[1:-5])):
    k = os.path.basename(f)[1:-5]; d = json.load(open(f))
    keep = {key: nd for key, nd in d['nodes'].items() if nd.get('result') == 'dense orbit' and nd.get('fullkey') not in seen}
    if not keep: continue
    for nd in keep.values(): seen.add(nd['fullkey'])
    tf = f'{outd}/L{k}_tree.json'; json.dump(dict(nodes=keep, roots=[]), open(tf, 'w'))
    o = f'{outd}/L{k}'
    subprocess.run([sys.executable, 'analysis/v2_69_stab.py', cert, block, tf, o], env=dict(os.environ, STAB_LEAF=k),
                   stdout=open(o + '.txt', 'w'), stderr=subprocess.STDOUT)
    if os.path.exists(o + '.json'):
        r = json.load(open(o + '.json')); tot.update(r['tally'])
        bad += [n for n in r['nodes'] if n.get('stab_components', 1) != 1 or 'result' in n]
    else: bad.append(('no output', k))
print(f"{block}: {dict(tot)}; problems {bad[:3]}")
