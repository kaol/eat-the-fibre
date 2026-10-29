"""v2.70: tally all certificates per quiver. Hard failures (any of them voids the certificate for that block):
gcert: a node without a certified point (candidates None), a failed check; uniform: an issue other than an integer
denominator (dimension, f ∈ I, radical membership failing outright, defect); cover: power and Rabinowitsch both failing.
Candidate primes: the union over blocks of gcert, uniform and cover candidates. Also cross-checks every direct char-p G1
check (it assumes the C1 basis of the node's block has no leading coefficient divisible by that prime).
Usage: python3 analysis/v2_70_tally.py [QUIVER ...]"""
import sys, os, json, glob, collections
QS = sys.argv[1:] or sorted(os.path.basename(d.rstrip('/')) for d in glob.glob('logs/v2_67/count/*/'))
SOFT = ('radical membership needs 1/N', '1 ∉ J (empty boundary)')
grand = collections.Counter(); out = {}
for Q in QS:
    cand = collections.Counter(); hard = []; nodes = collections.Counter(); direct = []
    for f in sorted(glob.glob(f'logs/v2_67/count/{Q}/*.json')):
        b = os.path.basename(f)[:-5]; base = f'logs/v2_70/{Q}/{b}'
        try: gc = json.load(open(base + '.gc.json')); u = json.load(open(base + '.u.json')); cv = json.load(open(base + '.cv.json'))
        except FileNotFoundError as e: hard.append((b, 'missing', str(e)[-30:])); continue
        lci = set()
        for r in u['results'].values():
            for c in r.get('lci') or []:
                if abs(int(c)) > 1: lci |= {str(p) for p in __import__('sympy').factorint(abs(int(c)))}
            for i in r.get('issues') or []:
                if not i.startswith(SOFT): hard.append((b, 'u', i))
        if u.get('problems'): hard.append((b, 'u problems', u['problems'][:2]))
        if cv.get('power_fail') and cv.get('rab_fail'): hard.append((b, 'cv', cv.get('rab_fail')))
        for n in gc['nodes']:
            nodes['points' if n['dim'] == 0 else 'dense'] += 1
            if n.get('candidates') is None: hard.append((b, 'gc', n['path'], n.get('checks') or [p['checks'] for p in n.get('points', [])][:1])); continue
            for p in n['candidates']: cand[str(p)] += 1
            for q, *_ in n.get('direct_G1', []) + n.get('direct', []):
                direct.append((b, n['path'], q))
                if str(q) in lci: hard.append((b, 'direct G1 at a C1 prime', n['path'], q))
        for p in u.get('candidates', []) + cv.get('candidates', []): cand[str(p)] += 1
    out[Q] = dict(nodes=dict(nodes), candidates=dict(sorted(cand.items())), hard=hard, direct=direct)
    grand.update(nodes)
    print(f"{Q}: nodes {dict(nodes)}; candidate primes {sorted(cand, key=int)} (node/block counts {dict(sorted(cand.items()))}); "
          f"direct char-p G1 checks {len(direct)} {sorted(set(q for _, _, q in direct))}; hard failures {len(hard)} {hard[:3]}")
json.dump(out, open('logs/v2_70/tally.json', 'w'), indent=1)
print(f"TOTAL nodes {dict(grand)}; primes {sorted(set(p for o in out.values() for p in o['candidates']), key=int)}; "
      f"hard failures {sum(len(o['hard']) for o in out.values())}")
