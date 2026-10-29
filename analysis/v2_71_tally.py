"""v2.71: tally per quiver. A block holds at a prime p ≥ 3 when
  * gcert (logs/v2_71): every node certified, p not among its candidates (after the direct checks);
  * C1–C4 (logs/v2_70 *.u.json, ℤ): p not a candidate, or (p = 3) the direct check *.u3.json has no issue;
  * C5 (logs/v2_70 *.cv.json): p not a candidate;
  * every direct check at q was made at a q prime to the block's C1 leading coefficients.
Prints, per quiver, the primes ≥ 3 that remain, the blocks failing at 3, and hard failures.
Usage: python3 analysis/v2_71_tally.py [QUIVER ...]"""
import sys, os, json, glob, collections
import sympy as sp
QS = sys.argv[1:] or sorted(os.path.basename(d.rstrip('/')) for d in glob.glob('logs/v2_67/count/*/'))
SOFT = ('radical membership needs 1/N', '1 ∉ J (empty boundary)')
summary = {}
for Q in QS:
    hard = []; prim = collections.Counter(); fail3 = []; direct = collections.Counter(); nodes = collections.Counter()
    for f in sorted(glob.glob(f'logs/v2_67/count/{Q}/*.json')):
        b = os.path.basename(f)[:-5]
        try:
            gc = json.load(open(f'logs/v2_71/{Q}/{b}.gc.json')); u = json.load(open(f'logs/v2_70/{Q}/{b}.u.json'))
            cv = json.load(open(f'logs/v2_70/{Q}/{b}.cv.json'))
        except FileNotFoundError as e: hard.append((b, 'missing', str(e)[-40:])); continue
        lci = set()
        for r in u['results'].values():
            for c in r.get('lci') or []:
                if abs(int(c)) > 1: lci |= set(sp.factorint(abs(int(c))))
            for i in r.get('issues') or []:
                if not i.startswith(SOFT): hard.append((b, 'u', i))
        bad = set()
        for n in gc['nodes']:
            nodes['points' if n['dim'] == 0 else 'dense'] += 1
            if n.get('candidates') is None: hard.append((b, 'gc', n['path'])); continue
            bad |= {int(c) for c in n['candidates']}
            for q, *_ in n.get('direct', []):
                direct[int(q)] += 1
                if int(q) in lci: hard.append((b, 'direct check at a C1 prime', n['path'], q))
        ub = {int(c) for c in u.get('candidates', [])}
        if 3 in ub and os.path.exists(f'logs/v2_71/{Q}/{b}.u3.json'):
            u3 = json.load(open(f'logs/v2_71/{Q}/{b}.u3.json'))
            iss = [i for r in u3['results'].values() for i in (r.get('issues') or []) if not i.startswith('radical membership needs 1/N')]
            if not iss and not u3.get('problems'): ub.discard(3)
            else: fail3.append((b, 'C1-C4 at 3', iss[:2] or u3.get('problems')))
        bad |= ub | {int(c) for c in cv.get('candidates', [])}
        if 3 in bad and not any(x[0] == b for x in fail3): fail3.append((b, 'node certificates at 3'))
        for p in bad: prim[p] += 1
    cs = f'logs/v2_71/charp/{Q}/summary_p3.json'
    okc = {r['block'] for r in json.load(open(cs)) if r['ok']} if os.path.exists(cs) else set()
    resolved3 = [x for x in fail3 if x[0] in okc]; fail3 = [x for x in fail3 if x[0] not in okc]
    if resolved3: prim[3] -= len({x[0] for x in resolved3})
    if prim.get(3, 0) <= 0: prim.pop(3, None)
    rest = sorted(p for p in prim if p >= 3)
    summary[Q] = dict(nodes=dict(nodes), primes=dict(prim), fail3=fail3, resolved3=[x[0] for x in resolved3], hard=hard, direct=dict(direct))
    print(f"{Q}: nodes {dict(nodes)}; failing primes ≥ 3: {rest} (blocks: { {p: prim[p] for p in rest} }); blocks failing at 3: "
          f"{len(fail3)} {fail3[:3]} (resolved by the char-3 route: {len(resolved3)}); direct checks {dict(sorted(direct.items()))}; hard failures {len(hard)} {hard[:3]}")
json.dump(summary, open('logs/v2_71/tally.json', 'w'), indent=1)
