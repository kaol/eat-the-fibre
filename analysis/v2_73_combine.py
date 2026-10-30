"""v2.73: combine the ℤ-certificates of a case-split block (A–E) into its candidate primes and its count.
Nodes are identified by full-coordinate key across the root trees; fixed point leaves are nodes too. Per point of a
dense node, the candidates are its certificate primes plus the primes separating it from every other node of the same
dimension (full coordinates); a prime survives if no point clears it (directly or by a direct check at a point whose
cross-separation does not contain it).
Usage: python3 analysis/v2_73_combine.py OUTDIR"""
import sys, os, json, glob, math, re
import sympy as sp
outd = sys.argv[1]
primes_of = lambda n: set(sp.factorint(abs(int(n)))) if abs(int(n)) > 1 else set()
roots = json.load(open(f'{outd}/roots.json'))['roots']
hard = []; block_primes = {}
def addb(tag, ps):
    for q in ps: block_primes.setdefault(int(q), set()).add(tag)
# leaves: A, B, C, P
fixed_points = []
for f in sorted(glob.glob(f'{outd}/leaf*.json')):
    L = json.load(open(f)); k = L['leaf']
    for q, tags in L['cand'].items(): addb('leaf' + str(k) + ':' + '/'.join(sorted(set(tags))), [q])
    ch = L['checks']
    if ch.get('B') not in (None, 'ok'): hard.append((k, 'B', ch['B']))
    if L['vs'] == 0 and L['ncomp']:
        if not (ch.get('P_G_fixes') and ch.get('P_W_fixes')):
            if not str(ch.get('C', '')).startswith('ok'): hard.append((k, 'point leaf neither fixed nor in a root'))
        else: fixed_points.append((k, L['point_full']))
    elif L['ncomp'] and not str(ch.get('C', '')).startswith('ok'): hard.append((k, 'C', ch.get('C')))
# trees
nodes = {}                                                  # fullkey → node record
for tf in sorted(glob.glob(f'{outd}/tree_L*_*.json')):
    if tf.endswith(('.gc.json', '.u.json', '.full.json')): continue
    base = tf[:-5]; T = json.load(open(tf))
    u = json.load(open(base + '.u.json')); gc = json.load(open(base + '.gc.json')); full = json.load(open(base + '.full.json'))
    for q in u.get('candidates', []): addb('C1-C4 ' + os.path.basename(base), [q])
    if u.get('problems'): hard.append((base, 'u problems', u['problems'][:2]))
    for r in u['results'].values():
        for i in r.get('issues') or []:
            if not i.startswith(('radical membership needs 1/N', '1 ∉ J (empty boundary)')): hard.append((base, 'u', i))
    bypath = {nd['path']: nd for nd in T['nodes'].values()}
    for n in gc['nodes']:
        fk = bypath[n['path']].get('fullkey')
        if fk is None: hard.append((base, 'no full key', n['path'])); continue
        rec = dict(src=os.path.basename(base), path=n['path'], dim=n['dim'], gc=n, full=full.get(n['path'], []))
        if fk in nodes and nodes[fk]['gc'].get('candidates') is not None and n.get('candidates') is None: continue
        nodes[fk] = rec
for k, pf in fixed_points:                                  # a fixed point leaf: key = its coordinates
    fk = ','.join(f'{n_}-({v})' for n_, v in sorted(pf.items()))
    nodes.setdefault('POINT:' + fk, dict(src=f'leaf{k}', path='point', dim=0, gc=dict(candidates=[]), full=[pf], point=True))
# point nodes inside trees: coordinates from their (linear) full key
for fk, rec in nodes.items():
    if rec['dim'] == 0 and not rec.get('point'):
        sol = {}
        for g in fk.split(','):
            e = sp.sympify(g.replace('^', '**')); v = list(e.free_symbols)
            if len(v) == 1 and sp.degree(e, v[0]) == 1: sol[str(v[0])] = str(sp.solve(e, v[0])[0])
        rec['full'] = [sol]
# the same point may be both a fixed point leaf and a point node of a tree: merge by coordinates
pts = [(fk, rec) for fk, rec in nodes.items() if rec['dim'] == 0]
for i, (fk1, r1) in enumerate(pts):
    for fk2, r2 in pts[i + 1:]:
        if fk1 not in nodes or fk2 not in nodes: continue
        a1, a2 = r1['full'][0], r2['full'][0]
        common = set(a1) & set(a2)
        if common and all(sp.Rational(a1[c]) == sp.Rational(a2[c]) for c in common) and len(common) >= min(len(a1), len(a2)):
            drop = fk2 if r2.get('point') else fk1         # keep the tree's node
            nodes.pop(drop)
# cross-separation
keyed = list(nodes.items()); evalcache = {}
def gens_of(fk):
    if fk.startswith('POINT:'): return [sp.sympify(t.replace('^', '**')) for t in fk[6:].split(',')]
    return [sp.sympify(g.replace('^', '**')) for g in fk.split(',') if g.strip() not in ('', '0')]
def sep(pt, fk2):
    g_ = 0; sub = {sp.Symbol(k): sp.Rational(v) for k, v in pt.items()}
    for g in gens_of(fk2):
        val = sp.Rational(g.subs(sub)); g_ = math.gcd(g_, abs(int(sp.fraction(val)[0])))
        if g_ == 1: break
    return primes_of(g_) if g_ else {'INSEP'}
summary = []
for fk, rec in keyed:
    others = [fk2 for fk2, r2 in keyed if fk2 != fk and r2['dim'] == rec['dim']]
    gcn = rec['gc']
    if gcn.get('candidates') is None: hard.append((rec['src'], 'node uncertified', rec['path'])); continue
    if rec['dim'] == 0:
        cross = set()
        for fk2 in others: cross |= sep(rec['full'][0], fk2)
        cand = {int(c) for c in gcn['candidates'] if str(c).isdigit()} | {c for c in cross if c != 'INSEP'}
        if 'INSEP' in cross: hard.append((rec['src'], 'point not separated', rec['path']))
    else:
        pts = gcn['points']; per = []
        for i, p in enumerate(pts):
            if not all(v for _, v in p['checks']) or i >= len(rec['full']) or rec['full'][i] is None: per.append(None); continue
            cross = set()
            for fk2 in others: cross |= sep(rec['full'][i], fk2)
            per.append(({int(q) for q in p['cand']}, {c for c in cross if c != 'INSEP'}, 'INSEP' in cross))
        good = [x for x in per if x is not None and not x[2]]
        if not good: hard.append((rec['src'], 'no point separated', rec['path'])); continue
        cand = set.intersection(*[a | b for a, b, _ in good])
        for q, pi, srcs in gcn.get('direct', []):          # a direct check at point pi clears q if q does not separate there
            if pi < len(per) and per[pi] is not None and q not in per[pi][1]: cand.discard(int(q))
    for q in cand: addb('node ' + rec['src'] + ':' + rec['path'], [q])
    summary.append((rec['src'], rec['path'], rec['dim'], sorted(cand)))
count = len(nodes)
out = dict(count=count, dense=sum(1 for r in nodes.values() if r['dim'] > 0), points=sum(1 for r in nodes.values() if r['dim'] == 0),
           primes={str(q): sorted(t) for q, t in sorted(block_primes.items())}, hard=hard, nodes=summary)
json.dump(out, open(f'{outd}/combined.json', 'w'), indent=1)
print(f"{outd}: count {count} ({out['dense']} dense, {out['points']} points); candidate primes {sorted(block_primes)}; hard failures {len(hard)} {hard[:4]}")
for q, t in sorted(block_primes.items()): print(f"  {q}: {sorted(t)[:6]}{' …' if len(t) > 6 else ''}")
