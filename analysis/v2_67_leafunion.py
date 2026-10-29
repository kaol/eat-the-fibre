"""v2.67: geometric orbit count on a case-split vertex block. Every component of the stratum is a leaf prime (v2.57),
and the orbit closures inside it are the nodes of its orbit stratification; leaf primes that are not H-stable are
skipped. Nodes from different leaves are identified by their reduced Gröbner basis in all stratum coordinates (the
leaf's prime plus the graph of its substitutions). Count = Σ over distinct nodes of their geometric components.
Usage: python3 analysis/v2_67_leafunion.py CERT BLOCK OUTDIR [J]"""
import sys, os, json, lzma, re, subprocess, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
cert, block, out = sys.argv[1:4]; J = int(sys.argv[4]) if len(sys.argv) > 4 else 2
os.makedirs(out, exist_ok=True)
op = lzma.open if cert.endswith('.xz') else open
leaves = {}
for l in op(cert, 'rt'):
    if f'"block": "{block}"' not in l: continue
    r = json.loads(l)
    if r['kind'] == 'stratum' and r['block'] == block: leaves[r['leaf']] = r
points = {k: r for k, r in leaves.items() if not r['vs'] and r['ncomp'] > 0}
todo = [k for k, r in sorted(leaves.items()) if r['vs'] and r['ncomp'] > 0]
run = {}
while todo or run:
    while todo and len(run) < J:
        k = todo.pop(0); f = os.path.join(out, f'L{k}')
        if os.path.exists(f + '.json'): continue
        run[k] = subprocess.Popen([sys.executable, os.path.join(HERE, 'v2_67_orbits.py'), cert, block, str(k), 'all', f, '--leafmode'],
                                  stdout=open(f + '.txt', 'w'), stderr=subprocess.STDOUT)
    time.sleep(1)
    for k in list(run):
        if run[k].poll() is not None: del run[k]
import fold_ideal as fi
Q = next(iter(leaves.values()))['quiver']
S = fi.stratum(Q, [[-1 if ch == 'm' else int(ch) for ch in row] for row in block.split(',')])
import numpy as np, sympy as sp
ctx = S['ctx']; d = ctx['d']; DLn = np.asarray(ctx['DL']); basis = ctx['basis']; arrows = ctx['arrows']; idx = {b: i for i, b in enumerate(basis)}
def fixed_point(r):
    """is the point leaf's fold fixed by H (all of Der_Δ, and every class weighting compatible at it)?"""
    gl = {str(g): g for g in S['R'].gens}; loc = {k_: sp.Symbol(k_) for k_ in gl}
    Pr = [row[:] for row in S['P']]
    for x, v in r['subs']:
        xv = gl[x]; vv = S['R'](sp.sympify(v, locals=loc)); Pr = [[e.compose(xv, vv) for e in row] for row in Pr]
    P = sp.Matrix(d, d, lambda a, b: sp.Rational(str(Pr[a][b])) if Pr[a][b] else 0)
    lead = next(rr for rr in leaves.values() if rr.get('DerD'))
    for D in lead['DerD']:
        Dm = sp.Matrix(d, d, [sp.Rational(x) for x in D])
        if not (Dm * P - P * Dm).is_zero_matrix: return False
    cl = lead['classes']; Ec = [sp.diag(*e) for e in lead['Ec']]
    def DLclass(c):
        D_ = np.zeros((d, d, d), dtype=np.int64)
        for l, b_ in enumerate(basis):
            if b_[0] == 'v': continue
            for kk, x in enumerate(b_[1]):
                if x not in c: continue
                pre = idx[('p', b_[1][:kk])] if kk > 0 else idx[('v', arrows[x][0])]
                suf = idx[('p', b_[1][kk + 1:])] if kk + 1 < len(b_[1]) else idx[('v', arrows[x][1])]
                D_[l, pre, suf] += 1
        return D_
    rows = []
    for c in cl:
        Dc = DLclass(c); col = []
        for l in range(d):
            lhs = sum((P[m, l] * sp.Matrix(Dc[m].tolist()) for m in range(d) if P[m, l]), sp.zeros(d, d))
            col += list(lhs - P * sp.Matrix(Dc[l].tolist()) * P.T)
        rows.append(col)
    Wn = sp.Matrix(rows).T.nullspace()                      # weightings compatible at P
    return all((sum((w[c] * Ec[c] for c in range(len(cl))), sp.zeros(d, d)) * P - P * sum((w[c] * Ec[c] for c in range(len(cl))), sp.zeros(d, d))).is_zero_matrix for w in Wn)
keys = {}; bad = []
for k, r in sorted(leaves.items()):
    if not r['ncomp']: continue
    if k in points and not fixed_point(r):                 # a single fold that H moves lies in an orbit counted elsewhere
        continue
    if k in points:                                         # a single fold: its key from the substitutions
        rel = [f"{x}-({v.replace('**', '^')})" for x, v in r['subs']]
        src = "\n".join([f"ring RF = 0,({','.join(S['names'])}),dp;", "option(redSB);", f"ideal G = std(ideal({','.join(rel)}));",
                         'print("FULL " + string(G));', 'quit;'])
        fn = os.path.join(out, f'P{k}.sing'); open(fn, 'w').write(src)
        o = subprocess.run(["Singular", "-q", "--no-rc", fn], capture_output=True, text=True).stdout
        keys.setdefault(re.search(r"FULL (.*)", o).group(1), (0, 1, f'L{k}')); continue
    f = os.path.join(out, f'L{k}.json')
    if not os.path.exists(f): bad.append((k, 'no output')); continue
    D = json.load(open(f))
    if D['stats']['nodense'] or D['stats']['unresolved']: bad.append((k, D['stats']))
    for nd in D['nodes'].values():
        if nd.get('result') in ('dense orbit', 'points'):
            if nd.get('geometric') is None: bad.append((k, 'geometric count unknown')); continue
            keys.setdefault(nd['fullkey'], (nd['dim'], nd['geometric'], f"L{k}:{nd['path']}"))
total = sum(g for _, g, _ in keys.values())
print(f"{block}: {len(leaves)} leaves, {len(keys)} distinct orbit closures, geometric orbits {total}; problems {bad}")
json.dump(dict(block=block, orbits=total, closures=len(keys), problems=bad), open(os.path.join(out, 'union.json'), 'w'))
