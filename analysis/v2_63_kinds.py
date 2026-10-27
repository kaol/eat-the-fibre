"""v2.63 (Program 3): the kind of every fold of Q₁ = Q₀ + ℓ (2→4) in characteristic 0, component by component.

kQ₁ = kQ₀ ⊕ J. For each irreducible component Z of each vertex-block stratum of Q₁ (certificates of v2.59, with the
integral case split for block 0000,0000,0000,0001), with P = Pr(vs) on Z:
  X ≡ 0 on Z  (the J-rows of the kQ₀-columns vanish: P(kQ₀) ⊆ kQ₀ on all of Z), or
  Y ≡ 0 on Z  (P(J) ⊆ J on all of Z),
decided by reducing the entries modulo a Gröbner basis of Z's prime over ℚ. "Kind" conditions are closed, so if
neither holds identically, the generic fold of Z is of the mixed kind. For components of the mixed kind the 2 × 2
minors of P are reduced too (rank ≤ 1 on all of Z).
Usage: python3 analysis/v2_63_kinds.py OUT_PREFIX CERT [CERT …]"""
import os, sys, json, lzma, subprocess, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE)); sys.path.insert(0, HERE)
import sympy as sp
import fold_ideal as fi
outp, certs = sys.argv[1], sys.argv[2:]
recs = []
for c in certs: recs += [json.loads(l) for l in (lzma.open(c, 'rt') if c.endswith('.xz') else open(c))]
groups = collections.OrderedDict()
for r in recs:
    g = groups.setdefault((r['block'], r['leaf']), {'stratum': None, 'comps': []})
    if r['kind'] == 'stratum': g['stratum'] = r
    else: g['comps'].append(r)
SS = {}; tally = collections.Counter(); rows = []
for (block, leaf), g in groups.items():
    st = g['stratum']; vs = st['vs']
    if block not in SS:
        V = [[-1 if ch == 'm' else int(ch) for ch in row] for row in block.split(',')]
        SS[block] = fi.stratum('pairs_in_short_tail', V)
    S = SS[block]; R = S['R']; b1 = S['ctx']['basis']; d = len(b1)
    gl = {str(x): x for x in R.gens}
    parse = lambda s_: R(sp.sympify(s_, locals={k: sp.Symbol(k) for k in gl}))
    Pr = [row[:] for row in S['P']]
    for x, v in st['subs']:
        xv, vv = gl[x], parse(v); Pr = [[e.compose(xv, vv) for e in row] for row in Pr]
    Jx = [i for i, b in enumerate(b1) if (b[0] == 'v' and b[1] == 3) or (b[0] == 'p' and 5 in b[1])]
    Qx = [i for i in range(d) if i not in Jx]
    ss = lambda e: str(e).replace('**', '^')
    X = [ss(Pr[r][c]) for r in Jx for c in Qx if Pr[r][c]]
    Y = [ss(Pr[r][c]) for r in Qx for c in Jx if Pr[r][c]]
    for c in g['comps']:
        if not vs:
            kind = ('X0' if not X else '') + ('Y0' if not Y else '')
            tally[kind or 'mixed (point)'] += 1; rows.append((block, leaf, c['comp'], kind)); continue
        gens = [x for x in c['gens'] if x.strip() not in ('0', '')] or ['0']
        minors = []
        src = [f"ring r = 0,({','.join(vs)}),dp;", f"ideal G = std(ideal({','.join(x.replace('**', '^') for x in gens)}));",
               f"ideal X = {','.join(X) or '0'}; ideal Y = {','.join(Y) or '0'};",
               "print(\"X0 \" + string(size(reduce(X, G)) == 0)); print(\"Y0 \" + string(size(reduce(Y, G)) == 0));",
               f"matrix PP[{d}][{d}] = " + ",".join(ss(Pr[i][j]) if Pr[i][j] else "0" for i in range(d) for j in range(d)) + ";",
               "ideal M2 = minor(PP, 2); print(\"RANK1 \" + string(size(reduce(M2, G)) == 0));", "quit;"]
        fn = f"{outp}.sing"; open(fn, 'w').write("\n".join(src))
        res = subprocess.run(["Singular", "-q", "--no-rc", fn], capture_output=True, text=True, timeout=3000).stdout
        assert "RANK1" in res and "? " not in res, res[-300:]
        x0 = "X0 1" in res; y0 = "Y0 1" in res; r1 = "RANK1 1" in res
        kind = 'split' if x0 and y0 else 'kQ-preserving' if x0 else 'J-preserving' if y0 else ('mixed, rank ≤ 1' if r1 else 'mixed, rank > 1')
        tally[kind] += 1; rows.append((block, leaf, c['comp'], kind))
    print(f"{block} leaf {leaf}: " + ", ".join(k for _, _, _, k in rows[-len(g['comps']):]), flush=True)
print("TALLY", dict(tally))
mixed = sorted({b for b, _, _, k in rows if k.startswith('mixed')})
print("blocks with mixed-kind components:", mixed)
