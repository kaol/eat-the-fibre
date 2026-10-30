"""v2.68: specialisation of an orbit tree (v2.67) to every prime.

For every node Y of the tree (a prime over ℚ, with the local group 𝔥_Y = Der_Δ + W(Y)) this computes, over ℤ
(Singular, ring integer; or over 𝔽_p with --p P for a direct check at one prime):
  C1  a strong Gröbner basis of I_Y (integral primitive generators of the ℚ-prime) and its leading coefficients;
  C2  J_Y: the closure of I_Y + (f̃) under the integral vector fields ξ̃ (each field scaled to integral coefficients;
      f̃ an e × e minor of the integral orbit matrix, e = dim Y), computed with ℤ-reductions, so every element of
      J_Y is a ℤ-combination of I_Y, f̃ and ξ̃-images of earlier elements;
  C3  K_Y = intersection over ℤ of the children's C1 ideals; for each generator ρ, the least k ≤ 8 with
      reduce(ρ^k, J_Y) = 0 over ℤ (then ρ^k ∈ J_Y ⊗ 𝔽_p for every p);
      a 0-dimensional boundary is the empty set: then 1 ∈ J_Y is needed (an integer N ∈ J_Y: primes of N);
  C4  for each local weighting w of Y (scaled integral): Σ_c w_c (Δ_c-defect) reduces to 0 modulo C1;
and at the roots
  C5  K_root = intersection over ℤ of the root components; ρ^k ∈ (stratum ideal) over ℤ.
Candidate bad primes: the primes of the C1 leading coefficients (flatness, hence H-stability and dimension of the
fibres), and, where an exact ℤ-membership of C3–C5 fails, the primes of the leading coefficients of the basis it was
tested against. Writes OUT.json and prints a summary.
Usage: python3 analysis/v2_68_uniform.py CERT BLOCK TREE.json OUT [--p P] [--only KEYINDEX,…]"""
import os, sys, json, re, time, argparse, math
HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser(); ap.add_argument('cert'); ap.add_argument('block'); ap.add_argument('tree'); ap.add_argument('out')
ap.add_argument('--p', type=int, default=0); ap.add_argument('--only', default=''); ap.add_argument('--tmo', type=int, default=20000)
U = ap.parse_args(); t0 = time.time()
# reuse the setup of v2_67_orbits.py (stratum, Pr, fields, orbit matrix, defects, minors)
sys.argv = ['v2_67_orbits.py', U.cert, U.block, os.environ.get('STAB_LEAF', '0'), 'all', U.out + '_o'] + (['--leafmode'] if os.environ.get('STAB_LEAF') else [])   # v2.73: a case-split leaf's chart
src = open(os.path.join(HERE, 'v2_67_orbits.py')).read(); src = src[:src.index("stats = {'primes'")]
G_ = {'__name__': 'v2_68', '__file__': os.path.join(HERE, 'v2_67_orbits.py')}
exec(compile(src, 'v2_67_orbits', 'exec'), G_)
sp = G_['sp']; R = G_['R']; vs = G_['vs']; nv = len(vs); ss = G_['ss']
T = json.load(open(U.tree)); NODES = T['nodes']; KEYS = list(NODES)
CH = 'integer' if not U.p else str(U.p)
DSHIFT = 1 if not U.p else 0                           # Krull dimension over ℤ is the fibre dimension + 1
def lcm(a, b): return a * b // math.gcd(a, b)
def integral(e):
    """a ring element scaled to integral primitive coefficients"""
    if not e: return e
    den = 1; num = 0
    for c in e.coeffs(): den = lcm(den, int(c.denominator))
    e2 = e * R(den)
    for c in e2.coeffs(): num = math.gcd(num, int(c.numerator))
    return e2 * R(sp.Rational(1, num)) if num > 1 else e2
def row_scale(row):
    den = 1
    for e in row:
        if e:
            for c in e.coeffs(): den = lcm(den, int(c.denominator))
    return den
def ints(gens): return [ss(integral(G_['parse'](g))) for g in gens if g.strip() not in ('', '0')] or ['0']
RINGZ = [f'LIB "primdec.lib";', f"ring Rz = {CH},({','.join(vs)}),dp;", "option(redSB);",
         "proc lcs(ideal G) { string s = \"\"; int i; for (i = 1; i <= size(G); i++) { if (G[i] != 0) { s = s + \" \" + string(leadcoef(G[i])); } } return(s); }"]
def run(lines, tag):
    return G_['singular'](lines, tag)
def primes_of(lcs):
    out = set()
    for c in lcs:
        c = abs(int(c))
        if c > 1: out |= set(sp.factorint(c))
    return out
results = {}; cand = set(); problems = []
sel = [int(x) for x in U.only.split(',') if x.strip()]
for n, key in enumerate(KEYS):
    if sel and n not in sel: continue
    nd = NODES[key]
    if nd.get('result') not in ('dense orbit', 'points'): problems.append((n, 'node result ' + str(nd.get('result')))); continue
    gens = nd['gens']; I = ints(gens); rec = dict(dim=nd['dim'])
    if nd['result'] == 'points':
        o = run(RINGZ + [f"ideal GI = std(ideal({','.join(I)}));", 'print("LCI" + lcs(GI));', 'print("DIM " + string(dim(GI)));', 'print("END"); quit;'], f'u{n}')
        lci = re.search(r"LCI(.*)", o).group(1).split(); rec.update(lci=lci, dimp=int(re.search(r"DIM (-?\d+)", o).group(1)) - DSHIFT)
        if not U.p: cand |= primes_of(lci)
        results[n] = rec; continue
    Wb = [[sp.Rational(x) for x in w] for w in nd.get('weights') or []]
    act = G_['active'](Wb); XI = G_['combine'](act, G_['XIB']); Mrows = G_['combine'](act, G_['MB'])
    sc = [lcm(row_scale(XI[a]), row_scale(Mrows[a])) for a in range(len(act))]
    XI = [[e * R(sc[a]) if e else e for e in XI[a]] for a in range(len(act))]
    Mrows = [[e * R(sc[a]) if e else e for e in Mrows[a]] for a in range(len(act))]
    dim_, iv, _ = G_['info'](gens, f'ui{n}')
    f, how = G_['choose_minor'](gens, iv, dim_, Mrows)
    if f is None: problems.append((n, 'no minor')); continue
    f = integral(f)
    kids = [NODES[c]['gens'] for c in nd.get('children', [])]
    # local weightings, integral
    defl = []
    for w in Wb:
        den = 1
        for x in w: den = lcm(den, int(sp.fraction(x)[1]))
        wi = [int(x * den) for x in w]
        defl.append([ss(sum((R(wi[c]) * G_['defects'][c][k] for c in range(G_['NCL']) if wi[c]), R(0))) for k in G_['nzd']])
    NH = len(XI)
    L = RINGZ + [f"matrix XI[{NH}][{nv}] = " + ",".join(ss(XI[a][b]) for a in range(NH) for b in range(nv)) + ";",
                 f"proc xi(poly h, int a) {{ poly g = 0; int b; for (b = 1; b <= {nv}; b++) {{ g = g + diff(h, var(b)) * XI[a, b]; }} return(g); }}",
                 f"ideal GI = std(ideal({','.join(I)}));", 'print("LCI" + lcs(GI));', 'print("DIM " + string(dim(GI)));',
                 f"poly f = {ss(f)};", 'if (reduce(f, GI) == 0) { print("FININ"); }',
                 "ideal G = std(GI + f); ideal NEW = f; int rounds = 0; int a; int i; int k; poly g;",
                 "while (size(NEW) > 0) { rounds++; ideal NX = 0; for (i = 1; i <= size(NEW); i++) { for (a = 1; a <= " + str(NH) +
                 "; a++) { g = reduce(xi(NEW[i], a), G); if (g != 0) { NX = NX + g; G = std(G + g); } } } NEW = simplify(NX, 2); kill NX; }",
                 'print("CLOSED " + string(rounds) + " DIMJ " + string(dim(G)));', 'print("LCJ" + lcs(G));']
    if kids:
        # K = intersection of the children over ℚ (fast), then an integral generating set. For p prime to the leading
        # coefficients of the strong ℤ-bases of K and of every child, K ⊗ ℤ_(p) and the children are saturated, so
        # ∏ C_i,p ⊆ K_p and V(K_p) ⊆ ∪ V(C_i,p).
        if len(kids) > 1:
            k0 = run([f'ring R0 = {"0" if not U.p else U.p},({",".join(vs)}),dp;', "option(redSB);"] +
                     [f"ideal C{i_} = std(ideal({','.join(ints(kg))}));" for i_, kg in enumerate(kids)] +
                     [f"ideal K0 = std(intersect({','.join(f'C{i_}' for i_ in range(len(kids)))}));", 'print("K " + string(K0));', 'print("END"); quit;'], f'k{n}')
            Kg = ints(re.search(r"^K (.*)$", k0, re.M).group(1).split(','))
        else:
            Kg = ints(kids[0])
        L += [f"ideal K = std(ideal({','.join(Kg)}));", 'print("LCK" + lcs(K));']
        for i_, kg in enumerate(kids):
            L += [f"ideal C{i_} = std(ideal({','.join(ints(kg))}));", f'print("LCC" + lcs(C{i_}));']
        L += ["int ok; int fails = 0; list bad;",
              "for (i = 1; i <= size(K); i++) { ok = 0; for (k = 1; k <= 3 && ok == 0; k++) { if (reduce(K[i]^k, G) == 0) { ok = k; } } "
              "if (ok == 0) { bad = insert(bad, i); } } print(\"PMFAIL \" + string(size(bad)) + \" OF \" + string(size(K)));",
              # Rabinowitsch for the rest: S = G + (1 - TT ρ); its integer generator N gives ρ ∈ rad(J ⊗ 𝔽_p) for p ∤ N
              f"ring Rt = {CH},({','.join(vs)},TT),dp; option(redSB); ideal Gt = imap(Rz, G); ideal Kt = imap(Rz, K); ideal S; int j;",
              'for (j = 1; j <= size(bad); j++) { S = std(Gt + ideal(1 - TT * Kt[bad[j]])); '
              'if (deg(S[1]) == 0) { print("RAB " + string(S[1])); } else { print("RABFAIL"); } }', "setring Rz;"]
    else:
        L += ['if (reduce(1, G) == 0) { print("EMPTY 1"); } else { print("EMPTY 0"); print("CONST " + string(G[1])); }']
    for j, dl in enumerate(defl):
        L += [f"ideal DW = {','.join(dl) or '0'};", f'if (size(reduce(DW, GI)) == 0) {{ print("DEF {j} ok"); }} else {{ print("DEF {j} fail"); }}']
    L += ['print("END"); quit;']
    o = run(L, f'u{n}')
    if o is None: problems.append((n, 'singular timeout/error')); print(f"node {n}: SINGULAR timeout/error", flush=True); continue
    lci = re.search(r"LCI(.*)", o).group(1).split(); lcj = re.search(r"LCJ(.*)", o).group(1).split()
    rec.update(lci=lci, dimp=int(re.search(r"DIM (-?\d+)", o).group(1)) - DSHIFT, dimJ=int(re.search(r"DIMJ (-?\d+)", o).group(1)) - DSHIFT,
               fin=("FININ" in o), defects=re.findall(r"DEF \d+ (\w+)", o))
    rec['lck'] = sum((l.split() for l in re.findall(r"LCK(.*)", o)), []); rec['lcc'] = sum((l.split() for l in re.findall(r"LCC(.*)", o)), [])
    m = re.search(r"PMFAIL (\d+) OF (\d+)", o)
    if m: rec['pm_fail'] = int(m.group(1)); rec['pm_total'] = int(m.group(2))
    rec['rab'] = re.findall(r"^RAB (-?\d+)", o, re.M); rec['rab_fail'] = o.count("RABFAIL")
    if "EMPTY" in o: rec['empty_ok'] = "EMPTY 1" in o
    if not U.p:
        cand |= primes_of(lci) | primes_of(rec['lck']) | primes_of(rec['lcc'])
        cand |= primes_of(rec['rab'])
        if rec.get('empty_ok') is False: cand |= primes_of(lcj)
        if 'fail' in rec['defects']: cand |= primes_of(lci)
    bad = []
    if rec['dimp'] != nd['dim']: bad.append('dim')
    if rec['fin']: bad.append('f in I')
    if rec.get('rab_fail'): bad.append(f"radical membership fails for {rec['rab_fail']} generators")
    if rec.get('rab') and any(abs(int(x)) != 1 for x in rec['rab']): bad.append(f"radical membership needs 1/N, N in {sorted(set(abs(int(x)) for x in rec['rab']))}")
    if rec.get('empty_ok') is False: bad.append('1 ∉ J (empty boundary): candidate primes from J')
    if 'fail' in rec['defects']: bad.append('defect')
    rec['issues'] = bad; results[n] = rec
    print(f"node {n}: dim {nd['dim']}, LC(I) primes {sorted(primes_of(lci))}, closure dim {rec['dimJ']}, "
          f"children {len(kids)}, {'issues ' + str(bad) if bad else 'all claims exact'}; t = {time.time() - t0:.0f}s", flush=True)
json.dump(dict(block=U.block, p=U.p, results=results, candidates=sorted(cand), problems=problems), open(U.out + '.json', 'w'))
print(f"RESULT nodes {len(results)}, problems {problems}, candidate primes {sorted(cand)}; t = {time.time() - t0:.0f}s")
