"""v2.68 claim C5: the stratum is covered by its char 0 components at every good prime.
K = intersection over ℚ of the root primes (integral generators); for every generator ρ: ρ^k ∈ (stratum ideal) over ℤ
for some k ≤ 3, else the Rabinowitsch ideal (stratum ideal) + (1 − Tρ) contains an integer N (then ρ ∈ rad(I_p) for
p ∤ N). With the roots and K saturated at p (leading coefficients prime to p), V(I_p) ⊆ V(K_p) ⊆ ∪ V(root_p).
Usage: python3 analysis/v2_68_cover.py CERT BLOCK TREE.json OUT"""
import os, sys, json, re, time, math
HERE = os.path.dirname(os.path.abspath(__file__))
cert, block, tree, out = sys.argv[1:5]; t0 = time.time()
sys.argv = ['v2_67_orbits.py', cert, block, '0', 'all', out + '_o']
src = open(os.path.join(HERE, 'v2_67_orbits.py')).read(); src = src[:src.index("stats = {'primes'")]
G_ = {'__name__': 'v2_68c', '__file__': os.path.join(HERE, 'v2_67_orbits.py')}
exec(compile(src, 'v2_67_orbits', 'exec'), G_)
sp = G_['sp']; R = G_['R']; vs = G_['vs']; ss = G_['ss']; st = G_['st']
T = json.load(open(tree)); roots = [T['nodes'][k]['gens'] for k in T['roots'] if k]
if any(k == '' for k in T['roots']):                        # a root is the zero ideal: it is the whole stratum, coverage is trivial
    res = dict(block=block, roots=len(T['roots']), K=0, lck_primes=[], power_fail=0, rab=[], rab_fail=0, candidates=[], trivial=True)
    json.dump(res, open(out + '.json', 'w')); print(f"COVER {res}; t = {time.time() - t0:.0f}s"); sys.exit(0)
def integral(e):
    den = 1; num = 0
    for c in e.coeffs(): den = den * int(c.denominator) // math.gcd(den, int(c.denominator))
    e2 = e * R(den)
    for c in e2.coeffs(): num = math.gcd(num, int(c.numerator))
    return e2 * R(sp.Rational(1, num)) if num > 1 else e2
ints = lambda gens: [ss(integral(G_['parse'](g))) for g in gens if g.strip() not in ('', '0')] or ['0']
o = G_['singular']([f"ring R0 = 0,({','.join(vs)}),dp;", "option(redSB);"] +
                   [f"ideal C{i} = std(ideal({','.join(ints(g))}));" for i, g in enumerate(roots)] +
                   [f"ideal K0 = std(intersect({','.join(f'C{i}' for i in range(len(roots)))}));" if len(roots) > 1 else "ideal K0 = C0;",
                    'print("K " + string(K0));', 'print("END"); quit;'], 'k')
K = ints(re.search(r"^K (.*)$", o, re.M).group(1).split(','))
rem = [r.replace('**', '^') for r in st['rem']]
o = G_['singular']([f"ring Rz = integer,({','.join(vs)}),dp;", "option(redSB);",
                    "proc lcs(ideal G) { string s = \"\"; int i; for (i = 1; i <= size(G); i++) { if (G[i] != 0) { s = s + \" \" + string(leadcoef(G[i])); } } return(s); }",
                    f"ideal K = std(ideal({','.join(K)}));", 'print("LCK" + lcs(K));',
                    f"ideal GS = std(ideal({','.join(rem)}));", 'print("NGS " + string(size(GS)));',
                    "int i; int k; int ok; list bad;",
                    "for (i = 1; i <= size(K); i++) { ok = 0; for (k = 1; k <= 3 && ok == 0; k++) { if (reduce(K[i]^k, GS) == 0) { ok = k; } } if (ok == 0) { bad = insert(bad, i); } }",
                    'print("PMFAIL " + string(size(bad)) + " OF " + string(size(K)));',
                    f"ring Rt = integer,({','.join(vs)},TT),dp; option(redSB); ideal Gt = imap(Rz, GS); ideal Kt = imap(Rz, K); ideal S; int j;",
                    'for (j = 1; j <= size(bad); j++) { S = std(Gt + ideal(1 - TT * Kt[bad[j]])); if (deg(S[1]) == 0) { print("RAB " + string(S[1])); } else { print("RABFAIL"); } }',
                    'print("END"); quit;'], 'cov')
lck = re.search(r"LCK(.*)", o).group(1).split(); rab = re.findall(r"^RAB (-?\d+)", o, re.M)
pm = re.search(r"PMFAIL (\d+) OF (\d+)", o)
cand = set()
for c in lck + rab:
    if abs(int(c)) > 1: cand |= set(sp.factorint(abs(int(c))))
res = dict(block=block, roots=len(roots), K=len(K), lck_primes=sorted(set().union(*[set(sp.factorint(abs(int(c)))) for c in lck if abs(int(c)) > 1]) if any(abs(int(c)) > 1 for c in lck) else []),
           power_fail=int(pm.group(1)), rab=sorted(set(abs(int(x)) for x in rab)), rab_fail=o.count("RABFAIL"), candidates=sorted(cand))
json.dump(res, open(out + '.json', 'w'))
print(f"COVER {res}; t = {time.time() - t0:.0f}s")
