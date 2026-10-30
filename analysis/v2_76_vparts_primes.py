"""v2.76: the components of the vertex-part ideal (v2_76_vertexpart.py output) with m0·m1 ≠ 0: saturate by m0·m1 and
take minimal associated primes over ℚ (Singular minAssGTZ).
Usage: python3 analysis/v2_76_vparts_primes.py VERTEXPART.json"""
import sys, json, subprocess, re
R = json.load(open(sys.argv[1])); gb = [g.replace('**', '^') for g in R['groebner']]
vs = sorted({v for g in R['groebner'] for v in re.findall(r'[xm]\d+(?:_\d+)?', g)})
L = ['LIB "primdec.lib"; LIB "elim.lib";', f"ring r = 0,({','.join(vs)}),dp;", f"ideal I = {','.join(gb)};",
     "def S_ = sat(I, ideal(m0*m1)); ideal J; if (typeof(S_) == \"list\") { J = S_[1]; } else { J = S_; }", "list L = minAssGTZ(J); int k;",
     'for (k = 1; k <= size(L); k++) { print("PRIME " + string(k) + ": " + string(std(L[k])) + " | dim " + string(dim(std(L[k])))); }', 'print("END"); quit;']
fn = sys.argv[1].replace('.json', '.sing'); open(fn, 'w').write("\n".join(L))
o = subprocess.run(['Singular', '-q', '--no-rc', fn], capture_output=True, text=True, timeout=3000).stdout
print(R['quiver'], 'variables', vs)
for l in o.splitlines():
    if l.startswith('PRIME'): print(' ', l)
