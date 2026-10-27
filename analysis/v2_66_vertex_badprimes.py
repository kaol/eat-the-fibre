"""v2.66: the vertex step at the primes dividing a non-unit leading coefficient.

v2_58_vertex_integral (mode int) gives a strong Gröbner basis over ℤ of the vertex ideal and of the ideal of its
points. Where all leading coefficients are 1, the points over every 𝔽̄_p are the reductions of the ℤ-points (v2.58).
For the primes p dividing some leading coefficient, decompose the vertex ideal over 𝔽_p directly (mode p=P,
minAssGTZ) and check: every associated prime is a point (linear), and the point set is the reduction mod p of the
char 0 point set, bijectively. Writes OUT_DIR/vertex_badprimes.json.
Usage: python3 analysis/v2_66_vertex_badprimes.py QUIVER OUT_DIR"""
import sys, os, re, json, subprocess, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
Q, out = sys.argv[1], sys.argv[2]
summ = json.load(open(os.path.join(out, 'summary.json'))); arrows = summ['vertex']['arrows']
txt = open(os.path.join(out, 'vertex_integral.txt')).read()
lcs = [int(x) for x in re.findall(r"^LC (-?\d+)", txt, re.M)]
rlc = re.search(r"^RLC (.*)", txt, re.M); lcs += [int(x) for x in rlc.group(1).split()] if rlc else []
primes = sorted({p for c in lcs if abs(c) != 1 for p in sp.factorint(abs(c))})
def pts(f, p):
    t = open(f).read(); res = []
    for bl in t.split("POINTS")[1].split("--")[:-1]:
        body = bl.split("\n", 1)[1]
        vals = {m.group(1): -int(m.group(2)) if m.group(2) else 0 for m in re.finditer(r"(k\d\d)([-+]\d+)?\s*[,\n]", body + "\n")}
        n_gens = len([x for x in body.replace('\n', '').split(',') if x.strip()])
        res.append((vals, n_gens))
    return res
nv = int(len(pts(os.path.join(out, 'vertex_char0.txt'), 0)[0][0]) ** 0.5)
p0 = [v for v, _ in pts(os.path.join(out, 'vertex_char0.txt'), 0)]
rep = {'non_unit_LC': sorted({c for c in lcs if abs(c) != 1}), 'primes': {}}
for p in primes:
    f = os.path.join(out, f'vertex_char{p}.txt')
    if not os.path.exists(f):
        with open(f, 'w') as fh:
            subprocess.run([sys.executable, 'analysis/v2_58_vertex_integral.py', Q, os.path.join(out, 'scratch', f'vb{p}'), f'p={p}', '--arrow', arrows],
                           stdout=fh, stderr=subprocess.STDOUT, cwd=ROOT, check=True)
    Pp = pts(f, p)
    linear = all(len(v) == nv * nv and n == nv * nv for v, n in Pp)
    red = [tuple(sorted((k, x % p) for k, x in v.items())) for v in p0]
    got = sorted(tuple(sorted((k, x % p) for k, x in v.items())) for v, _ in Pp)
    rep['primes'][p] = dict(zero_dim='PDIM 0' in open(f).read(), points=len(Pp), all_linear=linear,
                            reduction_injective=len(set(red)) == len(red), equal=sorted(red) == got)
json.dump(rep, open(os.path.join(out, 'vertex_badprimes.json'), 'w'))
print(Q, rep)
