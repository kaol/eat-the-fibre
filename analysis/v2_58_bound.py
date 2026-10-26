"""v2.58: the explicit census bound for a quiver from its certificates (Bézout counts).

For each vertex-block stratum X_V (n_V presolved coordinates, equations of degree ≤ D_V):
  #components of X_{V,p} ≤ D_V^{n_V}  (Heintz's Bézout inequality: deg V(f_1..f_m) ≤ D^n in A^n);
  each component Y lies in some Z_p, and the connected group H_Z has at most
  o_Z = 1 + Σ_{r < dim Z} D_{Z,r}^{n_V} orbits on Y (an r-dimensional orbit is dense in a component of V(J_{Z,r}),
  which has dimension ≤ r; D_{Z,r} = max degree of the generators of J_{Z,r} = gens_Z + minors).
Stabilisers: H̃ = G × T, G = ∏_classes {column-sum-1 m×m substitutions} (coordinates and an inverse of the
determinant), T ⊆ 𝔾_m^classes (coordinates and inverses); Stab(P) is cut out by Φ φ_w P = P Φ φ_w, of degree
≤ L·(1 + 1) in the coordinates (L = longest path), and the group equations; |π₀ Stab| ≤ D_s^{N}.
Bound: c(q) ≤ Σ_V D_V^{n_V} · max_Z o_Z · D_s^N. Printed as log10.
Usage: python3 analysis/v2_58_bound.py CERT.jsonl[.xz] QUIVER"""
import sys, json, lzma, math, collections, re
import sympy as sp
cert, q = sys.argv[1], sys.argv[2]
recs = [json.loads(l) for l in (lzma.open(cert, 'rt') if cert.endswith('.xz') else open(cert))]
groups = collections.OrderedDict()
for r in recs:
    g = groups.setdefault((r['block'], r['leaf']), {'stratum': None, 'comps': []})
    if r['kind'] == 'stratum': g['stratum'] = r
    else: g['comps'].append(r)
def deg(s_):
    # total degree of a polynomial string: max over terms of the sum of exponents (no parsing into sympy needed)
    best = 0
    for term in re.split(r'(?<=[^\^*eE(])[+-]', s_.replace(' ', '')):
        dgs = 0
        for m in re.finditer(r'([A-Za-z]\w*)(?:\*\*(\d+))?', term): dgs += int(m.group(2) or 1)
        best = max(best, dgs)
    return best
tot = 0.0; rows = []
for (b, leaf), g in groups.items():
    st = g['stratum']; n = len(st['vs'])
    if n == 0: tot += 1; rows.append((b, 0, 0, 0)); continue
    Dv = max([deg(e) for e in st['rem']] + [1])
    oz = 1
    for c in g['comps']:
        if not c['gens']: continue
        dg = max(deg(x) for x in c['gens'])
        s = 1 + sum(max([dg] + [deg(m['det']) for m in cut['minors']]) ** n for cut in c['cuts'])
        oz = max(oz, s)
    term = math.log10(Dv) * n + math.log10(oz)
    rows.append((b, n, Dv, round(term, 1))); tot += 10 ** term
for r_ in rows: print(r_)
print(f"Σ_V (components) × (orbits per component) ≤ 10^{math.log10(tot):.1f}")
