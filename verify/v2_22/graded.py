"""Orbits of G × 𝔾_m(grading) on the folds: φ_λ = λ^{|path|} satisfies Δφ_λ = λ(φ_λ⊗φ_λ)Δ, so P ↦ φ_λ P φ_λ⁻¹ preserves folds.
Prints counts under G and under G × F_q^*; with --sizes also the orbit-size distribution under G × F_q^*."""
import sys, collections, numpy as np
exec(open('gcensus.py').read().split("if __name__ == '__main__':")[0])
name, p = sys.argv[1], int(sys.argv[2])
res, folds = census(name, p, verbose=False)
nv, arrows = QUIVERS[name]; basis, M, DL = build(nv, arrows); d = len(basis)
lens = [0 if b[0] == 'v' else len(b[1]) for b in basis]
F = Field(p)
# rebuild G (same as census) via a small re-run of its group code
src = open('gcensus.py').read()
gcode = src[src.index("    classes = collections.defaultdict(list)"):src.index("    seen = set(); sizes = collections.Counter()")]
ns = dict(np=np, itertools=itertools, collections=collections, arrows=arrows, p=p, F=F, d=d, basis=basis)
exec("def _g():\n" + gcode + "    return G\n", ns); G = ns['_g']()
Gr = [(np.diag([pow(l, k, p) for k in lens]), np.diag([pow(l, (p - 1 - k % (p - 1)) % (p - 1), p) for k in lens])) for l in range(1, p)]
key = lambda A: (A % p).tobytes()
seen = set(); sizes = collections.Counter()
for P in folds:
    if key(P) in seen: continue
    o = set()
    for g, gi in G:
        Q = g @ P @ gi % p
        for e, ei in Gr: o.add(key(e @ Q @ ei))
    seen |= o; sizes[len(o)] += 1
print(f"{name} F_{p}: {res['N']} folds; G-orbits {res['orbits']}; (G × grading)-orbits {sum(sizes.values())}"
      + (f"; sizes {dict(sorted(sizes.items()))}" if '--sizes' in sys.argv else ""), flush=True)
