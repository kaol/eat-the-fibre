"""Do (H)-violations account for the positive fibres of looped 2-vertex quivers (v2.15)?  Compare, per quiver,
the vertex folds violating (H) with those having positive exact fibre (truncation L)."""
import sys
from hcheck import violates_H
from classify_loops import iso_classes
from cfibre import fibre
from extend import vertex_folds
L = int(sys.argv[1]); allok = True
for arrows in iso_classes(2):
    A = set(arrows); rows = []
    for P0 in vertex_folds(2):
        if not all(not [a for a in range(2) if P0[a, s]] or not [b for b in range(2) if P0[b, t]] or
                   all((a, b) in A for a in range(2) if P0[a, s] for b in range(2) if P0[b, t]) for s, t in arrows): continue
        R, J, n = fibre(2, arrows, P0, L); dim = J.dimension()
        rows.append((dim > 0, violates_H(2, arrows, P0, L + 1), [list(P0.column(v)) for v in range(2)], dim))
    ok = all(a == b for a, b, _, _ in rows); allok &= ok
    print(arrows, "positive fibre ⇔ (H) violated:", ok, [(c, d) for a, b, c, d in rows if a != b], flush=True)
print("all looped 2-vertex quivers:", allok)
