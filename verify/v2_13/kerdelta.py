"""dim ker Δ = n + Σ_{parallel classes} (m − 1) in degree 1 (+ higher-degree all-parallel products); = n exactly when no parallel arrows."""
from sage.all__sagemath_singular import matrix, QQ
from quiver import build
from collections import Counter
from gdim import Qs
ok = True
for name, (nv, arrows) in Qs.items():
    d, mult, Dl, allp = build(nv, arrows)
    Dm = matrix(QQ, d * d, d, lambda ab, l: int(Dl[l, ab // d, ab % d]))
    k = d - Dm.rank(); par = any(m > 1 for m in Counter(arrows).values())
    ok &= (k == nv) == (not par)
    print(f"{name}: dim ker Δ = {k}, vertices = {nv}, parallel arrows: {par}")
print("ker Δ = A_0 exactly when there are no parallel arrows:", ok)
