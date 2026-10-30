"""v2.75: the intrinsic cross-ratio φ of a fold P (sympy d×d matrix over any field of fractions) on a quiver with a
class C of parallel arrows i→j and two arrows x, x' out of j, each alone in its class.

  V = ⟨C⟩; P(V) ⊆ V is required; w spans ker(P|V) (one-dimensional is required);
  δ_y = P(w·y) for y = x, x', split as δ_y = (vertex part) + (V part) + (rest); rest = 0 is required;
  φ = cross-ratio of the lines (E ∩ ⟨vertices⟩, E ∩ V, kδ_x, kδ_x') in E = ⟨δ_x, δ_x'⟩
    = (v_x / v_x') / (u_x / u_x'), with v_y the vertex part and u_y the V part of δ_y, each pair proportional.
Under S ∈ H (vertices fixed, V → V, y ↦ t_y y): w ↦ S w (up to scale), δ_y ↦ t_y⁻¹ S δ_y (up to the same scale), so
v_y ↦ t_y⁻¹ v_y and u_y ↦ t_y⁻¹ S u_y: φ is unchanged. phi() returns (φ, None) or (None, reason)."""
import sympy as sp

def phi(P, basis, cls_arrows, x, x2, simplify=sp.cancel):
    d = len(basis); idx = {b: k for k, b in enumerate(basis)}
    Vi = [idx[('p', (a,))] for a in cls_arrows]
    vert = [k for k, b in enumerate(basis) if b[0] == 'v']
    for c in Vi:
        for r in range(d):
            if r not in Vi and simplify(P[r, c]) != 0: return None, 'P(V) not in V'
    PV = sp.Matrix([[P[r, c] for c in Vi] for r in Vi]).applyfunc(simplify)
    ns = PV.nullspace(simplify=simplify)
    if len(ns) != 1: return None, f'dim ker(P|V) = {len(ns)}'
    w = [simplify(e) for e in ns[0]]
    parts = []
    for y in (x, x2):
        vec = sp.zeros(d, 1)
        for k, a in enumerate(cls_arrows): vec[idx[('p', (a, y))]] += w[k]
        dl = (P * vec).applyfunc(simplify)
        rest = [r for r in range(d) if r not in Vi and r not in vert and dl[r] != 0]
        if rest: return None, f'δ has components outside vertices ⊕ V: {[basis[r] for r in rest]}'
        parts.append(([dl[r] for r in vert], [dl[r] for r in Vi]))
    def ratio(a, b):
        nz = [k for k in range(len(a)) if b[k] != 0]
        if not nz: return None
        rho = simplify(a[nz[0]] / b[nz[0]])
        return rho if all(simplify(a[k] - rho * b[k]) == 0 for k in range(len(a))) else None
    rv = ratio(parts[0][0], parts[1][0]); ru = ratio(parts[0][1], parts[1][1])
    if rv is None or ru is None or ru == 0: return None, f'not proportional / zero (vertex ratio {rv}, V ratio {ru})'
    return simplify(rv / ru), None
