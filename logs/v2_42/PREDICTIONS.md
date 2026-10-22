# v2.42 predictions (committed before any computation)

Quiver 1⇉2⇉3 + 1→3 (`two_pairs_short`). Orbits are taken under the full Δ-compatible group found in v2.28: G, the
local torus, and the shear U = {c ↦ c + s·δ_aδ_b}. Counts are 939 (𝔽₃) and 983 (𝔽₅) (logs/unipotent_merge_two_pairs_short.log).

1. The inner-curve search of v2.36 (`analysis/family_curves.py`), run against the merged orbit labels at q = 5, finds
   exactly 22 families: pairwise disjoint curves P_t = Inn(1 + t·x)·P₀·Inn⁻¹, each meeting ≥ 4 merged orbits.
   This is the slope (983 − 939)/(5 − 3) = 22.
2. Their integer closed forms, reduced mod 3, meet exactly 3 merged orbits each, pairwise disjoint. That would make
   939 − 3·22 = 983 − 5·22 = 873 orbits constant.
3. All families lie along x ∈ {δ_a, δ_b, δ_aδ_b}, as on 1⇉2⇉3. No direction outside the span of δ_a, δ_b and their
   product is needed.

Caveat recorded in advance: censuses over prime fields cannot see a family defined only over 𝔽_{p²} and paired with
its Frobenius conjugate. Such a pair contributes finitely many 𝔽_p-points. So "no other families" claims are about
prime-field counts only.
