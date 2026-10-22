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

## Outcomes of 1–3 (before the targeted search)

1. **Failed.** 35 curves meet ≥ 3 orbits at q = 5. Reducing all 35 mod 3 separates them:
   - 18 are genuine families: each meets q orbits at both primes, pairwise disjoint.
   - 17 are bounded curves: they meet 3 orbits at both primes.
   - The 18 explain 36 of the growth of 44, so 4 families are missing.
2. **Held for the 18 found** (they meet 3 orbits mod 3, pairwise disjoint). The constant count does not follow, because
   of the 4 missing families: the orbits they miss go 834 → 842.
3. **Failed as written.** Two of the 18 need x = δ_a − 2δ_aδ_b (reps 350 and 351, classes R and L). Whether that is
   equivalent to a pure direction is not checked.

The missing growth sits in exactly two signature groups, an L/R pair with image dim 7 and kernel dim 5. Each goes
23 → 27, i.e. two families each (`growth_3_5.txt`).

## Prediction 4 (committed before the targeted search)

4. Inner curves Inn(1 + t·x)·P₀ with x ∈ rad, *including the arrow c: 1→3 and the paths a_ib_j*, found by a tangent
   computation at liftable representatives of those two groups, give exactly the 4 missing families: 2 L and 2 R,
   each meeting q orbits at q = 3 and 5.
