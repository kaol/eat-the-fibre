# v2.45 predictions (committed before any computation)

Setting: 1⇉2⇉3 + 1→3, the two image-dimension-7 groups (L, R). For a fold P₀ with image A:
- V_A = {P₀ + N : N ∈ Hom(ker P₀, A), P₀ + N a fold}. It is the set of folds with image A; idempotence and image
  closure are automatic, so the conditions are kernel closure and Δ.
- Der_A = {D ∈ Der(kQ) : D(A) ⊆ A}. Its directions in V_A are [D, P₀] = −P₀D(1 − P₀).
- The local behaviour dimension is m(P₀) = dim T_{P₀}V_A − dim([Der_A, P₀] ∩ T_{P₀}V_A).

14. At every liftable representative of the growing classes, (b, e) = (1, 0) and (2, 0), of both groups at q = 5,
    m(P₀) = b. All non-derivation tangent directions of Folds are kernel tilts.
15. At the constant-class representatives of the same groups (b ≥ 1 there too, v2.42 be_targets), m(P₀) = b as well.
    So kernel tilts are present throughout the dim-7 groups, and what distinguishes the growing classes is not
    m(P₀). (Recorded as a genuine guess: the alternative is that m = 0 off the growing classes and b there is
    non-reducedness.)

## Outcomes of 14, 15 (`fixed_image_{L,R}.txt`)

- dim W = 35 and dim T_{P₀}V_A = 7 at every representative.
- **L group:** m = b at all representatives except two (1, 0) orbits, where m = 2 (dim Der_A = 11 there).
  - 14 **failed** (2 of 8 growing (1, 0) orbits).
  - 15 **held** in L (constant classes have m = b, up to 5).
- **R group: not interpretable.** The substring selection matched 47 orbits, more than the 27 of the R group at
  q = 5, so the output mixes signature groups. The R run must be redone with exact signature matching.
- **Reading:** m up to 5 at orbits whose count does not grow means the tangent dimensions overcount. V_A is singular
  or non-reduced at many points. First-order counts (b, m) are upper bounds, not dimensions.

## Prediction 16 (committed before the Singular run)

16. For P₀ = rep 5107 (L, class (1, 0)), V_A ⊂ 𝔸³⁵ (N = B_A·Z·R_K, Z a 7×5 matrix of unknowns) has Krull dimension 7
    over ℚ. That is the tangent dimension, so V_A is smooth of dimension 7 at P₀, and the local behaviour dimension is
    m = 1.
