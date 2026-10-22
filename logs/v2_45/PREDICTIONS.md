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

## Outcome of 16: **held** (`VA_5107.txt`)

V_A for rep 5107: 35 unknowns, 995 equations, max degree 2. Singular: dim 7, one minimal associated prime, of
dimension 7, containing Z = 0. The tangent dimension at Z = 0 is 7. So V_A is irreducible of dimension 7 and smooth at
P₀. Within it, the derivation directions [Der_A, P₀] ∩ T have dimension 6, so the behaviours with image A near P₀
form at least 1 dimension.

## Predictions 17, 18 (committed before the runs)

17. For rep 5097 (L, class (3, 0), constant), dim V_{A'} ≤ dim([Der_A', P₀] ∩ T_{P₀}V_{A'}) = 7 − 3 = 4, at least
    near P₀. So there are no behaviour moduli at a constant orbit, and its tangent count m = 3 is an overcount
    (singular or non-reduced point). Concretely: Singular's dimension of the component(s) of V_{A'} through Z = 0
    is ≤ 4.
18. With exact signature matching, the R group has 27 orbits at q = 5. Its growing classes (1, 0) and (2, 0) have
    m = 1 and m = 2 except at most 2 orbits, mirroring L.
