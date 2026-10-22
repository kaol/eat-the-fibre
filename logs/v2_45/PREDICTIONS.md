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
