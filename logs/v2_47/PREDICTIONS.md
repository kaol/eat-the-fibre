# v2.47 predictions (committed before any computation)

26. **The (2, 0) stratum is a divisor.** Take the line ℓ from rep 5103 ((b, e) = (2, 0), image A; the first family in
    `logs/v2_44/driver_L_20.txt`), with direction N_ℓ, and its own admissible weights W_ℓ (exact over ℚ[t]). Let
    𝔥′_ℓ = {D ∈ Der_Δ ⊕ Lie(T_{W_ℓ}) : D(A) ⊆ A}. At t = 2 and t = 3, exactly over ℚ:
    - rank [𝔥′_ℓ, P_t] = 5;
    - adding N_ℓ gives rank 6.

    So H′·ℓ ⊂ V_A is 6-dimensional, a divisor of the 7-dimensional V_A, with a one-parameter family of behaviours of
    its own. This would account for the second q in the L group's 2q growth.
27. **The mirror theorem for the R group.** For rep 6118 (R group, image B) and its first line in
    `logs/v2_44/driver_R_10.txt`:
    - V_B is irreducible of dimension 7 (Singular);
    - at t = 2 the slice has rank 6, and 7 with N, exactly over ℚ;
    - dim T_{P₂}V_B = 7.

    So the dominant-slice theorem of v2.46 holds for B as well.

## Outcomes of 26, 27

- 26 **failed** (`slice_5103.txt`). The line through the (2, 0) representative 5103 is generic in V_A:
  - W = span{(0, 1, 0), (1, 0, 1)} as for 5107, dim 𝔥′ = 6;
  - rank 6 (7 with N) at t = 2, 3, and 5 (6 with N) at t = 1.

  The (2, 0) points are special points on generic lines; the line does not lie in a divisor.
- 27 **held** (`VB_6118.txt`, `slice_6118.txt`, `exact_6118.txt`). For image B (R group):
  - V_B is irreducible of dimension 7 and dim T_{P₂}V_B = 7;
  - W = span{(1, 0, 0), (0, 1, 1)}, the mirror of L under swapping classes a and b;
  - rank 6, and 7 with N, at t = 1, 2, 3;
  - N ∉ [Der, P₂] (22 → 23);
  - the 6 R lines share image B.

  The dominant-slice theorem holds for B.

## Prediction 28 (committed before the Singular run)

28. The rank-drop locus D = {P ∈ V_A : rank [𝔥′, P] ≤ 5}, with 𝔥′ as for the 5107 line, has dimension 6, so it is
    a divisor of V_A. Computed as the incidence variety {(c, Z) : Z ∈ V_A, Σ cᵢ[Dᵢ, P(Z)] = 0} in the charts cₖ = 1,
    which has dimension dim D when the drop is exactly 1. The maximum over the charts is 6.

### Outcome of 28 (partial, `divisor_A.txt`)

- Incidence dimensions by chart c_k = 1, k = 0..5: 6, 6, 5, 6, 5, 5.
- The maximum, 6, matches the prediction only if generic fibres are points, i.e. the stabiliser jumps by exactly 1.
  To make it unconditional, the next run eliminates the c-variables (projection to V_A) in the dimension-6 charts and
  takes dim D directly. The prediction is unchanged: dim D = 6.
