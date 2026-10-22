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
