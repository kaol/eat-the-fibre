# v2.46 predictions (committed before any computation)

Brief: docs/briefs/chatgpt-v2.45-frontier.md (ChatGPT, 2026-10-22). Setting: rep 5107 (L group, q = 5), its exact
line P_t = P₀ + t·N (v2.44, `line_isFold`), A = im P₀ = im P_t, V_A irreducible of dimension 7 (v2.45).

Census symmetries that map folds to folds:
- Aut_Δ = G·U (Lie algebra Der_Δ, the derivations with ΔD = (D⊗1 + 1⊗D)Δ), by transport;
- class weightings φ_w at folds compatible with each class-split Δ_c (v2.23 Lemma).

Let 𝔥 = {D ∈ Der_Δ ⊕ span(E_a, E_b, E_c) : D(A) ⊆ A}, the Lie algebra of the A-stabiliser H.

21. P_t is compatible with each class-split Δ_c for all t (exactly over ℚ[t]). So H·P_t ⊆ Folds, and since H fixes A,
    H·P_t ⊆ V_A.
22. dim span{[D, P_t] : D ∈ 𝔥} = 6 at t = 1, exactly over ℚ, and N is independent of that span. So
    Φ(h, t) = h·P_t·h⁻¹ has differential rank 7 at (1, 1).
23. Consequence: Φ : H × 𝔸¹ → V_A is dominant (V_A irreducible, dim 7). With v2.44's finite intersection of the line
    with each Aut(kQ)-orbit, the behaviours represented on a dense subset of V_A are those of the line: a
    1-dimensional family. (This is the dominant-slice lemma of the brief, Task 3, with Task 2's ranks certified at
    one rational point.)
24. (Task 1) Among the 32 lines of v2.44 (`driver_L_20.txt`, `driver_R_10.txt`, `lines_5107.txt`), the images
    (exactly over ℚ) fall into exactly 2 distinct subspaces within the L lines and 2 within the R lines.

## Outcomes of 21, 22 (`slice_5107.txt`)

- 21 **failed as stated**. Over ℚ[t], P_t is compatible with Δ_b (class 2→3) but not with Δ_a (1→2) or Δ_c (1→3)
  separately. Compatibility is linear in the weights, so the admissible set W(P_t) is a subspace. It contains
  (1, 1, 1) and (0, 1, 0), hence (1, 0, 1), but not E_a or E_c alone. The 𝔥 used (dim 6) wrongly included E_a and E_c
  separately.
- 22 **failed as stated**: at t = 1 the rank is 5 (6 with N), so t = 1 is a special point. At t = 2 the rank is 6
  (7 with N).
- dim Der = 22, dim Der_Δ = 5, as before.

## Prediction 22′ (committed before the corrected run)

22′. Let W = {w : P_t is compatible with Δ_w for all t} (a subspace, computed exactly) and
     𝔥′ = {D ∈ Der_Δ ⊕ {Σ w_c E_c : w ∈ W} : D(A) ⊆ A}. Then at t = 2 (exactly over ℚ),
     rank [𝔥′, P_2] = 6 and adding N gives rank 7.

## Outcome of 22′: **held** (`slice_5107_W.txt`)

Exactly over ℚ:
- W = span{(0, 1, 0), (1, 0, 1)};
- dim 𝔥′ = 6;
- at t = 2 and t = 3, rank [𝔥′, P_t] = 6, and 7 with N. (t = 1 is special: 5 and 6.)

## Prediction 25 (committed before the run)

25. At P₂ = P₀ + 2N, exactly over ℚ:
    - dim T_{P₂}V_A = 7 (Jacobian of the 995 V_A equations at Z = 2·Z_N);
    - rank [Der(kQ), P₂] = 22, and N is not in its span.

    So m(P₂) = 1 exactly, and P₂ is a smooth point of V_A.

## Outcomes of 25, 24 (`exact_checks.txt`)

- 25 **held**, exactly over ℚ:
  - dim T_{P₂}V_A = 7 (995 equations, 35 unknowns);
  - dim Der = 22, rank [Der, P₂] = 22, and 23 with N.

  So P₂ is a smooth point and m(P₂) = 1.
- 24 **failed**. The 14 base representatives have exactly 2 distinct images: all 8 L-group lines share image A, and
  all 6 R-group lines share image B. So it is one fixed image per group, not two.

Hypothesis for the 2q growth per group with one image (to be pre-registered as a prediction next round): the (2, 0)
points of V_A form a divisor D on which H′ has 5-dimensional orbits. So D/H′ is one-dimensional as well, giving the
q orbits of class (2, 0) next to the q orbits of class (1, 0) from the generic stratum.
