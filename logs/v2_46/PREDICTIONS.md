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
