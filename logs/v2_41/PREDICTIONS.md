# v2.41 predictions (committed before writing the Lean proofs)

Quiver 1⇉2⇇3 (sources e₀, e₂; sink e₁). S is the source swap: e₀ ↔ e₂, a_i ↔ b_i. Already computed (sympy, before
this file): S·F1_t·S = F3_t, S·F2_t·S = F2_{1−t}, S·F4_t·S = F4_{1−t}.

1. Over an integral domain k, every multiplicative Φ on k(1⇉2⇇3) with a two-sided inverse satisfies
   Filt Φ ∨ Filt (S·Φ).
   - Proof plan: rad is preserved (b² = 0 and k reduced).
   - The sink kills rad on the left and the sources kill it on the right. This gives Φ₀₁ = Φ₂₁ = Φ₁₀ = Φ₁₂ = 0 and
     Φ₁₁ = 1.
   - The rest of the vertex block is a 2×2 matrix of idempotents with row sums 1 and orthogonal rows. Its columns
     are nonzero because rad² = 0 kills idempotents in rad. So it is the identity or the swap.
2. S is an algebra automorphism (Mult S, S·S = 1), and the three identities above hold over every commutative ring.
3. Over an integral domain, for t, t' ≠ 0: an automorphism intertwining F_i(t) and F_j(t') forces i = j or
   {i, j} = {1, 3}. F1_t and F3_t are intertwined by S. So, at nonzero parameters, the four families of Theorem 3 are
   at most three behaviours, with F1 ≅ F3. Exactly three would also need that F1_t and F1_t' are always intertwined,
   which this round does not prove.
