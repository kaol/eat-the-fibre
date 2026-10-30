# v2.75 predictions (committed before the runs)

Notation. 1⇛2, 2→3, 2→4 (K3_tgtfork): class arrows a0, a1, a2 (1→2), tails b = a3 (2→3), b' = a4 (2→4); V =
span(a0, a1, a2); σ: V → k the sum of coordinates (the class group G = {g ∈ GL(V): σg = σ}, times weightings).

Already observed before this file (exploratory, recorded so it is not mistaken for a prediction):
- Generic points of components 12 and 13 of block 1000,0100,0010,0001 have Im P = ⟨vertices, U, b, b', U·b, U·b'⟩
  (U ⊂ V a plane), ker(P|V) a line w ⊂ ker σ, and P(v·b) = P(v)·b + f(v)(λ_b c + m_b ℓ₀) with f the w-coordinate,
  ℓ₀ spanning U ∩ ker σ, and c = e₀ (component 12) or c = e₁ + e₂ + e₃ (component 13).
- The normal form P = P_{λ,m,λ',m'}: identity except P(a0) = a1, P(a0·b) = a1·b + λe₀ + m(a1 − a2),
  P(a0·b') = a1·b' + λ'e₀ + m'(a1 − a2), satisfies all fold equations identically in (λ, m, λ', m') (one run).

Predictions.

121. The variant P⁻ (vertex part λ(e₁ + e₂ + e₃) in place of λe₀) also satisfies the fold equations identically.
122. φ(P) := cross-ratio of the lines (E ∩ k⟨vertices⟩, E ∩ V, kδ_b, kδ_b') in E = ⟨δ_b, δ_b'⟩, δ_x = P(w·x) with
     w spanning ker(P|V), is H-invariant: φ(SPS⁻¹) = φ(P) as rational functions of a generic S ∈ H (symbolic g ∈ G,
     weightings s, t_b, t_b'). On the family φ = λm'/(λ'm). It agrees with v2.74's φ₁₂ (resp. φ₁₃) at the rational
     points of component 12 (resp. 13).
123. Hence, over every field k, the folds P_{1,1,1,t} (t ∈ k^×) are pairwise non-conjugate, as are the P⁻; the census
     of 1⇛2, 2→3, 2→4 has at least 2(q − 1) classes in block 1000,0100,0010,0001 for every q.
124. The obstruction at multiplicity 2 is ℓ₀ ∈ ker σ: on K3_tgtfork with ℓ₀ = a1 (∉ ker σ) and on 1⇉2, 2→3, 2→4
     with ℓ₀ = a1, the fold equations force m = m' = 0.
125. The family extends to every quiver with a class of ≥ 3 arrows i→j and two arrows b, b' out of j that are each
     alone in their class: P = the projection onto ⟨paths⟩ ∖ {a0, a0·b, a0·b'} along ⟨a0 − a1, k_b, k_b'⟩ (k_x =
     a0·x − a1·x − λ_x e_i − m_x(a1 − a2)) satisfies the fold equations identically on K4_tgtfork, K3_tgtfork3,
     K3_tgtfork_pre, K3_tgtfork_long, K3_tgtfork_in, K3_tgtfork_short. With 122's argument, all six have infinitely
     many orbits over every infinite field.
126. Census at q = 2 and q = 3 (torch engine):
     - flat (equal) for K3_tgtpair (1⇛2⇉3), pair_tgtfork3 (1⇉2 with three tails), K3_tgtinout (1⇛2→3, 4→2) and
       K3_srcout_after (1⇛2→4, 1→3): the family needs two non-parallel arrows out of j (or, dually, into i);
     - not flat for K3_tgtfork3 (three tails; the family has two independent cross-ratios).
