# v2.64 predictions (committed before any computation on the new quiver)

Program 3: test the leaf-extension mechanism on a second extension of Q₀, Q₂ = Q₀ + 3→4 (pairs_in_short_c4:
1⇉2⇇3, 1→3, 3→4; d = 13; unobstructed, the same two converging classes as Q₀ with joined ends).

77. The local-torus census count of Q₂ is the same at q = 3, 5 and 7.
78. The direct char 0 pipeline (v2.55/v2.57) certifies every component of every vertex block of Q₂.
79. As for Q₁: in char 0 every component of Q₂ is J- or kQ₀-preserving except components in the vertex blocks whose
    3 × 3 corner is not a Q₀ block, and those mixed components have rank ≤ 1.
80. The relative criterion over Q₀'s 366 classes (fibres of the J- and kQ₀-preserving kinds, λ = 0 and 1) passes on
    every positive-dimensional fibre component, and every λ = 0 fibre is finite (as Lemma A predicts for the J kind).

## Outcome

- **77 open**: 1773 at q = 3; q = 5, 7 queued as GPU jobs 18, 19 (the container's q = 5 run was not finished).
- **78 held**: 99/99 blocks certified (one with the integral case split after two out-of-memory presolves).
- **79 failed**: 43 mixed components (37 of rank > 1) in 12 blocks, 6 of them with a Q₀ corner.
- **80 half held**: all 205 positive-dimensional fibre components pass; J-preserving λ = 0 fibres finite; but the
  kQ₀-preserving λ = 0 fibres include 83 positive-dimensional components (the Lemma B idempotent families).
