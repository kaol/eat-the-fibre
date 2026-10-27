# v2.66 predictions (committed before the runs)

The whole per-quiver procedure (analysis/v2_66_prove.py: integral structure, vertex blocks over ℤ, certified char 0
pipeline, specialisation to every prime, Bézout bound) on the small unobstructed quivers of the census, those with one
parallel class (K2, K3, pair_after, pair_before, pair_srcout, pair_tgtin, pair_after2, pair_mid, pair_tgtfork,
pair_srcboth, K3_after, K3_srcout) and the three-plus-two class triple_pair_in_short (1⇛2⇇3 + 1→3).

83. On every one of these 13 quivers: the vertex ideal is zero-dimensional, every component passes the v2.54 criterion,
    every specialisation claim holds at every prime (no exceptional prime). So Conjecture 3ᶜ holds for all 13, uniformly
    in q and in every characteristic.
84. The number of vertex blocks equals the number of blocks seen by the census at q = 3 where that is known
    (pairs_in_short-type checks: v2.52, v2.56); for the one-class quivers it is not known in advance.

Added while the larger quivers run (K2 … pair_srcboth done, block sets agree with the q = 3 census reps of v2.60;
the q = 3 censuses of pair_before and pair_mid started but not finished):

85. Duality: path reversal r: kQ → kQ^op is an anti-isomorphism with Δ_op ∘ r = τ ∘ (r⊗r) ∘ Δ, so P ↦ rPr⁻¹ is a
    bijection Folds(Q) → Folds(Q^op) compatible with the census group and the weightings. pair_before = pair_after^op,
    so its q = 3 census has 71 local-torus classes (pair_after's); pair_tgtin = pair_srcout^op (both 101, already known).
86. The vertex blocks of pair_before and pair_mid at q = 3 are exactly the vertex-ideal points (23 and 100).
