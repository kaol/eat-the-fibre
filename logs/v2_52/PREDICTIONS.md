# v2.52 predictions (committed before any computation)

Context: ChatGPT's Program 2 (the joined-ends configuration). Theorem 3's modulus on the unjoined 1⇉2⇇3 is the vertex
block of P: the coefficients of vertices in P(e_v) (v2.27). v2.29 showed that its families fail on 1⇉2⇇3 + 1→3 in
the Δ condition of the arrow 1→3. Question: does the arrow kill the vertex-block modulus outright?

45. Over all folds (all fibres in the census checkpoints), with entries lifted symmetrically to integers:
    a. on 1⇉2⇇3 + 1→3 (pairs_in_short) the set of vertex blocks is one finite set of integer matrices, the same at
       q = 3, 5, 7;
    b. on 1⇉2⇇3 (pairs_in) the number of distinct vertex blocks grows with q (q = 3, 5, 7, 11), by Theorem 3.
46. On pairs_in_short, every vertex block has entries in {−1, 0, 1}.

## Added after 45–46 (which held), before the elimination below

47. A necessary system over ℚ: P_K (P restricted to K = ker Δ, 25 unknowns) is idempotent with image and kernel closed
    under products, and P(c) (10 unknowns) satisfies Δ(P c) = P(e₁) ⊗ P(e₃). Projected to the vertex block (9 unknowns),
    it has finitely many points, and they are exactly the 21 census blocks.
48. The same system without the arrow-c condition (which is all of P_K's conditions on the unjoined 1⇉2⇇3, whose K is
    the same 5-dimensional algebra) projects onto a positive-dimensional set of vertex blocks.

## Outcomes
45a, 45b, 46 held (entries are in {0, 1}). 47 failed as stated: P_K plus the c condition is 1-dimensional; with the
a, b conditions added, the elimination gives exactly the 21 census blocks. 48 held: dimension 2 without arrow
conditions, 1 with a, b.
