# v2.70 predictions (committed before the runs)

Goal: make Theorem G (v2.69, census = geometric orbits for p ≫ 0) effective: an explicit N_Q per quiver with
census = Σ geometric orbits for every prime p ∤ N_Q. Per dense node Y of every tree, certificates over ℤ, each giving
a finite set of candidate primes (computed at two independent rational points x of the dense orbit; a prime is a
candidate for Y only if it is one at both points):

  G1  orbit closure: a strong ℤ-Gröbner basis of the graph ideal of h ↦ h·x (H_Y = G × T_{W(Y)}), elimination order;
      the elimination ideal equals I_Y over ℚ; p ∤ (leading coefficients, denominators) ⇒ Y_p = closure(H_p x_p),
      geometrically irreducible.
  S   stabiliser: (i) the equations are affine-linear in g (already observed at all 6696 nodes), so Stab ∩ G is an
      open subset of an affine space in every characteristic; (ii) the torus part (elimination of g) is the lattice
      ideal of a saturated character lattice; (iii) strong ℤ-basis leading coefficients and the Smith invariants of
      the Jacobian at 1 ⇒ Stab_p smooth, connected, of dimension dim H − dim Y.
  W   the compatibility defect matrix at x: rank mod p = rank over ℚ (Smith invariants); W(Y) is a partition span
      (observed at all 6664 nodes), so the census's linear weighting moves are T_{W(Y)}(𝔽_p).
  C   the v2.68 claims C1–C5 (Lemma C at p) for every tree.

102. G1, S, W hold at every node over ℚ (elimination = I_Y; torus part saturated), and the candidate primes, after
     intersecting over the two points, lie in {2, 3, 5}.
103. The Lemma C claims C1–C5 hold on every tree with candidates in {2, 3}.
104. Hence N_Q has prime factors in {2, 3, 5} for every one of the 17 quivers: census = geometric orbits for all
     p ≥ 7 (expected: p ≥ 5).
