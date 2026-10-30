# v2.74 predictions (committed before the runs)

The full chain (v2.66 driver → v2.67 orbit trees → v2.70–72 certificates and char-p trees → v2.73 if a case split
appears) on unobstructed quivers beyond the census list: K₄, 1⇛2→3→4 (K3_after2), 1→2⇛3→4 (K3_mid),
1⇛2, 2→3, 2→4 (K3_tgtfork), 1⇛2, 1→3, 4→2 (K3_srcboth). (1→2⇛3 and 1⇛2, 3→2 are opposites of done quivers.)

115. K₄: the local-torus census is 48 = 14·4 − 8 (v2.18's Kronecker count) at q = 2 and 3, and Conjecture G holds
     with count 48 at every q.
116. For the four K3 quivers, the census at q = 2 equals the census at q = 3 (flat in q), and equals the geometric
     orbit count of the v2.67 trees, block by block.
117. The chain completes on all five without new methods: every node certified at every p (candidates {2, 3} resolved
     by direct checks or char-p trees); Conjecture G holds on all five at every q.

## Addendum (before the checks below; after the census and the identity-block pipeline)

Observed: census of 1⇛2, 2→3, 2→4 is 708 at q = 2 and 710 at q = 3; the difference sits in block
1000,0100,0010,0001 (111 vs 113); the v2.55 pipeline marks two of its 29 components (12, 13, dimension 7) "fail":
orbit rank ≤ 6 at every sampled point.

118. On components 12 and 13 the generic orbit rank (exact, by fraction-free elimination over the function field of the
     component, with 𝔥 = Der_Δ + all three class weightings) is 6 < 7: each carries a one-parameter family of orbits,
     so 1⇛2, 2→3, 2→4 has infinitely many orbits (a counterexample to Conjecture 3ᶜ at multiplicity 3; its
     multiplicity-2 version 1⇉2, 2→3, 2→4 is finite, 414).
119. Its census grows with q: the q = 5 count exceeds 710.
120. (Added while the q = 5 census runs, after certifying the invariants φ₁₂ = y0_7·y6_8/(y0_8·y6_7) and
     φ₁₃ = y3_7·y6_8/(y3_8·y6_7) on components 12 and 13.) The census of 1⇛2, 2→3, 2→4 is 706 + 2(q − 1): 714 at
     q = 5.
