# v2.76 predictions (committed before the runs)

Subject: the second mechanism on 1⇛2, 2→3→5, 2→4 (K3_tgtfork_long; census 3417 at q = 2, 3423 at q = 3, v2.75), and
Conjecture 3ᶜ′ (finite iff 3ᶜ's conditions and every class of size ≥ 3 has at most one arrow leaving its target and
at most one entering its source).

127. 1⇉2 with three tails (pair_tgtfork3) is flat: the q = 3 census equals 3714 (3ᶜ′: classes of size 2 never
     trigger).
128. The 6 extra classes of K3_tgtfork_long at q = 3 sit in the identity vertex block 10000,01000,00100,00010,00001
     (as for 1⇛2, 2→3, 2→4).
129. On that block the v2.55 pipeline finds a component whose exact generic orbit rank is below its dimension, and
     its generic points carry a one-parameter family whose invariant again pairs the two tails a3 and a4 (a
     cross-ratio of a vertex part and a V part), with the continuing tail's δ supported on paths through a5.
130. 1⇛2, 2→3, 2→4, 3→4 (K3_tgtfork_sq: both tails from j, neither pair to sinks) is not flat at q = 2, 3.

## Addendum (before any run of this round)

Idea: Theorem 5's proof uses the tails x, x′ only through (a) N(u·a₁·x) = u·(λe_i + mℓ₀), (b) x ends at a sink, (c)
S(x) ∝ x. All three make sense for two distinct *paths* π, π′ from j to sinks through single-class arrows.

131. The path version (N(u·a₁·π) = u·(λ_π e_i + m_π ℓ₀), N = 0 elsewhere) is a fold identically on
     1⇛2, 2→3→5, 2→4 (π = a3a5, π′ = a4), on 1⇛2, 2→3, 2→4, 3→4 (π = a4, π′ = a3a5, parallel paths) and on
     1⇛2→3, 3→4, 3→5 (π = a3a4, π′ = a3a5, sharing their first arrow), and φ is invariant on each. If so, the last one
     is infinite although only one arrow leaves the target of the class: Conjecture 3ᶜ′ fails, and the condition
     becomes "two distinct paths from j to sinks" (dually, from sources to i).
132. Hence 128 is refined: the 6 extra classes of K3_tgtfork_long include ≥ 2(q − 1) − ... from this family; the
     q = 2 → 3 jump of K3_after_fork's census is nonzero.
