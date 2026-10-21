# v2.33 pre-registration (committed before the behaviour census is run)

Behaviour of a fold P: the Aut(kQ)-orbit of the splitting (im P, ker P), i.e. P up to conjugation by *all* algebra
automorphisms of kQ (inner automorphisms and arrow substitutions with any column sums), not only the Δ-compatible
group G of the census. (P is determined by the subspaces im P, ker P, so behaviour must forget the embedding to be
coarser than P; Aut(kQ) is the coarsest ambient choice that keeps the algebra.)

Already observed (family computations over 𝔽₅, before the census):
- Theorem 2, Theorem 3 F2 and F4: every member of the family is inner-conjugate to every other.
- Theorem 4 (t ≠ 1) and Theorem 3 F1, F3 (t ∉ {0, −1}): single Aut(kQ)-orbits, but only via Inn ∘ Sub compositions.

Predictions for the behaviour census (1⇉2⇇3 and 1⇇2⇉3 over 𝔽₃, 𝔽₅, 𝔽₇; 1⇉2, 1→3 over 𝔽₃, 𝔽₅):
1. (Le Chat, sub-question 1) The number of behaviours is the same for every q tested on each quiver.
2. On 1⇉2⇇3 the q-dependence of the local-torus count 449 + 4(q−2) disappears entirely: the four Theorem 3 families
   are four behaviours with fibres of width q − 2 (counted in G-orbits: wider).
3. Module class (L, R, K, E) is constant on behaviours (it is Aut(kQ)-invariant by construction; a check, not a
   prediction).
4. Opposite duality: pairs_in and pairs_out have the same number of behaviours.
