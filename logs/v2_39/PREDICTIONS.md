# v2.39 predictions (committed before any computation)

Goal: a Lean theorem that the 19 families of 1⇉2⇉3 are 19 distinct behaviours (Aut(kQ)-orbits).

Filt(Φ) := (rows of vertices are unit rows: Φ v j = δ_vj) ∧ (Φ l c = 0 whenever len l < len c), with len the path
length of a basis element. This is ε∘Φ = ε together with Φ(rad^n) ⊆ rad^n. The plan: prove that the generators of
Aut(kQ) = Inn ⋊ Sub (v2.33; Aut(Q) acts within classes for this quiver) satisfy Filt, that Filt is closed under
products, and that no Filt pair Φ, Ψ = Φ⁻¹ intertwines two different families.

1. For every family and every vertex v with depth d_v ≥ 0, row v of P_t has an entry equal to the constant ±1 in some
   column of path length ≥ d_v. Then each ordered pair is refuted by one entry, using Φ when d_v(i) > d_v(j) and Ψ
   otherwise.
2. Inn(1 + x) for x ∈ rad, with x given by 8 free coordinates, has an explicit polynomial matrix that is Filt.
   subst ga gb is Filt for any ga, gb. Both are checked in Lean by `fin_cases`/`simp`.
3. The packaged theorem (Φ, Ψ Filt, mutually inverse, Φ·P^{(i)}_t = P^{(j)}_{t'}·Φ ⇒ i = j) builds, with every
   result on the standard axioms.

## Outcomes

1. **Held**: every nonzero vertex row of every family has a constant ±1 entry at its full depth.
2. **Held, and stronger than planned**:
   - Inn(1 + x) for x ∈ rad with all 8 (on 1⇉2⇇3: 4) coordinates free is proved Filt (`innM_filt`).
   - It is proved to be y ↦ u y u⁻¹ (`innM_conj`, `uU_uI`).
   - subst ga gb is Filt for any ga, gb (`subst_filt`).
3. **Held**: `fam_behaviours_distinct` builds, on the standard axioms.
   - Added (not predicted): the same for Theorem 3's F1–F4 on 1⇉2⇇3 with t, t' ≠ 0 (`thm3_behaviours_distinct`).
   - The ±1 witness of prediction 1 does not exist there: F1 and F3 have rows −t·(…), so the witness entry is ±t and
     t ≠ 0 is needed.
