# Discordianism: Open Questions Explored (v2.35): the nineteen families of 1⇉2⇉3, by behaviour type

*v2.28 asked for a classification of the moduli of 1⇉2⇉3: the local-torus count 678 + 19(q−1) has only two families
with closed forms (Theorems 2 and 4). v2.33 showed that families are wide behaviour fibres. So any invariant of algebra
automorphisms is constant on a family, and grouping orbits by such invariants locates the families without the
infeasible q⁸ inner-orbit search. This round does that, and checks the method on the converging quiver and on the
joined-ends quiver. Predictions for the second part were committed first (`logs/v2_35/PREDICTIONS.md`). v2.34's
statement of the missing lemma (L) is corrected (§4).*

---

## 0. Summary

| # | Statement | Status |
|---|---|---|
| 1 | On 1⇉2⇉3 the growth 19(q−1) sits in 9 signature groups, with slopes 9 and 3 (no class), 1 + 1 + 1 (L), 1 + 1 + 1 (R) and 1 (LRKE, Theorem 2). Every family has a semisimple image, k×k or k³ | computed, 𝔽₅ → 𝔽₇ (exploratory) |
| 2 | Every growing group's family members have their excess direction inside V = ⟨ad δ_aδ_b, ad δ_a, ad δ_b, E_a, E_b⟩: inner derivations by the canonical lines δ_a = a₁ − a₂, δ_b = b₁ − b₂, δ_aδ_b, plus the two class scalings | computed (exploratory) |
| 3 | 1⇉2⇇3: total slope 4, split 3 (no class) + 1 (L), at 𝔽₅ → 𝔽₇ and 𝔽₇ → 𝔽₁₁. Only 2 signature groups, because F1–F3 share one signature | pre-registered; holds in total, fails in group count |
| 4 | 1⇉2⇇3 + 1→3 (joined ends): no growing signature group (106 groups, constant over 𝔽₃ → 𝔽₅). This is Conjecture 3ᶜ at the level of behaviour types | pre-registered; holds |
| 5 | On 1⇉2⇇3 the growth needs ad δ_a and ad δ_b together, i.e. the δ_b − δ_a direction of v2.33 | pre-registered; holds |

## 1. Method

A *signature* of a fold P combines:
- its module class (L, R, K, E; v2.32);
- its vanishing pattern (which of the 18 typed words vanish; v2.32);
- isomorphism invariants of the algebras im P and ker P: dimension, powers, annihilators, centre.

The invariants that count points (idempotents, square-zero elements) are left out, because they are polynomials in q.
Every ingredient is preserved by conjugation with any algebra automorphism of kQ, so the signature is constant on
behaviours (v2.33), and hence on families. The script asserts that it is constant on each local-torus orbit.

`analysis/family_signatures.py` groups the local-torus orbits of two censuses by signature. A group whose size grows
with q contains families; its slope counts them, up to collisions of signatures.

## 2. The nineteen families of 1⇉2⇉3

| class | im P | ker P (dim, powers, ann_L, ann_R, centre) | local orbits, q = 5 → 7 | slope | excess carried by |
|---|---|---|---|---|---|
| none | k×k | (9, 9, 0, 0, 1) | 61 → 79 | 9 | ad δ_aδ_b with the torus (3); ad δ_a, ad δ_b with the torus, excess 2 (6) |
| none | k³ | (8, 8, 0, 0, 0) | 35 → 41 | 3 | one each: ad δ_a, ad δ_aδ_b, ad δ_b, each with the torus |
| L | k×k | (9, 9, 0, 0, 0) | 17 → 19 | 1 | ad δ_aδ_b, ad δ_a |
| L | k×k | (9, 9, 0, 5, 0) | 11 → 13 | 1 | ad δ_aδ_b, ad δ_b |
| L | k³ | (8, 8, 0, 4, 0) | 12 → 14 | 1 | ad δ_b |
| R | k×k | (9, 9, 0, 0, 0) | 17 → 19 | 1 | ad δ_aδ_b, ad δ_b |
| R | k×k | (9, 9, 5, 0, 0) | 11 → 13 | 1 | ad δ_aδ_b, ad δ_a |
| R | k³ | (8, 8, 4, 0, 0) | 12 → 14 | 1 | ad δ_a |
| LRKE | k×k | (9, 9, 6, 6, 4) | 11 → 13 | 1 | ad δ_aδ_b (Theorem 2) |

The slopes total 19, the coefficient of q in 678 + 19(q−1). (`logs/v2_35/signatures_two_pairs_5_7.txt`;
excess by motif from `analysis/v2_35/motif.py`.)

Reading the table:
- **Theorem 4** lies in the first row: its image is k×k, and its kernel invariants are (9, 9, 0, 0, 1), as computed in
  v2.33. **Theorem 2** is the last row.
- **Opposite duality** pairs L with R, swapping the annihilators and δ_a with δ_b. 1⇉2⇉3 is self-opposite, and the
  table is symmetric under this swap.
- **No family needs a derivation outside V.** Every family modulus is a combination of the class scalings with inner
  derivations by three canonical elements:
  - δ_a and δ_b, which exist for any class of size 2 and are killed by Δ;
  - δ_aδ_b, which exists only when the two classes are consecutive.

  The census group G rescales the lines δ_a and δ_b separately: a column-sum-1 substitution M sends δ_a to
  (M₁₁ − M₁₂)·δ_a. So only the span of each line is meaningful, never a fixed combination. A probe along the single
  direction ad(δ_b − δ_a) (`analysis/v2_35/probe_dir.py`) misses the growth for exactly this reason.
- **The point-level marker is not specific.** 296 of 754 local orbits over 𝔽₅ have an excess direction that extends to
  second order. Degenerate points that lead into larger orbits qualify too. But in every growing group the marked count
  grows exactly as fast as the group. (`logs/v2_35/excess_marked_two_pairs.txt`)

## 3. Out-of-sample checks

- **1⇉2⇇3.** Total slope 4, split 3 + 1, at both 𝔽₅ → 𝔽₇ and 𝔽₇ → 𝔽₁₁, as Theorem 3 requires. The prediction of four
  groups failed: F1, F2 and F3 have identical signatures. The invariants used here are too coarse to separate these
  behaviours.
- **1⇉2⇇3 + 1→3.** No growing group. All 106 signature groups have the same size over 𝔽₃ and 𝔽₅, so the joining arrow
  removes growth at the level of behaviour types, not only of orbits. This is weaker than Conjecture 3ᶜ in one
  direction (two primes) and stronger in another (it holds behaviour type by behaviour type).
- **Growth on 1⇉2⇇3 needs both lines.** The growing points need ad δ_a and ad δ_b together, the span containing
  v2.33's δ_b − δ_a.

## 4. Correction to v2.34

v2.34 §5 first stated the missing lemma (L) in pointwise form: T_e X_P = Lie H + Lie Stab(P) for every fold P. That is
e(P) = 0, which v2.34's own first-order table refutes at special points. The corrected statement, now in v2.34, is
generic: on each irreducible component of a twist variety X_P, the reduced dimension equals that of the double coset
through a generic point.

## 5. Next

1. **Closed forms for the seventeen unnamed families.** Each lies in a known signature group, with a known excess
   motif. Following v2.33, a family is an orbit of a curve φ_s = (torus) ∘ Inn(1 + s·x) with x in the motif's span.
   Lifting one member to ℚ and solving for φ_s, as in `analysis/family_behaviour.py`, should give the normal forms.
2. **A finer invariant** that separates F1, F2 and F3. The vertex block does, but it is not Aut(kQ)-invariant (v2.33).
   Candidates: the completion variety (v2.32), and invariants of the pair (im P, ker P) beyond the separate algebras.
3. **The mechanism of the joining arrow.** With c: u → w, the derivations a_i ↦ a_i + λ·c·b_j enlarge Der. Test
   whether at every point with δ-excess on 1⇉2⇇3 + 1→3, the excess is absorbed by these new directions or obstructed.

Files: `analysis/family_signatures.py`, `analysis/v2_35/`, `logs/v2_35/`.
