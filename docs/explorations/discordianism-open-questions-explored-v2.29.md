# Discordianism: Open Questions Explored (v2.29): Lean 4

*The fold results of v2.26–v2.28 rested on symbolic checks in sympy and on hand proofs. This round moves the checkable
part into Lean 4 with Mathlib (`lean/`, see `lean/README.md`). Everything below is kernel-checked and depends only on
Lean's standard axioms (`lean/EatTheFibre/Audit.lean`). Labels: **proved (Lean)** marks a statement checked by Lean.*

---

## 0. Summary

| # | Statement | Status |
|---|---|---|
| 1 | **Basis criterion.** For a finite-dimensional algebra given by structure constants, a linear `P` is a fold (idempotent, image and kernel closed under products, `ΔP = (P⊗P)Δ`) iff the conditions hold on basis vectors | **proved (Lean)**, `isFold_iff_basis`. This is the method of every symbolic and census check in the project |
| 2 | Conjugation by a Δ-compatible algebra automorphism maps folds to folds | **proved (Lean)**, `IsFold.conj` |
| 3 | Vertex-block lemma (v2.27): an automorphism fixing the vertices and preserving the radical leaves the vertex-to-vertex entries unchanged | **proved (Lean)**, `vertex_block`. In v2.27 it was proved by hand |
| 4 | Theorem 2 / Theorem 1 normal form on 1⇉2⇉3 is a fold for every t, over every commutative ring | **proved (Lean)** |
| 5 | Theorem 3's four families on 1⇉2⇇3 are folds for every parameter | **proved (Lean)** |
| 6 | Theorem 4's normal form on 1⇉2⇉3 and on 1⇉2⇉3 + 1→3 is a fold for every t | **proved (Lean)** |
| 7 | Controls: Theorem 2's family on 1⇉2⇉3 + 1→3, and Theorem 3's families on 1⇉2⇇3 + 1→3, fail as polynomial identities, in the Δ condition of the arrow 1→3 | **proved (Lean)**. This matches the hand proofs: the arrow between the outer ends breaks Δ-compatibility |
| 8 | Rigidity: for every invertible class substitution (any `ga, gb ∈ GL₂`), conjugation never maps a Theorem 3 member to a member with another parameter | **proved (Lean)**, `thm3_rigid_GT` |
| 9 | Rigidity: if a class substitution with column sums 1, scaled by class weights, intertwines Theorem 4's `P_t` and `P_t′`, then t = t′ (and the scalings on the image are 1) | **proved (Lean)**, `thm4_rigid_GT` |

## 1. How the checks are made trustworthy

The families have polynomial entries. Lean checks them by reflection:
- `Poly` is a computable integer polynomial type.
- `check c D P` computes the four basis conditions on `Poly`. The kernel evaluates it (`decide +kernel`), one row per
  declaration so that each check fits in about 5 GB.
- `check_sound` proves that `check = true` implies `IsFold` over every commutative ring at every parameter value.

No `native_decide`, so the compiler is not trusted. The first approach, entrywise `simp`/`ring`, took over 10 minutes on
the 7-dimensional case and ran out of memory on the 11-dimensional one. Reflection takes about 3 minutes per
11-dimensional algebra.

The concrete algebras are generated (`lean/gen/gen.py`) from `pcensus.build`, the census engine's own basis, structure
constants and Δ. Lean certifies the results relative to those tables, and the generator is the only link to the census.

## 2. What is not formalized

- The reduction in Theorems 3 and 4 from an arbitrary acyclic quiver to the subalgebra W (the two-class path algebra).
  It needs general path algebras of quivers, which Mathlib does not have in this form. The Lean results cover the base
  quivers, and the extension relies on the hand proof.
- That the census group G × T consists of exactly the matrices used in the rigidity theorems. Those theorems are stated
  for the explicit substitution matrices, which contain G × T.
- Census counts and Conjecture 3ᶜ.

## 3. Open questions

1. A general path-algebra layer (an acyclic quiver gives `CoordAlg` with paths as basis), proving the W-reduction once.
2. Generate the remaining quivers of v2.27/v2.28 (1⇉2⇇3 with tails, the triple classes, 1→2⇉3⇉4). They are larger, so
   the kernel checks need more splitting.
3. Formalize Theorem 4's invariant for the general member (x, y, z, w′) rather than the normal form.
