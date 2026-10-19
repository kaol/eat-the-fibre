# Discordianism: Open Questions Explored (v2.30): Theorems 3 and 4 in Lean on every acyclic quiver

*v2.29 kernel-checked Theorems 1–4 on their base quivers and left one gap: the reduction from an arbitrary acyclic
quiver to the two-class subalgebra. This round closes it. Theorems 3 and 4 are now proved in Lean on every finite
quiver with an acyclicity potential, together with their rigidity statements. All 22 audited results depend only on
Lean's standard axioms (`lean/EatTheFibre/Audit.lean`). Label **proved (Lean)** marks a statement checked by Lean.*

---

## 0. Summary

| # | Statement | Status |
|---|---|---|
| 1 | **Reduction to a sub-basis.** Let `S` be a set of basis indices whose complement spans a two-sided ideal. If `P` kills the basis vectors outside `S`, satisfies the fold conditions on basis vectors in `S`, and `P Δ(e_l) Pᵀ = 0` for `l ∉ S`, then `P` is a fold | **proved (Lean)**, `isFold_of_support` |
| 2 | **Pushing folds.** If a small algebra embeds in a big one (structure constants and Δ agree, the image is closed under products, Δ is local to it, and its complement spans an ideal), a fold of the small algebra extended by zero is a fold of the big one, given the Δ condition off the image | **proved (Lean)**, `isFold_push`; embeddings compose (`Embeds.trans`) |
| 3 | **Path algebras.** For a finite quiver, paths form a basis with the census conventions for products and arrow deletion. An acyclicity potential makes the paths finite. An injective quiver morphism embeds path algebras in the sense of 2 | **proved (Lean)**, `Path.finite`, `Hom.embeds` |
| 4 | **Theorem 4**: on every acyclic quiver with consecutive classes `a₁ ≠ a₂ : u → v`, `b₁ ≠ b₂ : v → w`, the normal form `P_t`, pushed into `k Q`, is a fold for every `t`, over every commutative ring | **proved (Lean)**, `thm4_general` |
| 5 | **Theorem 3**: on every acyclic quiver with converging classes `a₁ ≠ a₂ : u → v`, `b₁ ≠ b₂ : w → v`, `u ≠ w`, and no arrow with both ends in `{u, w}`, all four families push to folds | **proved (Lean)**, `thm3_general` |
| 6 | Rigidity on every acyclic quiver: an automorphism fixing vertices and preserving the radical never conjugates members of a Theorem 3 family with different parameters. For Theorem 4: an `h` that scales `δ_a, δ_b, δ_aδ_b` by `λ, μ, λμ` and intertwines `P_t, P_t′` forces `t = t′` | **proved (Lean)**, `thm3_general_rigid`, `thm4_general_rigid` |
| 7 | The census tables of 1⇉2⇉3 and 1⇉2⇇3 (generated from `pcensus.build`) are the path algebras of those quivers, as defined in Lean | **proved (Lean)** by `decide +kernel` (`W2_c`, `W2_D`, `W3_c`, `W3_D`) |
| 8 | Instances: Theorem 4 on 1⇉2⇉3 + 1→3, and Theorem 3 on 1⇉2⇇3 + 2→4 and on 4→1⇉2⇇3 | **proved (Lean)**, `Examples.lean` |

## 1. How the general proof goes

The two theorems share one argument, now formal.
1. The two classes give an injective quiver morphism φ from 1⇉2⇉3 (resp. 1⇉2⇇3) into Q. The vertices are distinct
   because the potential strictly increases along u → v → w, and for Theorem 3 because u ≠ w is assumed.
2. φ embeds the small path algebra into `k Q` (item 3). Item 7 identifies the small path algebra with the census
   table on which the families were already checked (v2.29).
3. `isFold_push` reduces the fold property in `k Q` to one condition: `P Δ(e_l) Pᵀ = 0` for every path `l` of `Q`.
   `conj_delta_eq_zero` reduces that to: Δ(e_l) has no component on two of P's column paths. Such a component would
   require an arrow from the end of one column path to the start of another:
   - Theorem 4: the columns are a₂ and b₂, so the arrow would run from v or w back to u or v, against the potential;
   - Theorem 3: the columns are e_u, e_w, a₂, b₂, so the arrow would run between u and w (excluded by hypothesis) or
     from v back to u or w (against the potential).

The hand proofs of v2.27/v2.28 argued the kernel closure by "paths outside W". Here that is the ideal condition of the
embedding, which `Hom.embeds` proves from prefixes and suffixes of image paths.

## 2. What is still not formalized

- Theorem 2 beyond 1⇉2⇉3. It needs `P(e_u), P(e_w) ≠ 0`, so its Δ condition involves vertex columns; the same route
  should work.
- That the class substitutions and weightings of an arbitrary quiver satisfy the hypotheses of `thm4_general_rigid`
  (it is shown for 1⇉2⇉3 itself, `subst_scalesV`). For Theorem 3, the vertex-fixing hypothesis is the natural one.
- Census counts and Conjecture 3ᶜ.

## 3. Open questions

1. Formalize Theorem 2 in general, and the families of v2.27's 18 lines on 1⇉2⇇3 with 2→4. `Examples.lean` already has
   four of them.
2. Use the reduction in the other direction: show that a fold of `k Q` supported on the image of φ restricts to a fold
   of the small algebra. That would turn the classification of the families of 1⇉2⇉3 (v2.26, 19 families) into
   statements about every quiver containing it.
