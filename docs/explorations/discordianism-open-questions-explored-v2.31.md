# Discordianism: Open Questions Explored (v2.31): the v2 axioms in Lean, and what they say about folds

*The earlier rounds left the v2 axioms (`docs/axioms/discordianism-axioms-v2.md`) as checked claims about models. This
round derives them in Lean 4 + Mathlib from the smallest hypotheses that suffice, model-checks the rest, and lists what
remains informal. It then applies the algebraic layer to the folds of the research project:*
- *folds are exactly the idempotent Rota–Baxter(−1) operators that commute with Δ;*
- *the dendriform layer transfers to folds; the Loday layer does not;*
- *v2's δ-D2 classification refines into four strata, keyed to how the kernel of the fold sits over its image.*

*The 66 audited results (22 earlier, 44 new) all depend only on Lean's standard axioms (one, `calendar`, on none; count corrected in v2.32). No `sorry`, no `native_decide`
(`lean/EatTheFibre/Audit.lean`).*

---

## 0. Summary

| # | Statement | Status |
|---|---|---|
| 1 | **Bridge.** On any coordinate algebra (associativity not needed), `P` is a fold ⟺ `P² = P`, `P` is Rota–Baxter of weight −1, and `ΔP = (P⊗P)Δ` | **proved (Lean)**, `isFold_iff_rb` |
| 2 | Every member of the Theorem 3 and Theorem 4 families, on every acyclic quiver, is an idempotent RB(−1) operator commuting with arrow deletion. This answers v2 open question 2 ("a fold that deconstructs") in the path-algebra setting | **proved (Lean)**, `thm3_isRB`, `thm4_isRB` |
| 3 | M4, the RB(−1) identity, the five Loday axioms δ-L, the Leibniz identity δ-L1 and E1 all follow from M3 alone (an idempotent ring endomorphism) | **proved (Lean)**, `Fold.lean` |
| 4 | δ-D (three dendriform axioms) and δ-D1 (star splitting, (A,⋆) ≅ im × ker) need only RB(−1). E1⋆ needs only idempotence | **proved (Lean)**, `dD_1..3`, `dD1_hom`, `E1star`, `fold_star` |
| 5 | **δ-L does not transfer to folds.** For Theorem 4's fold, Loday's A1 fails at (e₂, e₂, b₂), and `P(e₁a₂) ≠ Pe₁·Pa₂`. Both hold for every `t` over every nontrivial ring | **proved (Lean)**, `thm4_not_dialgebra`, `thm4_not_multiplicative` |
| 6 | **δ-D2 refined** (table in §3): Baxter point for every fold; line 1 + point 0 when ker·im ⊆ ker; v2's two lines + two points when ker is an im-bimodule (every idempotent endomorphism); the points thicken to lines when ker·im·ker = 0 | sufficiency **proved (Lean)**, `Completion.lean`; necessity generic, `analysis/completion_generic.py` |
| 7 | Census placement: Theorem 4 and Theorem 3 F1–F3 have only the Baxter point. Theorem 3 F4 has line 1 + point 0. Theorem 2 (an endomorphism) has four lines | computed, `analysis/completion_census.py`; Theorem 4's module failure **proved (Lean)**, `thm4_not_module` |
| 8 | 𝔼: L2 (𝔼 = S_fde as consequence relations), L3 (right variable-inclusion companion of LP), L4, the L5 reduction, Meta-1 (4 of 49) | **proved (Lean)**, `Logic.lean` |
| 9 | M2, N5, N6 witness, β4, B4 (Eckmann–Hilton, and B4 ⇒ "no neutral chaos"), B1 (cardinality and null-measure parts), B3 period table, §9 arithmetic | **proved (Lean)**, `Ops.lean` |

## 1. The bridge: folds are the "folds that deconstruct"

v2 open question 2 asked for an idempotent Rota–Baxter(−1) map on the grid algebra that is also an ε-bialgebra
morphism. The research project studied *folds*: idempotent `P` whose image and kernel are subalgebras and with
`ΔP = (P⊗P)Δ`. The two notions are the same. The whole proof is one bilinear identity, true for any bilinear product:

  `P x · y + x · P y − x y = P x · P y − (x − P x)(y − P y)`.

- *Fold ⇒ RB.* Apply `P`. The first term stays (the image is closed), and the second dies (the kernel is closed).
- *RB ⇒ fold.* Substitute `x, y ∈ im P` into the RB identity, which gives image closure. Substitute `x, y ∈ ker P`,
  which gives kernel closure.

`isFold_iff_rb` proves this for every `CoordAlg`, so it covers every census table and every path algebra. Combined
with v2.30's `thm3_general` and `thm4_general`, Lean now exhibits infinite families of such maps on every acyclic
quiver that has two consecutive (Theorem 4) or two converging (Theorem 3) pairs of parallel arrows. Axiom δ-ε2 ("the
fold f is not a coalgebra map") stays true of v2's f. The slogan "order does not respect deconstruction" does not
become a theorem: some orders do.

## 2. Which axioms survive the passage from f to a fold

v2 derives its algebra from M3: f is an idempotent algebra **endomorphism**. Project folds are not endomorphisms
(item 5). So each axiom is now tagged by the weakest hypothesis that proves it:

| axiom | hypothesis used in Lean | holds for project folds? |
|---|---|---|
| M4 splitting (`M4_split`) | idempotent endomorphism | yes, by definition of fold |
| M4 Rota–Baxter (`M4_rotaBaxter`) | idempotent endomorphism | yes (item 1) |
| δ-L A1–A5, δ-L1 Leibniz | idempotent endomorphism | **no**: A1 fails (11 basis triples for Theorem 4), as do A2, A4, A5 |
| E1 (Leibniz annihilators) | idempotent endomorphism | not applicable, since δ-L fails |
| δ-D dendriform (`dD_1..3`) | RB(−1) | yes |
| δ-D1 star splitting (`dD1_split`, `dD1_hom`, `fold_star`) | RB(−1); the split itself is a ring identity | yes: `P(x⋆y) = PxPy`, `(1−P)(x⋆y) = −(x−Px)(y−Py)` |
| E1⋆ (⋆-annihilator = Z(im) ⊕ Z(ker)) | idempotent only | yes |
| δ-D2 completion family | see §3 | only the Baxter point, generically |

E1 needed one correction. The left annihilator of the Leibniz bracket is the **centralizer** of im f (`E1_left`). It
equals im f exactly when im f is self-centralizing, as in the free model, where im f = k[a]. v2's wording ("exactly
im f") is a property of that model, not a consequence of M3.

## 3. δ-D2, derived and refined

Fix `x ≺ y = x P(y)`, and let `x ≻ y = b₁xy + b₂P(x)y + b₃xP(y) + b₄P(x)P(y)`. Write each variable as `p + n` with
`p ∈ im P` and `n ∈ ker P`. Then `P` of any two-piece product is determined by the class of the fold:

| fold class | completions (b₁, b₂, b₃, b₄) | Lean |
|---|---|---|
| any fold | Baxter point (−1, 1, 0, 0) | `baxter` |
| ker·im ⊆ ker (`LeftMod`) | line 1 (−1−s, 1+s, s, −s); point 0 | `line1`, `point0` |
| im·ker ⊆ ker only | Baxter point (`≺ = xP(y)` breaks the left/right symmetry) | (no new completion) |
| ker an im-bimodule | line 1; line 2 (−t, 1+t, t, −1−t); points 0 and (0, 1, −1, 0): **exactly v2's δ-D2** | `line2`, `pointB` |
| … and ker·im·ker = 0 | the points become lines (c, −c, −c, c) and (c, 1−c, −1−c, c) | `line0`, `lineB` |

- Every idempotent algebra endomorphism is in the bimodule row (`endo_bimod`). So v2's δ-D2 holds for *every* M3
  projection, not only in k⟨a,b⟩, where [R4] computed it.
- In the decomposition the lines read cleanly:
  - line 1: `x ≻ y = −n_x (q_y + (1+s) m_y)`;
  - line 2: `x ≻ y = (p_x − t n_x) m_y`.
- The D2 defect of line 1 on a general fold is exactly `s · n_x · P(n_y q_z)`. That is why ker·im ⊆ ker is the
  threshold.
- The Lean proofs are uniform: decompose, `simp` with the product rules of the class, close with `module`. Scalars
  range over any commutative ring `R`, with `A` an `R`-algebra.

**Necessity.** `analysis/completion_generic.py` treats words in the pieces as free and solves the dendriform equations
in each class. It returns exactly the rows above, so nothing else is a completion generically. Specific algebras can
have more, as the last row shows.

**Census placement** (`analysis/completion_census.py`, run on the generated Lean tables):

| fold | ker·im ⊆ ker failures | im·ker ⊆ ker failures | completions |
|---|---|---|---|
| Theorem 2 on 1⇉2⇉3 (endomorphism; ker·im·ker = 0) | 0 | 0 | four lines |
| Theorem 4 on 1⇉2⇉3 (symbolic t) | 5 | 3 | Baxter point |
| Theorem 3 F1, F2, F3 on 1⇉2⇇3 | 4, 16, 4 | 18, 10, 18 | Baxter point |
| Theorem 3 F4 on 1⇉2⇇3 | 0 | 18 | line 1 + point 0 |

The failure counts are basis pairs. `thm4_not_module` certifies one failure of each kind for Theorem 4 at every `t`.
F4 is the one family whose kernel is a right module over its image. This new invariant separates F4 from F1–F3,
and it is cheaper to test than the rigidity statements of v2.29–v2.30. ker·im·ker vanishes for Theorem 2's fold (0 of
1331 basis triples fail) but not for F4 (140 fail) or Theorem 4 (72 fail). That accounts for the four lines in the
first row and the unthickened point in the last.

**Correction to v2.** δ-D2's "exactly two lines and two points" is the generic answer for endomorphisms. When
ker·im·ker = 0, as for Theorem 2's fold on 1⇉2⇉3, the two points are lines.

## 4. 𝔼 in Lean

`Logic.lean` models formulas over ℕ-indexed atoms with ¬, ∧, ∨, as follows.
- 𝔼's values are subsets of {T, F, M} (nonempty for consequence), with the connective-wise lift of Bochvar's tables.
- S_fde is the four-valued matrix: LP on T, B, F, with N infectious.
- `collapse` (the classical part) is a matrix homomorphism `𝔼 → S_fde`, decided over all 8 × 8 subset pairs. It
  preserves and reflects designation and has a section, so `E_iff_Sfde` follows: **L2**.
- **L3** is proved for all finite Γ. The ⇒ direction uses the valuation "N on the missing atom, B elsewhere". The ⇐
  direction repairs N-atoms to B; designated premises have no N-atoms, so nothing they mention changes.
- **L4**: all-{M} kills every formula, and all-{T,F} designates every formula.
- **L5**: `¬SatCL γ ⟺ [γ] ⊨_𝔼 ⋁_{p∈Var γ}(p ∧ ¬p)`. This is the coNP-hardness reduction; membership needs no proof.
- **Meta-1**: `¬(μ∨ν) ∧ (μ∧ν)` is designated exactly when μ and ν both contain T and F. That is 4 of the 49
  assignments, by `decide`. "Glut" in v2 means a value containing both T and F, including {T,F,M}.

## 5. The operational layer

| axiom | Lean form |
|---|---|
| M2 | `x ++ y = [] → x = [] ∧ y = []` |
| N5 | a Mathlib `Coalgebra` with zero counit is trivial: `(ε ⊗ id)Δ = id` forces `x ⊗ 1 = 0` |
| N6 | `(ℤ/5, a − b)`: `ν a a = 0`, non-associative, no two-sided unit (`decide`) |
| β4 | `(a + b) − (a + c) = b − c` in any abelian group |
| B4 | Eckmann–Hilton (`EckmannHilton.mul`). `B4_no_neutral_chaos`: a non-associative chaos interchanging with a unital order has no unit. So N6's "no unit" is *derived* from B4 together with non-associativity |
| B1 | hand-off grids are countable and `#(ℕ → Bool) = 𝔠`; every singleton, and so the countable set of grid-named chaoses, is null under `infinitePi β` whenever `β{b} < 1` for both `b` (bound `q^n` on cylinders) |
| B3 | the minimal periods of the 256 words of length 8 are 2, 2, 6, 12, 30, 48, 82, 74 (`decide +kernel`) |
| §9 | `catalan 4 · 5! = 1680`, `5 · 73 = 365`, `73 % 5 = 3` |

## 6. Not formalized

- **Coinductive Nu side.** N1–N4 (finality, bisimulation, productivity, stochasticity); M1 and M5; δ-H (δ1–δ3).
  These need a coinductive chaos type, and Mathlib has no ready one for `Dist`-valued streams.
- **B1 topology.** Σ⁰₂-completeness via Wadge's lemma, and "not ω-regular". Mathlib has neither Wadge reducibility
  nor ω-automata. B1's non-identical product caveat (Σ(1−pₙ) < ∞) is also not covered: the Lean statement is i.i.d.
- **Unpinned definitions.** B2 entropy, the B3 hand-off cost claim (236 of 256), B5, and the measure layer C0–C4. The
  axioms do not fix definitions precise enough to state these.
- **Remaining logic and coalgebra items.** L3's PWK and B₃ neighbours; δ-ε1 and δ-ε3 (coassociativity on UT_n, no
  Apple element). The latter are finite computations and good candidates for the reflection machinery of v2.29.
- **δ-D2 necessity in Lean.** It is model-specific; the generic computation stands in for it.

## 7. Next

1. Prove the δ-D2 necessity in Lean for one free model (k⟨a,b⟩ truncated at length 3 suffices to separate the
   strata), using the `CoordAlg` reflection.
2. Use "ker is a right/left im-module" as a census invariant over F₁₁ alongside jobs 12–15. It is linear in P for
   fixed P and cheap to test.
3. δ-ε3 on UT₃ by reflection: no `r` with `Δ = Δ_r`.

Files: `lean/EatTheFibre/Discordia/{Fold,Completion,Bridge,Logic,Ops}.lean`, `analysis/completion_generic.py`,
`analysis/completion_census.py`.
