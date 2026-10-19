# Discordianism: Open Questions Explored (v2.28)

*Takes up v2.27's open questions 1 and 2. Question 2 (where the moduli of 1⇉2⇉3 + 1→3 come from) is answered by a
theorem. Every acyclic quiver with two consecutive parallel classes has a one-parameter family of folds that vanish on
every vertex, whatever other arrows are present (Theorem 4). With Theorem 3, this proves the "if" half of the post hoc
reading of v2.27. The "only if" half has passed every test that has finished, and the two remaining tests are handed
over.*

*Scripts: `analysis/theorem4.py`, `analysis/loc_blocks.py`, `analysis/unipotent_merge.py`, `analysis/fixedlocus.py` (new local-torus mode
`FL_LOCAL=1`, and `FL_SIZES`), `analysis/family.py` and `analysis/family_invariant.py` (both local-aware),
`tests/test_symmetry_plan.py`. Logs: `logs/theorem4.log`, `logs/fixedlocus_two_pairs_short_3_local_*.log`,
`logs/fixedlocus_pairs_in_3_local_exact.log`, `logs/results_container.jsonl`. Raw data: `logs/reps_two_pairs_short_{3,5}.pkl`
(G-orbit representatives and their local-torus components). Labels: **proved**, **verified (exact)**, **numerical**,
**conjectured**. Predictions were committed before their runs (commits f04e924, e8f33fa). Discordian glosses are
descriptive labels (L3).*

---

## 0. Summary

| # | Question | Answer |
|---|---|---|
| v2.27 Q2 | The moduli of 1⇉2⇉3 + 1→3 | **Theorem 4:** for consecutive classes A: u→v and B: v→w on any acyclic quiver, P_t(a₂) = (t − 1)e_u − δ_a + t·δ_aδ_b, P_t(b₂) = e_u + e_v − δ_b + δ_aδ_b, P_t = 0 elsewhere, is a fold for every t. t is an invariant of the image **[proved; symbolic on 6 quivers; invariance and distinct orbits checked over 𝔽₃, 𝔽₅]** |
| v2.27 Q1 | The post hoc reading | **"If" half proved** (Theorems 3 and 4). **"Only if" half:** P1 holds (366 over 𝔽₇ too). P2 and P4 are too large for the container and are handed over (`run_batch.py` jobs 12–15) |
| — | Where the 1⇉2⇉3 + 1→3 moduli sit | P3 holds: 21 vertex blocks over 𝔽₃ and 𝔽₅, so the moduli lie in the radical. Theorem 4 lies in the zero block; the full accounting is open |
| — | **Correction** | On 1⇉2⇉3 + 1→3, Aut_Δ° is larger than the census group G: c ↦ c + s·δ_aδ_b. With it, the orbit count is **939 over 𝔽₃ and 983 over 𝔽₅** (not 1668, 1840). It still grows. The converging quivers are unaffected |
| — | Method | The fixed-point-locus test (v2.26) used the whole class torus. That is wrong where weightings are not symmetries. **Local mode:** use the Lie algebra of the local torus at P, {E_w : P compatible with Δ_w}. On 1⇉2⇇3 it gives the same 39 candidates and removes all 36 anomalous rows |
| — | Engine | The census symmetry plan no longer stores H (10⁶ elements for a triple class over 𝔽₅). It uses generators of H restricted to K and BFS transport, and was checked identical to the old plan and to six full 𝔽₃ censuses |

---

## 1. Theorem 4: consecutive classes, no condition on the outer ends

Let Q be acyclic and A: u → v, B: v → w parallel classes with |A|, |B| ≥ 2. Write δ_a = a₁ − a₂, δ_b = b₁ − b₂,
δ_aδ_b = a₁b₁ − a₁b₂ − a₂b₁ + a₂b₂. For x, y, z, w′ ∈ k with x·w′ = y·z − 1 define P by

  P(a₂) = x·e_u − δ_a + y·δ_aδ_b,  P(b₂) = z·(e_u + e_v) − δ_b + w′·δ_aδ_b,  P = 0 on every other basis path.

**Theorem 4.** (i) P is a Δ-compatible fold. (ii) Its image I = span{P(a₂), P(b₂)} contains unique elements X, Y with
(δ_a, δ_b)-coordinates (−1, 0) and (0, −1). Put t(P) = (δ_aδ_b-coefficient of X)·(e_u-coefficient of Y); for P as above,
t = yz. Then t is preserved by G, by every class weighting, and by Aut_Δ° on the quivers checked. On 1⇉2⇉3 + 1→3,
Aut_Δ° also contains c ↦ c + s·δ_aδ_b (§2), which fixes V pointwise. (iii) The normal
forms P_t (x = t − 1, y = t, z = w′ = 1) lie in pairwise distinct orbits. So Folds(kQ) has infinitely many local-torus
orbits over every infinite field.

*Proof.*
- **Reduction to W.** Let W be the span of e_u, e_v, e_w, A, B and the paths a_i b_j: the path algebra of u ⇉ v ⇉ w.
  P maps into V = span{e_u, e_v, e_w, δ_a, δ_b, δ_aδ_b} ⊆ W.
- **Inside W.** Idempotency and closure of the image are computations in V. Kernel closure follows as in Theorem 3:
  ker P = (span of the paths outside W) ⊕ (ker P ∩ W). Products of two paths outside W, and of such a path with an
  element of W, lie outside W. So what remains is (ker P ∩ W)² ⊆ ker P, which is a computation in W.
- **Δ-compatibility.** V ⊆ ker Δ, so ΔP = 0. Since P kills every vertex, (P⊗P)Δ(b) ≠ 0 needs b = π·c·π′ with π, π′
  ∈ {a₂, b₂}. That needs an arrow c out of v or w into u or v, which closes a cycle. Acyclicity excludes it, so no
  hypothesis on the outer ends is needed. This is the difference from Theorem 2, whose family has P(e_u), P(e_w) ≠ 0,
  and an arrow u → w then gives (P⊗P)Δ(c) = P(e_u)⊗P(e_w) ∉ Im Δ.
- **Invariant.** Every h ∈ G and every class weighting fixes the vertices and scales δ_a, δ_b and δ_aδ_b by λ, μ and λμ.
  So h(I) has X′ = λ⁻¹h(X), Y′ = μ⁻¹h(Y), and the product of the two coefficients does not change. ∎

`theorem4.py`:
1. The general member over ℚ(x, y, z) (w′ = (yz − 1)/x) and the normal form over ℚ(t) are folds on six quivers:
   - 1⇉2⇉3, 1⇉2⇛3, 1→2⇉3⇉4 and 1⇉2⇉3→4;
   - **1⇉2⇉3 + 1→3**, where Theorem 2 fails;
   - 1⇉2⇉3 + 1→3, 1→4→3, 2→4, a dense control.
2. Over 𝔽₃ and 𝔽₅:
   - t_of(P_t) = t;
   - every conjugate of a member under G × (all class weightings) that is a fold keeps t (864 and 128 000 conjugates
     on 1⇉2⇉3 + 1→3);
   - on 1⇉2⇉3 + 1→3, the p members lie in p distinct local-torus orbits (from the `--save-reps` dumps).

Theorem 4 makes Theorem 2 a special case for finiteness, though not as a family: the two families differ (Theorem 2's
moves e_u and e_w, Theorem 4's kills every vertex). Both exist on 1⇉2⇉3; how they sit among the 19 families of v2.26 is open (§7, Q2).

**How it was found.** Over 𝔽₃, `fixedlocus.py` in local mode found 1469 orbits that lift to ℚ (199 do not) and 288
candidates. The numerical stage was restricted to the growing sizes 36, 72, 144, one candidate per type. Singular on the
first FAMILY candidate gave one component through P, of dimension 3 > dim N·P = 2, with a single nonlinear generator.
`family_invariant.py` parametrized it and verified the fold equations symbolically. The local torus at P is the full
3-torus. The invariant came out as 1/(y₄₃(y₁₀ + 1)), and the normal form above is that parametrization with the orbit
directions fixed.

**Gloss.** *Two twin steps in a row: a split can drop every place, keep only the second twins, and let them lean on the
twin-difference product. How much they lean, times how far the second pair reaches back to the start, is the modulus.
Short cuts between the ends do not matter, because nothing at the ends is kept.*

---

## 2. 1⇉2⇉3 + 1→3: where the growth sits, and a correction

**Correction to v2.27: the census group is too small on this quiver.** `theorem3.lie_aut_delta` gives dim Lie(Aut_Δ) = 5
here, against dim Lie G = 4. The extra Δ-derivation is c ↦ δ_aδ_b. The arrow c: 1 → 3 runs parallel to δ_aδ_b ∈ ker Δ, so
φ_s: c ↦ c + s·δ_aδ_b is a Δ-compatible algebra automorphism. G and the class weightings normalize the group U = {φ_s}.
- **Affected:** only quivers where consecutive classes have a shortcut arrow (also the dense control).
- **Not affected:** every converging quiver tested (1⇉2⇇3 + 1→3, 1⇛2⇇3 + 1→3, 1⇉2⇇3 + 1→3 + 2→4, 1⇉2⇇3) and
  1⇉2→3⇉4, where dim Lie(Aut_Δ) = dim Lie G.

Merging local-torus orbits along U (`analysis/unipotent_merge.py`, from the `--save-reps` dumps):

| | 𝔽₃ | 𝔽₅ |
|---|---|---|
| local-torus orbits under G (census, v2.27) | 1668 | 1840 |
| **orbits under G, U and the local torus (Aut_Δ° level)** | **939** | **983** |
| vertex blocks | 21 | 21 |
| blocks whose orbit count grows (G-level) | — | 13 of 21 |

The count still grows, as Theorem 4 requires: U fixes V pointwise, so t is a U-invariant too. So v2.27's refutation of
3ʰ by this quiver stands, now by a theorem rather than by counts. Below, block-level numbers are G-level.

The two largest increases are in the block P(e₁) = 0, P(e₂) = e₁ + e₂, P(e₃) = e₃ and its mirror (224 → 274 each). In
the zero block the growth is in the types q²(q − 1)² (109 → 127) and q²(q − 1)³ (86 → 98). Theorem 4's family lies in
the zero block (it kills every vertex). The other blocks need other families, plausibly Theorem 4 composed with vertex
idempotents. They are not yet classified. Numerical family tests in local mode (one candidate per type):
NUMERIC_PLACEHOLDER

---

## 3. The finiteness criterion

| configuration (two classes of size ≥ 2) | status |
|---|---|
| consecutive (any outer ends) | infinitely many orbits: **Theorem 4** (and Theorem 2 when the ends are non-adjacent) |
| converging or diverging, outer ends not joined | infinitely many: **Theorem 3** |
| converging or diverging, outer ends joined | 1⇉2⇇3 + 1→3: 366 over 𝔽₃, 𝔽₅, 𝔽₇ (P1). 1⇛2⇇3 + 1→3 (P2) and 1⇉2⇇3 + 1→3 + 2→4 (P4): handed over |
| vertex-disjoint | 1⇉2→3⇉4: 1158 over 𝔽₃, 𝔽₅ (v2.26) |
| one class | no growth in any quiver tested (Conjecture 3ᵖ) |

The diverging case with joined ends is the converging case seen through op-duality (and a swap of 1 and 3), so it
needs no separate test.

**Conjecture 3ᶜ** (the "only if" half). If no two parallel classes of size ≥ 2 are consecutive, and no converging or
diverging pair has unjoined outer ends, then Folds(kQ) has finitely many local-torus orbits. The "if" direction is
Theorems 3 and 4.

This is recorded with the record of v2.27 in mind: two conjectures failed there. 3ᶜ is the post hoc reading, it has
passed one test since (P1 at 𝔽₇), and P2 and P4 are its next tests.

---

## 4. Engine and method notes

- **Symmetry plan.** The old plan stored all of H = G × grading and an R × |H_K| table. The new one enumerates H once and
  keeps one element per restriction to K, together with N = ker(H → H|_K). It finds a generating set of H_K greedily and
  computes orbits of torus representatives by BFS over generator edges. Each edge is φ_l·g, where the grading element
  φ_l returns the image to its torus representative.
  - Sources: identical to the old plan.
  - Transport elements may differ from the old plan, but each is checked to carry the source to its representative.
    Any such element gives the same fibre.
  - Tests: `tests/test_symmetry_plan.py` on 6 quivers, and full 𝔽₃ censuses of 1⇛2⇇3 + 1→3, 1⇉2⇉3 + 1→3 and
    1⇉2⇇3 with 2→4, all with JSON identical except for the wall time.
  - 1⇛2⇇3 + 1→3 over 𝔽₅ now gets through the plan (5080 H-orbits from 572 944 torus orbits), but is OOM-killed while
    holding its fibres.
- **`--save-reps FILE`** writes the G-orbit representatives, their sizes and their local-torus components, so that
  invariants such as the vertex block can be computed per local orbit (`loc_blocks.py`).
- **Fixed loci, local mode.** The torus part of Lie H at P is {E_w : P is compatible with Δ_w = Σ_c w_c Δ_c}, a linear
  condition on w: the Lie algebra of W(P) (v2.24). The local-torus orbits are exactly the (G × T)-orbits intersected with Folds, so the orbit labels are
  unchanged.

---

## 5. Predictions and outcomes

The predictions as committed (verbatim from the draft):

*Predictions recorded before the censuses and analyses run (this file is completed when they finish).*

The post hoc reading of v2.27 §5: Folds(kQ) has infinitely many local-torus orbits iff Q has consecutive parallel
classes (size ≥ 2), or converging/diverging ones whose outer ends are not joined by an arrow.

- **P1.** 1⇉2⇇3 + 1→3 over 𝔽₇: 366 local-torus orbits (as over 𝔽₃ and 𝔽₅).
- **P2.** 1⇛2⇇3 + 1→3 (a larger converging class, outer ends adjacent): the local-torus count is the same over 𝔽₃ and 𝔽₅.
- **P3.** 1⇉2⇉3 + 1→3: its moduli are in the radical, as for 1⇉2⇉3, so the number of vertex blocks is the same over 𝔽₃
  and 𝔽₅.

**Outcomes so far (appended in order):**

- **P1 holds.** 1⇉2⇇3 + 1→3 over 𝔽₇: N = 502 306, 651 G-orbits, 392 under G × T, **366** under the local torus.
- **P2: not yet decided.** Over 𝔽₃: 752 local-torus orbits (N = 316 612). The 𝔽₅ run exceeds the container: first the
  symmetry plan (H has ~10⁶ elements; the engine is now rewritten, see `tests/test_symmetry_plan.py`), then memory
  (OOM at 6.1 GB with ~5.7·10⁵ fibres). Handed over.
- **P3 holds.** 1⇉2⇉3 + 1→3 has 21 vertex blocks over 𝔽₃ and over 𝔽₅ (`analysis/loc_blocks.py`). Its local-torus
  growth (1668 → 1840) sits inside 13 of them, in orbit types of size q²(q − 1)^k, k = 2, 3, 4.

**P4 (recorded before the census).** 1⇉2⇇3 + 1→3 + 2→4 (converging classes, adjacent outer ends, a tail at v): the
local-torus count is the same over 𝔽₃ and 𝔽₅.

**Final status:**
- P1 holds.
- P3 holds.
- P2 and P4 are open: `run_batch.py` jobs 12 (P2), 13–14 (P4) and 15 (1⇉2⇇3 + 1→3 over 𝔽₁₁, predicted 366).
  - The 𝔽₃ census for P4 has 18 basis paths. In the container it did 4096 of 17 748 fibres in 30 minutes, so it was
    stopped.

---

## 6. Proposed amendments (v2.27 → v2.28)

| amendment | content |
|---|---|
| **Q1²⁸ Theorem 4** | Consecutive parallel classes on any acyclic quiver ⇒ a one-parameter family vanishing on vertices, with an image invariant |
| **Q2²⁸ Criterion, "if" half** | Consecutive classes, or converging/diverging classes with unjoined ends ⇒ infinitely many local-torus orbits (Theorems 3, 4) |
| **Q3²⁸ Conjecture 3ᶜ** | The converse (§3), with Conjecture 3ᵖ as its one-class case |
| **Q4²⁸ Method** | Fixed-point loci use the local torus at P (`FL_LOCAL=1`); full-torus verdicts are void where weightings are not symmetries |

---

## 7. Open questions

1. Conjecture 3ᶜ: P2 and P4 (handed over), then a proof. The converging case with joined ends is the one place where the
   kernel fibration must forbid moduli. Why does the arrow u → w kill Theorem 3's families but not Theorem 4's?
   Theorem 3 needs P(e_u) and P(e_w) nonzero together, and Theorem 4 needs neither. A proof could look for a family
   that kills e_u and e_w on converging classes and show none exists.
2. Classify the families of 1⇉2⇉3 + 1→3 (13 growing blocks), and of 1⇉2⇉3 (19, v2.26) in terms of Theorems 2 and 4.
3. Prove the 18 lines of 1⇉2⇇3 with 2→4 (v2.27).
4. Carried: Conjecture 3ᵖ; unobstructed violations; Kronecker plus a loop; zigzag; dihedral; ρ.

---

## 8. Sources

- W. Decker, G.-M. Greuel, G. Pfister and H. Schönemann, *Singular* (via passagemath-singular): minimal associated primes.
- M. Aguiar, *On the associative analog of Lie bialgebras*, J. Algebra 244 (2001).

*"Two twin steps in a row: drop every place, keep the second twins, and let them lean on the product of the
differences. Short cuts between the ends change nothing, because nothing at the ends is kept."*
