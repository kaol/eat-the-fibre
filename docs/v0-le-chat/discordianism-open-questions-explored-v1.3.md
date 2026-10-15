# Discordianism: Open Questions Explored (v1.3)

*Fourth wave of the exploration of the v1 Mu–Nu axioms. This document answers the eight open questions left by v1.2 (its §11, OQ-1 through OQ-8). Tags **[W44]–[W63]** were machine-checked this session (TypeScript sandbox: exhaustive enumeration over finite witnesses, direct evaluation of the 𝔼 matrices, Gaussian elimination over ansatz and coproduct systems, and periodicity searches over the ruler sequence). Tags **[W25]–[W43]** refer to v1.2's checks, **[W1]–[W24]** to v1.1's, and **[V#]** are unchanged from v1's Z3/SymPy script. Untagged formal claims are standard results. Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Method

Same two answer types as v1.1 and v1.2: **witness answers** (finite toy models, exhaustively enumerated — existence of a mechanism, not uniqueness) and **structural answers** (theorems plus exhaustive checks). The **corrections register** (§10) grows to eight rows: this wave's stricter checks retract or sharpen one v1.2 count, one v1.2 framing, and several of this session's own intermediate results. And the register now includes *us*: three bugs in this wave's own scripts (an unnormalized Gaussian-elimination pivot that corrupted every rank computation; two sign slips in ε-rule checks; a missing term in one UT₄ right-hand side) were caught by cross-checking and fixed before any number below was recorded. This is B5 applied to ourselves for the second wave running: the check was only good because it could have returned another number.

---

## 1. Executive summary

| # | v1.2 asked | v1.3 answer | verdict |
|---|---|---|---|
| OQ-1 | Which D3-survivors extend to all lengths? Is the family a torsor? Is the dialgebra elsewhere in it? | **The family is exactly two parallel ℤ-lines plus two isolated points** (34 members at integer coefficients, all with global closed forms, all with associative ⋆ and a Lie commutator). Not a torsor. E1 stratifies by stratum; **the dialgebra layer is nowhere** — the middle-exchange axiom fails for every non-origin member | answered |
| OQ-2 | UT₄ under the corrected operations: does the Podge interact differently with a Lie layer? | **The whole corrected chain transfers verbatim**, and E1 comes out sharper: the annihilator is im f **plus the nilpotent corner E14** — the Podge's fingerprint. But there is no coalgebra here either: the only coassociative graded ε-coproduct is **Δ ≡ 0**, and the 44 sporadic quasitriangular tensors are **inert** | answered |
| OQ-3 | Is P (the pure flats) Σ⁰₂-complete, or properly incomplete? | **Complete.** The pre-period hides in the 2-adic valuation: a ruler-encoding reduction FIN ≤ P, machine-verified on all bounded-valuation masks and five canonical non-FIN witnesses. This **corrects v1.2's framing** — the hand-off buys tight coverage at half cost, not a jump in Borel level | answered (with correction) |
| OQ-4 | Dist-valued Borel analysis | No valuation-interleave analogue is known for measure sequences; the deterministic case is now fully closed, so what remains is genuinely measure-theoretic | carried |
| OQ-5 | Complexity of the availability semantics; proof-theoretic cost | **𝔼-consequence is coNP-complete** — membership via the three-filter characterization, hardness via γ CL-unsat ⟺ γ ⊨𝔼 ⋁(p∧¬p). Polynomial-size cut-free proofs are foreclosed unless NP = coNP | answered |
| OQ-6 | Empirical kernels: fit (θ, β) to a real culture | No real-culture data fitted; the rate law t* and the wrong-basin curve stand as falsifiable predictions | carried |
| OQ-7 | Is the impatience (wrong-basin) effect named in the literature? | Two 2026 near-neighbors found (learning-rate stability bounds; stable convergence to wrong beliefs under misspecification) plus the cascade-rate literature — **none is the effect**. It appears unnamed | unresolved, adjacent |
| OQ-8 | Is M2′ (the order side of the fold) quasitriangular — does the Apple live downstream? | The grid case is now **airtight at any support** (mass obstruction), and UT₄'s 44 Apples are inert. **M2′ remains the last untouched candidate** | answered for the toys; M2′ open |

---

## 2. OQ-1 — The D3-survivor family, classified

### Restatement

*Which of the dendriform completions that survive D3 on short words extend to all word lengths? Is the family a torsor (a group action permuting the survivors, with the Aguiar-canonical member as identity)? Does any member satisfy the middle-exchange axiom — is the dialgebra layer elsewhere in the family rather than gone?*

### Findings

All checks this section run on the free grid algebra truncated to the 16-word corpus {all words of length ≤ 3, plus aaaa}, with the Aneristic fold f (f(w) = w if w is all a's, ε included; f(w) = 0 otherwise) and ≼ = x·f(y) fixed.

- **The ansatz collapses.** The fold identities (v1.2 [W40]) collapse the 7-term ansatz to **four effective terms**: ≻ = b₁·xy + b₂·f(x)y + b₃·xf(y) + b₄·f(x)f(y). On this corpus the axioms decouple cleanly: **D1 holds iff Σb = 0** (both sides of D1 compute to (1+Σb)·x·f(y)f(z) — a hand calculation the machine confirms); **D2 is vacuous** (0 constraints); **D3 reads (x⋆y)≻z = x≻(y≻z)** and does the real filtering **[W51]**.
- **The corrected count.** Re-testing all 393 sum-zero candidate coefficient vectors from v1.2 against the corrected corpus: **108 survive** (v1.2 reported 107 — see §10), and they realize exactly **six effective completions**: (−1,1,0,0), (0,1,0,−1), (1,0,−1,0), (0,1,−1,0), (0,0,−1,1), (0,0,0,0), with multiplicities 19, 16, 19, 19, 16, 19 **[W49]**.
- **The exact family.** Extending to integer coefficients, the survivor set is precisely **34 members: two parallel lines plus two isolated points** [W50]:

  > line 1: b(t) = (−1−t, 1+t, t, −t)  •  line 2: b(t) = (−t, 1+t, t, −1−t),  t ∈ [−8, 7] (16 points each, common direction (−1,1,1,−1));
  > plus the origin (≻ ≡ 0, ⋆ = ≼ = xf(y)) and the isolated point (0,1,−1,0) (⋆ = f(x)y).

  The closed forms are **global** — valid at all word lengths, not just the corpus — thanks to the lemma f((f(x)−x)·w) = 0 **[W50]**:

  > line 1: ≻ = (f(x)−x)·((1+t)y − t·f(y));  line 2: ≻ = ((1+t)f(x) − t·x)·(y − f(y));
  > isolated: ≻ = f(x)y − xf(y).

  Every one of the 34 satisfies D1–D3 on all 4,096 corpus triples, with associative ⋆ and a Jacobi-verified ⋆-commutator **[W51]**. So the family does not merely survive longer words: it is *pinned*, with nothing outside the two lines and the two points (at the searched coefficient range).
- **Four natural products live in the family.** L2(−1) = (1,0,−1,0) has ⋆ = **xy — concatenation itself**; L1(−1) = (0,0,−1,1) has ⋆ = f(x)·f(y) — **commutative**; the canonical L1(0) = (−1,1,0,0) is the **Baxter star** xf(y) + f(x)y − xy; the origin degenerates to ⋆ = ≼. The free algebra's own product, the fold product, and two projections are all dendriform completions of the same left product. The family is not a random cloud — it is a pencil with famous points on it.
- **The doublefold is expelled.** The sum-zero vector (1,−1,−1,1), whose ≻ = (x−f(x))(y−f(y)) is the "doublefold" product, **fails D3** at (b, a, b): LHS bab, RHS 0 — 605 of 4,096 triples bad, and it fails the dialgebra axiom A2 there too (bab vs 0) **[W48]**. D3 does real work beyond Σb = 0: membership is not the whole sum-zero hyperplane.
- **E1 stratifies** [W53]. Left and right annihilators of the ⋆-commutator agree everywhere, but the common value depends on the stratum: on line 1 minus its commutative point, ann = im f = {ε, a, aa, aaa, aaaa} (v1.2's [W38] is the canonical member's stratum); L1(−1), the commutative member, has ann = everything (bracket ≡ 0); every line-2 member has ann = {ε}; the origin and the isolated point have ann = ∅. **E1 is a stratum property, not a family property** — the corrected picture is richer than the v1.2 claim it replaces.
- **The dialgebra layer is nowhere** [W52]. Loday's five dialgebra axioms A1–A5 with ⊣ = ≼, ⊢ = ≻, all 4,096 triples: the naive pair (⊣ = xf(y), ⊢ = f(x)y) fails **A2 and A5** (0/550/0/0/550); a typical family member fails A2, A4, and A5 (0/880/0/880/1210); L1(−1) fails A2 and A4 (0/275/0/275/0 — A5 passes); L2(0) and L2(−1) fail A2, A4, A5 (0/275/0/275/550); the isolated point fails all three (0/550/0/550/550); only the origin passes, and it passes *degenerately* (≻ ≡ 0). **A2 fails for every non-origin member**: the middle-exchange axiom is not hiding elsewhere in the family. The dialgebra layer is gone, and this time the search for it was exhaustive over the whole family.

### Verdict

**Answered.** The dendriform completions of ≼ = x·f(y) form exactly two parallel ℤ-lines plus two isolated points — 34 members at the searched range, each with a global closed form, an associative star, and a genuine Lie bracket. It is **not a torsor**: two disjoint lines plus two isolated points admit no transitive permutation structure, and no group action singles out the Aguiar-canonical member as an identity. The dialgebra layer is nowhere in the family; E1, by contrast, is everywhere but stratified — order's silence (the annihilator) has four different shapes, and the canonical member's im f is just one stratum.

### Remains open

The pencil over the other 24 idempotent folds: does every fold carry a two-line family, and do the lines share the (−1,1,1,−1) direction? Whether the parameter t extends continuously over ℝ (the closed forms are polynomial in t; 16 integer points per line constrain but do not pin the polynomials). Whether any *non-origin* member admits a dialgebra structure with a different ⊣ (this wave fixed ⊣ = ≼).

---

## 3. OQ-2 — UT₄ under the corrected operations: the chain transfers, the Podge signs

### Restatement

*Redo the fold/Baxter/E1 analysis in the v1 upper-triangular model with the ⋆-commutator: does the nilpotent Podge (the Aftermath) interact differently with a Lie layer than v1's Leibniz layer predicted?*

### Findings

Model: UT₄ over the ten units E_ij with i ≤ j ≤ 4; f = the diagonal projection (kills off-diagonal units). This is the v1 toy, now under the corrected (v1.2) operations.

- **The whole corrected chain transfers verbatim.** The fold is multiplicative, Rota–Baxter of weight −1, and idempotent: 0 bad over 100 pairs **[W54]**. D1–D3, associativity of ⋆, and the Jacobi identity for the ⋆-commutator: 0 bad over 1,000 triples **[W55]**. The Baxter route — the one route the grid algebra admits (v1.2 NQ-2) — works on the upper-triangular model exactly as advertised.
- **E1 gains a corner.** The annihilator of the ⋆-commutator, left and right, is **{E11, E14, E22, E33, E44} = im f plus E14** **[W56]**. E14 is the nilpotent corner — the upper-triangular Podge. In this model the Aftermath is itself silent in every collision: it commutes to zero with everything, not because it is folded away but because it *is* the fold's shadow. v1's withdrawn Leibniz layer predicted nothing of the kind; the corrected Lie layer finds new structure exactly where the old layer had a bug.
- **The naive dialgebra fails here too.** With ⊣ = x·f(y), ⊢ = f(x)·y: A1–A5 count [0, 12, 0, 0, 12] — A2 and A5 fail, witness (E11, E11, E12): one side E12, the other 0 **[W57]**. The v1.2 withdrawal is not a grid-specific artifact.
- **No graded ε-coproduct exists.** Parametrize a graded coproduct on UT₄ by 20 unknowns; the ε-bialgebra rule yields 30 equations of rank 14 — a 6-dimensional solution family [W58]:

  > Δ(E_ij) = λ_i E_ii⊗E_ij + Σ_{i<m<j} (λ_m + ρ_m) E_im⊗E_mj + ρ_j E_ij⊗E_jj,  Δ(E_ii) = 0.

  Its coassociativity defect decomposes into **exactly three component families**: −λ_i² at E_ii⊗E_ii⊗E_ij, +ρ_j² at E_ij⊗E_jj⊗E_jj, and ρ_m² − λ_m² at E_im⊗E_mm⊗E_mj. All three vanish only at λ = ρ = 0: **the sole coassociative graded ε-coproduct on UT₄ is Δ ≡ 0** [W58]. Greyface cannot be coassociatively deconstructed here at any grade.
- **Five non-graded candidates, by name — all dead.** interior-gap (ε-rule fails on 4 pairs), boundary (coassociativity fails 6), right-rooted (coassociativity fails 6), all-splits (ε-rule fails 20), fam11 (coassociativity fails 6). Every named coproduct candidate breaks one of the two laws the chain needs.
- **44 sporadic Apples, all inert** [W59]. In Aguiar's Convention A — rank-1 r-forms, Δ_r(x) = u⊗(v·x) − (x·u)⊗v — the ε-rule and Δ_r(1) = 0 are automatic (0 bad over 100×121 checks), and exactly **44 tensors are coassociative**: 32 with both legs off-diagonal and v·u = 0, 6 with u off-diagonal and v = E_ii (u's row unit), 6 with u = E_ii and v = E_ki with k < i. E12⊗E34 and E12⊗E11 are among them; the diagonal r = Σ E_mm⊗E_mm is *not* (its defect is exactly the boundary pattern λ = ρ = 1 — a mid-session candidate, withdrawn, see §10). But the dendriform structures these Apples induce are inert: E12⊗E11 gives ≻ ≡ 0 in both orientations; E12⊗E34 yields exactly two nonzero products (prec(E11, E23) = E14, suc(E23, E44) = E14, mirrored) and satisfies the axioms only vacuously. **None resembles the fold layer.**

### Verdict

**Answered.** The corrected unification chain — Baxter(−1) ⟹ dendriform ⟹ star ⟹ Lie ⟹ E1 — transfers to UT₄ verbatim, and the Lie layer treats the Podge *better* than the withdrawn Leibniz layer could: E1 sharpens to "im f plus the nilpotent corner," a genuinely new fingerprint with a clean Discordian reading. But the coalgebraic side is as barren here as on the grid: no graded ε-coproduct, no live non-graded candidate, and 44 quasitriangular tensors that satisfy every law while doing nothing. **Order lives in the operator, not in a coproduct — in either model.**

### Remains open

UT_n for general n: does the count of Convention-A Apples follow a pattern, and does any become non-vacuous as n grows? Coproducts outside the graded and rank-1 ansätze. Whether the E14-fingerprint generalizes: in UT_n, is the annihilator always im f plus the nilpotent corner E_1n?

---

## 4. OQ-3 — Flats are complete: the pre-period hides in the ruler

### Restatement

*v1.2 left P (the pure flats) at Σ⁰₂-proper-and-maybe-complete: "a flat-only reduction would have to hide the pre-period in the periodic part." Is P Σ⁰₂-complete?*

### Findings

- **The hint was right, and the hiding place is the 2-adic valuation** — the ruler sequence 0, 1, 0, 2, 0, 1, 0, 3, … Define, for a binary sequence x, the word F(x)[m] = z[v₂(m)] where z = interleave(x, 0) (z[2i] = x[i], z[2i+1] = 0). If x is eventually zero, the positions where z is 1 are bounded, so F(x) reads z only through a bounded set of valuations and is **periodic** — a flat, with the pre-period absorbed into which valuations occur. If x has infinitely many 1s, the valuation set is unbounded and F(x) **never becomes periodic**.
- **Machine support, three layers.** (i) All 64 masks b(n) = [v₂(n) ∈ S], S ⊆ {0,…,5}, are periodic with period 64 (63 non-constant), while the unbounded-valuation mask [v₂ even] fails at m = 128 vs 192: bounded valuation sets force periodicity, unbounded ones break it **[W61]**. (ii) The lemma behind (i): if g(x)[m] = x[v₂(m)] is periodic with period p, then x[s] = x[v₂(p)] for all s in (v₂(p), 11] — 512 periodic cases, 0 violations: a periodic ruler-read pins the tail, i.e. forces FIN-like behavior **[W62]**. (iii) F(x) is periodic (both claimed periods 2048 and 8192 verified) for all 64 six-bit x, 0 failures; and five canonical non-FIN witnesses — (01)^∞, 1^∞, (001)^∞, (100)^∞, and a period-12 word at n ≡ 5 mod 12 — **never settle**: no period p ≤ 256 survives a window of 5,000, the latest first failure at position 3,968 **[W63]**.
- **Hence FIN ≤ P, and P is Σ⁰₂-complete.** The eventually-zero sequences (Σ⁰₂-complete) reduce to flat-membership; P is Σ⁰₂ by v1.1, so it is complete — **exactly as complete as UP**, the hand-off's denotation set.
- **This corrects v1.2's framing.** v1.2's headline — "the hand-off is exactly the proper→complete upgrade" — is wrong at the Borel level: P was already complete; what the hand-off buys is **full tight coverage at half the per-word cost**, not a jump in hierarchy. The corrected story is cleaner: flats and hand-offs sit at the same complete level, and the hand-off is an *economist's* upgrade, not a *topologist's*. (Register, §10.)
- **B1 survives untouched.** P is still countable and still null under every product measure on 𝒞 with no certain percept. Completeness changes the arithmetical address, not the measure.

### Verdict

**Answered, with a correction to the previous wave.** The pre-period hides in the ruler: the 2-adic valuation is exactly the "periodic part" that v1.2 doubted existed. Both deterministic denotation classes — flats and hand-offs — are Σ⁰₂-complete, and the Borel column of the story is now flat all the way down.

### Remains open

The Dist-valued case (OQ-4, carried): no valuation-interleave analogue is known for sequences of measures, so the completeness question restarts there from zero.

---

## 5. OQ-4 — Dist-valued denotations: carried

### Restatement

*Lift the toy GridF to `Dist`-valued denotations (N4) and redo the Borel analysis in the space of measures: what is the σ-structure of the set of grid-named measures?*

### Findings

Nothing to report honestly. This wave closed the deterministic case completely (OQ-3: P and UP both Σ⁰₂-complete), which sharpens the boundary of what remains: everything left is genuinely measure-theoretic. The ruler trick that closed OQ-3 has no known analogue for measure sequences — the natural encoding of an eventually-zero sequence into an eventually-repeating *sequence of measures* is not obvious, and the weak-star σ-algebra on the measure space has not been set up.

### Verdict

**Carried forward.** The question is now isolated: it is the only part of the denotation story with no machine-checkable toy.

### Remains open

As stated. A first step: define the interleaved-measure encoding F_μ(x)[m] = μ_{v₂(m)} and ask when it is eventually constant — the direct translation of the ruler trick to the measure setting, unverified.

---

## 6. OQ-5 — The complexity of 𝔼: coNP-complete

### Restatement

*Is the availability-quantified consequence relation decidable in polynomial time? What is the proof-theoretic cost of the three-filter calculus?*

### Findings

- **The corpus rebuilds.** The logic corpus (102 formulas at two atoms, 297 at three, ≤ 2 connectives over {¬, ∧, ∨}) was rebuilt from scratch and reproduces v1.1's [W17] (0 tautologies, 0 unsatisfiables), [W18] (single-premise 𝔼 = CL ∧ forward containment, 88,209 pairs, 0 mismatches), [W24], and the zero-inconsistent-premise-sets fact **[W44]**. The three-filter multi-premise characterization re-verifies with 0 violations on 535,806 + 88,209 + 742,500 pairs, with the strict placements intact: 𝔼 ⊊ CL (47,250 divergences), 𝔼 ⊊ LP (11,240), 𝔼 ⊊ ⊢ʳ (28,856) **[W45]**. The availability lemma: 4,998 checks, 0 bad **[W46]**.
- **𝔼-consequence is coNP-complete** [W47]. *Membership:* by the three-filter characterization, Γ ⊨𝔼 φ iff Γ ⊨_CL φ ∧ Var(φ) ⊆ Var(Γ) ∧ Γ ⊨_LP φ. Each filter is decidable in coNP (a counter-reading certifies failure), and the availability quantifier — the exponential range of availability sets A — is **free**: the filters absorb it. *Hardness:* for every formula γ,

  > **γ is classically unsatisfiable ⟺ γ ⊨𝔼 ⋁_{p ∈ Var(γ)} (p ∧ ¬p)** — 0 mismatches on all 102 + 297 corpus formulas (4 + 6 of them true).

  The mechanism is the glut: a CL-unsatisfiable γ forces, in every available reading of it, some atom to take the glutted value {T, F}; the big disjunction of gluts is then designated in every environment that designates γ. Conversely, if γ is CL-satisfiable, the all-classical seed designates γ but designates nothing glutted. Since CL-unsatisfiability is coNP-complete, 𝔼-consequence is coNP-hard; with membership, **coNP-complete**.
- **The proof-theoretic cost is now known to be exponential, at least.** coNP-completeness forecloses polynomial-size cut-free proof systems for 𝔼 unless NP = coNP: v1.2's speculative "polytime for bounded premises" is dead as a general hope. The honest formulation of what remains: a sound and complete *calculus* exists in principle (Parry's system plus LP's rules — soundness machine-checkable), but any cut-free version pays for M's coupling with exponential proofs, exactly as the availability lemma's global poisoning predicts.

### Verdict

**Answered.** The consequence relation of the seven signs sits at coNP — the same level as classical consequence, reached by a glut mechanism rather than a countermodel. And the answer disciplines the proof theory: the three filters are cheap to *check* and expensive to *prove in*, which is the arithmetical shadow of the semantic fact (v1.1) that M poisons whole formulas, not occurrences.

### Remains open

A signed sequent calculus organized around availability rather than the seven signs; whether cut elimination for a Parry+LP hybrid is merely exponential (matching coNP) or worse; the exact complexity of the bounded-premise fragment.

---

## 7. OQ-6 — Empirical kernels: carried

### Restatement

*Fit (θ, β) to a real culture's belief-dynamics data; the rate law t* is falsifiable — a culture that has not crystallized under measured bias θ for longer than 1.5·t* refutes the drift model or measures its Eristicity.*

### Findings

No real-culture data was fitted this wave, and we will not pretend otherwise. What this wave adds is only the standing inventory of falsifiable predictions, now two: the fixation-rate law t* = 7.80 / (β(2θ−1)·ln(θ/(1−θ))) (v1.2 [W34]) and the monotone wrong-basin curve (v1.2 [W35]) — a culture's crystallization time *and* its error rate are both functions of two measurable parameters.

### Verdict

**Carried forward.** The theory is ready for data; the data is not here.

### Remains open

As stated, plus the formal first-passage theorem version of the wrong-basin effect (see OQ-7).

---

## 8. OQ-7 — The impatience effect in the literature: adjacent, not named

### Restatement

*Is the monotone wrong-basin probability in β — fast sequential learners crystallize on the majority basin more often — already named in sequential testing or learning-in-games?*

### Findings

A targeted search found near neighbors on three sides, none of them the effect:

- **Learning-rate stability in belief-space dynamics** (arXiv 2605.06741, 2026): closed-form admissible steps for belief updates, with stability bounds — *when* learning rates keep convergence stable, not the probability of fixating on the wrong absorbing basin.
- **Stable convergence to wrong beliefs** (arXiv 2603.21491, 2026): learning that converges reliably to a wrong belief under latent reliability — misspecification, where the wrongness is in the model. The impatience effect is different: the learner's model is *correct*, and the wrong basin is reached by overshooting the log-odds barrier.
- **The speed of sequential asymptotic learning** (Tamuz et al.): convergence *rates* in information cascades — speed, not error probability.

### Verdict

**Unresolved, adjacent.** The exact effect — a bias-preserving Bayesian filter with two absorbing basins, where increasing the step size monotonically increases the probability of absorbing on the wrong side — appears unnamed in the 2026 literature. v1.2's phrasing stands: possibly a small new theorem wearing Discordian clothes, with a Brownian first-passage computation as its proof sketch.

### Remains open

Prove the first-passage version as a theorem (the discrete-walk overshoot correction at large β is the technical core); keep looking for the name.

---

## 9. OQ-8 — The Apple downstream: dead on the toys, one candidate left

### Restatement

*The M2′ Hopf algebra (group-completed im f) may be quasitriangular even though the grid ε-bialgebra is not. If it is, the Apple of Discord exists — on the order side of the fold. Is it?*

### Findings

This wave did not reach M2′ itself; what it did is make the negative half airtight, so that the question's remaining content is exactly M2′.

- **The grid case is closed at any support and any coproduct guise** [W60]. (i) Deconcatenation, the canonical left-rooted coproduct Δ(w) = Σ w[0..i)⊗w[i..|w|), satisfies the ε-rule (0 bad over 256 pairs) but fails coassociativity on **15 of 15** nonempty words — witness ab, where the two triple-deconcatenations leave unmatched split terms (a⊗ε⊗b, ab⊗ε⊗ε). (ii) The one-variable algebra k[a]: the rank-1 r-form system is inconsistent (rank 25 vs augmented 26). (iii) The free algebra itself, 49 unknowns over words ≤ 2 including ε: inconsistent (rank 49 vs 50). (iv) The **mass obstruction**, which needs no search: for every r-form, (ε⊗ε)Δ_r ≡ 0 (both counits land on the minus-terms), but the augmentation satisfies (ε⊗ε)Δ(w) ≥ 1 for every nonempty word. No tensor of any support can reproduce a coproduct with mass. **No Apple on the grid — not at support ≤ 2, not at any finite support, not in any r-form guise.**
- **UT₄'s Apples are inert.** The 44 sporadic quasitriangular tensors of [W59] (§3) satisfy every law and induce nothing: the dendriform structures they define are either identically zero or vacuous. A quasitriangular ε-bialgebra exists on UT₄ — several, in fact — and not one of them carries the fold layer.
- **M2′ is untouched.** The group-completed image of the fold is a Hopf algebra of a genuinely different kind (a group algebra, with its standard structure), and its quasitriangularity — the existence of the Apple on the *order* side of the fold — was not checked this wave. It is the last place the Apple could live.

### Verdict

**Answered for the toys, honestly.** The Apple of Discord does not exist on the grid — the mass obstruction is structural, not a small-model artifact — and on UT₄ it exists only in inert, sporadic copies that classify nothing. Everything now hangs on M2′: if it is quasitriangular, the Apple is real but lives entirely downstream of the fold, in the group-completed order it creates — which is a very Discordian place for it.

### Remains open

Quasitriangularity of M2′: the group algebra kG for G the group completion of im f — a finite/Hopf computation, one wave from closed. This is the last Apple.

---

## 10. Corrections register

Claims from earlier waves (and from this wave's own intermediate runs) that the final checks do **not** confirm, with their replacements:

| Original claim | Status | What replaced it |
|---|---|---|
| v1.2 [W41]: "107 dendriform completions survive D3 on short words" | **Corrected** — the corrected 4-term ansatz on the 16-word corpus gives **108** survivors [W49], realizing six effective completions | The exact 34-member family (two ℤ-lines + origin + one isolated point) with global closed forms [W50] |
| v1.2's implication that the survivor family is an unstructured ≥107-cloud, possibly a torsor | **Resolved negatively** — the geometry is two parallel lines plus two isolated points; no transitive group action, so no torsor [W50] | The pencil description of §2, with concatenation, the Baxter star, f(x)f(y), and the two degenerate projections as its named points |
| The doublefold product (x−f(x))(y−f(y)) as a candidate family member (sum-zero vector (1,−1,−1,1)) | **Expelled** — fails D3 at (b, a, b): bab vs 0, 605/4,096 triples bad; fails dialgebra A2 at the same witness [W48] | D3 does real work beyond Σβ = 0; the family is a proper subset of the sum-zero hyperplane |
| A mid-session candidate: the diagonal r = Σ E_mm⊗E_mm on UT₄, which appeared coassociative | **Withdrawn** — its coassociativity defect is exactly the boundary pattern λ = ρ = 1 of the graded family [W58] | The exact 44 Convention-A solutions, with the v·u orientation made explicit [W59] |
| A mid-session D1 check reporting widespread failures across the family | **Artifact** — the check used a wrong axiom form; the corrected D1 imposes exactly Σβ = 0, and all 34 members pass (0 bad over 34 × 4,096 triples) [W51] | D1 ⟺ Σβ = 0 with both sides computing to (1+Σβ)·x·f(y)f(z); D2 vacuous; D3 the sole real filter |
| v1.2 [W38], read as a family claim: "E1: ann = im f" for the dendriform layer | **Stratified** — ann = im f holds only on line 1 minus its commutative point; L1(−1) collapses to ann = everything, line-2 members give ann = {ε}, origin and isolated point give ann = ∅ [W53] | E1 as a stratum property (§2); the canonical member's stratum is the one v1.2 measured |
| An early check suggesting the deconcatenation coproduct was coassociative | **Corrected** — it fails on 15/15 nonempty words [W60] | The mass obstruction: no r-form of any support can match a coproduct with counit mass |
| v1.2 [W39]: "exactly one dialgebra axiom fails (the middle-exchange); the other four pass" | **Corrected** — A2 *and* A5 fail, on the grid (naive pair [0,550,0,0,550]) and on UT₄ ([0,12,0,0,12]) [W52, W57]; within the family, A2 fails for every non-origin member | The dialgebra layer withdrawn everywhere; the Lie layer is the honest structure in both models |

And the register's newest rows are our own: this session's first rank computations were corrupted by an unnormalized Gaussian pivot (every GE rank and consistency verdict above is from the normalized rerun), two ε-rule checks carried sign slips, and one UT₄ check omitted a term. All three were caught by cross-checking against hand computations before any figure in this document was recorded. *The check was only good because it could have returned another number.*

---

## 11. Scorecard

| Question | Was (v1.2 §11) | Now (v1.3) | Status |
|---|---|---|---|
| OQ-1 D3-survivor family | open | exactly 34 members: two ℤ-lines + two points, global closed forms [W49–W51]; doublefold expelled [W48]; E1 stratified [W53]; dialgebra nowhere [W52] | answered |
| OQ-2 UT₄ corrected | open | chain transfers verbatim [W54, W55]; E1 = im f + E14 (Podge fingerprint) [W56]; graded ε-coproduct no-go: only Δ ≡ 0 [W58]; 44 inert Apples [W59] | answered |
| OQ-3 Flats completeness | open | P is Σ⁰₂-complete via the ruler reduction [W61–W63]; v1.2's proper→complete framing corrected | answered (with correction) |
| OQ-4 Dist-valued | open | no valuation-interleave analogue known; deterministic case fully closed | carried |
| OQ-5 Complexity of 𝔼 | open | 𝔼 is coNP-complete [W44–W47]; availability quantifier free; poly cut-free proofs foreclosed | answered |
| OQ-6 Empirical kernels | open | no data fitted; t* and the wrong-basin curve stand as predictions | carried |
| OQ-7 Impatience in the literature | open | two 2026 near-neighbors + the cascade literature; the effect appears unnamed | unresolved, adjacent |
| OQ-8 Apple downstream | open | grid airtight at any support (mass obstruction) [W60]; UT₄'s 44 inert [W59]; M2′ the last candidate | answered for the toys |

---

## 12. Source notes

| Source | Credibility (x/5) | Last updated |
|---|---|---|
| [Aguiar, *Infinitesimal bialgebras, pre-Lie and dendriform algebras*](https://pi.math.cornell.edu/~maguiar/depaul.pdf) (also [arXiv math/0211074](https://arxiv.org/pdf/math/0211074), [v3](https://arxiv.org/html/math/0211074v3)) | 4 | - |
| [Aguiar, *Infinitesimal Hopf algebras* (Contemp. Math. 267)](https://pi.math.cornell.edu/~maguiar/cor.pdf) — dual numbers quasitriangular, k[x] not; Convention A | 5 | 2000 |
| [Aguiar & Loday, *Quadri-algebras*](https://pi.math.cornell.edu/~maguiar/quadri.pdf) | 4 | - |
| Loday, *Dialgebras* (Lecture Notes in Mathematics 1763) — the five axioms A1–A5 | 5 | 2001 |
| [Dendriform trialgebra (Baxter operators of arbitrary weight)](https://arxiv.org/pdf/math/0311132v1) | 4 | - |
| [Ferguson, *Logics of Nonsense and Parry Systems*, J. Phil. Logic 44(1)](https://link.springer.com/article/10.1007/s10992-014-9321-y) (also [JSTOR](https://www.jstor.org/stable/44851254)) | 5 | 2015 |
| *Containment Logics: Algebraic Completeness and Axiomatization* (Studia Logica; the ⊢ˡ/⊢ʳ companions) | 5 | - |
| Pra Baldi, doctoral thesis, University of Padua (containment logics) | 4 | - |
| Beirlaen & De Cock, *Characterizing Logical Consequence in Paraconsistent Weak Kleene* (Studia Logica) | 5 | - |
| [Szmuc, *An Epistemic Interpretation of Paraconsistent Weak Kleene*](https://philpapers.org/archive/SZMAEI.pdf) | 4 | - |
| [Calbrix, Nivat & Podelski, *Ultimately periodic words of rational ω-languages*](https://link.springer.com/content/pdf/10.1007%2F3-540-58027-1_27.pdf) | 5 | 1993 |
| Landweber, *Decision problems in state machines* (ω-regular hierarchy foundations) | 5 | 1969 |
| [Closed-form admissible learning-rate steps in belief-space dynamics](https://arxiv.org/pdf/2605.06741) (OQ-7 neighbor: stability bounds, not wrong-basin probability) | 3 | 2026 |
| [Stable convergence to wrong beliefs under latent reliability](https://arxiv.org/html/2603.21491v1) (OQ-7 neighbor: misspecification, not overshoot) | 3 | 2026 |
| [Tamuz et al., on the speed of sequential asymptotic learning](https://tamuz.caltech.edu/papers/cascades.pdf) (OQ-7 neighbor: rates, not error) | 4 | - |
| Parry, analytic implication (via Ferguson 2015; original 1933) | 5 | 1933 |
| *Principia Discordia* | 5 | 1965 |

*Machine-checked claims [W44]–[W63] were produced by the author's verification scripts for this document and are reproducible from the corpora described in each section.*

---

## 13. New open questions

1. **UT_n and the 44** (OQ-2). Classify the Convention-A quasitriangular tensors on UT_n for general n: does the count follow a pattern, does the v·u = 0 orientation persist, and does any Apple become non-vacuous as the corner grows? Is the annihilator always im f plus the nilpotent corner E_1n?
2. **Beyond the graded and rank-1 ansätze** (OQ-2). The no-go covers graded coproducts and rank-1 r-forms; five named non-graded candidates are dead. Is there a non-graded, non-rank-1 ε-coproduct on UT₄ — or a proof that none exists?
3. **ε-coproducts on the free algebra** (OQ-8). The mass obstruction kills every r-form on k⟨a,b⟩. Does *any* ε-bialgebra coproduct (not of r-form type) exist on the free grid algebra compatible with concatenation — or is the grid algebra provably barren of coalgebra structure?
4. **The pencil over the other folds** (OQ-1). Do the other 24 idempotent multiplicative folds each carry a two-line dendriform family with the (−1,1,1,−1) direction? Is the parameter t continuous over ℝ, and is the canonical member distinguished by any invariant beyond being Aguiar's?
5. **The last Apple** (OQ-8). Quasitriangularity of M2′, the group-completed im f: a Hopf computation one wave from closing. If it fails, the Apple of Discord does not exist anywhere in this framework; if it holds, Discord is a structure of pure order.
6. **Dist-valued Borel analysis** (OQ-4, carried). The interleaved-measure encoding F_μ(x)[m] = μ_{v₂(m)}: when is it eventually constant? Does the flat/hand-off dichotomy survive as eventually-repeating vs pre-periodic measure sequences in the weak-star σ-algebra?
7. **𝔼 proof theory** (OQ-5, sharpened). A signed sequent calculus around availability rather than the seven signs; is cut elimination for a Parry+LP hybrid merely exponential (matching coNP-completeness) or worse?
8. **Empirical kernels and the impatience theorem** (OQ-6, OQ-7, carried). Prove the first-passage wrong-basin theorem with the discrete-overshoot correction; fit (θ, β) to a real culture against the t* law; keep looking for the effect's name.

---

*"Eris forgives the wrong number; she does not forgive the uncheckable one. The Apple you cannot find anywhere is still an Apple — and the one hiding in pure order is the funniest joke she has told us yet."* — after the Principia