# Discordianism: Open Questions Explored (v2.2)

*Answers the open questions of `discordianism-open-questions-explored-v2.1.md` §8.*

*Tags **[Un]** are machine checks in `verify_v2_2.py` (10/10 pass). Tags [Tn] refer to `verify_v2_1.py`. Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Summary

| # | Question (v2.1 §8) | Answer | Kind |
|---|---|---|---|
| 1 | Non-finite-state sources | **The Borel class never climbs.** Every countable grid language names a Σ⁰₂-complete set of measures. Expressivity reappears in the *equivalence problem*: polynomial time for finite-state sources **[U5]**, equivalent to arithmetic-circuit identity testing for visibly pushdown sources, and Π⁰₁-complete for computable sources | theorem |
| 2 | Folds that are not endomorphisms | **Bigraded case classified: four operators [U1].** Besides 0, id and the identification fold, there is a new **framed fold** f□, which is not multiplicative **[U2, U3]**. On UT_n the monomial deconstructing folds are exactly the 2ⁿ **vertex supports** **[U4]** | theorem + classification |
| 3 | A label-free calculus; interpolation | **No-go.** No leaf-local axiom set is complete with the invertible LP rules **[U6]**. Craig and uniform interpolation hold, trivially, because of containment **[U7]** | theorem |
| 4 | Data | **Not fitted here; there is a published anchor.** The v2.1 culture model is exactly a random-field threshold model **[U10]**. Michard & Bouchaud fitted that class to birth rates, mobile-phone adoption and applause endings | literature |
| 5 | The noisy bureaucracy | **Lock-in survives noise; the 2-cycle does not.** The efficient policy's second eigenvalue tends to 1 (metastable lock-in) **[U8]**. The fair policy's most negative eigenvalue tends to −½ (anti-persistence, no period) **[U9]**. v2.1's 2-cycle came from tie-breaking by name | theorem (exact chain) + correction |

---

## 1. Expressive sources: the Borel class is blind to expressivity

**Theorem 1.** Let a grid language have countably many grids, and let it name at least the Dirac laws of all flats. Then the set N of named measures is Σ⁰₂-complete in P(𝒞), however expressive the language is: pushdown, counters, Turing-complete.

*Proof.*
1. N is countable, hence Σ⁰₂.
2. The set W of words whose Dirac law is named is countable and contains the flats, so it is dense in 2^ω and not G_δ. By Wadge's lemma, W is Σ⁰₂-complete in the zero-dimensional space 2^ω, via some reduction g.
3. Since δ_w ∈ N ⟺ w ∈ W, the map δ∘g reduces FIN to N.

v2.1 Theorem 1 was the finite-state instance. Its pumping lemma was needed only to identify W, not to bound the class.

**Where expressivity does show up: the equivalence problem.** Asking whether two grids name the same chaos is the Mu-lossiness question, i.e. N2's bisimulation equality.

| sources | deciding ⟦g⟧ = ⟦h⟧ |
|---|---|
| finite-state (n₁, n₂ states) | decidable from prefixes of length ≤ n₁+n₂−1. Checked on 5,349 inequivalent random pairs; every one disagrees within that bound **[U5]**. Polynomial time (Tzeng); randomized NC (Kiefer et al.) |
| probabilistic visibly pushdown | logspace-equivalent to arithmetic circuit identity testing (Kiefer, Murawski, Ouaknine, Wachter & Worrell). Membership in deterministic P is open |
| general probabilistic pushdown | studied by Forejt, Jančar, Kiefer & Worrell (Inf. Comput. 2014); I did not confirm its exact status |
| computable (Turing-complete `Loop`) | **Π⁰₁-complete** (proof below) |

*Proof of the last row.*
- *Membership:* two computable laws are equal iff every cylinder probability agrees, and equality of computable reals is Π⁰₁.
- *Hardness:* let grid g_e emit 0^ω, except a 1 at position t if machine e halts at step t. Then ⟦g_e⟧ = δ_{0^ω} ⟺ e never halts.

**Gloss.** *A Pope can always write a new grid; whether two Popes wrote the same one becomes undecidable once grids can loop.*

**Real parameters with expressive loops.** The map from parameters to laws can be discontinuous: termination probabilities jump. The image is then only analytic (Σ¹₁). Whether N is Borel in that case is **open**.

---

## 2. Folds beyond endomorphisms

**Setting.** Idempotent RB(−1) operators are exactly the splittings A = A₁ ⊕ A₂ into two subalgebras, where P projects onto A₁ along A₂. Commuting with letter deletion means:
- Δ(A₁) ⊆ A₁⊗A₁;
- A₂ is a coideal: Δ(A₂) ⊆ A₂⊗A + A⊗A₂.

Work in the letters a and u = a − b, so that Δa = 1⊗1, Δu = 0, and Δ deletes one a.

**Reduction lemma.** For any such P:
- Δ(P uʲ) = (P⊗P)Δ(uʲ) = 0, so **P preserves k[u]**.
- Δ is injective off k[u], so P is determined by two pieces of data: its restriction to k[u], and a k[u]-valued correction on each word.

**Theorem 2 (bigraded classification).** Call P bigraded if it preserves a-degree and u-degree. Exactly four bigraded deconstructing folds exist **[U1]**:

| P(1) | P on k[u] | fold | Hodge A₁ | Podge A₂ |
|---|---|---|---|---|
| 1 | id | identity | A | 0 |
| 0 | 0 | zero | 0 | A |
| 1 | 0 | **identification** f₁ (b ↦ a) | k[a] | (u) = (a − b) |
| 0 | id | **framed** f□ | span u⁺(a u⁺)\* | everything else, including **Void** |

*Proof.*
1. On each 1-dimensional block, P is 0 or the identity.
2. {j : uʲ ∈ A₁} and its complement in ℕ⁺ are both closed under addition, so one of them is empty. Hence P is all-or-nothing on k[u].
3. Every other block is then forced by the recursion Δ(Pw) = Σ P(left) ⊗ P(right).

The four cases all exist. Idempotence and the Rota–Baxter law were checked on words of length ≤ 5 **[U1]**.

**The framed fold f□** is new **[U2, U3]**:
- **Hodge** is spanned by the words in which every a is isolated and flanked by u on both sides: u, uau, uuauu, …
- **Podge** is spanned by the rest: words that begin or end in a, or contain aa, and the empty word.
- Both are subalgebras. Hodge is closed under Δ; Podge is a coideal (checked to length 7).
- f□ is **not multiplicative**: f□(uau) = uau, but f□(u)f□(a)f□(u) = 0.
- f□(1) = 0: **the Void is chaos.**

f□ supports the dendriform half (δ-D): its star is still pq − nm. It does **not** support the dialgebra half (δ-L), which needs a multiplicative fold.

**Gloss.** *f□'s order consists of events framed by difference on both sides. Recall that u = a − b is exactly what deconstruction cannot see (ker Δ = k[u]). So order is the visible, isolated event held in place by invisible differences.*

**UT_n with arrow deletion.** The monomial deconstructing splittings are exactly the 2ⁿ **vertex supports** (verified for n = 2..5) **[U4]**:
- Hodge = paths that never leave a vertex set V′;
- Podge = paths that visit some vertex outside V′.

These are non-unital algebra endomorphisms: quotient by the outside vertices, then include back. *Order is territory: the paths that never leave home.*

**Remains open.**
- *Non-bigraded folds on k⟨a,b⟩.* By the reduction lemma, these hinge on the subalgebra splittings of k[u]. Conjecture: the only splittings are 0 ⊕ k[u], k[u] ⊕ 0, k ⊕ (u−c) and (u−c) ⊕ k. If true, every deconstructing fold is a k[u]-twisted deformation of the four above. (The endomorphism family of v2.1 Theorem 2 is the k ⊕ (u−c) branch.)
- *Non-monomial folds on UT_n.*

---

## 3. The 𝔼 calculus: why the label is necessary, and interpolation

**Theorem 3 (locality no-go).** No set of leaf axioms makes the invertible LP rules complete for 𝔼 **[U6]**.

*Proof.* The valid sequent p∨q ⇒ p∨q has leaves {p ⇒ p,q ; q ⇒ p,q}. The invalid p ⇒ p∨q has the leaf {p ⇒ p,q}, which is a subset. Any leaf predicate that accepts the first proof accepts the second.

So every complete calculus for 𝔼 must do one of the following:
- **thread state**: G_𝔼's availability label, i.e. framework response (b);
- **use non-invertible rules with side conditions**, the route of Szmuc's GPAI_fde for Parry, which restructures the global invariant;
- **carry ghost formulas** recording the variables of discarded disjuncts, which is the label distributed into the sequent.

*Containment is irreducibly a global L1 invariant of proofs.*

**Interpolation.**
- *Craig:* if Γ ⊨_𝔼 φ, then Var(φ) ⊆ Var(Γ), so φ is its own interpolant (867 valid random instances) **[U7]**. Every containment logic interpolates trivially.
- *Uniform interpolation:* also holds. For V ⊆ Var(γ), let ι be the conjunction of γ's consequences over V. This is a finite conjunction, because the 4-valued matrix has finitely many designation sets over V. Then ∧-introduction, ∧-elimination and transitivity give γ ⊨ φ ⟺ ι ⊨ φ for every φ with Var(φ) ∩ Var(γ) ⊆ V.

*The interesting proof theory of 𝔼 is the locality question, not interpolation.*

---

## 4. Data: the empirical anchor

The v2.1 Glauber culture, P(adopt) = ½(1 + tanh(β(Jm + h))), is **exactly** a zero-temperature random-field threshold model. Agent i adopts iff φᵢ + Jm + h > 0, with propensities φᵢ ~ Logistic(0, 1/(2β)) **[U10]**. Hence the critical point K = βJ = 1 is the statement that **J_c scales with the spread of propensities**.

Michard & Bouchaud (EPJ B 47, 2005) fitted exactly this model class to three datasets:

| dataset | regime found |
|---|---|
| European birth rates, 1950–2000 | smooth: J < J_c |
| Mobile-phone adoption, 1994–2003 | smooth: J < J_c |
| Applause endings (17 recordings) | mostly persistent, but some end discontinuously, read as J > J_c |

The peak-slope exponents (≈ −0.6 to −0.7) are close to the mean-field near-critical value −2/3. The authors call the agreement noisy.

**Reading for C5.** Real cultures appear to sit **near but mostly below** K_c: strong social amplification without lock-in. Small, tightly coupled groups (an audience) can cross into the regime where wrong consensus is metastable.

This is a published fit of the same model class, not a fit performed here. Fitting (K, H) to a belief-dynamics dataset and testing the spinodal H_s(K) remains open.

---

## 5. The noisy bureaucracy

**Model.** As in v2.1 Q6, but restructuring sees only the empirical counts of a window of M forms, and ties are broken uniformly. The transition matrix over the six structures is computed exactly from the multinomial.

**Efficient policy** (fast lane to the most-filed form) **[U8]**:

| M | 2 | 5 | 10 | 20 | 40 |
|---|---|---|---|---|---|
| λ₂ | 0.40 | 0.65 | 0.83 | 0.954 | 0.996 |

Lock-in becomes **metastable switching** between the three self-fulfilling orders. The residence time ~1/(1−λ₂) grows with the window, which is the bureaucratic analogue of v2.1's Kramers law.

**Fair policy** (fast lane to the least-filed form) **[U9]**:

| M | 2 | 5 | 10 | 20 | 40 |
|---|---|---|---|---|---|
| λ_min | −0.32 | −0.41 | −0.47 | −0.495 | −0.4998 |

λ_min tends to **−½, not −1**. The probability that the fast form persists falls to ~10⁻⁴. With fair tie-breaking, the fast lane leaves its form and jumps *uniformly* to one of the other two. That is a random walk on the complete graph K₃, with spectrum {1, −½, −½}: perpetual reorganization **without a period**.

**Correction to v2.1 [T16].** Its "every orbit is a 2-cycle" is true only for its deterministic, name-based tie-break. The robust statement: *the fair diagonal has no fixed point, and the fast lane is anti-persistent (eigenvalue −½)*. The efficient policy's lock-in, by contrast, is robust: eigenvalue → 1.

**Gloss.** *Aneristic order is an eigenvalue near +1; Eristic reorganization is an eigenvalue near −½. Neither is visible in any single office.*

---

## 6. Proposed amendments (v2.1 → v2.2)

| amendment | content |
|---|---|
| **B1⁺⁺ Expressivity** | Every countable grid language names a Σ⁰₂-complete set of measures. Equality of grids is P / ACIT-equivalent / Π⁰₁-complete for finite-state / visibly pushdown / computable sources. Mu-lossiness is the axis that expressivity moves |
| **M3″ Fold catalogue** | Bigraded deconstructing folds: {0, id, f₁, f□}. On path algebras: vertex supports. f□ puts Void in Podge |
| **L6′ Locality** | Completeness of 𝔼 forbids leaf-local axioms with invertible rules. Craig and uniform interpolation are trivial |
| **C5′ Calibration** | The culture model is a logistic random-field threshold model. The empirical regime (Michard–Bouchaud) is near-critical, mostly below K_c |
| **Γ1′ Robust Greyface** | Under noise, lock-in survives (λ₂ → 1) and fair reorganization becomes anti-persistent (λ → −½), not periodic |

---

## 7. New open questions

1. **Splittings of k[u].** Are the four listed the only subalgebra splittings of k[u]? That would complete the classification of deconstructing folds on k⟨a,b⟩.
2. **Non-monomial folds on UT_n**, and deconstructing folds on general path algebras: is it always vertex supports up to twisting?
3. **Borel or not.** For real-parameter grids with loops, is the set of named measures Borel, or properly analytic?
4. **Fit.** Find a belief-dynamics dataset with measurable coupling and evidence rate, and test H_s(K).
5. **Larger bureaucracies.** For k forms, the fair policy's walk on K_k has eigenvalue −1/(k−1), so anti-persistence weakens with size. Does any policy produce a genuine period (an eigenvalue → −1) robustly under noise?

---

## 8. Sources

- Michard & Bouchaud, *Theory of collective opinion shifts: from smooth trends to abrupt swings*, Eur. Phys. J. B 47 (2005): https://arxiv.org/pdf/cond-mat/0504079
- Kiefer, Murawski, Ouaknine, Wachter & Worrell, *On the Complexity of the Equivalence Problem for Probabilistic Automata* (2011): https://arxiv.org/abs/1112.4644
- Forejt, Jančar, Kiefer & Worrell, *Language equivalence of probabilistic pushdown automata*, Inf. Comput. (2014): https://dl.acm.org/doi/10.1016/j.ic.2014.04.003
- Tzeng, *A polynomial-time algorithm for the equivalence of probabilistic automata*, SIAM J. Comput. 21 (1992).
- Szmuc, *A simple logical matrix and sequent calculus for Parry's logic of Analytic Implication*: https://philarchive.org/archive/SZMASL
- Kechris, *Classical Descriptive Set Theory* (1995), Thm 22.10.

*"The fair bureaucrat never leaves the fast lane where it was, and never knows where it will go; Eris calls this progress."*
