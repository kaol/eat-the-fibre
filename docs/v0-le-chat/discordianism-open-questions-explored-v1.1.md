# Discordianism: Open Questions Explored (v1.1)

*Companion to the v1 Mu–Nu axioms document. Every formal claim tagged **[W#]** was machine-checked this session (TypeScript sandbox: exhaustive enumeration over finite witnesses and direct evaluation of the 𝔼 matrices). Claims tagged **[V#]** are unchanged from v1's `verify_discordian.py` (Z3, SymPy, NumPy; 18/18). Untagged formal claims are standard results. Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Method

Each of v1's five open questions (§10) is answered in the same format: a precise restatement, findings with [W#] tags, a verdict against what the question conjectured, and what remains genuinely open. Two kinds of answer appear:

- **Witness answers** (Q1, Q2, Q3): finite toy models, exhaustively enumerated. They show a mechanism exists; they do not prove it is the only mechanism.
- **Structural answers** (Q4, Q5): theorems plus exhaustive checks, which pin down the shape of the result.

---

## 1. Executive summary

| # | v1 asked | v1.1 answer | verdict |
|---|---|---|---|
| Q1 | Does the Eristic/Aneristic ratio flow while the grid set stays fixed? | **Yes.** The ratio is a measure, not a structure; the E/A distinction is channel bias, confirming §0's thesis that Hodge and Podge are both folds | answered (toy model) |
| Q2 | Do social dynamics have Arrow fibers? | **Yes, witness constructed.** β3 (counter-dependent observation) is exactly the obstruction to global Arrow; failure occurs only at institution transitions | answered (witness) |
| Q3 | Is grid-definable Chaos exactly a Borel level? | **Yes: exactly Σ⁰₂-proper**, given β1's density refinement. No internal Wagner-style hierarchy; ⟦𝒢⟧ is a non-ω-regular set of ultimately periodic words | answered |
| Q4 | Does the δ-L bracket arise from an Aguiar infinitesimal bialgebra on the grid algebra? | **Yes, in precise form**: grid ε-bialgebra + Rota–Baxter operator of weight −1 ⟹ dendriform ⟹ Loday dialgebra, with M2′ (Hodge Hopf) as the group completion of im f. The Sacred Chao, algebraically | answered; converses open |
| Q5 | Is 𝔼 axiomatizable, and is it PWK on the M-free fragment? | **The conjecture is precisely inverted.** M-free fragment = LP (Priest); the PWK-flavor comes from M and points forward, the dual of PWK's direction. 𝔼 is incomparable with PWK and B₃; finitely axiomatizable by signed sequent calculus, not by a non-vacuous Hilbert system | answered; multi-premise characterization open |

---

## 2. Q1 — Measure: does the ratio flow on a fixed grid?

### Restatement

*Which leaf-data measure on 𝒢 does a culture realize, and does the Eristic/Aneristic ratio flow (measure dynamics) while the grid set is fixed?*

### Setup

A fixed universe of 256 eight-bit "grids". The fold: a word with minimal period ≤ 3 is **Aneristic** (10 words), otherwise **Eristic** (246). A culture is a probability kernel on this fixed set, updated by Bayesian multiplicative updates (learning rate β = 0.8, 15 steps) under four observation channels: neutral, Aneristic-biased, Eristic-biased, and five-rewarding (a word counts only if it has exactly five set bits — the five-fingered hand as a bit pattern; there are C(8,5) = 56 such words at baseline).

### Findings

- **Neutral channel: nothing flows.** P(Aneristic) stays at 0.0391 and entropy at 5.545 bits (= log₂ 256), constant across all 15 steps **[W1]**. This is B5's "a five counts only if it could have returned another number", exhibited as a channel property: an unbiased observer finds no pattern.
- **Aneristic-biased channel: crystallization.** P(Aneristic) → 1.0 and entropy → 2.305 bits = log₂ 10 **[W2]**: the kernel collapses onto exactly the 10 order-words. The **Aneristic Illusion is this collapse**: reading the posterior as a fact about 𝒞 instead of about the channel.
- **Eristic-biased channel: dispersion.** P(Aneristic) → 0 and entropy → 5.505 bits = log₂ 246 **[W3]**. Symmetric collapse in the opposite direction — the Eristic Illusion is the same defect with the opposite fold label.
- **Five-rewarding channel: the Law of Fives as channel bias.** The frequency of popcount-5 words → 1.0 and entropy → 4.026 bits = log₂ 56 **[W4]**. "The more you look, the more fives you find" is a fixed point of a five-selective channel, not a property of the world.

### Verdict

The ratio flows while the grid set is fixed: **the Eristic/Aneristic distinction is carried entirely by the measure, not by the structure**. This confirms v1 §0 (Hodge and Podge are both Mu, differing only by fold label) with a dynamics witness. Greyface is the error of freezing the kernel at its crystallized limit and calling it the territory.

### Remains open

Empirical kernel identification for real cultures (which channels is a culture actually running?); continuous leaf-data measures on 𝒢 rather than bit patterns; the phase structure of crystallization (the P(Aneristic) trajectory suggests a sigmoidal collapse — is there a sharp threshold in β, i.e. a genuine phase transition between Eristic and Aneristic regimes?).

---

## 3. Q2 — Real systems: do Arrow fibers exist?

### Restatement

*Social dynamics fail the Arrow test globally. Do they have Arrow fibers (fix institutions → behaviour is Arrow; fix behaviour → institutions are Arrow)? If so, the framework applies per fiber, and Greyface is the error of ignoring the diagonal.*

### Setup

A minimal social witness. State = (counter c ∈ ℤ/4, institution i ∈ {0,1}, behaviour b ∈ {0,1}). Three endomaps:

- **F** (routine): increment the counter.
- **H** (institutional change): flip the institution *iff* the counter reads 3 — the counter-dependence of β3, the bureaucrat whose observation depends on how many forms were filed.
- **G** (fibre behaviour): b ↦ b ⊕ i ⊕ 1.

### Findings

- **Within a fiber, Arrow holds.** Fixing the institution, all permutations of {F, G} give the same composite **[W5]**: with institutions frozen, behaviour is Arrow.
- **Globally, Arrow fails, exactly at the seams.** Permutations of {F, H, G} diverge: F∘H∘G and F∘G∘H reach (3, 1, 0) and (3, 1, 1) respectively **[W6]**. Every divergence is localized at an institution transition, where the fixed-institution reasoning no longer applies — these transition points are precisely the λ-defects of §4.

### Verdict

Arrow fibers exist in a concrete witness, and the obstruction to global Arrow is exactly β3 (counter-dependent observation). The diagnosis is sharpened: **Greyface is the error of promoting fiber-wise Arrow commutation to the diagonal** — of assuming that because behaviour is Arrow at fixed institutions and institutions behave at fixed behaviour, the whole social map commutes.

### Remains open

Larger witnesses (does fiber-Arrow genericity survive noise, i.e. N4's Dist-valued steps?); the converse — is every β3-defect decomposable into fibers; classification of the transition λ-defects by their B4/β2 regime.

---

## 4. Q3 — Topology: is grid-definable Chaos exactly a Borel level?

### Restatement

*𝒞 with the prefix topology is a Stone space. Which subsets are ⟦−⟧-definable? Compare ω-regular languages (Büchi): is "grid-definable Chaos" exactly a Borel level?*

### Setup

A toy GridF over binary percepts with the flat denotation: each grid denotes an ultimately periodic word (a flat). All grids up to size 9 were enumerated — 199,827 grids — and their denotations checked against depth-8 prefix approximations of 𝒞 = 2^ω.

### Findings

- **Coverage is exact and periodic.** The depth-8 prefixes covered by ⟦𝒢⟧ are exactly the words of minimal period ≤ 5: 52 covered prefixes against 52 words of minimal period ≤ 5, an exact match **[W7]**. (This "5" is an artifact of the size-9 cutoff — B5 applied to itself — not a Law of Fives.)
- **Minimal grid size for observing a flat u is 2·minper(u) − 1**, proved as a bound and verified on every word **[W8]**.
- **Mu lossiness, quantified.** 199,827 grids denote only **63** distinct chaoses. The flat "0" is named by 10,249 different grids; the worst case, "00", by 14,072 **[W9]**. The non-injectivity of ⟦−⟧ (§1 of v1) is not marginal — the average grid shares its chaos with thousands of rivals.
- **Period statistics.** The minimal-period distribution over the 256 depth-8 words is {1: 2, 2: 2, 3: 6, 4: 12, 5: 30, 6: 48, 7: 82, 8: 74} **[W10]** — most words are aperiodic at depth 8, consistent with N5's inexhaustibility: almost every chaos has no small grid.

### Theorems (standard)

- ⟦𝒢⟧ is countable, hence **F_σ = Σ⁰₂**.
- With β1's hand-off interface, the grid points are **dense** in 𝒞; a dense countable set is not G_δ (Baire), so ⟦𝒢⟧ is **exactly Σ⁰₂-proper**.
- The σ-algebra generated by the point images is flat: σ(⟦g⟧ : g ∈ 𝒢) = {S, S ∪ (𝒞 ∖ ⟦𝒢⟧) : S a countable or co-countable subset of ⟦𝒢⟧}. Within the definable sets there is **no Wagner-style internal hierarchy** — contrast the Wadge–Wagner hierarchy of ω-regular languages (Landweber; Selivanov; Carton–Perrin), which lives at these same low Borel levels but with rich internal structure. Grid-definability collapses it.
- **⟦𝒢⟧ is not ω-regular.** Every singleton flat is ω-regular, but ⟦𝒢⟧ is countable and infinite, whereas every infinite ω-regular language is uncountable. By Calbrix–Nivat–Podelski, an ω-regular language is determined by its ultimately periodic words; ⟦𝒢⟧ is a set of ultimately periodic words that is not ω-regular — so the theorem's converse fails sharply, and ⟦𝒢⟧ is a canonical witness.
- **B1 is robust.** ⟦𝒢⟧ is null under *every* product measure on 𝒞 in which no single percept has probability 1 — order covers a null set, whichever coin you flip.

### Verdict

**Yes: grid-definable Chaos is exactly Σ⁰₂-proper** (given β1 density), with a trivial internal σ-structure. The Büchi comparison resolves cleanly: grid-definable ⊂ ultimately-periodic-named, incomparable with ω-regular, and strictly coarser than both.

### Remains open

Dist-valued denotations (N4) — the toy model used deterministic flats; what Borel complexity does a `Dist`-valued ⟦−⟧ reach? Bounded-variability words (the realistic image of β1's hand-off). Effectivity: is ⟦−⟧ computable, and does ⟦𝒢⟧ become Σ⁰₂-*complete* (not merely Σ⁰₂) in the full, non-toy GridF?

---

## 5. Q4 — Unification: does δ-L arise from an infinitesimal bialgebra?

### Restatement

*Aguiar's infinitesimal bialgebras satisfy a Leibniz-type rule Δ(ab) = a₁⊗a₂b + ab₁⊗b₂. Does the δ-L bracket arise from such a structure on the grid algebra, unifying M2′ (Hodge Hopf) and δ-L?*

This was the highest-risk question, and it closes completely.

### Findings

- **The free grid algebra is an infinitesimal bialgebra.** On k⟨a,b⟩ with the left-rooted deconcatenation coproduct, Aguiar's rule Δ(ab) = a₁⊗a₂b + ab₁⊗b₂ holds on all 961 word-pairs (words up to length 4) **[W11]**.
- **The hand-off is mandatory.** The symmetric strict version of the rule fails, and it fails by exactly the term u⊗v **[W12]** — the same hand-off interface that β1 demands between sequential grid components. The ε-bialgebra law is not decoration: **β1 and Aguiar's rule are the same defect, read from two sides.**
- **A second witness.** UT₄ (upper triangular matrices, the v1 model) with the gap coproduct is also an infinitesimal bialgebra **[W13]**: the unification is not an artifact of the free model.
- **The Aneristic collapse is a Rota–Baxter operator.** The idempotent multiplicative f (all-a words ↦ a^len, else 0; the diagonal on UTₙ) is a **Rota–Baxter operator of weight −1** and an infinitesimal-bialgebra morphism, in both models **[W14]**.
- **Loday's axioms follow.** With x ⊣ y = x·f(y) and x ⊢ y = f(x)·y, all five δ-L axioms hold on the free model **[W15]** (re-deriving v1's [V11] witness inside the unified structure).
- **E1 survives unification.** The left annihilator of the bracket [x, y] = x f(y) − f(y) x is exactly im f = span{ε, a, aa, aaa, aaaa}; the right annihilator is kε ⊕ ker f **[W16]** — v1's [V13], now verified in the free model, where it is furthest from trivial.

### The unification chain

```
grid algebra: ε-bialgebra + Rota–Baxter f (weight −1)
        │  Aguiar: Baxter ⟹ dendriform
        ▼
   dendriform algebra
        │  Loday: dendriform ⟹ dialgebra
        ▼
   Loday dialgebra (δ-L), right-Leibniz bracket [x,y] = xf(y) − f(y)x
        │  im f, group-completed
        ▼
   kG Hopf algebra (M2′, the Hodge Hopf quotient)
```

M2′ and δ-L are two faces of one object: an infinitesimal bialgebra with a Baxter projection. The grid side and the collision algebra are not parallel axiom groups; the second is manufactured from the first. *The Sacred Chao is the pair (ε-bialgebra, Baxter operator), with the denotation map between.*

**Bonus on the Aftermath.** The free algebra is a domain, so its Podge (ker f, the non-Hodge ideal) is **never nilpotent**: chaos composed with itself never vanishes. In UT₄ the Podge *is* nilpotent (v1's Aftermath). Both satisfy E1, so E1 does not decide the Aftermath: **the Aftermath is a model choice, not a theorem** — some cosmologies compose to zero, some do not, and the axioms are silent.

### Verdict

**Yes, in precise form.** The δ-L bracket arises from an infinitesimal bialgebra structure on the grid algebra via a weight −1 Baxter operator, and this same operator's image group-completes to the M2′ Hopf algebra. Q4's conjecture is confirmed with both of Aguiar's routes (Baxter ⟹ dendriform) on the table.

### Remains open

The converse: which infinitesimal bialgebras admit an idempotent Baxter(−1) morphism, and is the Loday dialgebra uniquely recoverable from the (ε-bialgebra, Baxter) pair? The quasitriangular refinement: Aguiar's second route (quasitriangular ε-bialgebras ⟹ dendriform) — is the grid coproduct quasitriangular, i.e. is there an r-matrix on the grid algebra, and does it classify the β2 fused nodes?

---

## 6. Q5 — 𝔼: axiomatizability, and the true face of the seven values

### Restatement

*Is the lifted-Bochvar 𝔼 (seven values: nonempty subsets of {T, F, M}; M infectious; designated = contains T) axiomatizable by a finite Hilbert system, and is its consequence relation that of a known paraconsistent logic — v1 conjectured PWK on the M-free fragment?*

### Setup

Exhaustive corpora: all formulas with ≤ 2 connectives over 2 atoms (102 formulas) and over 3 atoms (297); all 88,209 single-premise pairs; all 7,956 multi-premise sets (up to 2 premises).

### Findings

- **Theorem set is empty, and so is the counter-theorem set.** There are no tautologies (the all-{M} seed kills every formula) and dually no unsatisfiable formulas (every formula is designated at the all-{T,F} seed) **[W17]**. A "finite Hilbert axiomatization of the theorems" is therefore *vacuous* — the theorems are ∅. The honest question is a finite rule-based calculus for consequence, which exists (below).
- **Single-premise consequence is exactly characterized.** Verified on all 88,209 pairs **[W18]**:

  > γ ⊨_𝔼 φ ⟺ **Var(φ) ⊆ Var(γ)** and γ ⊨_CL φ.

  A classical core with **forward** variable-sharing (the conclusion's variables must occur in the premise). There is no vacuous branch: no premise is ever rejected wholesale. This is the formal dual of PWK's **reverse** variable inclusion (Beirlaen–De Cock's characterization: PWK single-premise consequence holds iff φ is a classical tautology, or γ ⊨_CL φ with Var(γ) ⊆ Var(φ)) **[W19]** — reproduced exactly by the check.
- **The M-free fragment is LP, not PWK.** Restricted to the values {T}, {F}, {T,F}, 𝔼's consequence coincides exactly with Priest's LP **[W20]**: gluts behave pair-independently. **v1's conjecture is precisely inverted**: the PWK-flavor comes from M, and its variable-blocking points forward — the dual direction of PWK — while gluts alone give plain LP.
- **𝔼 is its own logic, incomparable with both neighbours.** 𝔼 vs PWK: the separator is p ⊨ q∨¬q (PWK's implosion, true; 𝔼, false); 13,668 of 88,209 single-premise pairs diverge, 1,240 of 7,956 multi-premise **[W21]**. 𝔼 vs B₃ (Bochvar's own three-valued logic): 1,668 single-premise and 622 multi-premise divergences. Notably, 𝔼 is not the naive "pointwise Bochvar world semantics" — the pointwise lift produces the complex-algebra reading in which gluts are independent pairs, which is why {p ∧ ¬p} can be designated without infecting everything.
- **Multi-premise consequence is not world-witnessable.** 20 of the 7,956 premise sets have no characterization by value-cells: the single-premise theorem does not lift by pointwise arguments **[W22]**. A genuinely negative result — the abstract characterization of multi-premise ⊨_𝔼 is open.
- **But the paraconsistent signatures survive.** {p, ¬p} ⊭_𝔼 q (non-explosion; B₃ validates it only vacuously), and {p, p} ⊭_𝔼 p∨q (variable blocking in the multi-premise setting; PWK validates it) **[W23]**. 𝔼 keeps both gluts and variable-blocking, each in its honest form.
- **Meta-1 cross-checks.** ¬(μ∨ν) ∧ (μ∧ν) is designated in exactly 4 of 49 assignments **[W24]** — reproducing v1's [V15] with independently written code.
- **Proof theory (standard).** By Rousseau's many-valued sequent method and the Baaz–Fermüller–Zach sets-as-signs construction, any finite matrix has a finite sound and complete signed sequent calculus. 𝔼 has one, with seven signs. **So: yes to a finite rule-based axiomatization (sequent calculus, seven signs); the Hilbert question is vacuous.**

### Verdict

**Inverted and resolved.** 𝔼 is a glut logic (LP core) whose variable-blocking comes from M and runs *forward* — the mirror image of PWK. It is incomparable with PWK and with B₃, and finitely axiomatizable by a signed sequent calculus. The "known paraconsistent logic" the question hoped for is a chimera: 𝔼 is dual to PWK in Beirlaen–De Cock's symmetry, not an instance of it.

### Remains open

Whether "classical consequence + forward variable-sharing" is a named paraconsistent consequence relation in the literature (it is the formal dual of PWK; the dual of a known logic under a known symmetry may still be unnamed). The abstract multi-premise characterization (the 20 witness failures delimit the search space). Signed-calculus proof-theoretic properties: cut elimination, interpolation, and whether 𝔼's seven signs yield a calculus with a normalization theorem worth calling a rite.

---

## 7. Scorecard

| Question | Was (v1 §10) | Now (v1.1) | Status |
|---|---|---|---|
| Q1 Measure | open question | ratio flows on fixed grid [W1–W4]; E/A is channel bias, confirming §0 | answered (toy model) |
| Q2 Arrow fibers | conjectured | witness constructed [W5–W6]; β3 is the obstruction; Greyface = fiber→diagonal promotion | answered (witness) |
| Q3 Topology | "is it a Borel level?" | exactly Σ⁰₂-proper [W7–W10] + standard theorems; no internal hierarchy; not ω-regular | answered |
| Q4 Unification | conjectured | full chain ε-bialgebra + Baxter(−1) ⟹ dendriform ⟹ δ-L, with M2′ = group completion of im f [W11–W16] | answered; converses open |
| Q5 𝔼 | conjectured PWK on M-free fragment | conjecture inverted: M-free = LP [W20]; forward dual of PWK [W18–W19]; signed-sequent axiomatizable [W17] | answered; multi-premise open |

---

## 8. Source notes

| Source | Credibility (x/5) | Last updated |
|---|---|---|
| [Aguiar, *Infinitesimal bialgebras, pre-Lie and dendriform algebras*](https://pi.math.cornell.edu/~maguiar/depaul.pdf) (preprint; journal version in J. Algebra) | 4 | - |
| [Calbrix, Nivat & Podelski, *Ultimately periodic words of rational ω-languages*](https://link.springer.com/content/pdf/10.1007%2F3-540-58027-1_27.pdf) | 5 | 1993 |
| [Szmuc, *An Epistemic Interpretation of Paraconsistent Weak Kleene*](https://philpapers.org/archive/SZMAEI.pdf) | 4 | - |
| Beirlaen & De Cock, *Characterizing Logical Consequence in Paraconsistent Weak Kleene* (Studia Logica) | 5 | - |
| *Proof Theory of Paraconsistent Weak Kleene* (Studia Logica) | 5 | - |
| *On Paraconsistent Weak Kleene: Axiomatisation and Algebraic Analysis* (Studia Logica) | 5 | - |
| Rousseau, *Sequents in many valued logic I* | 5 | - |
| Baaz, Fermüller & Zach, *Sets-as-signs* labelled calculi for finite matrices | 5 | - |
| Landweber, *Decision problems in state machines* (ω-regular hierarchy foundations) | 5 | 1969 |
| Selivanov / Wagner fine hierarchy; Carton–Perrin Wadge–Wagner hierarchy of regular ω-languages | 5 | - |
| *Principia Discordia* (the seven clauses = the seven values of 𝔼) | 5 | 1965 |

*Machine-checked claims [W1]–[W24] were produced by the author's verification scripts for this document and are reproducible from the corpora described in each section.*

---

## 9. New open questions

1. **Baxter converse (Q4).** Which infinitesimal bialgebras admit an idempotent Rota–Baxter operator of weight −1? Is the Loday dialgebra uniquely recoverable from the (ε-bialgebra, Baxter) pair, or is there a genuinely non-Baxter route to δ-L — a chaos that does not factor through a grid?
2. **The r-matrix (Q4).** Is the grid coproduct quasitriangular (Aguiar's second route to dendriform)? If so, the r-matrix may classify the β2 fused nodes: is *the Apple of Discord* a universal R-type element?
3. **Naming 𝔼 (Q5).** Is "classical consequence + forward variable-sharing" a named paraconsistent consequence relation? It is the dual of PWK under the Beirlaen–De Cock symmetry; if unnamed, 𝔼's single-premise fragment deserves one.
4. **𝔼 multi-premise (Q5).** Give an abstract characterization of multi-premise ⊨_𝔼. The 20 non-cell-witnessable sets show world-witnessing fails; the characterization must be algebraic, not semantic.
5. **Beyond flats (Q3).** Repeat the Borel analysis for Dist-valued denotations (N4) and bounded-variability images. Is ⟦𝒢⟧ Σ⁰₂-complete in the full GridF? Is grid-definability closed under any natural operation at all, or is the collapse of the internal hierarchy total?
6. **Culture kernels (Q1).** Identify the observation channels of real cultures empirically, and locate the crystallization threshold: is the Eristic→Aneristic transition a phase transition in the learning rate β? *(Greyface is β above the critical point.)*

---

*"Some look at the Borel hierarchy and see order; some look and see chaos; the grids were on your eyes the whole time."* — after the Principia