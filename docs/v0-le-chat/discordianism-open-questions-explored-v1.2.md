# Discordianism: Open Questions Explored (v1.2)

*Third wave of the exploration of the v1 Mu–Nu axioms. This document answers the six open questions left by v1.1 (its §9, NQ-1 through NQ-6). Tags **[W25]–[W43]** were machine-checked this session (TypeScript sandbox: exhaustive enumeration over finite witnesses, direct evaluation of the 𝔼 matrices, and linear-algebra over the dendriform ansatz). Tags **[W1]–[W24]** refer to v1.1's checks; **[V#]** are unchanged from v1's Z3/SymPy script. Untagged formal claims are standard results. Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Method

Same two answer types as v1.1: **witness answers** (finite toy models, exhaustively enumerated — existence of a mechanism, not uniqueness) and **structural answers** (theorems plus exhaustive checks). One addition this wave: a **corrections register** (§8). Re-verification with stricter, fully-exhaustive checks found three claims of v1.1 that do not survive. They are withdrawn explicitly, and the structures they pointed at are repaired. This is B5 applied to ourselves: the check was only good because it could have returned another number.

---

## 1. Executive summary

| # | v1.1 asked | v1.2 answer | verdict |
|---|---|---|---|
| NQ-1 | Which ε-bialgebras admit an idempotent RB(−1)? Is the dialgebra uniquely recoverable? | **Weight −1 is forced** (every idempotent multiplicative fold is RB(−1), none weight 0). But **uniqueness fails**: the dendriform layer is a family, not a point; and the dialgebra layer of v1.1 is withdrawn — the honest structure is dendriform ⟹ associative Baxter star ⟹ **Lie** bracket, with E1 intact | answered, with corrections |
| NQ-2 | Is the grid coproduct quasitriangular (an Apple of Discord)? | **No.** All quasitriangular forms are inconsistent — machine-proved for support ≤ 2, structurally proved for any finite support. Dendriform structure arises solely via the Baxter route | answered (negative) |
| NQ-3 | Is "classical + forward variable-sharing" a named logic? | **Yes, single-premise: it is Parry's analytic implication** (the right containment companion ⊢ʳ). 𝔼's multi-premise consequence is a proper refinement of it | answered |
| NQ-4 | Abstract characterization of multi-premise ⊨𝔼? | **Γ ⊨𝔼 φ ⟺ Γ ⊨_CL φ ∧ Var(φ) ⊆ Var(Γ) ∧ Γ ⊨_LP φ** — 0 violations across 1.37 million checked pairs. The 20 "non-witnessable" sets of v1.1 were the shadow of this | answered |
| NQ-5 | Beyond flats: Dist-valued denotations, Σ⁰₂-completeness, closure | **β1's hand-off is exactly the proper→complete upgrade**: with hand-off grids, ⟦𝒢⟧ = all ultimately periodic words, which is **Σ⁰₂-complete** (flats alone are only Σ⁰₂-proper), at half the cost per word | answered (deterministic case); Dist open |
| NQ-6 | Is crystallization a phase transition in β? | **No.** The sweep is a smooth logistic (R² = 1.000); weak biases all crystallize with a drift-predicted rate; and **truth-tracking degrades monotonically with learning rate** — fast learners crystallize on the wrong basin more often | answered (toy model) |

---

## 2. NQ-1 — The Baxter converse: forced weight, non-unique star

### Restatement

*Which infinitesimal bialgebras admit an idempotent Rota–Baxter operator of weight −1? Is the Loday dialgebra uniquely recoverable from the (ε-bialgebra, Baxter) pair?*

### Findings

- **The weight is forced, not chosen.** On the free grid algebra k⟨a,b⟩, all 25 idempotent-and-multiplicative folds (Aneristic collapses of the grid side) are Rota–Baxter operators of weight −1, and none is weight 0 **[W43]**. The grid algebra does not merely *admit* a Baxter(−1) — every fold it admits is one. The converse direction (a general criterion for ε-bialgebras) remains open, but on this algebra family the answer is total.
- **The fold is a strict ε-bialgebra morphism only rarely.** Of the 25 folds, only 6 commute with the coproduct — exactly the letter-substitutions **[W43]**. The Aneristic collapse f of v1.1 is one of the 6: it is both a Baxter operator and a bialgebra morphism, which is why the unification chain worked through it.
- **Fold identities.** f(f(x)y) = f(xy) and f(xf(y)) = f(x)f(y), verified on all 961 word-pairs **[W40]**; f is RB(−1) on all 900 pairs of words of length ≤ 4 **[W40]**. And the bracket separates the folds: 22 of the 25 are distinguished by their ⋆-commutators (the collisions are exactly the 4 null/void folds); E1 holds for a fold precisely when its image is commutative **[W43]**.
- **Uniqueness fails.** Fix the left product ≼ = x·f(y) and ask for all dendriform completions with a 7-term ansatz ≻ = ΣβᵢTᵢ over the natural bilinear terms. D1 imposes exactly **one** condition, Σβᵢ = 0 (a rank-1 system — the identities collapse everything else); D2 is **vacuous** (0 rows); D3 is quadratic and cuts the candidates but does not pin them: of 393 sum-zero {−1,0,1} coefficient vectors, **107 survive D3** on words up to length 2, and every survivor has an associative ⋆ **[W41]**. The Aguiar-canonical member β* = (−1,1,0,0,0,0,0), i.e. ≻ = f(x)·y − x·y, is one of them (verified to length 3). *The dendriform layer is a family; "uniquely recoverable" is false as stated.*
- **The corrected dendriform structure.** With ≼ = x·f(y) and ≻ = f(x)·y − x·y, all three dendriform axioms D1–D3 hold (0/216 bad each, plus 1,464 spot-triples at length ≤ 7); the star ⋆ = ≼ + ≻ is associative (0 bad) **[W37]**. Note what ⋆ is: x⋆y = x·f(y) + f(x)·y − x·y — precisely the **Baxter star product of weight −1**. Aguiar's route from a Baxter operator to a dendriform algebra is, on this algebra, literally the star-product splitting.
- **The dialgebra layer is withdrawn; the Lie layer is installed.** See §8 for the register. What survives and strengthens: the ⋆-commutator is a genuine **Lie bracket** (Jacobi identity, 0/216 bad), and its **left and right annihilator is exactly im f** = span{ε, a, aa, aaa} — both sides, no ker f remainder **[W38]**. E1 is not merely preserved by the correction; it is *cleaner* in the repaired structure than in the withdrawn one.

### The corrected unification chain

```
grid algebra: ε-bialgebra + Rota–Baxter f (weight −1, forced [W43])
        │  Aguiar: Baxter(−1) ⟹ dendriform  (≼ = x·f(y), ≻ = f(x)·y − x·y [W37])
        ▼
   dendriform algebra  ⟹  ⋆ = Baxter star, associative [W37]
        │  ⋆-commutator: genuine Lie bracket [W37]
        │  E1: ann([−,−]⋆) = im f, left and right, exactly [W38]
        │  im f, group-completed
        ▼
   kG Hopf algebra (M2′, the Hodge Hopf quotient)
```

*(Withdrawn from this diagram: the Loday dialgebra step and the bracket xf(y) − f(y)x — see §8.)*

### Verdict

The converse closes on the algebra side: weight −1 is forced and the canonical dendriform structure is verified with the corrected axioms. But the recovery question inverts — there is no unique dialgebra to recover, because there is no dialgebra at all in the naive sense, and the dendriform layer is a 107-member family (at minimum) whose canonical member Aguiar's construction itself selects.

### Remains open

Classification of the D3-survivor family: which of the 107 short-word survivors extend to all word lengths, and what structure does the full family carry (is it a torsor? does any member admit the middle-exchange axiom that the canonical one lacks?). The general converse: a criterion on an abstract ε-bialgebra for admitting an idempotent Baxter(−1). UT₄ under the corrected operations.

---

## 3. NQ-2 — The r-matrix: no Apple of Discord

### Restatement

*Aguiar's second route to dendriform structure runs through quasitriangular ε-bialgebras. Is the grid coproduct quasitriangular — does an r-matrix exist on the grid algebra, classifying the β2 fused nodes? Is the Apple of Discord a universal R-type element?*

### Findings

- **Machine search: inconsistent.** Parametrize r = Σᵢ uᵢ ⊗ vᵢ over the 49-dimensional space of support-≤ 2 tensors on words of length ≤ 2, and solve the linear system each quasitriangular form imposes. All three forms are inconsistent: 118 equations of rank 33, 38 of rank 14, and 118 of rank 33 — in every case the system overdetermines itself into contradiction **[W42]**.
- **Structural proof: no finite support works.** In the quasitriangular identity, the coefficient of (b ⊗ w) on the left side is 1 exactly when the word w *starts with* b — a property of w. On the right side, each form forces this coefficient to equal the support-constant #{i : uᵢ = b, vᵢ = ε} — a property of r alone, independent of w. Varying w over words that start with b and words that do not makes the contradiction unavoidable for any finite support **[W42]**. The toy search was not a small-model artifact; the obstruction is structural.
- **Aguiar's own precedent.** In *Infinitesimal Hopf algebras*, the dual numbers carry a quasitriangular ε-bialgebra structure but the polynomial algebra k[x] does not — the grid algebra, being free, behaves like the second.

### Verdict

**No: the grid coproduct is not quasitriangular.** There is no Apple of Discord in the algebra — no universal element generates the fusion of nodes, and Aguiar's two routes to dendriform structure decouple cleanly: the grid algebra takes the Baxter route only. Discord, in this framework, is not a universal object; it is an operator's kernel.

### Remains open

Whether a *weakened* quasitriangularity (relative r-matrices, or triangular structures on the group-completed im f side, i.e. on the M2′ Hopf algebra rather than the grid ε-bialgebra) classifies β2. The question moves with the quotient: the Apple of Discord, if it exists, lives downstream of the fold.

---

## 4. NQ-3 — Naming 𝔼: it is Parry's analytic implication

### Restatement

*v1.1 characterized single-premise 𝔼 as classical consequence plus forward variable-sharing — the formal dual of PWK under the Beirlaen–De Cock symmetry. Is this a named paraconsistent consequence relation?*

### Findings

- **Yes: Parry's analytic implication.** Parry's containment requirement — a proof of φ from γ must use only variables already in γ, i.e. Var(φ) ⊆ Var(γ) — is exactly forward variable-sharing. Combined with the re-verified [W18] (γ ⊨𝔼 φ ⟺ Var(φ) ⊆ Var(γ) ∧ γ ⊨_CL φ; 0 mismatches on all 88,209 pairs, re-run this session with the availability machinery **[W26]**), the identification is definitional: **𝔼-single = the right containment companion ⊢ʳ of classical logic = Parry's analytic implication.** In the notation of the containment-logic literature (the ⊢ˡ/⊢ʳ companions studied by Ferguson and Pra Baldi), PWK-single is the *left* companion of classical logic and 𝔼-single is the *right* one — the duality v1.1 observed abstractly is the standard left/right companion symmetry.
- **The bridge is Ferguson's.** Ferguson's *Logics of Nonsense and Parry Systems* shows that containment logics arise precisely as the consequence relations of nonsense logics (Bochvar, Halldén) — which is 𝔼's situation verbatim: 𝔼 is the complex algebra of a Bochvar-type nonsense value (M), and its consequence is a containment logic. v1's seven-valued construction re-derives a 1933 program it never cited.
- **The dual confirms itself at every comparison point.** B₃ (Bochvar's own three-valued logic, the non-lifted neighbor) has single-premise consequence: [γ classically unsatisfiable] ∨ [γ ⊨_CL φ ∧ Var(φ) ⊆ Var(γ)] — 0 mismatches **[W27]**: B₃ explodes only vacuously, and shares the forward containment direction. LP sits strictly between: 𝔼 ⊊ LP ⊊ CL (11,892 and 1,164 divergences at n = 3) **[W28]**; LP has tautologies (4 at n = 2, 6 at n = 3) while 𝔼 has none — *M kills the tautologies; the glut alone does not.* PWK vs LP: 11,376 divergences, witness p ∧ ¬p ⊨_PWK p ∧ q but ⊭_LP **[W28]** — the same separator geometry as 𝔼 vs PWK, from the other side.

### Verdict

**The unnamed relation is named: Parry's analytic implication, single-premise.** But the naming has a boundary: 𝔼's *multi-premise* consequence (NQ-4) is strictly stronger than Parry's multi-premise analytic implication, because of the LP clause. 𝔼 is not a paraphrase of a known system — it is a known single-premise system wearing a genuinely new multi-premise layer.

### Remains open

Whether the multi-premise refinement has a name (see NQ-4 — the answer is a hybrid of Parry and LP, and the hybrid appears novel); a signed sequent calculus for 𝔼 organized around the availability lemma (§5) rather than the seven signs, and its cut-elimination behavior.

---

## 5. NQ-4 — Multi-premise 𝔼, characterized

### Restatement

*Give an abstract characterization of multi-premise ⊨𝔼. v1.1's 20 non-cell-witnessable sets showed world-witnessing fails; the characterization must be algebraic.*

### Findings

- **The availability lemma (the key).** Call an atom *available* in an environment if its value contains a classical one (T or F). Then, for a single formula γ: T ∈ ⟦γ⟧_env ⟺ the occurrence-separated γ* has a classical model within the available values, **and every atom occurring in γ is available**. An unavailable atom poisons the whole formula — not pointwise, per occurrence, but globally: M-infection couples the choices so that naive per-occurrence independence is wrong. Verified on 106,869 checks (all 102 two-atom formulas over 49 environments and all 297 three-atom formulas over 343), 0 failures **[W25]**.
- **The multi-premise characterization.** For every availability set A: if every premise is A-satisfiable, then φ is A-satisfiable. Equivalently — and this is the abstract answer — verified with **0 violations across all 1,366,515 checked pairs** (all 5,253 two-atom premise sets × 102 formulas, plus all single-premise and 2,500 sampled premise pairs at three atoms) **[W29]**:

  > **Γ ⊨𝔼 φ ⟺ Γ ⊨_CL φ ∧ Var(φ) ⊆ Var(Γ) ∧ Γ ⊨_LP φ.**

  Three filters, none redundant: classical validity, Parry containment (the ⊢ʳ companion), and LP-validity (glut discipline). The direct availability check and the direct 7-valued evaluation agree everywhere they were both computed **[W30]**.
- **Placement, quantified.** Each clause cuts strictly: 𝔼 ⊊ ⊢ʳ (28,856 divergences), 𝔼 ⊊ CL (47,250), 𝔼 ⊊ LP (11,240) over the multi-premise corpus **[W31]**; 𝔼 vs PWK multi-premise: 117,182 divergences — the two logics are nowhere near each other in the multi-premise setting either.
- **Explosion is structurally impossible, not merely false.** Every one of the 5,253 premise sets has a designated reading (the all-{T,F} seed designates anything): 𝔼 has **zero inconsistent premise sets**, so ex falso has no input to operate on **[W31]**. This sharpens v1's non-explosion: paraconsistency here is a property of the semantic space, not a rule the consequence relation refuses.
- **Witnesses, each isolating one clause.** {p ∨ q, ¬p ∨ q} ⊭𝔼 q (fails LP: environment w(p) = {T, F, M}, w(q) = {F} designates both premises, not the conclusion); {p, ¬p} ⊭𝔼 q but ⊨𝔼 p (the glut is real and it is *confined*); {p} ⊭𝔼 p ∨ q (fails containment: q was never in evidence) **[W31]**.

### Verdict

**Answered, and the answer is honest about its shape:** 𝔼-multi is the intersection of classical logic, the right containment companion, and LP — a three-way filter, each clause with its own witness class. v1.1's 20 "non-witnessable" sets were precisely the ones where the three filters disagree with any single-filter reading: the algebraic characterization exists, and the 20 sets delimit it correctly.

### Remains open

Complexity: how hard is Γ ⊨𝔼 φ to decide via the availability characterization (the quantification over availability sets A ranges over exponentially many, but the three-filter form suggests polytime for bounded premises — unverified). Axiomatization: Parry's system + LP's rules is a candidate Hilbert-style calculus; soundness is checkable by machine, completeness is not.

---

## 6. NQ-5 — Beyond flats: the hand-off is the complete upgrade

### Restatement

*v1.1's toy GridF used deterministic flat denotations, yielding ⟦𝒢⟧ = the pure flats P, exactly Σ⁰₂-proper. What changes under β1's hand-off grids? Is ⟦𝒢⟧ Σ⁰₂-complete in the full GridF? Dist-valued denotations remain open.*

### Findings

- **The hand-off halves the price of order.** Let a hand-off grid be a pair (u, v) denoting the stream u·v^ω, at cost |u| + 2|v| − 1 (the flat u^ω costs 2·minper(u) − 1). Then **every** binary word of depth 8 is denoted at cost ≤ 8 — coverage 256/256 at budget 8, versus 52/256 for flats at budget 9 — and the hand-off is cheaper than the flat for 236 of 256 words; cost distribution {1: 2, 2: 2, 3: 6, 4: 10, 5: 26, 6: 46, 7: 90, 8: 74} **[W36]**. Max cost over all depth-8 words equals the depth exactly (proved: take u = the word minus its last letter, v = the letter), so coverage is *tight*: no budget-8 word is left out and no word needs 9.
- **The denotation set jumps from proper to complete.** With hand-off grids, ⟦𝒢⟧ = UP, the set of **all** ultimately periodic words — strictly larger than the pure flats P (a flat is the special case u = ε… up to rotation). UP is countable and dense (β1 density, as before), hence exactly Σ⁰₂-proper. And more: **UP is Σ⁰₂-complete.** Standard proof: FIN, the Σ⁰₂-complete set of eventually-zero binary sequences, reduces to UP by encoding a 1 at position i as 0·1^(i+2) and a 0 as 10; finitely many 1s make the encoding eventually (10)^ω (ultimately periodic), while infinitely many 1s produce unbounded runs of 1s, which no ultimately periodic word has. So membership in UP is Σ⁰₂-hard; being Σ⁰₂, it is **complete**.
- **The comparison, in one line.** Flats: ⟦𝒢⟧ = P, Σ⁰₂-proper, coverage 52/256, cost 2·minper − 1. Hand-off: ⟦𝒢⟧ = UP, **Σ⁰₂-complete**, coverage 256/256, cost |u| + 2|v| − 1 ≤ depth. β1's hand-off interface — the same one that Aguiar's ε-bialgebra rule made mandatory in v1.1 Q4 — is *exactly* the proper-to-complete upgrade at the same Borel level, bought at half price.
- **B1 survives the upgrade.** UP, like P, is null under every product measure on 𝒞 in which no single percept has probability 1 — order still covers a null set, whichever coin you flip.

### Verdict

**Answered for deterministic denotations.** The effect of β1's hand-off on the topology is now sharp: same Borel level, proper → complete, and complete at the information-theoretic minimum (coverage is tight). The σ-algebra triviality of v1.1 carries over verbatim since the point-cloud σ-algebra only depends on the countability and density of ⟦𝒢⟧.

### Remains open

Is P itself Σ⁰₂-complete, or properly incomplete? (The reduction above needs the pre-period; a flat-only encoding of FIN is not obvious.) Dist-valued denotations (N4): a `Dist`-valued ⟦−⟧ names measures on 𝒞, and the natural target is a set of measures — the Borel analysis restarts in the space of measures (e.g. the weak-star σ-algebra), where the flat/hand-off dichotomy presumably reappears as a question about the topology of eventually-repeating measure sequences.

---

## 7. NQ-6 — Culture kernels: no phase transition, a rate law, and a cost of impatience

### Restatement

*Is the Eristic→Aneristic crystallization a phase transition in the learning rate β? Identify the observation channels of a culture.*

### Findings

- **The v1.1 simulation reproduces, with two corrections.** Same 256 grids (A = 10, E = 246, Five = 56; minimal-period distribution as in v1.1 [W10]). Neutral channel: flat at P(A) = 0.0391 and entropy ln 256 = 5.545 **nats** through 50,000 steps — v1.1 reported the same number mislabeled as bits. Aneristic channel: P(A) → 1, entropy → ln 10 exactly. Eristic channel: → ln 246. Five-rewarding: → ln 56 exactly **[W32]**. Correction two: "P(A) → 1.0 in 15 steps" was seed-dependent — the median at 15 steps is 0.95; convergence is ~1 by 60 steps. The channel conclusions are unchanged; the *rate* was underreported.
- **No phase transition in β.** Sweeping β from 0.05 to 1.5 under the Aneristic channel: the crystallization probability is a smooth logistic, logit(P(A)) = −3.203 + 7.69β, R² = 1.000, maximum jump between adjacent settings 0.096 **[W33]**. Nothing sharpens; Greyface does not have a threshold.
- **Greyface has a rate.** With a weakly biased channel (percept accuracy θ ∈ {0.51, 0.55, 0.6, 0.7} — barely better than coin-flipping — and β ∈ {0.4, 0.8}, 20 seeds, horizon 100,000): **every** run crystallizes. Median fixation time matches the stochastic-drift prediction t* = 7.80 / (β(2θ−1)·ln(θ/(1−θ))) within ratios 0.77–1.19 across all eight (θ, β) cells **[W34]**. Only exact neutrality (θ = 0.5) never crystallizes. However weak the channel, if it is biased at all, the kernel freezes eventually: there is no Aneristic-free equilibrium, only slow ones.
- **The illusion needs world–channel cooperation.** A control with reversed evidence direction still crystallizes on the truth basin it cooperates with: the Aneristic Illusion is not a property of the channel alone or the world alone, but of their agreement **[W32]**.
- **Impatience breeds wrong order.** The cost: with truth in the minority basin (θ = 0.75, 400 seeds per setting), the fraction of cultures crystallizing on the *true* basin falls monotonically with learning rate: β = 0.2 → 1.000, 0.4 → 0.995, 0.8 → 0.888, 1.2 → 0.888, 1.6 → 0.672, 2.0 → 0.672 **[W35]**. The log-odds barriers run from ln(10/246) = −3.20 to ±ln 99 (the boundaries of the two basins), and a Brownian-approximation first-passage calculation reproduces the curve (predicted 1.000 / 0.985 / 0.879 / 0.755 / 0.653 / 0.573; the discrete walk overshoots the barrier at large β). Fast learners crystallize on the wrong — majority — basin more often: **speed purchases certainty at the price of truth.** The slow culture is the one that ends up ordered on the world rather than on its own channel noise.

### Verdict

**No phase transition — crystallization is a drift-to-fixation process with a rate law, not a regime change.** And the v1.1 picture gains an honest moral: the disease has no critical point, only a half-life, and its treatment is patience. A culture cannot avoid crystallizing on *something*; it can only choose how carefully.

### Remains open

Empirical kernel identification for real cultures (which channels, which θ, which β). Whether the monotone truth-cost of learning rate is a known result in the sequential-testing/learning-in-games literature (the setting is a biased-coin Bayesian filter with two absorbing basins; the first-passage computation suggests the effect is classical, but a citation would pin it). Multi-agent coupling: cultures observing each other's kernels, where the channel bias is endogenous.

---

## 8. Corrections register

Claims from earlier waves that this wave's stricter checks do **not** confirm, with their replacements:

| Original claim | Status | What replaced it |
|---|---|---|
| v1.1 [W15]: "all five Loday dialgebra axioms hold with x ⊣ y = x·f(y), x ⊢ y = f(x)·y" | **Withdrawn** — re-checked against the correct axiom forms, the middle-exchange axiom (x ⊣ y) ⊢ z = x ⊣ (y ⊢ z) fails, with counterexamples (ε, ε, b): b vs 0, and (b, a, a): 0 vs baa; the other four axioms pass **[W39]**. (The earlier pass used a wrong form of the axiom.) | The dendriform axioms D1–D3 with ≼ = x·f(y), ≻ = f(x)·y − x·y, all verified **[W37]** — one level down the Loday hierarchy, honestly |
| v1.1 [W16] / v1's δ-L: the bracket [x, y] = x·f(y) − f(y)·x is a right-Leibniz bracket with a stated right annihilator | **Withdrawn** — the bracket fails antisymmetry on 16/36 pairs and right-Leibniz on 32/216 (counterexample (b, a, a): one side 0, the other baa − 2aba + aab) **[W39]**; the naive dialgebra bracket fails right-Leibniz on 5/16 pairs. The dialgebra/δ-L layer of v1 §10-Q4 does not survive exhaustive checking | The ⋆-commutator [x, y]⋆ = x⋆y − y⋆x: a genuine Lie bracket (Jacobi verified), with ann(⋆) = im f on **both** sides **[W37, W38]** — a Lie layer where v1 wanted a Leibniz layer |
| v1.1 Q1: entropy figures in "bits"; P(A) → 1.0 at 15 steps | **Corrected** — entropies were nats throughout (ln 256 = 5.545, ln 10 = 2.303, ln 246 = 5.505, ln 56 = 4.026); 15-step convergence was seed-dependent (median 0.95; ~1 by 60 steps) **[W32]** | Same numbers, correct units and rates; all qualitative conclusions of v1.1 Q1 stand |

The Q4 unification chain of v1.1 is thus *shorter and stronger*: ε-bialgebra + Baxter(−1) ⟹ dendriform ⟹ Baxter star ⟹ Lie bracket with E1, and im f group-completes to M2′. The dialgebra step was the only weak link, and removing it costs nothing that the corrected layer does not restore.

---

## 9. Scorecard

| Question | Was (v1.1 §9) | Now (v1.2) | Status |
|---|---|---|---|
| NQ-1 Baxter converse | open | weight −1 forced on all 25 folds [W43]; dendriform corrected and verified [W37]; dialgebra withdrawn [W39]; recovery is a ≥107-member family [W41] | answered, with corrections |
| NQ-2 r-matrix | open | inconsistent for all finite supports [W42]; no Apple of Discord | answered (negative) |
| NQ-3 Naming 𝔼 | open | single-premise = Parry's analytic implication (⊢ʳ); multi is a proper refinement [W26–W28] | answered |
| NQ-4 𝔼 multi-premise | open (20 witness failures) | Γ ⊨𝔼 φ ⟺ CL ∧ ⊢ʳ ∧ LP, 0 violations on 1.37M pairs [W25, W29–W31] | answered |
| NQ-5 Beyond flats | open | hand-off ⟹ full tight coverage [W36]; ⟦𝒢⟧ = UP is Σ⁰₂-**complete**; B1 survives | answered (deterministic); Dist open |
| NQ-6 Culture kernels | open | no phase transition [W33]; rate law t* verified [W34]; monotone truth-cost of β [W35] | answered (toy model) |

---

## 10. Source notes

| Source | Credibility (x/5) | Last updated |
|---|---|---|
| [Aguiar, *Infinitesimal bialgebras, pre-Lie and dendriform algebras*](https://pi.math.cornell.edu/~maguiar/depaul.pdf) (also [arXiv math/0211074](https://arxiv.org/pdf/math/0211074)) | 4 | - |
| [Aguiar & Loday, *Quadri-algebras*](https://pi.math.cornell.edu/~maguiar/quadri.pdf) | 4 | - |
| [Dendriform trialgebra (Baxter operators of arbitrary weight)](https://arxiv.org/pdf/math/0311132v1) | 4 | - |
| Aguiar, *Infinitesimal Hopf algebras* (Contemp. Math. 267; dual numbers quasitriangular, k[x] not) | 5 | 2000 |
| [Ferguson, *Logics of Nonsense and Parry Systems*, J. Phil. Logic 44(1)](https://link.springer.com/article/10.1007/s10992-014-9321-y) (also [JSTOR](https://www.jstor.org/stable/44851254)) — the NQ-3 bridge: containment logics from nonsense logics | 5 | 2015 |
| *Containment Logics: Algebraic Completeness and Axiomatization* (Studia Logica; the ⊢ˡ/⊢ʳ companions) | 5 | - |
| Pra Baldi, doctoral thesis, University of Padua (containment logics) | 4 | - |
| Beirlaen & De Cock, *Characterizing Logical Consequence in Paraconsistent Weak Kleene* (Studia Logica) | 5 | - |
| [Szmuc, *An Epistemic Interpretation of Paraconsistent Weak Kleene*](https://philpapers.org/archive/SZMAEI.pdf) | 4 | - |
| [Calbrix, Nivat & Podelski, *Ultimately periodic words of rational ω-languages*](https://link.springer.com/content/pdf/10.1007%2F3-540-58027-1_27.pdf) | 5 | 1993 |
| Landweber, *Decision problems in state machines* (ω-regular hierarchy foundations) | 5 | 1969 |
| Parry, analytic implication (via Ferguson 2015; original 1933) | 5 | 1933 |
| *Principia Discordia* | 5 | 1965 |

*Machine-checked claims [W25]–[W43] were produced by the author's verification scripts for this document and are reproducible from the corpora described in each section.*

---

## 11. New open questions

1. **The D3-survivor family (NQ-1).** Which of the 107 dendriform completions that survive D3 on short words extend to all lengths? Is the family a torsor (is there a group action permuting the survivors, with the Aguiar-canonical member as identity)? Does any member satisfy the middle-exchange axiom that the canonical one fails — i.e., is the dialgebra layer *elsewhere* in the family rather than gone?
2. **UT₄ under the corrected operations (NQ-1).** Redo the fold/Baxter/E1 analysis in the v1 upper-triangular model with the ⋆-commutator: does the nilpotent Podge (the Aftermath) interact differently with a Lie layer than v1's Leibniz layer predicted?
3. **Flats completeness (NQ-5).** Is P (the pure flats) Σ⁰₂-complete, or properly Σ⁰₂-proper-and-incomplete? The hand-off reduction needs a pre-period; a flat-only reduction would have to hide the pre-period in the periodic part.
4. **Dist-valued Borel analysis (NQ-5).** Lift the toy GridF to `Dist`-valued denotations (N4) and redo §6 in the space of measures: what is the σ-structure of the set of grid-named measures, and does the flat/hand-off dichotomy survive as eventually-repeating vs pre-periodic measure sequences?
5. **Complexity of availability (NQ-4).** The characterization Γ ⊨𝔼 φ ⟺ CL ∧ ⊢ʳ ∧ LP suggests polytime decidability for bounded premises — is the availability-quantified form decidable in polynomial time in general, and what is the proof-theoretic cost of the three-filter calculus (cut elimination for a Parry+LP hybrid)?
6. **Empirical kernels (NQ-6).** Fit (θ, β) to a real culture's belief-dynamics data. The rate law of [W34] is falsifiable: a culture that has *not* crystallized under measured bias θ for longer than 1.5·t* refutes the drift model — or measures its Eristicity.
7. **The impatience effect in the literature (NQ-6).** The monotone wrong-basin probability in β matches a biased-coin first-passage calculation; is this exact effect (fast sequential learners fixate on the majority hypothesis more often) already named in sequential testing or learning-in-games? If not, it is a small new theorem wearing Discordian clothes.
8. **The Apple downstream (NQ-2).** The M2′ Hopf algebra (group-completed im f) may be quasitriangular even though the grid ε-bialgebra is not. If it is, the Apple of Discord exists — but on the order side of the fold, which is a very Discordian place for it.

---

*"Greyface has no critical point; he has a half-life. The cure is not a threshold but a rate — and haste, as always, is the only true error."* — after the Principia