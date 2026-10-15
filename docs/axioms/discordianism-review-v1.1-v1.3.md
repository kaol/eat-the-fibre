# Review of Le Chat's v1.1–v1.3 (open questions of the v1 Mu–Nu axioms)

*Tags **[Rn]** refer to `audit_algebra.py` (11/11 pass) and **[Sn]** to `audit_logic_topology.py` (14/14 pass). Both scripts use SymPy and exhaustive enumeration. Tags [Wn] are Le Chat's; [Vn] are v1's.*

---

## 0. Bottom line

Le Chat's numerical reproductions are almost all correct: every table I re-ran matches. The failures are in **what the checks were checking**. There are four systematic errors, and each one produced a "result":

1. **Wrong axiom forms.** v1.2's corrections register withdrew v1's Loday dialgebra [V11] and right-Leibniz bracket [V12]. Both withdrawals rest on identities that are not Loday's axioms and not the right-Leibniz identity. **Both v1 results stand** [R1–R3, R7]. v1.3 then compounded the error ("A2 and A5 fail") using the same wrong forms.
2. **A coproduct that is not an ε-bialgebra coproduct.** The whole ε-bialgebra half of the "unification chain" ([W11], [W14], the Apple/r-matrix program, the "mass obstruction") runs on deconcatenation, which is not coassociative [R8]. The genuine ε-bialgebra structure (Aguiar's letter or arrow deletion) exists on both models [R8, R10]. The fold f is not a morphism of it [R9].
3. **A corpus too small to separate the logics.** [W18] (single-premise 𝔼 = CL + containment) fails at three connectives [S1, S2]. That sinks the "Parry's analytic implication" identification (NQ-3). The correct statement: **𝔼 = LP ∩ right containment**, and 𝔼 *is* Deutsch's S_fde [S3–S5, S14].
4. **Tautologies presented as findings.** Several results hold by construction:
   - the logistic in β with R² = 1;
   - "weight −1 is forced";
   - fiber commutativity of maps acting on disjoint coordinates;
   - "proper → complete" at Σ⁰₂;
   - the "unnamed" impatience effect, which is gambler's ruin [S11–S13].

What genuinely survives: the dendriform classification, the coNP-completeness of 𝔼, the Σ⁰₂-completeness of both P and UP, the ⋆-commutator's annihilator, and every reproduced table.

---

## 1. Verdict table

| Claim | Where | Verdict | Evidence |
|---|---|---|---|
| (xf(y), f(x)y) fails a Loday axiom; v1 [V11] withdrawn | v1.2 [W39], v1.3 [W52, W57] | **Wrong.** All five Loday axioms hold on k⟨a,b⟩ (4096 triples) and on UT₄. The "middle exchange" (x⊣y)⊢z = x⊣(y⊢z) is not an axiom; Le Chat's counterexample reproduces against that non-axiom | [R1, R2, R7] |
| [x,y] = xf(y) − f(y)x is not right-Leibniz; v1 [V12] withdrawn | v1.2 [W39] | **Wrong.** Right Leibniz holds on 0/4096 failing triples. The witness (b,a,a) violates only the *left* identity, which v1 never claimed. Failure of antisymmetry is the point of a Leibniz algebra, and v1 stated [x,x] ≢ 0 | [R3] |
| Baxter(−1) ⟹ dendriform ⟹ dialgebra | v1.1 Q4 | **Wrong step.** Dendriform does not imply dialgebra; the two operads are Koszul dual (Dias^! = Dend). The same f gives both structures *in parallel*: dialgebra → Leibniz bracket, dendriform → associative ⋆ → Lie bracket | — |
| k⟨a,b⟩ with deconcatenation is an ε-bialgebra | v1.1 [W11] | **False.** It satisfies the ε-rule but is not coassociative (15/15 words; v1.3 [W60] found this but did not retract [W11]) | [R8] |
| Free algebra may be "barren of coalgebra" (new OQ-3) | v1.3 §13 | **Answered: no.** Letter deletion is coassociative and satisfies the ε-rule; it is Aguiar's own example for path algebras | [R8] |
| f is an ε-bialgebra morphism; M2′ via that morphism | v1.1 [W14], v1.2 [W43] | **Only for the non-coassociative coproduct.** For letter deletion f fails: Δf(ab) = 0 ≠ a⊗ε | [R9] |
| UT₄: only Δ ≡ 0 is a coassociative graded ε-coproduct | v1.3 [W58] | **Ansatz artifact.** UT₄ = kA₄ (path algebra) has arrow deletion, Δ(E_ij) = Σ_{i≤m<j} E_im⊗E_{m+1,j}: nonzero, coassociative, ε-rule. The ansatz only covered degree-0 splittings. (v1.1 [W13]'s "gap coproduct" is probably exactly this and was right; v1.3 contradicts it without registering) | [R10] |
| No quasitriangular structure (grid; UT₄'s 44 inert Apples) | v1.2 [W42], v1.3 [W59–W60] | **Conclusion right, argument on the wrong coproduct.** Clean obstruction: Δ has path-length degree −1, while Δ_r(x) = x·r − r·x has degree ≥ 0. Confirmed: no r on UT₃ | [R11] |
| "The last Apple": is M2′ = kG quasitriangular? | v1.3 OQ-8 | **Two notions conflated.** As a Hopf algebra (Drinfeld), kG is cocommutative, so R = 1⊗1 works trivially. In Aguiar's sense kG is not an ε-bialgebra at all: Δ(gh) = gh⊗gh ≠ g⊗gh + gh⊗h | — |
| Dendriform completions: two ℤ-lines + two points, "34 members" | v1.3 [W49–W51] | **Right shape, wrong cardinality.** Exact variety over any field: two affine lines (t ∈ k, verified symbolically) plus two points. "34 / t ∈ [−8,7]" and "is t continuous over ℝ?" are search-range artifacts | [R4] |
| E1 for the ⋆-commutator: ann = im f; UT₄ gets "im f + E14" | v1.2 [W38], v1.3 [W56] | **Right, and it has a one-line reason.** The Baxter star is x⋆y = pq − nm (with p = f(x), n = x − f(x)), so (A,⋆) ≅ im f × ker f and ann = Z(im f) ⊕ Z(ker f). For UT_n this is diag ⊕ k·E_1n, the centre of the strictly upper-triangular nilradical, for all n | [R5, R6] |
| "Weight −1 is forced on all 25 folds" | v1.2 [W43] | **Tautology.** f(f(x)y + xf(y) − xy) = f(x)f(y) for any idempotent algebra endomorphism. Also "25" is unspecified: k⟨a,b⟩ has infinitely many idempotent endomorphisms (a↦a, b↦aⁿ) | — |
| Single-premise 𝔼 = CL + Var(φ) ⊆ Var(γ) | v1.1 [W18], v1.3 [W44] | **False outside the ≤2-connective corpus.** ((p∨q)∧¬p) ⊨_CL q with variable inclusion, but ⊭_𝔼 q. Also (p∧¬p)∧q ⊭_𝔼 q∧¬q | [S1, S2] |
| 𝔼-single = Parry's analytic implication | v1.2 NQ-3 | **Falls with [W18].** Le Chat's own [W27] shows CL + containment (with explosion) is Bochvar's B₃, from which 𝔼 diverges by 1,668 pairs | [S1] |
| Multi-premise: CL ∧ containment ∧ LP, "none redundant" | v1.2 [W29] | **Right relation, wrong claim of irredundancy.** LP ⊆ CL, so the CL filter never binds. **𝔼 = LP ∩ right containment** (exhaustive single-premise up to 3 connectives, 1,236,544 pairs; 4,000 random multi-premise) | [S3, S4] |
| 𝔼 is "its own logic" | v1.1 Q5 | **It is Deutsch's S_fde.** The 7 values collapse to {T, F, B, N} with N infectious, LP/FDE tables on T,B,F, and D = {T,B}. Consequence relations coincide on all tested pairs. Since LP has no antitheorems, 𝔼 is also LP's containment (right variable-inclusion) companion | [S5, S14] |
| 𝔼 is coNP-complete | v1.3 [W47] | **Correct.** The reduction γ CL-unsat ⟺ γ ⊨_𝔼 ⋁(p∧¬p) holds on all ≤3-connective formulas. Membership is immediate for any finite matrix. "No poly-size cut-free proofs unless NP = coNP" holds for *every* proof system (Cook–Reckhow), not only cut-free ones | [S6] |
| Meta-1 designated 4/49 | v1.1 [W24] | Correct | [S7] |
| Minimal-period and hand-off cost tables; "236 cheaper" | [W7, W10, W36] | Correct | [S8, S9] |
| ⟦𝒢⟧ not ω-regular "because every infinite ω-regular language is uncountable" | v1.1 Q3 | **False reason, true conclusion.** a\*b^ω is ω-regular and countably infinite. Correct argument: a countable ω-regular language is a finite union of U·v^ω (U regular), so it has finitely many periodic tails; P has unboundedly many | [S10] |
| Hand-off = "proper → complete" upgrade | v1.2 NQ-5 | **Wrong** (v1.3 corrected). But v1.3's ruler reduction, though correct, was unnecessary. By Wadge's lemma (Kechris 22.10), a Σ⁰₂ subset of 2^ω that is not Π⁰₂ is Σ⁰₂-complete. Every countable dense set qualifies by Baire, so P and UP were complete as soon as v1.1 showed "not G_δ" | — |
| B1 robust under "every product measure with no certain percept" | v1.1 Q3 | **Too strong.** It holds for i.i.d. (or marginals bounded away from 0 and 1). A product with Σ(1−p_n) < ∞ puts positive mass on 1^ω, a pure flat | — |
| No phase transition in β (R² = 1.000) | v1.2 [W33] | **True by construction.** Intercept −3.203 = ln(10/246); exponential weights give an exact logistic in β. A single-agent filter cannot have a phase transition; that needs a thermodynamic limit | [S12] |
| Rate law t* = 7.80/(…) | v1.2 [W34] | Correct but elementary: 7.80 = ln 99 − ln(10/246), distance over drift | [S12] |
| Impatience effect "appears unnamed" | v1.2 [W35], v1.3 OQ-7 | **It is gambler's ruin.** The exact discrete formula (1−ρ^{k₁})/(1−ρ^{k₁+k₂}) with k = ⌈distance/step⌉ matches every Monte Carlo point (0.889 vs 0.888, 0.669 vs 0.672). The monotonicity in step size is Dubins–Savage *timid play is optimal in a favourable game* | [S11] |
| Arrow fibers witness | v1.1 [W5–W6] | **Reproduces, but it is a vocabulary collision.** F and G act on disjoint coordinates, so they commute trivially. The Arrow condition is about structure being independent of data, not about order-independence of endomaps (framework calibration item 9) | [S13] |
| [W8], [W9] (grid size 2·minper − 1; 199,827 grids → 63 chaoses) | v1.1 | **Unverifiable.** The toy GridF grammar and cost model are never specified | — |

---

## 2. What changes in v1

Nothing in v1 needs withdrawing. The audits add three facts worth folding in:

- **δ-L is half of a Koszul pair.** The idempotent fold f yields:
  - the Loday dialgebra (xf(y), f(x)y), whose bracket is right-Leibniz (v1);
  - via its Rota–Baxter(−1) property, a dendriform pair whose star is x⋆y = f(x)f(y) − (x−f(x))(y−f(y)), i.e. the product algebra im f × ker f.

  Dias and Dend are Koszul dual, so Hodge/Podge appear twice: as Leibniz collision and as Lie collision.
- **E1 has a ⋆-twin.** ann([−,−]⋆) = Z(im f) ⊕ Z(ker f). For UT_n this is diag ⊕ k·E_1n: the apex of the Aftermath is the one chaotic element inert under ⋆-collision.
- **𝔼 = S_fde.** Deutsch (1977), studied by Fitting and Oller; Szmuc's dS_fde is the dual with infectious B. Consequence is LP ∩ containment. Meta-1's glut reading is unchanged.

---

## 3. Citations

- **"Beirlaen & De Cock, *Characterizing Logical Consequence in PWK*".** I could not find a paper by those authors with that title. The PWK characterization it cites is, I believe, Ciuni & Carrara (2016); I have not confirmed that either.
- **The two 2026 arXiv neighbours exist.** [2605.06741](https://arxiv.org/abs/2605.06741) (Li & Li, admissible learning-rate steps) and [2603.21491](https://arxiv.org/abs/2603.21491) (Zhang, "Stable but Wrong"). Le Chat's description of the second ("latent reliability") does not match its abstract (a persistent biased update signal).
- **Confirmed to exist:** Bonzio & Pra Baldi, [*Containment logics*](https://link.springer.com/article/10.1007/s11225-020-09930-1), and Ferguson, [*Logics of Nonsense and Parry Systems*](https://www.researchgate.net/publication/271952110_Logics_of_Nonsense_and_Parry_Systems).
- **S_fde tables:** from Petrukhin, [*Normalisation for Some Infectious Logics*](https://arxiv.org/pdf/2204.06722v1) (EPTCS 358). The tables reached me through a garbled fetch, but they agree with the computation [S14].

---

## 4. Open questions, re-stated after the audit

1. **Dist-valued denotations (OQ-4).** Still open, and still the only genuinely open topology question.
2. **Empirical kernels (OQ-6).** Open, but the predictions to test are now closed forms: gambler's ruin and drift.
3. **Fold and genuine coproduct.** Is there an idempotent Rota–Baxter(−1) map compatible with letter deletion, i.e. a fold that *is* an ε-bialgebra morphism? f is not [R9]. This replaces the "Apple" program, which collapses under the degree obstruction [R11].
4. **Proof theory of S_fde.** A cut-free calculus organised around availability. The literature on S_fde and dS_fde (Petrukhin; Szmuc; Ferguson) is the starting point, rather than Parry.
5. **The Koszul pair.** Do the Leibniz (Dias) and Lie (Dend star) brackets from the same f determine each other? Is either recoverable from the pair (im f, ker f) alone? The second is, by [R5].

*"Eris forgives the wrong number; she is less patient with the right number checked against the wrong axiom."*
