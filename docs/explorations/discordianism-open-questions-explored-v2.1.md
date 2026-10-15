# Discordianism: Open Questions Explored (v2.1)

*Answers the six open questions of `discordianism-axioms-v2.md` §10.*

*Tags **[Tn]** are machine checks in `verify_v2_1.py` (16/16 pass). Tags [V], [R], [S], [N] refer to the v2 scripts. Untagged formal claims are standard results with references, or proofs given inline.*

*Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Summary

| # | Question | Answer | Kind |
|---|---|---|---|
| Q1 | Borel class of grid-named *measures* | **Σ⁰₂-complete and meagre** in P(𝒞), for rational or real parameters. Key lemma: a Dirac law produced by a finite-state source sits on an ultimately periodic word **[T1, T2]** | theorem |
| Q2 | A fold that respects deconstruction? | **Yes, and exactly classified.** The endomorphism folds commuting with letter deletion are a ↦ a+g(u), b ↦ a+g(u)−c with u = a−b and g(c) = 0 **[T3–T5]**. The *exclusion* fold b ↦ 0 never qualifies; the *identification* fold b ↦ a does | theorem |
| Q3 | Does one half of the Koszul pair determine the other? | **No.** They share the left operation. δ-L is blind to Podge·Podge, δ-D to Hodge acting on Podge. Together they recover the product: xy = x⊢y − x≻y **[T6, T7]**. The Leibniz bracket recovers im f when im f is self-centralizing **[T8]** | theorem + counterexamples |
| Q4 | Proof theory of 𝔼 = S_fde | **G_𝔼**: a cut-free invertible sequent calculus for LP plus a threaded availability label **[T9, T11]**. Containment is a *global* invariant of the proof tree and cannot be checked leaf-locally **[T10]**. coNP-completeness (v2 L5) rules out polynomial-size proofs in any system unless NP = coNP | calculus |
| Q5 | Real cultures; where could a phase transition appear? | **In coupled cultures.** Mean-field consensus appears at K = βJ = 1. Wrong consensus is metastable below the spinodal H_s(K), and its lifetime grows with culture size **[T12–T14]**. No real data fitted | model |
| Q6 | Arrow fibers, properly posed | **A bureaucracy witness.** With the structure fixed, a fold over the static tree is sound **[T15]**. The diagonal either locks in (3 self-fulfilling fixed points) or never settles (only 2-cycles), depending on the restructuring policy **[T16]**. Greyface = promoting fiber soundness to the diagonal | witness |

---

## 1. Q1 — Dist-valued denotations

**Restatement.** Let grids denote *measures* on 𝒞 = 2^ω (N4). What is the descriptive complexity of the set N ⊆ P(𝒞) of grid-named measures, in the weak-* topology?

**Setup.** Read a Dist-valued grid as a finite-state stochastic source: a hidden Markov model with Moore emissions. Its parameters are rational (countably many grids) or real (a "coin" constructor with any bias). P(𝒞) is compact metrizable but **not zero-dimensional**. Wadge's lemma (used for B1 in v2) therefore does not apply, and completeness needs an explicit reduction.

**Lemma (Dirac pumping).** If an n-state source has output law δ_w, then w = u·v^ω with |u| + |v| ≤ n.

*Proof.* Follow any positive-probability state path. Some state s recurs at times i < j ≤ n. The output law from s is a fixed distribution, so w[i..] and w[j..] are both draws from δ, hence equal.

Exhaustive check over n ≤ 3 with all probabilities in {0, ½, 1}: 6,604 Dirac sources to depth 10, 5,500 of them with genuinely random state paths, 0 violations **[T1]**.

**Theorem 1.** For rational or real parameters, N is **Σ⁰₂-complete and meagre**.

1. **Σ⁰₂.**
   - *Rational:* N is countable.
   - *Real:* for each n, the map from the compact parameter cube to P(𝒞) is weak-* continuous, so its image is compact. N is a countable union of compacta.
2. **Hardness.** x ↦ δ_{F(x)}, where F encodes FIN into ultimately periodic words (1 ↦ 0·1^{i+2}, 0 ↦ 10), is continuous because the Dirac embedding 𝒞 → P(𝒞) is a homeomorphism onto its image. By the lemma, δ_{F(x)} ∈ N ⟺ F(x) is ultimately periodic ⟺ x ∈ FIN **[T2]**.
3. **Meagre.** Every compact piece has dense complement, so N is a countable union of nowhere-dense sets. The dense complement:
   - *Atoms lie on periodic words.* An atom of a finite-state law lies on an ultimately periodic word. Proof: conditioning on longer prefixes drives the atom's conditional mass to 1. So infinitely often some state's law gives mass > ½ to the current suffix. A recurring state then has two atoms of mass > ½, which must coincide.
   - *Density.* The mixtures (1−ε)μ + εδ_w with w non-periodic therefore lie outside N, and they are dense.

**Verdict.** B1 lifts to measures:
- *Topologically,* almost every chaos-law is grid-less, in the Baire sense rather than the measure sense.
- *Arithmetically,* N sits exactly at Σ⁰₂-complete.

The flat/hand-off distinction stays economic (v2 B3), not topological.

**Remains open.** Grids whose `Dist` is not finite-state, e.g. a `Loop` over unbounded counters, which give pushdown-like sources. Is N then still Σ⁰₂, or does it climb?

---

## 2. Q2 — A fold that deconstructs

**Restatement.** Is there an idempotent Rota–Baxter(−1) fold on the grid algebra that is also a morphism for the ε-coproduct Δ (letter deletion)? v2's f (b ↦ 0) is not one [R9].

**Key observation.** Change variables to a and u = a − b. Then Δ(a) = 1⊗1 and Δ(u) = 0. By the ε-rule, Δ deletes one occurrence of a, and the term from the first a is injective on words containing a. Hence:

**ker Δ = k[u]**

Deconstruction cannot see pure difference. Confirmed in degree ≤ 2 **[T3]**.

**Theorem 2.** Let φ be a unital, idempotent algebra endomorphism of k⟨a,b⟩ with Δφ = (φ⊗φ)Δ. Then either φ = id, or

  φ(u) = c (a constant), φ(a) = a + g(u), φ(b) = a + g(u) − c, with g ∈ k[u] and g(c) = 0.

*Proof.*
1. φ(1) = 1, so Δφ(ℓ) = 1⊗1 for each letter ℓ. Therefore φ(ℓ) ∈ a + ker Δ = a + k[u].
2. φ(u) = φ(a) − φ(b) ∈ k[u], so φ acts on k[u] by u ↦ m(u). Idempotence gives m∘m = m. Comparing degrees, deg m ∈ {0, 1}.
3. If deg m = 1, write m(u) = αu + β. Then α² = α and αβ + β = β, so m = id. Then idempotence on a forces g = 0, i.e. φ = id.
4. If m is constant, m = c, and idempotence on a forces g(c) = 0.
5. Conversely, every such φ commutes with Δ: check on generators, then use the ε-rule together with multiplicativity.

Exhaustive search (images of degree ≤ 2, coefficients in {−1, 0, 1}) finds exactly the 14 predicted folds, in both directions **[T4]**.

**Corollary.** Every deconstructing fold has:
- **Hodge** im φ = k[a + g(u)] ≅ k[z]: commutative and generated by one element;
- **Podge** ker φ = (a − b − c), the ideal of *differences*.

The exclusion fold b ↦ 0 is never a coalgebra map: Δ(0) = 0 ≠ 1⊗1. The identification fold b ↦ a is one **[T5]**. Both have the same Hodge, k[a]. They differ in the Podge: (b) versus (a − b).

**Verdict.** *An order that respects deconstruction must conflate, not exclude.*

| fold | action | Podge | respects deconstruction? |
|---|---|---|---|
| Greyface's | b ↦ 0: ignores what is not order | (b) | no |
| the mystic's | b ↦ a, "All is One" | (a − b) | yes |

The price of the mystic's fold is that its order can say only one word. v2's axioms (E1, E1⋆, δ-L, δ-D) hold for both folds, since they need only an idempotent algebra endomorphism.

**Remains open.**
- Idempotent RB(−1) operators that are *not* endomorphisms. These are splittings A = A₁ ⊕ A₂ into two subalgebras, with A₂ not necessarily an ideal.
- The analogue on UT_n with arrow deletion.

---

## 3. Q3 — The Koszul pair

**Restatement.** Does δ-L (dialgebra) determine δ-D (dendriform), or the converse? Does the Leibniz bracket determine (im f, ker f)?

**Block analysis.** Write x = p + n and y = q + m with p, q ∈ im f and n, m ∈ ker f **[T7]**:

| operation | blocks it sees |
|---|---|
| ⊣ = ≺ = x f(y) | pq + nq *(shared by both halves)* |
| ⊢ = f(x) y | pq + pm *(δ-L only)* |
| ≻ = f(x)y − xy | −nq − nm *(δ-D only)* |

- **δ-L is blind to nm** (Podge·Podge).
- **δ-D is blind to pm** (Hodge acting on Podge from the left).

**Theorem 3.** Neither half determines the other.
- UT₂ and k × k[ε] have identical δ-D and different δ-L: they differ only in pm.
- k[x]/(x³) and the square-zero extension k ⊕ k² have identical δ-L and different δ-D: they differ only in nm.

Together the halves recover the product, xy = x⊢y − x≻y, and f(y) = 1 ⊣ y when A is unital **[T6]**.

**Bracket recovery.**
- The left annihilator of the Leibniz bracket is the centralizer of im f. It equals im f exactly when im f is self-centralizing, as in the free model and UT_n.
- The right annihilator is ker f ⊕ (im f ∩ Z(A)). On UT₃ this is ker f ⊕ k·1 **[T8]**.

**Verdict.** The two collision algebras are complementary *views* of one Peirce decomposition, each missing one block. *Leibniz collision sees order acting on chaos; Lie collision sees chaos acting on chaos; only the pair sees the world.*

---

## 4. Q4 — Proof theory of 𝔼 = S_fde

**Restatement.** Give a calculus for 𝔼, organised around availability.

**G_LP.** Work with negation-normal-form, multi-conclusion sequents Γ ⇒ Δ.
- *Rules:* the invertible classical rules for ∧ and ∨.
- *Initial sequents:* identity ℓ ⇒ ℓ and excluded middle ⇒ p, ¬p.
- *No explosion:* p, ¬p ⇒ is not an initial sequent.

**G_𝔼.** Sequents carry a label V, the *available atoms*, fixed to Var(Γ) at the root. A leaf closes only if it is a G_LP initial sequent *and* every atom on its right is in V.

**Theorem 4.**
- G_LP is sound and complete for LP.
- G_𝔼 is sound and complete for 𝔼.
- Cut is admissible in G_𝔼: from Γ ⊢ φ and Γ, φ ⊢ ψ infer Γ ⊢ ψ. The label stays monotone because Var(φ) ⊆ Var(Γ).

Evidence: all 113,424 single-premise pairs (premise ≤ 3 connectives, conclusion ≤ 2), 3,000 random multi-premise instances, and 2,912 random cut instances, with 0 failures **[T9, T11]**.

*Soundness of the label.* The label is exactly the availability lemma made syntactic. An assignment that designates Γ makes every atom of V available, so an unavailable right atom has no witness.

**Containment is global (L1).** Without the label, a leaf-local check ("right atoms ⊆ left atoms of the leaf") rejects the valid p∨q ⇒ p∨q: its leaf p ⇒ p, q fails the check **[T10]**.

In framework terms, variable inclusion is a **global L1 invariant** of the proof tree. G_𝔼 takes response (b) of the framework's Mu-quality check: it threads auxiliary state (the label) rather than restructuring. Szmuc's GPAI_fde for Parry's logic takes the other route: side conditions on non-invertible rules (restricted weakening and ¬R).

**Complexity.** By v2 L5 and Cook–Reckhow, no proof system for 𝔼 has polynomial-size proofs unless NP = coNP. G_𝔼's invertible proof trees are exponential in the worst case, inherited from LP.

**Relatives.**

| logic | definition | relation to 𝔼 |
|---|---|---|
| Parry's PAI_fde | classical + variable inclusion (Szmuc) | validates disjunctive syllogism, which 𝔼 rejects [S1] |
| Petrukhin's natural deduction for S_fde | normalising | the natural-deduction counterpart of G_𝔼 |

**Remains open.** A calculus with *local* rules and no label. This would mean restructuring the global invariant into per-rule side conditions, as GPAI_fde does for Parry. Interpolation for 𝔼.

---

## 5. Q5 — Coupled cultures: where the phase transition lives

**Restatement.** Fit (θ, β) to real cultures. Find where multi-agent coupling makes a genuine phase transition.

**Data.** None fitted. This remains the one empirical gap, and the predictions below are what data would test.

**Model.** N agents. Each agent's crystallized basin is s_i = ±1, with +1 the true basin. Updates follow Glauber dynamics:

  P(s_i = +1) = ½(1 + tanh(K·m₋ᵢ + H))

where:
- m₋ᵢ is the mean of the other agents;
- K = βJ is learning rate × social coupling;
- H = β·h, with h the evidence rate of the world channel (v2 §7, C3).

**Theorem 5** (Curie–Weiss mean field). The fixed points of m = tanh(Km + H) behave as follows:
1. **Critical coupling K = 1.** At H = 0 there is consensus without evidence, m ≠ 0, iff K > 1.
2. **Spinodal.** For K > 1, a *wrong-consensus* branch (m < 0 despite H > 0) exists iff

   H < H_s(K) = K√(1−1/K) − artanh√(1−1/K)

   For example H_s(1.5) = 0.2075 and H_s(2) = 0.5328 **[T12]**.

**Simulation** (N = 400) **[T13]**:
- |m| ≈ 0.03 at K = 0.5 and 0.09 at K = 0.8.
- |m| ≈ 0.72 at K = 1.3 and 0.89 at K = 1.6.
- At K = 1.5, an all-wrong culture stays wrong for 200 sweeps at H = ½H_s, and flips at H = 1.5H_s.

**Lifetime.** The median escape time from wrong consensus grows with N: 23.5, 25.6 and 123 sweeps at N = 16, 32, 64 **[T14]**. This is consistent with Kramers escape, τ ~ exp(N·ΔF).

**Verdict.** The single-agent no-go (v2 C2) and the many-agent transition are two sides of one picture:
- **Greyface needs a crowd.** Below K = 1 a culture tracks evidence. Above it, the culture crystallizes *without* evidence: the Aneristic Illusion becomes spontaneous symmetry breaking.
- **Wrong order is stable when evidence is weak.** For H < H_s(K), a culture that crystallized on the wrong basin stays there, for a time that grows exponentially with its size.
- **v2 C4, impatience, reappears.** Raising β raises both K and H. Fast learners in a crowd cross K = 1 first.

**Remains open.**
- Data. Candidate observables are the consensus order parameter m against measured coupling (network density × imitation rate) and evidence rate.
- Beyond mean field: on sparse networks, does the transition persist (it does for the Ising model on lattices of dimension d ≥ 2), and how do community structures shift H_s?

---

## 6. Q6 — Arrow fibers: the bureaucracy

**Restatement.** Find a system whose structure adapts to data globally, yet which is an Arrow on each fiber.

**Model.** Three form types X, Y, Z are routed to three lanes (fast/mid/slow, latencies 1/3/6).
- A *structure* is an assignment of forms to lanes: a GridF tree `Par[Seq(lane, form)…]`, six in all.
- *Traffic* responds to structure: citizens favour the form currently in the fast lane (q = 0.6).
- *Restructuring* maps traffic to a new structure.

**Fiber 1: structure fixed.** The run is an Arrow. A fold over the static tree, expected latency Σ π(form)·latency(lane), predicts the simulated mean latency for every structure and traffic mix **[T15]**. The Mu–Nu homomorphism holds on this fiber.

**Fiber 2: traffic fixed.** Restructuring is a pure function of the measure. It is Arrow-like in the other direction.

**The diagonal.** Φ = restructure ∘ traffic. Its behaviour depends on the policy **[T16]**:

| policy | rule | diagonal behaviour |
|---|---|---|
| *efficient* | fast lane to the most-filed form | 3 fixed points, one per favoured form: **self-fulfilling order**. Which one the bureaucracy ends in depends on its history (crystallization again, cf. Q5 with K > 1) |
| *fair* | fast lane to the least-filed form | **no fixed point**: every orbit is a 2-cycle. The bureaucracy reorganizes forever |

**Theorem 6** (the Curse of Greyface, formalized). Arrow soundness on every fiber implies nothing about the diagonal. It is compatible both with many fixed points (lock-in) and with none (perpetual reorganization). The global system is monadic: the wiring at t + 1 is bound to the data at t. Mu–Nu analysis is valid per fiber; the diagonal needs a fixed-point analysis of Φ, which no fiber-wise fold can perform.

*Greyface examines each office, finds every office well-designed, and declares the bureaucracy stable.* This is precisely the framework's "monadic systems with Arrow fibers".

**Remains open.** Stochastic Φ, with noisy traffic, which turns the orbits into a Markov chain over structures. Does the "fair" 2-cycle survive noise, or become a two-state stationary distribution (Eristic order)? Larger structure spaces, where Φ can have long cycles or chaotic orbits.

---

## 7. Proposed amendments to the axioms (v2 → v2.1)

**M3′ Two folds.** Distinguish:
- the **exclusion fold** f₀ (b ↦ 0), as in v2;
- the **identification fold** f₁ (b ↦ a).

Both satisfy M3–M4, δ-L, δ-D, E1 and E1⋆. Only f₁, and exactly the family of Theorem 2, is an ε-coalgebra map. The **ε-kernel axiom**: ker Δ = k[a − b]; deconstruction is blind to pure difference.

**δ-K Koszul visibility.** δ-L and δ-D share ≺ = ⊣. Each is blind to one Peirce block (nm, respectively pm), and xy = x⊢y − x≻y.

**B1⁺ Lifted duality.** The grid-named measures form a Σ⁰₂-complete, meagre subset of P(𝒞).

**L6 Calculus.** 𝔼 has the cut-free, availability-labelled sequent calculus G_𝔼. Containment is a global L1 invariant of proofs.

**C5 Coupled cultures.**
- K_c = 1;
- wrong-consensus metastability iff H < H_s(K);
- lifetime exponential in N.

**Γ1 Greyface.** Fiber-wise Arrow soundness does not determine the diagonal. Lock-in and perpetual reorganization are both consistent with it.

---

## 8. New open questions

1. **Non-finite-state sources (Q1).** Which Borel classes do grid-named measures reach once `Loop` has unbounded counters?
2. **Non-endomorphic folds (Q2).** Classify idempotent RB(−1) operators on k⟨a,b⟩, i.e. subalgebra splittings, that commute with letter deletion. Then do the same for UT_n with arrow deletion.
3. **Local G_𝔼 (Q4).** A label-free calculus with local side conditions. Does 𝔼 have Craig interpolation?
4. **Data (Q5).** Measure (K, H) in any real belief-dynamics dataset and test the spinodal prediction.
5. **Noisy bureaucracy (Q6).** Does perpetual reorganization survive noise, and what is the stationary law over structures?

---

## 9. Sources

- Szmuc, *A simple logical matrix and sequent calculus for Parry's logic of Analytic Implication* (PhilArchive preprint): https://philarchive.org/archive/SZMASL
- Petrukhin, *Normalisation for Some Infectious Logics and Their Relatives*, EPTCS 358 (2022): https://arxiv.org/pdf/2204.06722v1
- Kechris, *Classical Descriptive Set Theory* (1995), §17 (weak-* topology on P(X)) and §22.
- Ellis, *Entropy, Large Deviations, and Statistical Mechanics* (1985) — Curie–Weiss mean field.
- Aguiar, *Infinitesimal Hopf algebras*, Contemp. Math. 267 (2000).
- Loday, *Dialgebras*, LNM 1763 (2001).

*"Greyface examined every office and found each one sound; the Bureaucracy has been reorganizing ever since."*
