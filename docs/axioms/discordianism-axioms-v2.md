# Discordianism: Mu–Nu Axioms (v2)

*Supersedes v1. It integrates the results of Le Chat's exploration waves v1.1–v1.3 that survived independent re-checking, and corrects the ones that did not (§11).*

*Every formal claim tagged **[Vn]**, **[Rn]**, **[Sn]** or **[Nn]** is machine-checked. The tags refer to four scripts, 49/49 passing:*

| tags | script |
|---|---|
| [V] | `verify_discordian.py` (Z3, SymPy, NumPy) |
| [R] | `audit_algebra.py` |
| [S] | `audit_logic_topology.py` |
| [N] | `verify_v2_new.py` |

*Untagged formal claims are standard results with a reference. Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Thesis

Under the Mu–Nu reading, μ and ν are not "order" and "disorder". μ is the **initial algebra** of a functor and ν its **final coalgebra**. Principia's metaphysics is:

> Order and disorder are both grids laid over an underlying Chaos.

This metaphysics becomes the type structure:

- **Chaos is the Nu.** It is productive, observed only by unfolding, and equal only up to bisimulation.
- **Hodge and Podge are both Mu.** The Aneristic and Eristic grids are finite descriptions. They differ by a fold label, not by type.
- **The Sacred Chao** is the denotation ⟦−⟧ from grids into Chaos.

v2 adds two layers to v1:

1. **The fold f.** One idempotent algebra endomorphism, the Aneristic projection, generates both collision algebras. It splits the grid algebra into Hodge ⊕ Podge, and the two collision algebras are Koszul-dual readings of that split (§5).
2. **A measure layer.** A culture is a measure on a *fixed* set of grids. The Aneristic and Eristic Illusions are crystallizations of that measure (§7).

---

## 1. Ontology

```
GridF atom r  =  Atom atom | Seq r r | Par r r | Neg r | Void      -- closed, 5 constructors
𝒢  =  μ GridF                                                      -- grids (Hodge ∪ Podge)
ObsF X        =  Obs → Dist (Percept × X)                           -- Auto Dist
𝒞  =  ν ObsF                                                       -- Chaos
⟦−⟧ : 𝒢 → 𝒞                                                         -- denotation (δ1)
regularity : 𝒢 → Claim,  Claim = Aneristic | Eristic | Both | Neither   -- a cata; L3 label
culture    : Dist 𝒢                                                  -- leaf-data measure (§7)
```

**Properties of ⟦−⟧.**
- 𝒢 is countable and 𝒞 has cardinality 2^ℵ₀.
- ⟦−⟧ is neither surjective nor injective (Mu lossiness).
- Its image is Σ⁰₂-complete and null (B1).

**Named errors.**
- **Aneristic Illusion:** an interpretive-L3 failure. A fold over a grid is read as a fact about 𝒞.
- **Eristic Illusion:** the same failure with the opposite label. Calling something "random" also uses a grid.
- **Curse of Greyface:** imposing the Arrow condition (structure independent of data) on a domain whose structure adapts to data. Chaos is monadic globally, so Mu–Nu analysis applies per fiber at most.

---

## 2. Mu axioms — the grid side

**M1 Monoid.** `Seq` is associative with unit `Void`. 𝒢 is initial: every GridF-algebra receives exactly one cata. Every analysis of a doctrine is a fold.

**M2 Freeness (no opposites).** In 𝒢, `Seq x y = Void ⇒ x = y = Void`. A law has no inverse in the free grid.

**M2′ Hodge Hopf quotient.** v0's M1–M5 axiomatize exactly the group Hopf algebra kG: Δg = g⊗g, ε(g) = 1, S(g) = g⁻¹ **[V10]**.
- These axioms hold only after group completion, which identifies grids that cancel.
- Reading Δg = g⊗g: dogma copies exactly. *Greyface's quotient.*

**M3 The fold.** Let A be the grid algebra: linear combinations of grids, multiplied by `Seq`. The Aneristic projection f : A → A is an idempotent algebra endomorphism. In the free model k⟨a,b⟩, f(w) = w if w ∈ a\*, and 0 otherwise. Gridding twice is gridding once.

**M4 Splitting.** f splits the grid algebra as A = im f ⊕ ker f.
- **Hodge = im f**, a subalgebra.
- **Podge = ker f**, a two-sided ideal.

Any idempotent algebra endomorphism is automatically a Rota–Baxter operator of weight −1:

f(x)f(y) = f(f(x)y + x f(y) − x y)

This is a one-line identity, so the weight is not an empirical finding.

**M5 Closure and seriality.** GridF is closed, so every code-walker is exhaustive. Grids serialize; Chaos does not.

---

## 3. Nu axioms — the Chaos side (unchanged from v1)

**N1 Finality.** Every ObsF-coalgebra has a unique map into 𝒞.

**N2 Bisimulation is equality.** Two chaoses are equal iff no finite experiment separates them.

**N3 Productivity.** Every step emits a percept. A process that stalls is not chaos but void.

**N4 Stochasticity.** Steps land in `Dist`. No deterministic operation raises entropy above the sum of its inputs' entropies **[V8]**.

**N5 Inexhaustibility.** For every n there are distinct chaoses that agree to depth n. Chaos has a step map, not a counit: a counit ε ≡ 0 forces 𝒞 = 0 **[V2]**.

**N6 The chaos operation.** Where a binary ν coexists with order μ:
- ν(a, a) = e;
- ν is non-associative for *some* triple;
- ν has **no unit** (forced by B4).

Witness: (ℤ/5, +, −) **[V1, V7, V9]**.

---

## 4. β — Bialgebraic axioms: the distributive law between grid and Chaos (unchanged)

Setup: a GSOS law λ : GridF(Id × ObsF) ⇒ ObsF ∘ GridF\*. Where λ exists, ⟦−⟧ is compositional (Turi–Plotkin). **Discord is the defect of λ.**

**β1 Sequential congruence.** ⟦Seq g h⟧ ≈ ⟦g⟧ >>> ⟦h⟧, with a small interface handed off between the two.

**β2 Parallel regimes.** Each `Par` node is in exactly one of three regimes:
- *segregated*: ⟦Par g h⟧ ≈ ⟦g⟧ \*\*\* ⟦h⟧;
- *fused*: a designed λ-failure, i.e. emergence;
- *coupled*: h ≈ φ g.

Undecided nodes are defects. *The Apple of Discord* is the canonical fused node. It is a λ-defect, not an element of any algebra (see δ-ε).

**β3 Seed purity (L2).** A counter-dependent observer is not compositional **[V18]**. *Greyface's bureaucrat.*

**β4 Equivariance.** ν(μ(a,b), μ(a,c)) = ν(b,c) **[V9]**.

---

## 5. δ — Dialgebraic axioms: two collision algebras from one fold

f generates two structures. They are the two halves of a Koszul-dual pair of operads (Dias^! = Dend; Loday 2001). Both are axioms of v2. Neither implies the other.

### δ-L Loday dialgebra: the Leibniz collision (v1, confirmed)

Define x ⊣ y = x·f(y) and x ⊢ y = f(x)·y. They satisfy Loday's five axioms:

```
A1  (x⊣y)⊣z = x⊣(y⊣z)        A4  (x⊣y)⊢z = (x⊢y)⊢z
A2  x⊣(y⊣z) = x⊣(y⊢z)        A5  (x⊢y)⊢z = x⊢(y⊢z)
A3  (x⊢y)⊣z = x⊢(y⊣z)
```

All five hold on UT₃ symbolically and on k⟨a,b⟩ and UT₄ exhaustively **[V11, R1, R7]**. The identity (x⊣y)⊢z = x⊣(y⊢z) is *not* among them **[R2]**.

**δ-L1 Leibniz bracket.** [x,y] = x⊣y − y⊢x = x f(y) − f(y) x satisfies the *right* Leibniz identity:

[[x,y],z] = [[x,z],y] + [x,[y,z]]

with 0 failures, and [x,x] ≢ 0. The left identity fails, as it should for a non-Lie Leibniz algebra **[V12, R3]**.

### δ-D Dendriform: the Lie collision (new, from Le Chat v1.2–v1.3)

Aguiar's construction from a Rota–Baxter operator of weight −1 (M4) gives:

x ≺ y = x·f(y),  x ≻ y = f(x)·y − x·y

These satisfy the three dendriform axioms. Their star is the Baxter product:

x ⋆ y = x f(y) + f(x) y − x y

**δ-D1 Star splitting.** Write p = f(x), n = x − f(x), q = f(y), m = y − f(y). Then

x ⋆ y = pq − nm  **[R5]**

So (A, ⋆) ≅ im f × ker f as algebras, via n ↦ −n. Under ⋆, Hodge and Podge stop interacting: the star is the Sacred Chao drawn as a direct product.

**δ-D2 The completion family.** Fix ≺ = x f(y) and let ≻ range over b₁xy + b₂f(x)y + b₃xf(y) + b₄f(x)f(y). The dendriform completions are exactly two affine lines over k plus two isolated points **[R4]**:

| completion | b(t) | ⋆ |
|---|---|---|
| line 1 | (−1−t, 1+t, t, −t) | contains the Baxter star (t = 0) and f(x)f(y) (t = −1) |
| line 2 | (−t, 1+t, t, −1−t) | contains concatenation xy (t = −1) |
| point | (0, 0, 0, 0) | x f(y) |
| point | (0, 1, −1, 0) | f(x) y |

- The classification is due to Le Chat (v1.3). The exact variety over a general field is from **[R4]**.
- The "doublefold" product (x−f(x))(y−f(y)) is not a completion **[R4]**.

### δ-ε The ε-bialgebra substrate (new)

**δ-ε1 Existence.** The grid algebra carries Aguiar's infinitesimal bialgebra structure:
- on k⟨a,b⟩, Δ = letter deletion;
- on UT_n = kA_n, Δ = arrow deletion: Δ(E_ij) = Σ_{i≤m<j} E_im ⊗ E_{m+1,j}.

Both coproducts are coassociative and satisfy the ε-rule Δ(ab) = a₁⊗a₂b + ab₁⊗b₂ **[R8, R10]**. Deconcatenation is *not* coassociative **[R8]**.

**δ-ε2 The fold is not a coalgebra map.** For x = ab: Δf(x) = 0, but (f⊗f)Δ(x) = a⊗ε **[R9]**. Order does not respect deconstruction.

**δ-ε3 No Apple element.** The ε-coproduct has path-length degree −1. Every principal coproduct Δ_r(x) = x·r − r·x has degree ≥ 0. Hence the grid ε-bialgebra is not quasitriangular. Checked on UT₃: no r exists **[R11]**.

*The Apple of Discord is a λ-defect (β2), not a universal R-element.* Group Hopf algebras are not ε-bialgebras **[N5]**, so asking whether "M2′ is quasitriangular" conflates Drinfeld's notion with Aguiar's. In Drinfeld's sense kG is cocommutative, so R = 1⊗1 works trivially.

### δ-H Hagino: maps F X → G X, F ≠ G (unchanged)

**δ1 Ritual realization.** `enact = unfoldIter realize : Free GridF Slot → 𝒞`. Every slot must be filled within its deadline. *Operation Mindfuck* chooses `realize` so that β2 fuses where the audience expects segregation.

**δ2 Anticipation.** A believer is a biKleisli arrow `Cofree GridF Seen → Dist Next`. Open ambiguity must be rooted at the open span. *(This result is inherited from the music-axioms verification and was not re-checked here.)*

**δ3 Every Discordian is a Pope.** Doctrine is Church-encoded:

`Doctrine = ∀a. (P a a → a) → a`

No interpretation algebra is privileged.

---

## 6. Eris and the bounds

**E1 Principle of Eris, Leibniz form (v1).** For the δ-L bracket:
- the left annihilator is exactly im f;
- ker f lies in the right annihilator **[V13]**.

> **x collides with something ⇔ x ∉ im f.**

Only what carries chaos collides.

**E1⋆ Principle of Eris, Lie form (new).** By δ-D1, the annihilator of the ⋆-commutator is

ann [−,−]⋆ = Z(im f) ⊕ Z(ker f)  **[R5]**

| model | annihilator | evidence |
|---|---|---|
| free model | im f, since Z(ker f) = 0 | [R5] |
| UT_n | diag ⊕ k·E_1n, where E_1n is the centre of the nilradical (checked n = 3, 4, 5) | [R6] |

The apex of the Aftermath, E_1n, is the one chaotic element that collides with nothing.

The ⋆-annihilator depends on the stratum of the completion family **[N6]**:
- line 1 (generic member): im f;
- the commutative member f(x)f(y): everything;
- line 2: the scalars;
- the two points: 0.

(Classification due to Le Chat v1.3.)

**B1 Duality.**
- *Countability:* |𝒢| = ℵ₀ < 2^ℵ₀ = |𝒞|.
- *Measure:* ⟦𝒢⟧ is null under every i.i.d. non-degenerate product measure on 𝒞. Order covers a null set. This fails for non-identical products with Σ(1−pₙ) < ∞, which put mass on a pure flat.
- *Topology:* for flats and for hand-off grids alike, ⟦𝒢⟧ is countable and dense. So it is Σ⁰₂ and, by Baire, not G_δ. By Wadge's lemma (Kechris 22.10) it is therefore **Σ⁰₂-complete**. Le Chat's ruler reduction FIN ≤ P (v1.3) is a correct but unnecessary second proof.
- *Not ω-regular:* a countable ω-regular language is a finite union U₁v₁^ω ∪ … with finitely many periodic tails, while ⟦𝒢⟧ has unboundedly many. (The v1.1 argument "infinite ω-regular ⇒ uncountable" is false: a\*b^ω **[S10]**.)

**B2 Entropy.** H(μ(a,b)) ≤ H(a)+H(b). Strict increase under ν requires stochasticity (N4) **[V8]**.

**B3 Economy of grids.** At depth 8:
- minimal periods follow {1:2, 2:2, 3:6, 4:12, 5:30, 6:48, 7:82, 8:74};
- a hand-off grid (u, v) ↦ u·v^ω, at cost |u|+2|v|−1, names every word at cost ≤ 8, and beats the cheapest flat on 236 of 256 words **[S8, S9]**.

The hand-off is an economic upgrade over flats, not a topological one (B1).

**B4 Eckmann–Hilton exclusion.** Interchange plus units for both μ and ν forces μ = ν. So if order and chaos interchange, there is no neutral chaos **[V7]**.

**B5 Law of Fives (overproduction bound).** A five counts only if the check could have returned another number. §7 exhibits the Law as a channel fixed point.

---

## 7. The measure layer: cultures (new, from Le Chat v1.1–v1.2)

**C0 Setup.** A culture is a measure on a fixed grid set, updated by an observation channel. In the toy model:
- the grid set is the 256 eight-bit words;
- *Aneristic* means minimal period ≤ 3 (10 words);
- updates are exponential weights.

**C1 Fixed grid, flowing measure.** **[N3]**
- A neutral channel leaves the measure uniform: entropy ln 256 nats.
- A biased channel crystallizes the measure onto exactly its selected set: entropy → ln 10 (Aneristic), ln 246 (Eristic), ln 56 (five-bit words).

So the Eristic/Aneristic distinction is carried by the measure, not by the structure. *The Aneristic and Eristic Illusions are the two crystallizations.* The Law of Fives is the five-selective channel's fixed point. (Le Chat's v1.1 reported these entropies in bits; they are nats.)

**C2 No critical point.** In the deterministic update, P(Aneristic) is exactly a logistic in the learning rate β, with intercept ln(10/246) **[S12]**. A single-agent filter has no thermodynamic limit, so a phase transition is impossible by construction. *Greyface has a half-life, not a threshold.*

**C3 Rate law.** Under a channel of accuracy θ, the median fixation time is

t* ≈ (ln 99 − ln(10/246)) / (β(2θ−1) ln(θ/(1−θ)))

i.e. distance over drift (Le Chat v1.2 [W34]; constant checked in **[S12]**).

**C4 Impatience.** With truth in the minority basin, the probability of crystallizing on the wrong basin is gambler's ruin. On the lattice it is

P = ρ^{k₁}(1−ρ^{k₂}) / (1−ρ^{k₁+k₂}),  ρ = (1−θ)/θ,  k = ⌈distance/step⌉

This matches Le Chat's Monte Carlo at every point **[S11]**.
- In the diffusion limit it increases strictly with step size **[N4]**: haste reaches the wrong basin.
- On the lattice it is *not* monotone in general. P falls in k₁ and rises in k₂, and both shrink as the step grows. The largest reversal found is 0.073 **[N4]**.
- The principle is classical: timid play is optimal in a favourable game (Dubins–Savage), not a new theorem.

---

## 8. 𝔼, the seven-valued logic

**L1 Values.** The non-empty subsets of {T, F, M}, i.e. the seven clauses of the Principia epigraph **[V14]**. Connectives are the pointwise lift of Bochvar's tables (M infectious). A value is designated iff it contains T.

**L2 Collapse to four.** Designation depends only on the classical part of each atom's value **[S5]**. So 𝔼 is the four-valued matrix on {T, F, B, N}:
- LP tables on T, B, F;
- N infectious;
- designated {T, B}.

This is **Deutsch's S_fde**: the consequence relations coincide on every tested pair **[S14]**. (S_fde tables as given in Petrukhin 2022.)

**L3 Consequence.**

> **Γ ⊨_𝔼 φ ⟺ Γ ⊨_LP φ ∧ Var(φ) ⊆ Var(Γ)**  **[S3, S4]**

So 𝔼 is the containment (right variable-inclusion) companion of LP. LP has no antitheorems, so the companion needs no exceptional clause. Its neighbours, verified against their matrices:

| logic | consequence | evidence |
|---|---|---|
| PWK | ∃Γ′ ⊆ Γ. Var(Γ′) ⊆ Var(φ) ∧ Γ′ ⊨_CL φ — the left companion of CL (Ciuni & Carrara 2016) | [N1] |
| B₃ (Bochvar) | Γ is CL-unsatisfiable, or Γ ⊨_CL φ ∧ Var(φ) ⊆ Var(Γ) — the right companion of CL | [N2] |
| 𝔼 = S_fde | Γ ⊨_LP φ ∧ Var(φ) ⊆ Var(Γ) — the right companion of LP | [S3, S4] |

The rejected characterization "CL + containment" fails on disjunctive syllogism: (p∨q)∧¬p ⊭_𝔼 q **[S1, S2]**.

**L4 Theorems and antitheorems.** 𝔼 has no tautologies and no unsatisfiable premise sets. The all-{M} assignment kills every formula; the all-{T,F} assignment designates every formula.

**L5 Complexity.** 𝔼-consequence is coNP-complete (Le Chat v1.3).
- *Membership:* any finite matrix logic is in coNP.
- *Hardness:* γ is CL-unsatisfiable ⟺ γ ⊨_𝔼 ⋁_{p∈Var γ}(p∧¬p) **[S6]**.

By Cook–Reckhow, no proof system has polynomial-size proofs unless NP = coNP.

**Meta-1 Principle of Absurdity.** ¬(μ∨ν) ∧ (μ∧ν) is designated in exactly 4 of 49 assignments: those where both μ and ν are gluts **[V15, S7]**.

**Meta-2 Principle of Eris.** E1 and E1⋆.

**Meta-3 Principle of Popes.** δ3.

---

## 9. Worked examples

**The five-fingered hand.**
- μ(five fingers) has one reading.
- ν(five fingers) has 14 bracketings × 5! orders = 1680 readings **[V16]**. The divisibility of 1680 by 5 is a B5 artifact.
- In δ-L the anatomical hand is in im f and collides with nothing. A gesture's ker f component is what collides.
- In δ-D the gesture's ⋆-product with any hand is pq − nm. Hand and gesture compose only componentwise.

**The calendar.**
- 5 seasons × 73 days = 365, and 73 ≡ 3 (mod 5) **[V17]**.
- St. Tib's Day is the escape-hatch constructor, celebrated rather than mined.

---

## 10. Open questions

1. **Dist-valued denotations.** What are the Borel class and Wadge degree of the set of grid-named *measures* (weak-star topology)? This is the only open topology question.
2. **A fold that deconstructs.** Is there an idempotent Rota–Baxter(−1) map on the grid algebra that is also an ε-bialgebra morphism for letter deletion? f is not **[R9]**. If none exists, "order does not respect deconstruction" becomes a theorem.
3. **The Koszul pair.** Is δ-L recoverable from δ-D, or conversely, without f? δ-D1 shows ⋆ alone determines (im f, ker f) up to the sign twist. Does the Leibniz bracket?
4. **Proof theory of 𝔼 = S_fde.** Start from the S_fde / dS_fde literature (Deutsch; Fitting; Oller; Szmuc; Petrukhin 2022). Look for a calculus organised around availability.
5. **Real cultures.** Fit (θ, β) to belief-dynamics data. The predictions are closed forms (C3, C4). Multi-agent coupling, where the channel bias is endogenous, is where a genuine phase transition could appear.
6. **Arrow fibers, properly posed.** Find a social system whose global structure adapts to data but which yields an Arrow sub-system when either side is held fixed. Le Chat's v1.1 witness tested commutativity of maps on disjoint coordinates, which is a different property (framework calibration item 9).

---

## 11. Provenance

### From v1

All of v1's axioms stand. Le Chat's v1.2 withdrawals of [V11] and [V12] used non-axioms and are reversed **[R1–R3, R7]**.

### Adopted from Le Chat v1.1–v1.3

- the dendriform completion classification (δ-D2);
- the ⋆-annihilator, including the E₁₄ observation and the strata (E1⋆, generalized to UT_n, strata re-checked in [N6]);
- coNP-completeness (L5);
- Σ⁰₂-completeness (B1; simplified proof);
- the cost and period tables (B3);
- the measure layer (§7; units and claims corrected).

### Rejected

| claim | why | evidence |
|---|---|---|
| deconcatenation ε-bialgebra | not coassociative | [R8] |
| "the fold is an ε-morphism" | false for the genuine coproduct | [R9] |
| "UT₄: only Δ ≡ 0" | artifact of a degree-0 ansatz | [R10] |
| single-premise 𝔼 = CL + containment, and the Parry identification | corpus too small | [S1–S3] |
| "none of the three filters is redundant" | the CL filter is redundant | [S3, S4] |
| "infinite ω-regular ⇒ uncountable" | counterexample a\*b^ω | [S10] |
| "34 members" | search-range artifact | [R4] |
| R² = 1 as evidence | true by construction | [S12] |
| the Arrow-fiber witness | tests commutativity, not the Arrow condition | [S13] |

### Corrections to my own review (2026-10-10)

1. **Impatience monotonicity.** The review called the effect monotone in step size by Dubins–Savage. It is monotone only in the diffusion limit; on the lattice it can reverse **[N4]**.
2. **Bochvar and the right companion of CL.** The review used this as an argument without a confirmed source. It is now computationally verified **[N2]**. A literature statement of it is still unconfirmed here.

### Citation correction (supplied by Le Chat, checked here)

The PWK characterization is Ciuni & Carrara, *Characterizing Logical Consequence in Paraconsistent Weak Kleene*, in Felline, Ledda, Paoli & Rossanese (eds.), *New Developments in Logic and Philosophy of Science*.

| element | status |
|---|---|
| authors, title, book title, editors | confirmed |
| College Publications 2016, pp. 165–176 | not independently confirmed |
| the characterization itself | matches the PWK matrix on 17,136 pairs **[N1]** |

Le Chat's correction note pairs its second link label with the wrong URL. The labelled "PWK axiomatisation and algebraic analysis" (Bonzio, Gil-Férez, Paoli & Peruzzi, Studia Logica 2017) points to the Ciuni–Carrara page instead. The "extended version" it mentions is Ciuni & Carrara, *Semantical analysis of weak Kleene logics*, J. Applied Non-Classical Logics 29(1):1–36, 2019.

---

## 12. Sources

- Aguiar, *Infinitesimal bialgebras, pre-Lie and dendriform algebras*, arXiv math/0211074; *Infinitesimal Hopf algebras*, Contemp. Math. 267 (2000).
- Loday, *Dialgebras*, LNM 1763 (2001).
- Ciuni & Carrara, *Characterizing Logical Consequence in PWK* (in Felline et al., eds., *New Developments in Logic and Philosophy of Science*); *Semantical analysis of weak Kleene logics*, JANCL 29(1), 2019.
- Bonzio & Pra Baldi, *Containment logics: algebraic completeness and axiomatization*, Studia Logica (2021).
- Bonzio, Gil-Férez, Paoli & Peruzzi, *On Paraconsistent Weak Kleene Logic*, Studia Logica (2017).
- Petrukhin, *Normalisation for Some Infectious Logics and Their Relatives*, EPTCS 358 (2022) — S_fde and dS_fde tables.
- Kechris, *Classical Descriptive Set Theory* (1995), Thm 22.10 (Wadge).
- Calbrix, Nivat & Podelski, *Ultimately periodic words of rational ω-languages* (1993).
- Dubins & Savage, *How to Gamble If You Must* (1965).
- Turi & Plotkin, *Towards a mathematical operational semantics* (1997).
- *Principia Discordia* (1965).

---

*"All affirmations are true in some sense, false in some sense, meaningless in some sense, true and false in some sense, true and meaningless in some sense, false and meaningless in some sense, and true and false and meaningless in some sense."* — **Principia Discordia** (seven clauses, seven values; four of which suffice to decide anything)
