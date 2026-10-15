# Discordianism: Mu–Nu Axioms (v1)

*Refines v0 (Le Chat draft). Every formal claim tagged **[Vn]** is machine-checked by `verify_discordian.py` (Z3, SymPy, NumPy; 18/18 pass). Untagged formal claims are standard results. Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Diagnosis of v0

v0 placed order and chaos on two carriers with two binary operations and glued them with brackets. Under the Mu–Nu reading this is a vocabulary collision (framework §Cross-domain, item 9): μ and ν are not "order" and "disorder" but the **initial algebra** and **final coalgebra** of a functor. And Principia's own metaphysics refuses v0's split: order *and* disorder are both grids laid over an underlying Chaos. Hence:

- **Chaos is the Nu.** The territory: productive, observed only by unfolding, equal only up to bisimulation.
- **Hodge and Podge are both Mu.** The Aneristic grid (order) and the Eristic grid (disorder, "it's all random") are finite descriptions; the difference between them is a fold label, not a type.
- **The Sacred Chao** (pentagon and apple in one circle) is the pair (grids, Chaos) with the denotation map between them.

Five v0 axiom groups were inconsistent or ill-typed (§9). Every repair below is the smallest change that makes the group satisfiable, with a witness model.

"Dialgebra" in v0 conflated two notions; v1 uses both, labelled: **δ-L** (Loday: two products ⊣, ⊢, five axioms, induced Leibniz bracket) and **δ-H** (Hagino: maps F X → G X with F ≠ G).

---

## 1. Ontology

```
GridF atom r  =  Atom atom | Seq r r | Par r r | Neg r | Void      -- closed, 5 constructors
𝒢  =  μ GridF                                                      -- grids (Hodge ∪ Podge)
ObsF X        =  Obs → Dist (Percept × X)                           -- Auto Dist
𝒞  =  ν ObsF                                                       -- Chaos
⟦−⟧ : 𝒢 → 𝒞                                                         -- denotation (§5, δ1)
regularity : 𝒢 → Claim,  Claim = Aneristic | Eristic | Both | Neither   -- a cata; L3 label
```

- 𝒢 is countable (finite trees over countable atoms); 𝒞 has cardinality 2^ℵ₀ once |Percept| ≥ 2. So ⟦−⟧ is **not surjective**, and its image is null under the fair-coin measure on 𝒞: almost all Chaos has no grid. ⟦−⟧ is also **not injective** (Mu lossiness): many grids name the same Chaos.
- **Aneristic Illusion** = interpretive L3 failure: reading a fold over a grid as a fact about 𝒞.
- **Eristic Illusion** = the same failure with the opposite label. Calling something "random" is also a grid.
- **Curse of Greyface** = imposing the Arrow condition (structure independent of data) on a domain whose structure adapts to data. Chaos is monadic globally; Mu–Nu analysis holds only per fiber (framework: *monadic systems with Arrow fibers*).

---

## 2. Mu axioms — the grid side

**M1 Monoid.** `Seq` is associative and `Void` is its unit. 𝒢 is initial: for every GridF-algebra there is exactly one cata out of 𝒢. Every analysis of a doctrine is a fold. *(v0 M1–M2, kept.)*

**M2 Freeness (no opposites).** In 𝒢, `Seq x y = Void ⇒ x = y = Void`. A law has no inverse in the free grid. *(Replaces v0 M3.)*

**M2′ Hodge Hopf quotient.** v0's M1–M5 are consistent, and they axiomatize exactly the group Hopf algebra kG: Δg = g⊗g, ε(g) = 1, antipode S(g) = g⁻¹ **[V10]**. They hold on the group completion of 𝒢, which is a quotient: grids that cancel are identified. Read Δg = g⊗g as "order is copyable": dogma duplicates exactly. Read ε(g) = 1 as "observed order is one thing". The price of opposites is a quotient, so v0 M3 is admissible only post-quotient. *Greyface's quotient.*

**M3 Idempotent grid.** The Aneristic projection is an algebra endomorphism f with f∘f = f **[V11]**. Gridding twice is gridding once.

**M4 Closure.** GridF is closed, so every code-walker over grids is exhaustive. New kinds of grid need a new constructor, never an escape hatch (but see St. Tib's Day, §8.2).

**M5 Seriality.** Grids serialize; Chaos does not. Whatever can be written down is Mu.

---

## 3. Nu axioms — the Chaos side

**N1 Finality.** Every ObsF-coalgebra (S, γ) has a unique map S → 𝒞. Chaos is what any process becomes when observed forever.

**N2 Bisimulation is equality.** Two chaoses are the same iff no finite experiment separates them. This is the coinduction principle.

**N3 Productivity.** Every step emits a percept. A process that stalls is not chaos but void (`fix id`).

**N4 Stochasticity.** Chaos is effectful: the step lands in `Dist`. No deterministic operation raises entropy above its inputs' sum **[V8]**, so v0's "chaos increases entropy" needs its own randomness source.

**N5 Inexhaustibility.** For every n there are c ≠ c′ ∈ 𝒞 that agree to approximation depth n. No finite observation exhausts Chaos. *(Replaces v0 N5. A counit ε ≡ 0 together with the counit law forces 𝒞 = 0 **[V2]**. Chaos has a step map, not a counit.)*

**N6 The chaos operation.** Where a binary ν is wanted on a carrier that also has an order μ, the consistent form is:
- ν(a, a) = e, so self-confrontation yields the norm. This is v0's ∅ read as the unit e. v0 N2 ∧ N3 cannot both hold under any reading **[V1]**.
- ν is non-associative for some triple (∃, not v0's ∀).
- ν has **no unit** (forced by B4).

Witness model: (ℤ/5, μ = +, ν = −) **[V7, V9]**.

---

## 4. β — Bialgebraic axioms: the distributive law between grid and Chaos

**Setup.** An abstract GSOS law is λ : GridF(Id × ObsF) ⇒ ObsF ∘ GridF\*. A bialgebra (X, α, β) satisfies β∘α = ObsF α♯ ∘ λ_X ∘ GridF⟨id, β⟩. Where λ exists at a constructor, ⟦−⟧ is compositional there: bisimilar parts give bisimilar wholes (Turi–Plotkin). **Discord is the defect of λ.**

**β1 Sequential congruence.** ⟦Seq g h⟧ ≈ ⟦g⟧ >>> ⟦h⟧, with a small hand-off interface between the two. This is the TArrow homomorphism, `static (f >>> g) = static f <> static g`, read in reverse.

**β2 Parallel regimes.** Each `Par` node is designed in one of three regimes:
- *segregated*: ⟦Par g h⟧ ≈ ⟦g⟧ \*\*\* ⟦h⟧;
- *fused*: a designed λ-failure, where the emergent percept is not a product;
- *coupled*: h ≈ φ g for an observable φ.

An undecided node is a defect. *The Apple of Discord* is the canonical fused node. Three segregated claimants plus one inscription ("to the fairest") yield an observation that does not factor through the parts. Emergence and discord are the same λ-defect.

**β3 Seed purity (L2).** An observer of Chaos must be a function of its seed. A counter-dependent observer makes the whole observe differently from its parts **[V18]**. This is the call-counter L2 violation of the Anno corner, and also the replay-desync bug. *The bureaucrat's observation depends on how many forms were filed before yours.*

**β4 Equivariance.** A grid transformation the system relies on must lift to a bisimulation-equivariance. Concretely, order acting uniformly is invisible to chaos: ν(μ(a,b), μ(a,c)) = ν(b,c) **[V9]**. *(Replaces v0 Rule 1, μ(a, ν(b,c)) = ν(μ(a,b), μ(a,c)), which fails in the witness model **[V9]**.)*

---

## 5. δ — Dialgebraic axioms

### δ-L (Loday): the collision algebra

A Loday dialgebra (D, ⊣, ⊢) satisfies exactly five axioms:

```
L1  (x⊣y)⊣z = x⊣(y⊣z)        L4  (x⊣y)⊢z = x⊢(y⊢z)
L2  (x⊣y)⊣z = x⊣(y⊢z)        L5  (x⊢y)⊢z = x⊢(y⊢z)
L3  (x⊢y)⊣z = x⊢(y⊣z)                     — middle associativity
```

**Model.** Take an associative algebra A with an idempotent endomorphism f (M3). Set x⊣y = x·f(y) ("chaotic left meets gridded right") and x⊢y = f(x)·y. All five axioms hold **[V11]**. The canonical instance is A = UTₙ (upper-triangular matrices) with f the diagonal part. This is a Wedderburn–Malcev split A = S ⊕ N:
- S = im f is semisimple: the **Hodge**;
- N = ker f is the nilpotent radical: the **Podge**.

N is nilpotent (Nⁿ = 0): chaos composed with itself enough times vanishes. *The Aftermath.*

**δ-L1 Leibniz bracket.** Define [x, y] = x⊣y − y⊢x = x f(y) − f(y) x. It satisfies **right Leibniz**:

[[x,y],z] = [[x,z],y] + [x,[y,z]]

and [x, x] ≢ 0, so the algebra is genuinely Leibniz, not Lie **[V12]**. *(Replaces v0 D1/D2. Those were one equation stated twice **[V3]**, with a sign that contradicts standard Leibniz unless [y,[x,z]] = 0 **[V4]**.)*

**δ-L2 Middle associativity.** L3 is the honest form of v0's "Emergent Symmetry". *(v0 Rule 3, the full interchange law, is constrained by B4.)*

### δ-H (Hagino): maps F X → G X, F ≠ G

**δ1 Ritual realization.** A rite is a script with holes, `Free GridF Slot`. Enacting it is a dialgebra into Chaos: `enact = unfoldIter realize : Free GridF Slot → 𝒞`. Its productivity precondition is that every slot is filled within its deadline. *Operation Mindfuck* is a rite whose `realize` is chosen to make β2 fuse where the audience expects segregation.

**δ2 Anticipation.** A believer predicts as a biKleisli arrow `Cofree GridF Seen → Dist Next`. Such arrows compose only through a distributive law W∘M ⇒ M∘W. The traversal law fails duplication exactly when a non-root node is uncertain. *(This result is inherited from the music-axioms verification and was not re-checked here.)* Open ambiguity must be rooted at the open span.

**δ3 Every Discordian is a Pope.** Doctrine is a mixed-variance fixpoint, `P : Believer^op × Doctrine → Doctrine`, so a naive procata over it is ill-typed. The sound form is Church-encoded:

`Doctrine = ∀a. (P a a → a) → a`

A doctrine exists only as the collection of its interpretations, and each believer is an algebra. No algebra is privileged.

---

## 6. Eris and the bounds

**E1 Principle of Eris (refined Meta-2).** In the δ-L model:
- the left annihilator {x : ∀y [x,y] = 0} is **exactly** im f;
- ker f lies in the right annihilator **[V13]**.

> **∀x. x ∉ im f ⇔ ∃y. [x, y] ≠ 0.**

Only what carries chaos collides. Pure order is inert, and pure chaos cannot be struck. *(v0 Meta-2 fails at x = 0 in every bilinear bracket **[V5]**, and in this model at every element of S.)*

**B1 Duality (repaired v0 Bound 1).** v0's μ(a,b) with b ∈ Nu is ill-typed, since μ acts on Mu × Mu. The true duality is cardinal and measure-theoretic: |𝒢| = ℵ₀ < 2^ℵ₀ = |𝒞|, and ⟦𝒢⟧ is null (§1). Order and chaos do not cancel; order covers a null set.

**B2 Entropy.** H(μ(a,b)) ≤ H(a) + H(b) holds for every deterministic μ. H(ν(a,b)) ≥ H(a) + H(b) requires ν to be injective with a ⟂ b (equality), or stochastic (strict). The random search never exceeded the sum (max gap −0.0009) **[V8]**.

**B3 Finite generation.** 𝒢 is finitely generated over its atoms, and 𝒞 is infinite in depth (N5). *(v0 Bound 3, kept.)*

**B4 Eckmann–Hilton exclusion.** Interchange, ν(μ(a,b), μ(c,d)) = μ(ν(a,c), ν(b,d)), plus a unit for each operation forces μ = ν, commutative and associative. This contradicts N6 **[V7]**. So if order and chaos interchange, **there is no neutral chaos**. (ℤ/5, +, −) satisfies interchange with a unitless, non-associative ν **[V7]**.

**B5 Law of Fives (overproduction bound).** "The more you look, the more fives you find" is the framework's own calibration item 9 (vocabulary-driven correspondence), sanctified. **A five counts only if it survives a check that could have returned another number.** In this document:

| five | verdict |
|---|---|
| Loday's five axioms | structural: the variety needs all five |
| 1680 gestures ÷ 5 (§8.1) | artifact: 5! divides every count with ≥ 5 labelled parts |
| 5 seasons × 73 = 365 **[V17]** | cultural: chosen, not found |

---

## 7. Meta-axioms in 𝔼, the seven-valued logic

**Truth values.** The non-empty subsets of {T, F, M} (true, false, meaningless). There are exactly seven **[V14]**, and they are precisely the seven clauses of the Principia epigraph below.

**Connectives.** Pointwise lift of Bochvar's logic, with M infectious: v(φ∘ψ) = {a∘b : a ∈ vφ, b ∈ vψ}. A value is designated iff it contains T.

**Meta-1 Principle of Absurdity.** The formula ¬(μ∨ν) ∧ (μ∧ν) is classically unsatisfiable **[V6]**. In 𝔼 it is designated in exactly 4 of the 49 assignments: those where both μ and ν contain {T, F}. In 𝔼 the axiom therefore *means* "both principles are gluts". Ex falso fails, and the Boolean fragment stays classical **[V15]**.

**Meta-2 Principle of Eris.** This is E1.

**Meta-3 Principle of Popes.** This is δ3: no interpretation algebra is privileged.

---

## 8. Worked examples

### 8.1 The five-fingered hand (repaired)

v0 applied binary μ and ν to five arguments.
- By M1, μ(thumb, …, pinky) has **one** reading: the hand.
- ν is non-associative and non-commutative, so the gesture has **14** bracketings (Catalan₄) × 120 orders = **1680** readings **[V16]**.

v0's [μ(hand), ν(gesture)] bracketed two unary applications. The repair works in δ-L, with f the anatomical projection:
- The anatomical hand lies in im f, so it collides with nothing (E1).
- A gesture has a nonzero ker f component, so some context y gives [gesture, y] ≠ 0.

The salute means something only through its chaotic part.

### 8.2 The Discordian calendar

The calendar has 5 seasons × 73 days = 365 **[V17]**: Chaos, Discord, Confusion, Bureaucracy, The Aftermath. Since 73 ≡ 3 (mod 5), seasons and the five-day week drift against each other. The calendar is a grid with two incommensurable folds.

*St. Tib's Day* falls outside every season and every weekday. It is the calendar's **escape-hatch constructor**: a leaf the closed functor cannot otherwise place. The morphism-vocabulary pathology is here celebrated rather than mined.

---

## 9. Audit: v0 → v1

| v0 | verdict | evidence | v1 |
|---|---|---|---|
| M1–M2 | sound | — | M1 |
| M3 inverses | sound only post-quotient | [V10] | M2 (free) / M2′ (Hopf quotient) |
| M4–M5 | sound for group-likes | [V10] | M2′ |
| N1 ∀ non-assoc | too strong; witness model has associative triples (c = e) | [V7] | N6 (∃) |
| N2 ∧ N3 | **inconsistent** | [V1] | N6: ν(a,a) = e |
| N4 non-coassoc | vacuous once N5 holds | [V2] | N1–N5 (step map, no counit) |
| N5 ε ≡ ∅ | **forces 𝒞 = 0** | [V2] | N5 inexhaustibility |
| D1, D2 | duplicate; non-standard sign | [V3, V4] | δ-L1 right Leibniz |
| D3 | unmotivated | — | δ-L2 middle associativity |
| Bound 1 | ill-typed | — | B1 |
| Bound 2 (ν) | false for deterministic ν | [V8] | B2 + N4 |
| Rule 1 | fails in the witness model | [V9] | β4 equivariance |
| Rule 2 μ(a,ν(a,b)) = ε(b) | ill-typed (ε lands in scalars) | — | dropped; N6 carries its intent |
| Rule 3 interchange | consistent only if ν is unitless | [V7] | B4 |
| Meta-1 | classical contradiction | [V6] | 𝔼, glut reading [V15] |
| Meta-2 | fails at 0 | [V5] | E1 [V13] |
| Five-fingered hand | arity error | [V16] | §8.1 |

---

## 10. Open questions (sharpened)

1. **Measure.** N4 already makes Chaos Dist-valued. Which leaf-data measure on 𝒢 does a culture realize, and does the Eristic/Aneristic ratio flow (measure dynamics) while the grid is fixed?
2. **Real systems.** Social dynamics fail the Arrow test globally. Do they have Arrow fibers (fix institutions → behaviour is Arrow; fix behaviour → institutions are Arrow)? If so, the framework applies per fiber, and Greyface is the error of ignoring the diagonal.
3. **Topology.** 𝒞 with the prefix topology is a Stone space. Which subsets are ⟦−⟧-definable? Compare ω-regular languages (Büchi): is "grid-definable Chaos" exactly a Borel level?
4. **Hodge Hopf ↔ Leibniz.** Aguiar's infinitesimal bialgebras satisfy a Leibniz-type rule Δ(ab) = a₁⊗a₂b + ab₁⊗b₂. Does the δ-L bracket arise from such a structure on the grid algebra, unifying M2′ and δ-L?
5. **𝔼 proof theory.** Is the lifted-Bochvar 𝔼 axiomatizable by a finite Hilbert system, and is its consequence relation that of a known paraconsistent logic (e.g. Paraconsistent Weak Kleene on the M-free fragment)?

---

*"All affirmations are true in some sense, false in some sense, meaningless in some sense, true and false in some sense, true and meaningless in some sense, false and meaningless in some sense, and true and false and meaningless in some sense."* — **Principia Discordia** (seven clauses = the seven values of 𝔼)
