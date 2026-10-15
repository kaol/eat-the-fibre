# Discordianism: Open Questions Explored (v2.12)

*Answers questions from `discordianism-open-questions-explored-v2.11.md` §4: parameter complexity (Q4), quiver moduli (Q3, carried since v2.4) and the counting asymptotics (Q2).*

*Tags **[Sn]** are machine checks in `verify_v2_12.py` (folder `v2_12/`, zipped as `verify_v2_12.zip`). Labels: **proved**, **verified** (exact computation) or **numerical** (floating-point evidence on samples). Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Summary

| # | Question | Answer |
|---|---|---|
| v2.11 Q4 | Which parameters name a given μ ∈ S? | **Exactly the reals that compute μ: Name(μ) = cone(μ)** (proved; toy checks **[S1]**). The deadline function v2.11 suspected is always computable from μ itself, so it adds nothing. **Name(μ) can lack a least Turing degree** (proved via Miller's non-total continuous degrees and an explicit homeomorphism **[S5]**) |
| v2.4 → v2.11 Q3 | Quiver moduli | **Rigid modulo automorphisms** (numerical, **[S2, S3]**). On the sampled finite acyclic path algebras, every fold's local moduli dimension equals the dimension of its orbit under Δ-compatible automorphisms. v2.4's "continuous Kronecker moduli" were these orbits. **Retraction** of the moduli reading of v2.4 |
| v2.11 Q2 | Decay constant ρ of the compatible fraction | **Not determined by 26 terms** **[S4]**. Corrected fits give 0.39–0.41 and drift with cutoff and ansatz; v2.11's ≈ 0.42 (raw ratios) overstates it |
| — | Dihedral = C_L ∪ C_R; framed folds | Carried |

---

## 1. Parameter complexity: Name(μ) = cone(μ)

**Setting** (v2.10–v2.11). S = ⋃ B(U, ε^L) is the language-wide named set. For μ ∈ S let **Name(μ)** be the set of reals a such that some program of the coin language, with parameters computable from a, has law μ. Write μ ≤_T a when a computes μ, i.e. computes every μ(w) as a real, uniformly in w.

**Theorem 1 (proved; checked [S1]).** For μ ∈ S, with positive set U and lower bound ε^L,

  **Name(μ) = {a : μ ≤_T a}.**

The proof goes through an intermediate characterisation. Call d : U × ℕ → ℕ a **deadline** for μ if, whenever w ∈ U and μ(wb′) ≥ 2⁻ᵏ μ(w), the word wb′ is enumerated into U by stage d(w, k). This notion does not depend on the enumeration: changing enumerations changes d by a computable reindexing.

**Step 1: Name(μ) = {a : μ ≤_T a and a computes a deadline for μ}.**
- **Necessity.** Suppose a program with a-computable parameters names μ. Then μ ≤_T a: lower approximations of μ(w) come from exploring coin histories, and productivity makes the sum over a node's children converge to μ(w), so upper approximations follow. For a deadline, explore histories at w until the unexplored mass is below 2⁻ᵏ⁻¹ μ(w). Any child wb′ of mass at least 2⁻ᵏ μ(w) must by then have been emitted by an explored history, which certifies wb′ ∈ U. Return that stage.
- **Sufficiency.** This is v2.11's give-up-coin program, with the encoded real deadline T*_w replaced by an a-computable schedule. The give-up hazard stays at (δ/2)·2⁻ᵗ⁻¹ until d(w, k) and is then raised so that the survival probability halves, piece by piece, on each successive band 2⁻ᵏ⁻¹ ≤ c_w < 2⁻ᵏ. A branch with c_w ≥ 2⁻ᵏ is verified by stage d(w, k), while the survival ρ is still at least 2c_w. So the emission ratio c_w/ρ ≤ ½ is a legitimate probability, and P(emit b′) = c_w exactly. A branch outside U is never verified and gives up a.s.

**Step 2: a deadline is always computable from μ.** Given μ, define

  d(w, k) = the stage at which wb′ appears in U, if an approximation of μ(wb′)/μ(w) within 2⁻ᵏ⁻² is ≥ 3·2⁻ᵏ⁻²; and 0 otherwise.

- In the first case the true ratio is ≥ 2⁻ᵏ⁻¹ > 0, so wb′ ∈ U and the search halts.
- In the second case the true ratio is < 2⁻ᵏ, so the deadline condition is vacuous.

Hence d ≤_T μ, and Step 1 collapses to the theorem. ∎

**Examples.**
- **Halting-½ measure** (conditional ½ at every halting node). Here μ ≡_T 0′, since μ(w1) > 0 decides halting. So Name(μ) = cone(0′): exactly the reals that compute the halting problem.
- **Fair-coin halting measure** (v2.10's geometric-budget program). μ is computable, so Name(μ) = cone(0), i.e. every real. The computable schedule d(n, k) = k works directly **[S1]**.
- **Decidable support.** Name(μ) = cone(μ), as before.

**Machine checks [S1]** (exact rational arithmetic on a toy halting table):
- a deadline extracted from a naming program certifies every heavy branch;
- every dominating deadline reproduces all target masses exactly, while a non-dominating computable d misses the target on 11 halting nodes;
- the deadline computed from μ alone (by the threshold rule above) certifies every heavy branch and reproduces every branch mass, for both the halting-type and the fair-coin measures.

v2.11 §4 guessed that the relevant degree was that of the support's verification-time function. That guess is **superseded**: verification time is needed only on heavy branches, and μ already knows which branches are heavy.

### No least degree

**Theorem 2 (proved; homeomorphism checked [S5]).** There is μ ∈ S, with full (decidable) support, such that Name(μ) has **no least Turing degree**.

*Proof.* Let B = B(full tree, ¼): measures whose conditional probabilities all lie in [¼, ¾].
- The map x ↦ μ_x(w) = ∏_{v ⊏ w} (x_v or 1 − x_v) is a computable homeomorphism from the Hilbert-cube box [¼, ¾]^{2^{<ω}} onto B.
- Its inverse μ ↦ μ(w1)/μ(w) is computable, with the explicit modulus 2η·4^{|w|}/(1 − η·4^{|w|}), since μ(w) ≥ 4^{−|w|} **[S5]**.

So a computes μ_x iff a computes x. J. S. Miller (2004) showed that some points x of the Hilbert cube have **non-total continuous degree**: the set of reals computing x has no least Turing degree. For such x, Name(μ_x) = {a : x ≤_T a} has no least element, and B ⊆ S by v2.10. ∎

This fits a known pattern: probability measures are a standard home of non-total continuous degrees. Levin's *neutral* measures (every sequence is random relative to them) have non-Turing continuous degree (Day–Miller).

**Gloss.** *A Pope can name a world exactly when the Pope already knows the world. Telling the Church how long to wait is no extra burden, because knowing how much a path weighs already tells you whether waiting for it is worthwhile. Yet some worlds have no least-informed Pope: every Pope who knows them knows strictly more than some other Pope who also does.*

---

## 2. Quiver moduli: rigid modulo automorphisms

**Setting** (v2.4). Folds on the path algebra kQ of a finite acyclic quiver, with Δ the arrow-deletion coproduct. v2.4 found continuous families of folds on the Kronecker quiver and read them as moduli.

**Engine** (`quiver.py`). For a quiver it builds the path basis, the multiplication and Δ tensors, the fold residual, and the Lie algebra of Δ-compatible derivations. It then samples folds and compares two dimensions at each fold.
- The **local dimension** of the solution variety is the number of tangent directions whose displacement under an ε-step scales like ε. A direction counts only if its singular values at ε = 10⁻⁴ and 10⁻⁵ both exceed 0.2 and their ratio lies in (¼, 4).
- The **orbit dimension** is the rank of the infinitesimal action of the Δ-compatible automorphisms.

The ε-scaling test fixes a false positive: at ε = 10⁻³, curvature inflated the count.

**Result (numerical [S2]).**

| quiver | dim kQ | dim G | folds sampled (verifier; extended run) | max (local − orbit) |
|---|---|---|---|---|
| A₂ | 3 | 0 | 20 | 0 |
| A₃ | 6 | 0 | 23 | 0 |
| 1 → 2 ← 3 | 5 | 0 | 25 | 0 |
| Kronecker 2 | 4 | 2 | 25; 59 | 0 |
| Kronecker 3 | 5 | 6 | 25; 60 | 0 |
| 1 ⇉ 2 → 3 | 8 | 2 | —; 58 (`TRIALS=60 python3 quiver.py '1⇉2→3'`) | 0 |

With dim G = 0 the folds are isolated outright. (Sample counts are below TRIALS because some random starts fail to converge.)

For the m-Kronecker quiver, G acts on the arrow span by matrices with **row sums 1**, so dim G = m² − m: 2 for m = 2 and 6 for m = 3.

**Obstructed tangents [S3].** At Kronecker folds whose first-order tangent space is 3-dimensional, only 2 directions are genuine (displacement ∝ ε). The third displacement sits at a constant noise floor, so that direction is obstructed at second order. The 2 genuine directions are exactly the G-orbit.

**Conjecture 3.** Every Δ-compatible fold on a finite acyclic path algebra is rigid modulo Δ-compatible automorphisms: the fold moduli space is discrete.

**Retraction.** v2.4's "continuous Kronecker moduli" are automorphism orbits, not moduli.

**Contrast.** The free algebra k⟨a, u⟩ is the path algebra of a quiver with **loops**, and there the moduli are genuinely continuous: the identification theorem's R ∪ L (v2.8), and the one-sided dihedral families. Acyclicity kills paths of unbounded length, and with them the room for C-type series corrections (v2.9 Theorem 1), which are what carry the moduli. A proof should run through that mechanism: in a finite acyclic path algebra the cofree corrections terminate after finitely many cuts, and the remaining finite-dimensional fold equations should be shown to have no unobstructed deformations transverse to G.

**Gloss.** *Where history cannot loop, the split between order and chaos has no room to drift; what looks like drift is only the same split seen from a renamed vantage point. Moduli need loops: chaos needs room to repeat itself.*

---

## 3. Counting asymptotics: not determined

The terms a(1…26) of v2.11 (in `terms26.txt`) give the compatible fraction r(n) = a(n)/(2(n+2)^{n−1}).

**Result [S4].** Two ansätze for log r(n) were solved exactly on the last 3 (resp. 4) terms before a cutoff:
- log r(n) = A + n log ρ + α log n;
- the same plus β/n.

Their estimates of ρ drift with the cutoff, and the two ansätze disagree at the largest cutoff:

| ansatz | ρ estimates, increasing cutoff |
|---|---|
| A + n log ρ + α log n (cutoffs 22, 24, 26) | 0.4088, 0.4100, 0.4107 |
| … + β/n (cutoffs 22, 24, 26) | 0.3921, 0.4015, 0.4074 |

**Amendment.** v2.11's "≈ 0.42" came from raw successive ratios, which carry the n^β drift. The corrected estimate lies in roughly 0.39–0.41, and **ρ is not identified**. A closed form must come from the tree-function EGF: the composition structure of the DP of v2.11 Theorem 2 should give a singularity analysis at the branch point of T(z) = z e^{T(z)}.

**Gloss.** *The survivors thin out by about two fifths per step, but where exactly the thinning settles, twenty-six generations do not say.*

---

## 4. Proposed amendments (v2.11 → v2.12)

| amendment | content |
|---|---|
| **B1⁷ Names** | For μ ∈ S, the parameters that name μ are exactly the reals computing μ. Some μ, even with full support, have no least such degree (non-total continuous degrees) |
| **Q1¹² Quivers** | On finite acyclic path algebras, folds are rigid modulo Δ-compatible automorphisms (numerical; conjectured in general). Moduli need loops |
| **Retracted** | The v2.4 reading of Kronecker fold families as moduli; v2.11's ρ ≈ 0.42 (replaced by "0.39–0.41, unidentified") |

---

## 5. Open questions

1. **Dihedral = C_L ∪ C_R** (carried), and the all-order completion of the first-order family (e₀, F(M)) of v2.10.
2. **Quiver rigidity proof** (Conjecture 3), via the terminating cofree corrections. Then: quivers with loops, where the free-algebra moduli suggest a classification by **cycles**: is the moduli dimension determined by the loop structure?
3. **ρ** via singularity analysis of the tree-function EGF.
4. **Framed folds** (non-unital Hodge) in a C-type framework (carried).
5. **Degrees of named laws with given support.** Every continuous degree is the degree of some full-support μ ∈ S: by Theorem 2's homeomorphism and the computable rescaling [0,1]^ω ≅ [¼,¾]^ω. For a c.e. support U the picture varies: the halting-½ measure computes 0′, yet on the same halting tree the conditionals c_n = ε^{T(n)} (T = halting time, 0 if divergent) give a *computable* μ. Which continuous degrees occur for a given U, and is that set always upward closed?

---

## 6. Sources

- J. S. Miller, *Degrees of unsolvability of continuous functions*, J. Symbolic Logic 69 (2004), 555–584 (non-total continuous degrees in the Hilbert cube).
- A. R. Day and J. S. Miller, *Randomness for non-computable measures*, Trans. Amer. Math. Soc. 365 (2013) (neutral measures; measures of non-total continuous degree).
- M. S. Keane and G. L. O'Brien, *A Bernoulli factory*, ACM Trans. Model. Comput. Simul. 4 (1994); Ş. Nacu and Y. Peres, *Fast simulation of new coins from old*, Ann. Appl. Probab. 15 (2005).
- A. S. Kechris, *Classical Descriptive Set Theory*, Springer, 1995.
- H. Rogers, *Theory of Recursive Functions and Effective Computability*, MIT Press, 1987.
- P. Flajolet and R. Sedgewick, *Analytic Combinatorics*, CUP, 2009, §VII (tree-function singularity analysis).

*"To name a world, know it; some worlds have no simplest knower. Where history cannot loop, order and chaos split rigidly."*
