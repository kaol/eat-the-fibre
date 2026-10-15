# Discordianism: Open Questions Explored (v2.10)

*Answers questions from `discordianism-open-questions-explored-v2.9.md` §4. It finally takes up **language-level Borelness**, open since v2.4.*

*Tags **[Qn]** are machine checks in `verify_v2_10.py` (folder `v2_10/`, zipped as `verify_v2_10.zip`). Labels: **proved**, **verified** (exact computation in a stated window) or **conjectured**. Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Summary

| # | Question | Answer |
|---|---|---|
| v2.9 Q2–Q3 | Full one-sided dihedral family; why roots? | **First-order classification** (proved, checked **[Q1, Q2]**). The edge corrections ψ(h) = C(h·a·q₀) solve Mψ − ψM = Φ((1−p)(q₀·))ψ(1), where M = p∘(q₀·) on A₀/Q₀ ≅ k·1̄ ⊕ k[M]·ū. The commuting solutions are indexed by a constant e₀ and a series F(M). **Roots:** root deviation is the case F = 1, e₀ = −r; it solves the equation because ev_r kills Q₀, i.e. r is a root. It is also the **only** exact fold supported on one-a words **[Q3]** |
| v2.4 Q5 | Is the language-wide named set N Borel? | **Yes on decidable supports, and Σ⁰₂-complete there** (proved; key lemma checked **[Q4, Q5]**). Globally, **N lies inside an explicit Σ⁰₂ set S**, with N ∩ Dec = S ∩ Dec. Per-program non-Borel sets wash out at the language level. What remains open is measures whose support is c.e. but undecidable |
| v2.9 Q1 | Prove dihedral = C_L ∪ C_R | Carried |
| — | Counting sequence, quiver moduli, framed folds | Carried |

---

## 1. The dihedral one-sided family at first order

Let M = p∘(q₀·) act on H₀ ≅ A₀/Q₀, and Φ(Σ βₙ q₀ⁿ) = Σ βₙ Mⁿ.

**Theorem 1 (proved; checked [Q1, Q2]).**
1. **Module structure.** As a k[M]-module, A₀/Q₀ ≅ k·1̄ ⊕ k[M]·ū. Here M kills the class of 1, because q₀ ∈ Q₀, and acts freely on the class of u: A₀ is free of rank 2 over k[q₀], with basis 1, u.
2. **The equation.** The first-order condition on ψ(h) = C(h·a·q₀) is v2.6's Q₀-bimodule condition rewritten:

   **M ψ(h) − ψ(M h) = Φ((1−p)(q₀h)) · ψ(1).**

   Its solutions extend the choice of ψ(1) by the **commuting family**: ψ(1̄) = 0 and ψ(ū) = e₀·1̄ + F(M)·ū, for an arbitrary constant e₀ and series F. Then ψ(f(M)ū) = f(0)e₀·1̄ + f(M)F(M)·ū.
3. **Agreement with the engine.** In the window, the solution space of this equation coincides with the all-a-degree engine's one-sided tangent: 5 = 5 dimensions **[Q2]**.

**Why roots (proved).** Root deviation ψ_r(x) = x − ev_r(x)·1̄ is the commuting solution with F = 1 and e₀ = −r. To satisfy the equation it needs ev_r(p(q₀h)) = 0 for all h, i.e. ev_r ∘ M = 0. Since p(q₀h) ≡ q₀h modulo Q₀, that holds exactly when ev_r kills Q₀, that is, **when r is a root of (u−b)² − s**.

**Exactness (verified [Q3]).** With C supported on one-a words only, the exact folds in this family are precisely the root deviations. Every direction involving M (F ≠ const) needs corrections on words with more a's; v2.9 showed generic one-sided directions do integrate.

**Gloss.** *Chaos can lean on an event's edge in infinitely many ways to first order, one for each power series in its own action. Only the lean measured against a point where chaos vanishes is exact as it stands; every other lean has to recruit longer histories.*

---

## 2. Language-level Borelness

**Setting** (v2.3–v2.4). Programs of a Turing-complete coin language have real parameters p ∈ [0,1]^k (coins of bias pᵢ) and fair coins. N is the set of laws μ_p over all programs and all productive p. Per program, the named set can be non-Borel (v2.4 [Y5]).

**Key lemma (proved; checked [Q4]).** Fix a program with p ∈ (0,1)^k. Boundary parameters reduce to other programs with deterministic coins.
- **(a) Support.** Every finite raw-coin history has positive probability. So the **positive set** U = {w : μ_p(w) > 0} is the set of words emitted by some finite history. It is independent of p and c.e.
- **(b) Lower bound.** μ_p(w) ≥ ε^{L(w)}, where ε = min(½, minᵢ min(pᵢ, 1−pᵢ)) and L(w) is the length of the shortest emitting history, a partial computable function, total on U.

Checked with exact rational arithmetic on random productive finite-state coin programs, for all words of length ≤ 4 and three values of p.

**Definitions.**
- For a tree U and bounds δ = (δ_v), let **B(U, δ)** = {μ : μ(v) = 0 off U, and μ(wb) ≥ δ_{wb}·μ(w) at every branching node w}. These are closed conditions, so B(U, δ) is **closed**. Every μ in it has positive set exactly U.
- Let **S** = ⋃ B(U, ε^{L}), over c.e. trees U, partial computable L and rational ε: a countable union of closed sets, hence **Σ⁰₂**.
- Let **Dec** be the set of measures whose positive set is decidable.

**Theorem 2.**
1. **Upper bound.** N ⊆ S. By the lemma, μ_p(wb) ≥ ε^{L(wb)} ≥ ε^{L(wb)}μ_p(w).
2. **Lower bound on decidable supports.** For every decidable tree U and computable rational δ, **B(U, δ) ⊆ N**. Construction:
   - **Parameter.** A single parameter x ∈ [¼, ¾] is mapped by a computable continuous surjection onto the product of the boxes [δ_{w1}, 1 − δ_{w0}], one coordinate per branching node. Use the iterated Hilbert curve x ↦ (H₁x, H₁H₂x, H₁H₂H₂x, …) **[Q5]**, rescaled.
   - **Emission.** At a branching node w (decidable), the program emits 1 using a **Bernoulli factory** for the coordinate f_w(x). This is continuous and bounded away from 0 and 1, so it is exactly simulable with a.s. termination (Keane–O'Brien; constructively Nacu–Peres). Non-branching nodes are deterministic.
   - **Coverage.** The law has conditional probabilities f_w(x), and as x varies these cover the whole box product.
3. **Conclusion.** **N ∩ Dec = S ∩ Dec = ⋃_{U decidable, δ computable} B(U, δ)**. This is a **Σ⁰₂** set, in particular Borel.
4. **Completeness.** N ∩ Dec is **Σ⁰₂-complete**. A named Dirac measure δ_x has c.e. positive set {prefixes of x}, so x is computable, and conversely every computable x is named. So N ∩ {Diracs} is the countable dense set {δ_x : x computable} in the perfect closed set of Diracs, which is not G_δ by Baire. Hence N ∩ Dec is not Π⁰₂, and a Σ⁰₂ set that is not Π⁰₂ is Σ⁰₂-complete by Wadge.

**What this settles.** The per-program non-Borel sets of v2.4 consist of Bernoulli(p₁) laws: full support, which is decidable, with conditional probabilities in boxes. They sit inside B(full tree, δ) ⊆ N, so **the language-wide union washes them out**, exactly as v2.4's gloss guessed.

**What remains open.** Measures with c.e. but *undecidable* support. Such measures are named. Example: at step n, sample a geometric budget s, run machine n for s steps, and emit a fair bit if it halted (else 0). The positive set contains w1 iff machine n halts. For these, the lower-bound program would have to decide branching, which it cannot. So the question becomes: is N ∩ M_U = S ∩ M_U for every undecidable c.e. U? If yes, then N = S, and N is Σ⁰₂.

**Gloss.** *Each Pope may name a world no Borel description captures, but the Church's collective naming is Σ⁰₂ wherever it can tell which paths are open. The only residual mystery lives where even the Church cannot decide which paths exist.*

---

## 3. Proposed amendments (v2.9 → v2.10)

| amendment | content |
|---|---|
| **B1⁵ Borel** | The language-wide named set lies in an explicit Σ⁰₂ set S (c.e. support plus computable exponential lower bounds). On decidable supports it equals S and is Σ⁰₂-complete. Open only for undecidable supports |
| **M3¹⁰ Dihedral one-sided** | First order: Mψ − ψM = Φ(·)ψ(1). The commuting solutions are (e₀, F(M)); root deviations are F = 1, e₀ = −r with ev_r(Q₀) = 0 |

---

## 4. Open questions

1. **Undecidable supports.** Is every μ ∈ S with undecidable c.e. support named? Equivalently, is N = S?
2. **Dihedral = C_L ∪ C_R** (carried from v2.9), and the all-order completion of the commuting first-order family (e₀, F(M)).
3. Carried: the counting sequence; quiver moduli; framed folds in a C-type framework.

---

## 5. Sources

- M. S. Keane and G. L. O'Brien, *A Bernoulli factory*, ACM Trans. Model. Comput. Simul. 4 (1994).
- Ş. Nacu and Y. Peres, *Fast simulation of new coins from old*, Ann. Appl. Probab. 15 (2005).
- A. S. Kechris, *Classical Descriptive Set Theory*, Springer, 1995: Hahn–Mazurkiewicz (§4), Baire category, Wadge's lemma (§21–22).
- H. Rogers, *Theory of Recursive Functions and Effective Computability*, MIT Press, 1987: c.e. sets.

*"The Church names Σ⁰₂ where it can see the paths; beyond that, it is still praying."*
