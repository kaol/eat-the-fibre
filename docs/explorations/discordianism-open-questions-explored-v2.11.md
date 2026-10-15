# Discordianism: Open Questions Explored (v2.11)

*Answers questions from `discordianism-open-questions-explored-v2.10.md` §4, plus the counting sequence carried since v2.5.*

*Tags **[Rn]** are machine checks in `verify_v2_11.py` (folder `v2_11/`, zipped as `verify_v2_11.zip`). Labels: **proved**, **verified** (exact computation) or **conjectured**. Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Summary

| # | Question | Answer |
|---|---|---|
| v2.10 Q1 | Undecidable supports: is N = S? | **Yes. The language-wide named set is N = S, a Σ⁰₂-complete set, so Borel** (proved; construction checked on a toy halting table **[R4, R5]**). The trick is to **encode verification deadlines in a real parameter**. Real parameters are essential here: with computable parameters N is countable, and halting-type measures are not named |
| v2.5 → v2.10 | The counting sequence 2, 6, 23, 100, … | **A polynomial-time formula** that agrees with exhaustive enumeration for n ≤ 10 **[R1]**. The sequence now runs to n = 26 **[R2]**. No small holonomic recurrence fits **[R3]**. The compatible fraction decays by a factor ≈ 0.42 per step |
| v2.9–v2.10 | Dihedral = C_L ∪ C_R | Carried |
| — | Quiver moduli; framed folds | Carried |

---

## 1. Language-level Borelness: settled

**Recall** (v2.10). S = ⋃ B(U, ε^L): c.e. trees U, partial computable L total on U, rational ε, where

  B(U, ε^L) = {μ : μ = 0 off U, μ(v) ≥ ε^{L(v)} on U},

a closed set. The upper bound **N ⊆ S** holds because every finite coin history has positive probability **[R5]**. The lower bound was proved only for decidable U.

**Theorem 1 (proved).** **N = S.** Hence N is Σ⁰₂, and Σ⁰₂-complete by the Dirac argument of v2.10.

*Lower bound for arbitrary c.e. U.* Let μ ∈ B(U, ε^L). At a branching node w the conditional probability satisfies c_w ∈ [ε^{L(w1)}, 1 − ε^{L(w0)}], as in v2.10. The program π* works as follows.
1. **Find a known child.** Search U for some child b of w. This halts, since w ∈ U and U has no dead ends, and uses no coins.
2. **Verify the other child b′ under a budget.** Run the c.e. search for wb′ ∈ U. Before step t, flip a **give-up coin** with bias q_{w,t}. If it comes up heads, emit b.
   - The bias q_{w,t} is a Bernoulli-factory coin whose value depends continuously on a real parameter p₂. That parameter encodes the deadline T*_w: the verification time of wb′ if wb′ ∈ U, and 0 otherwise.
   - Before the deadline, q_{w,t} = (δ/2)·2^{−t−1}, tiny but positive. From the deadline on, q_{w,t} = ½. Here δ = ε^{L(wb)} is computable once b is found.
3. **On verification at step T.** The survival probability is ρ_T = ∏_{t<T}(1 − q_{w,t}) ≥ 1 − δ/2 > c_w. Emit b′ with probability c_w/ρ_T, otherwise emit b. Then P(emit b′) = c_w exactly.
4. **If wb′ ∉ U.** Verification never succeeds, the give-up coin fires geometrically after the encoded deadline, and the program emits b with probability 1.

**Productivity.** Each emission finishes a.s.

**Positive set.** It is exactly U, because b′ is emitted only after verification.

**Parameters.** The conditional probabilities come from a second real x, through the computable Hilbert-cube surjection of v2.10 [Q5], rescaled at each node once the bounds δ are known. Both reals are interleaved into **one** parameter y by a robust digit code: base-4 digits {1, 2}, so y never sits on a cell boundary. Every function the program needs is then continuous near y and bounded away from 0 and 1, so the univariate Keane–O'Brien factories suffice. ∎

The toy check uses a halting table with 8 halting and 4 non-halting nodes and exact rational arithmetic **[R4]**:
- with encoded deadlines, P(emit 1) equals the target exactly on halting nodes and is 0 elsewhere, and the remaining mass is below 2⁻¹⁰⁰;
- any fixed computable schedule, here q_t = 1/(t+2), is capped at ρ_T = 1/(T+1) and misses every target.

**Why real parameters are essential (proved).**
- With computable parameters there are only countably many named laws.
- More sharply, suppose a computable-parameter program named a measure on the halting tree with conditional probability ½ at every halting node. Then μ(w), μ(w0) and μ(w1) would be computable, by lower approximation plus productivity. Halting would then be decided by waiting for μ(w1) > 0 or μ(w0) > ¾μ(w).

The non-computable real carries exactly the information that verification needs: *when to stop waiting*.

**Gloss.** *The Church can name every world whose open paths it can eventually confirm, provided a Pope has whispered to it, in an uncomputable number, how long each confirmation may take. What results is Σ⁰₂: Borel, and no more complex than what a finite-state grid already names.*

---

## 2. The counting sequence: a polynomial-time formula

**Theorem 2 (proved; agrees with enumeration [R1]).** For the leaf-type rule of v2.6 Theorem 3, let:
- **E** = Podge leaves, **HL** = Hodge leaves, **P1** = Podge vertices with one leaf child, **B** = the rest;
- e, h, m = |E|, |HL|, |P1|;
- B_H, B_P = Hodge- and Podge-coloured B-vertices, with b_H, b_P their numbers.

The number of coloured forests with a prescribed type word, with each P1 matched to a designated HL child, is

  N(e, h, m, b_H, b_P) = Σ_{z_H, z_P} G(b_H, b_P, z_H, z_P) · Att₁(e+m, b_H, z_H) · Att₂(h−m, b_P, z_P),

where:
- **G** counts bipartite rooted forests on the B-vertices with exactly z_H Hodge and z_P Podge vertices childless inside B. It is computed by inclusion–exclusion from f(x, y) = (x+y+1)(x+1)^{y−1}(y+1)^{x−1}, the number of spanning trees of K_{1,x,y};
- **Att₁** attaches the E-vertices and P1-units to B_H ∪ {root}, with each of the z_H childless B_H vertices getting at least one;
- **Att₂** attaches the remaining HL-vertices to B_P ∪ {root}, with each of the z_P childless B_P vertices getting at least two (one leaf child alone would make it P1).

a(n) sums N over type words that obey the adjacency rule, with weight C(b, b_H). The rule "P1-runs have consecutive children" becomes a left-to-right DP. It marks HL-blocks and matches each P1-run with a marked block of the same length at the moment the second of the two appears, which counts every matching exactly once.

**Values** (n = 1…25; n ≤ 10 agree with exhaustive enumeration, the n = 10 run covering 2·12⁹ ≈ 10¹⁰ forests):

2, 6, 23, 100, 563, 3734, 28777, 246164, 2395445, 25665390, 303286713, 3899364820, 54534356615, 820099501690, 13236989251335, 227927850738224, …, a(25) = 343416730684520842761581941, a(26) = 9706864926713557355624856470.

**Shape.**
- No linear recurrence with polynomial coefficients of small order and degree fits a(n) or a(n)/n! **[R3]**. That is expected for sequences built on (n+2)^{n−1}, whose exponential generating functions involve the tree function and are not D-finite.
- The compatible fraction a(n)/2(n+2)^{n−1} falls by a factor near 0.42 per step: successive ratios run 0.427 … 0.422 at n = 25, still drifting. A closed constant is not identified.

**Gloss.** *Most ways of splitting a chain into order and chaos break some arrow; the survivors thin out geometrically, by a factor of about 0.42 per extra vertex.*

---

## 3. Proposed amendments (v2.10 → v2.11)

| amendment | content |
|---|---|
| **B1⁶ Borel** | For the Turing-complete coin language with real parameters, the named set is **N = S, Σ⁰₂-complete**. The Borel class never climbs, from finite-state grids to Turing-complete ones (v2.2 → v2.11) |
| **P2¹¹ Counting** | Polynomial-time formula for arrow-compatible splittings of path algebras; terms to n = 25 |

---

## 4. Open questions

1. **Dihedral = C_L ∪ C_R** (carried), and the all-order completion of the first-order family (e₀, F(M)) of v2.10.
2. **Asymptotics of the compatible fraction.** Identify the decay constant ≈ 0.42 via the tree-function EGF.
3. **Quiver moduli** (carried since v2.4); **framed folds** in a C-type framework.
4. **Parameter complexity.** Which Turing degrees of parameters suffice to name a given μ ∈ S? The deadline function T* suggests the degree of the support's verification-time function.

---

## 5. Sources

- M. S. Keane and G. L. O'Brien, *A Bernoulli factory*, ACM Trans. Model. Comput. Simul. 4 (1994); Ş. Nacu and Y. Peres, Ann. Appl. Probab. 15 (2005).
- A. S. Kechris, *Classical Descriptive Set Theory*, Springer, 1995 (Σ⁰₂, Wadge, Baire).
- R. P. Stanley, *Enumerative Combinatorics* Vol. 2, §5.3, and the spanning-tree count of complete multipartite graphs: n^{k−2}∏(n − nᵢ)^{nᵢ−1} (e.g. R. Lewis, *The number of spanning trees of a complete multipartite graph*, Discrete Math. 197/198 (1999)).

*"An uncomputable number tells the Church how long to wait; with it, every confirmable world is named, and the named worlds stay Σ⁰₂."*
