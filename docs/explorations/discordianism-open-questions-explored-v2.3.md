# Discordianism: Open Questions Explored (v2.3)

*Answers the open questions of `discordianism-open-questions-explored-v2.2.md` §7.*

*Tags **[Xn]** are machine checks in `verify_v2_3.py` (9/9 pass); [Un] and [Tn] refer to earlier rounds. Each claim is labelled **proved**, **verified** (exhaustive or symbolic up to the stated bound) or **conjectured**. Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Summary

| # | Question (v2.2 §7) | Answer |
|---|---|---|
| 1 | Is k[u] split only in the four obvious ways? | **No — v2.2's conjecture is false.** **Dihedral splittings** exist: k[u] = k[(u−a)²] ⊕ ((u−b)²−s)·k[(u−b)²] for any a ≠ b **[X2]**. Every splitting of k[u] extends run-wise to a fold that respects deconstruction **[X4]**. The endomorphism folds of v2.1 are conjugates of these extensions **[X5]** |
| 2 | Folds on UT_n beyond the monomial ones | **Exhaustive for n = 2, 3:** the folds correspond exactly to the splittings of the vertex algebra kⁿ that satisfy *arrow compatibility*. That is 6 of 8 for n = 2 and 23 of 50 for n = 3 **[X6]**. Necessity is proved for all n |
| 3 | Borel, or properly analytic? | **Borel (Σ⁰₂-complete)** for grids that are guarded or pushdown/recursive. **All of P(𝒞)** with exact oracle access to reals. **Open** only for unguarded Turing-complete grids with coin parameters, where the set is analytic **[X7]** |
| 4 | Fit to data | Carried. The published anchor is still Michard–Bouchaud (v2.2 §4) |
| 5 | Can any policy make a bureaucracy keep time? | **Not a symmetric one.** Every S_k-symmetric policy has a real spectrum inside [−1/(k−1), 1] **[X8]**. Breaking the symmetry to a cyclic order makes a data-driven clock, with eigenvalues at e^{±2πi/k} **[X9]**. *The Discordian calendar is such a symmetry-breaking clock* |

---

## 1. Splittings of k[u], and the deconstructing folds they generate

**Correction to v2.2 §2.** v2.2 conjectured that k[u] splits only as 0 ⊕ k[u], k[u] ⊕ 0, k ⊕ (u−c) and (u−c) ⊕ k. That is false.

**Lemma 1 (unit lemma; proved).** In any splitting k[u] = B₁ ⊕ B₂ into subalgebras, 1 lies in one of the summands.

*Proof.* Write 1 = e₁ + e₂ with eᵢ ∈ Bᵢ.
- If e₁ = t is nonconstant, then t − t² = (1−t) − (1−t)² is a nonzero element of t·k[t] ⊆ B₁ and of (1−t)·k[1−t] ⊆ B₂ **[X1]**. That contradicts B₁ ∩ B₂ = 0.
- If e₁ is a constant c ∉ {0, 1}, then 1 ∈ B₁ ∩ B₂, again a contradiction.

**Proposition 2 (structure; proved).** Take 1 ∈ B₁ and let C = k ⊕ B₂ (the unitization). Then:
- B₁ ∩ C = k and B₁ + C = k[u];
- B₂ is the kernel of a character of C.

So a splitting is the same thing as a complementary pair of unital subalgebras plus a point. If B₁ = k, then C = k[u] and B₂ = (u − c): this is the character splitting.

**Theorem 3 (dihedral splittings; verified).** For a ≠ b and any s,

  k[u] = k[(u−a)²] ⊕ ((u−b)² − s)·k[(u−b)²]

Checked by basis computation in every degree ≤ 2m, for several (a, b, s) **[X2]**.

*Why it works.* Both summands are fixed algebras of reflections, u ↦ 2a − u and u ↦ 2b − u. Together the two reflections generate the infinite dihedral group, whose invariants are only the constants: a polynomial invariant under translation by 2(b − a) is constant. A dimension count in each degree, (m+1) + m = 2m + 1, closes the sum. The simplest instance is k[u²] ⊕ (u² + u)·k[u² + u].

**Conjecture 4 (exhaustiveness).** Every splitting of k[u] is trivial, a character splitting, or a dihedral splitting.

*Partly proved.* A unital subalgebra ≠ k has finite codimension in some k[f], via Lüroth's theorem and finiteness of the normalization. So a nontrivial splitting needs two such k[f] and k[g] whose sum is all of k[u]. Computationally, that happens only for two quadratics with distinct centres; cancellations up to degree 30 were allowed **[X3]**. The missing step is a "density" bound that controls cancellation between k[f] and k[g]. This is a Ritt-type problem.

**Theorem 5 (run-wise extension; verified, with a proof sketch).** Fix a splitting of k[u] with projection P. Define P̂ on words over {a, u} by applying P to every maximal u-run and keeping the letters a:

  P̂(u^{i₀} a u^{i₁} ⋯ a u^{i_k}) = P(u^{i₀}) a P(u^{i₁}) ⋯ a P(u^{i_k})

Then P̂ is an idempotent RB(−1) fold that commutes with letter deletion. This was checked for the character, framed, dihedral and swapped-dihedral splittings **[X4]**.

*Proof sketch.*
- Letter deletion cuts at an a, so it never merges runs: Δ commutes with P̂ by construction.
- Concatenation merges exactly one pair of runs. The image of P̂ (all runs in B₁) and its kernel (some run in B₂) are both closed under concatenation.

**Theorem 6 (conjugation; proved and verified).** The Δ-compatible automorphisms of k⟨a,u⟩ are exactly

  G = {u ↦ αu + β, a ↦ a + h(u)}

*Proof.* Δθ(u) = 0 forces θ(u) ∈ k[u], and invertibility makes it affine. Δθ(a) = 1⊗1 forces θ(a) ∈ a + k[u].

Conjugating a fold by an element of G gives a fold. v2.1's endomorphism folds φ_{c,g} are exactly θ_g ∘ P̂_char ∘ θ_g⁻¹ **[X5]**.

**Conjecture 7 (classification).** Every fold on k⟨a,b⟩ that commutes with letter deletion is a G-conjugate of the run-wise extension of a splitting of k[u]. It holds for every fold found so far: all bigraded ones [U1] and all endomorphism ones [T4, X5]. Together with Conjecture 4, it would give the complete catalogue:

| fold | splitting of k[u] it extends | Hodge | Podge |
|---|---|---|---|
| identity, zero | trivial | everything / nothing | nothing / everything |
| identification / endomorphism folds | character | k[z], generated by one element | the ideal of differences |
| framed fold | character, with roles swapped (1 ∈ B₂) | framed events | contains Void |
| **dihedral folds** | dihedral | runs symmetric about a | runs symmetric about b |

**Gloss.** *Two mirrors suffice to split the world. What looks ordered about one centre and what looks ordered about another share only constants, and between them they account for everything. Hodge and Podge are not order and disorder but two incompatible symmetries.*

---

## 2. Path algebras: vertex splittings and arrow compatibility

**Reduction (proved).** On UT_n with arrow deletion, ker Δ is the diagonal D ≅ kⁿ (the vertex algebra). So every fold preserves D and restricts to a splitting of kⁿ into two subalgebras.

Subalgebras of kⁿ are spans of indicator functions of disjoint vertex sets, so these splittings are combinatorial. There are 8 for n = 2 and 50 for n = 3.

**Necessity of arrow compatibility (proved for all n).** For each arrow i → i+1, Δ(P α) = P(eᵢ) ⊗ P(eᵢ₊₁). This tensor has rank one and must lie in Δ(UT_n) ∩ (D⊗D) = span{e_k ⊗ e_{k+1}}. Hence at least one of the following holds:
- P(eᵢ) = 0;
- P(eᵢ₊₁) = 0;
- (P(eᵢ), P(eᵢ₊₁)) ∝ (e_k, e_{k+1}) for some arrow k.

**Theorem 8 (verified for n = 2, 3).** A splitting of the vertex algebra extends to a fold on UT_n **iff** it is arrow-compatible. This was established by exhaustive polynomial solving over all splittings: 6 of 8 for n = 2, 23 of 50 for n = 3 **[X6]**.

The six folds on UT₂ are:
- the four vertex supports (v2.2);
- the two **character folds**, P(x) = x_vv·1, which project onto the unit along a character. These are the path-algebra analogue of the character splittings of k[u].

**Correction to my own exploration.** A first pass reported 18 folds on UT₃. Its enumeration of splittings was incomplete: an ordering bug generated only 34 of the 50. The criterion matched in both runs.

**Gloss.** *On a path algebra, order is chosen vertex by vertex; an arrow tolerates only an order that does not tear it.*

---

## 3. The Borel question

**Setting.** Grids with real coin parameters and possibly unbounded loops. N is the set of named measures.

**Theorem 9 (proved).**
1. **Exact oracle access to reals** (grids can read the bits of a parameter): **N = P(𝒞).** Every measure is computable relative to some real, and fair coins together with an oracle for the cylinder probabilities sample it exactly. Expressivity saturates and the answer is trivially Borel.
2. **Guarded productive grids** (a bounded number of coin flips per emitted percept): every cylinder probability is a polynomial in the parameters. The law-map is therefore continuous, N is σ-compact, and N is **Σ⁰₂-complete** (hardness as in v2.1).
3. **Unguarded productive grids.** Each cylinder probability is a sum of nonnegative polynomials, hence lower-semicontinuous. On the productive region its complement, 1 minus the other cylinders of the same length, is lsc as well, so the law-map is **continuous on the productive region**. That region is G_δ.
4. **Pushdown / recursive Markov chain grids.** Termination and cylinder probabilities are least fixed points of polynomial systems. They are first-order definable over ℝ, hence **semialgebraic** in the parameters; for example, a biased walk returns with probability min(1, (1−p)/p) **[X7]**. The productive region is then semialgebraic, hence σ-compact, and N is **Σ⁰₂**.

**Open.** Unguarded Turing-complete grids with coin access. N is analytic, being the continuous image of a G_δ set. Whether some program's productive region and law-map produce a non-Borel image is open. Any example would need a productive region that is not σ-compact, such as the irrationals.

**Gloss.** *Productivity (N3) is what keeps the grid-named world Borel. Order that keeps talking stays describable; order allowed to deliberate forever might not be.*

---

## 4. Bureaucracies that keep time

**Theorem 10 (proved).** Let k forms share the lanes, with traffic depending only on which form holds the fast lane. Then the "fast form" process is an exact Markov chain (it is lumpable). If the restructuring policy is **S_k-equivariant** (it treats forms symmetrically), then:
- the chain's matrix commutes with every permutation matrix, so it lies in span{I, J};
- its spectrum is {1, 1 − kβ}, all real;
- stochasticity forces 1 − kβ ≥ −1/(k−1).

**No symmetric bureaucracy can keep time.** The most anti-persistent behaviour available is v2.2's walk on K_k. Checked for k = 3, 4, 5 under three symmetric policies **[X8]**:

| policy | minimum eigenvalue |
|---|---|
| most-filed | ≥ 0.87 (lock-in) |
| least-filed | reaches −0.478 for k = 3, against the bound −1/2 |
| second-most-filed | within the bound for every k |

**Theorem 11 (verified).** Break the symmetry to a cyclic order: give the fast lane to the *successor* of the most-filed form. The chain is then circulant. Its complex eigenvalues sit exactly at e^{±2πi/k}, and their modulus rises toward 1 as the observation window grows **[X9]**:

| k | modulus at M = 4, 8, 12 |
|---|---|
| 3 | 0.60, 0.78, 0.87 |
| 4 | 0.70, 0.87, 0.94 |
| 5 | 0.74, 0.90 (M = 4, 8) |

This is a **data-driven clock**.

**Gloss.** *A clock needs an arbitrary order on things that are otherwise alike. The Discordian calendar's five seasons (Chaos → Discord → Confusion → Bureaucracy → The Aftermath) are exactly such a C₅ order. A data-independent rotation is the purest Mu, a clock with no traffic at all. Greyface's symmetric bureaucracy can lock in or churn, but it cannot keep time. The Law of Fives is the symmetry breaking that lets order have a period.*

---

## 5. Proposed amendments (v2.2 → v2.3)

**M3‴ Fold catalogue.**
- Splittings of k[u] include the dihedral family. Every splitting extends run-wise to a deconstructing fold.
- G = {u ↦ αu + β, a ↦ a + h(u)} acts by conjugation.
- Conjectured: a fold is precisely a G-orbit of a run-wise extension, and the splittings of k[u] are only trivial, character and dihedral.

**P1 Path algebras.** Folds on UT_n correspond to arrow-compatible splittings of the vertex algebra. Necessity is proved; sufficiency is verified for n ≤ 3.

**B1‴ Borel.** Productivity plus guardedness, or pushdown structure, keeps N at Σ⁰₂-complete. Oracle access gives N = P(𝒞).

**Γ2 Clocks.** S_k-equivariant dynamics have real spectra bounded below by −1/(k−1). Periodic order needs a broken symmetry (C_k), as in the calendar.

---

## 6. New open questions

1. **Density lemma.** Prove that k[f] + k[g] = k[u] forces deg f = deg g = 2. This would settle Conjecture 4. Literature on Ritt decompositions and on polynomials s(f) − r(g) of small degree (Fried, Zieve) is the likely route.
2. **Conjecture 7.** Are all deconstructing folds G-conjugates of run-wise extensions? Equivalently: are the k[u]-valued corrections always coboundaries?
3. **Arrow compatibility for n ≥ 4**, and for quivers with cycles. There the free algebra (one vertex, two loops) and UT_n (a line) should be two ends of one theorem.
4. **A non-Borel grid world.** Construct an unguarded coin-access program whose productive region is the irrationals, or prove that productivity forces σ-compactness.
5. **Clocks under coupling.** Combine v2.1's coupled culture with v2.3's cyclic policy. Does social coupling sharpen the clock (|λ| → 1 faster), or lock it?

---

## 7. Sources

- Michard & Bouchaud (2005), *Theory of collective opinion shifts*: https://arxiv.org/pdf/cond-mat/0504079
- Etessami & Yannakakis, *Recursive Markov chains, stochastic grammars, and monotone systems of nonlinear equations*, J. ACM 56 (2009) — termination probabilities as least fixed points.
- Kechris, *Classical Descriptive Set Theory* (1995), §18 (Lusin–Souslin), §35 (Arsenin–Kunugui).
- Aguiar, *Infinitesimal Hopf algebras* (2000); Gubarev, *Rota–Baxter operators on unital algebras*, arXiv:1805.00723 — RB operators of nonzero weight as subalgebra splittings.

*"Two mirrors and a calendar: everything splits, and only the arbitrary keeps time."*
