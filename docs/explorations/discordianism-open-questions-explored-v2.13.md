# Discordianism: Open Questions Explored (v2.13)

*Answers questions from `discordianism-open-questions-explored-v2.12.md` §5: degrees of named laws with a given support (Q5) and quiver rigidity (Q2).*

*Tags **[Tn]** are machine checks in `verify_v2_13.py` (folder `v2_13/`, zipped as `verify_v2_13.zip`). Labels: **proved**, **verified (exact)** (exact computer algebra over ℚ, here Gröbner bases and primary decomposition in Singular) or **numerical**. Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

*The default run takes minutes. `SLOW=1` adds the 3-Kronecker orbit test, the larger fibres, the numerical cross-check and the A₄–A₆ certificates. The logs of those slow runs are included (`verify_slow.log`, `lin_big.log`). Install the dependencies with `pip install -r requirements.txt` (sympy, numpy, scipy, mpmath, passagemath-singular).*

---

## 0. Summary

| # | Question | Answer |
|---|---|---|
| v2.12 Q5 | Which degrees do named laws with a given c.e. support U have? | **A trichotomy, with no dependence on how complicated U is** (proved; constructions checked **[T1, T2]**). Let b(U) be the number of branching nodes. Then the degrees are **{0}** if b(U) = 0, **all Turing degrees** if 0 < b(U) < ∞, and **all continuous degrees** if b(U) = ∞. In particular, every c.e. support, the halting tree included, carries a **computable** law in S |
| v2.12 Q2 | Quiver rigidity | **Vertex reduction theorem** (proved). Without parallel arrows, ker Δ = A₀, every fold preserves the vertex algebra, its vertex part obeys an **arrow rule**, and every vertex fold obeying the rule has a canonical fold extension P₀^{⊗•}. **The extension is unique**, by an exact linear certificate made possible by a conicity lemma, on 12 quivers (**[T5, T10]**), so the folds there are **finite and counted by the arrow rule**: A₄ has 100, D₄ has 182, the diamond has 88. **The counting sequence 2, 6, 23, 100, 563, … of v2.5–v2.11 is the number of folds of kAₙ** (equality verified for n ≤ 6, lower bound proved for all n). With parallel arrows, the 2- and 3-Kronecker fold varieties are **finite unions of Aut_Δ-orbits** (exact **[T4]**), and dim Aut_Δ = Σ m_ij(m_ij − 1) **[T7]** |
| — | Dihedral = C_L ∪ C_R; ρ; framed folds | Carried |

---

## 1. Degrees of named laws with a given support

**Setting.** S = ⋃ B(U, ε^L) (v2.10), and Name(μ) = cone(μ) (v2.12). For a c.e. tree U without dead ends, let

  **Deg(U)** = {continuous degree of μ : μ ∈ S, positive set of μ = U}.

A node w ∈ U is **branching** if both w0 and w1 lie in U. Let b(U) ∈ {0, 1, …, ∞} be the number of branching nodes.

**Theorem 1 (proved; checked [T1, T2]).**
- If b(U) = 0, then Deg(U) = {0}.
- If 0 < b(U) < ∞, then Deg(U) = all Turing (total) degrees.
- If b(U) = ∞, then Deg(U) = all continuous degrees.

Consequently, the possible sets Name(μ) for laws supported on U are exactly the cones above members of Deg(U). How hard it is to decide U plays no role.

*Proof.*

**(a) A computable law on every U: stage weights.** Fix an enumeration of U. Give each node v ∈ U the weight 2^{−s(v)}, where s(v) is the stage at which v enters U, and set

  c_w = weight(w1) / (weight(w0) + weight(w1)).

Then:
- **The S-bound holds.** c_w ≥ 2^{−s(w1)−1}, so μ ∈ B(U, ½^{L}) with L(w1) = s(w1) + 1, a partial computable function that is total on U.
- **The positive set is U.** Every conditional at a node of U is positive exactly on children in U.
- **μ is computable.** At stage s, the approximate weights are off by at most 2^{−s}. Once the first child of w has appeared (at stage m), the error in the ratio is at most 2^{m−s}. This bound is derived from |ab′ − a′b| ≤ 2^{−s}(a′ + b′) together with a + b ≥ 2^{−m}. Taking the product of the clipped lower bounds along w gives lower approximations of μ(w), so μ(w) is uniformly lower semicomputable. For a probability measure on Cantor space, uniform lower semicomputability is equivalent to computability, because μ(w) = 1 − Σ_{|v| = |w|, v ≠ w} μ(v).

**[T1]** checks this on random toy c.e. trees with exact rationals, using enumeration stages up to 40 per edge:
- the positive set is exactly U;
- the S-bound holds;
- the stagewise lower bounds never exceed μ;
- at stage 1000 they capture all but less than 2⁻⁵⁰⁰ of the mass at depth 9.

The halting tree therefore carries a computable law (for instance the stage-weight law, whose conditional at a halting node is about 2^{−T(n)}), alongside the halting-½ law, which has degree 0′.

**(b) Every continuous degree when b(U) = ∞.** The branching nodes form an infinite c.e. set. List them in discovery order as d₀, d₁, …. Take the computable μ from (a) and a point y ∈ [0,1]^ω, and define μ′ by changing the conditionals only at the listed nodes:
- at d_{2i}, put the code c′ = ¼ + y_i/2 ∈ [¼, ¾];
- at d_{2j+1}, put the displaced original conditional c_{d_j}.

Every original conditional survives somewhere: the one at d_j moves to d_{2j+1}, the one at d_{2j+1} moves to d_{4j+3}, and so on along each chain. Consequences:
- μ′ ≤_T μ ⊕ y, directly from the construction.
- y ≤_T μ′: read the even nodes.
- μ ≤_T μ′: every conditional at a node of U is recovered, and μ(w) is lower semicomputable from those conditionals once all prefixes of w have been enumerated into U.
- μ′ ∈ S: take ε′ = min(ε, ¼) and set L′(d_{2j+1}b) = L(d_j b) + 2.
- The positive set is still U.

So μ′ ≡ y. Every continuous degree is the degree of some point y ∈ [0,1]^ω (Miller 2004). **[T2]** checks the decoding exactly: 64 designated nodes, y and every displaced conditional recovered, all other nodes untouched.

**(c) Finitely many branching nodes.** μ is determined by the finitely many reals c_w at the branching nodes, together with the finite, non-uniform information of which nodes those are. A finite tuple of reals has a total degree, and every total degree occurs: put any real r ∈ (0,1) at one branching node and ½ at the others, with ε small enough.

**(d) No branching nodes.** U is a single path, c.e. without dead ends, hence computable, and μ is the Dirac measure on it. ∎

**Gloss.** *The shape of a world's open paths limits how secret the world can be, but not how simple it can be. Every support, however undecidable, admits a world any Pope can name. If the paths fork infinitely often, a world can hide any secret at all, including secrets that no simplest Pope holds. With finitely many forks it can hide only ordinary, total secrets; with none it hides nothing.*

---

## 2. Quiver folds: exact rigidity, and what the counting sequence counts

**Decoupling lemma (proved).** A linear map P on a path algebra A is a Δ-compatible fold iff A = H ⊕ Q with H = im P and Q = ker P, where:
- H is a subalgebra and a subcoalgebra (Δ H ⊂ H ⊗ H);
- Q is a subalgebra and a coideal (Δ Q ⊂ Q ⊗ A + A ⊗ Q).

*Proof.* Expand v2.12's residual. P² = P. P(xy) = xy on H × H. P(xy) = 0 on Q × Q. And ΔP = (P ⊗ P)Δ, which says Δh ∈ H ⊗ H on H and (P ⊗ P)Δq = 0 on Q. ∎

So **Folds = {(H, Q) ∈ X_H × X_Q : H ∩ Q = 0}**, an open subset of a product of two projective varieties.

### 2.1 Quivers without parallel arrows: reduction to the vertices

Let Q be finite and acyclic with no parallel arrows. Write A₀ = k^{Q₀} for the vertex algebra (functions on vertices).

**Theorem 2 (proved).**
1. **ker Δ = A₀.** Δ preserves path length (Aₘ → ⊕_{i+j=m−1} Aᵢ ⊗ Aⱼ). Deleting the i-th arrow of a path gives (prefix) ⊗ (suffix), and with no parallel arrows that term determines the path. So no cancellation is possible in degree ≥ 1. With parallel arrows, ker Δ also contains the differences of parallel arrows and their products **[T6]**.
2. **Every fold preserves A₀.** Δ e_v = 0, so Δ P(e_v) = (P ⊗ P)Δ e_v = 0, hence P(e_v) ∈ ker Δ = A₀. So P₀ = P|A₀ is a fold of k^{Q₀}, i.e. a Cayley forest (v2.5): P₀(e_v) = ±1_{S_v}.
3. **Arrow rule (necessary).** For an arrow α : s → t, ΔP(α) = P₀(e_s) ⊗ P₀(e_t) ∈ A₀ ⊗ A₀. Comparing degrees forces P(α) ∈ A₁ ⊕ A₀ and P₀(e_s) ⊗ P₀(e_t) ∈ Δ(A₁) = span{e_a ⊗ e_b : a → b}. So **supp P₀(e_s) × supp P₀(e_t) ⊂ Q₁**.
4. **Canonical extension (sufficient).** Embed A in the vertex-sequence algebra V = ⊕ₘ A₀^{⊗(m+1)}:
   - a path v₀ → ⋯ → vₘ maps to e_{v₀} ⊗ ⋯ ⊗ e_{vₘ};
   - the product glues the last factor of the left element to the first factor of the right one, using multiplication in A₀;
   - Δ is deconcatenation between factors.

   Then **T = ⊕ₘ P₀^{⊗(m+1)}** is a Δ-compatible fold of V:
   - T is idempotent and commutes with deconcatenation;
   - im T = ⊕ H₀^{⊗(m+1)} is closed under gluing, since H₀ is a subalgebra;
   - ker T = ⊕ Σᵢ A₀^{⊗i} ⊗ Q₀ ⊗ A₀^{⊗(m−i)} is closed under gluing. In a product of two kernel elements, at least one Q₀-slot survives as a slot of the product. The only exception is when both Q₀-slots are the glued ones, and then the glued slot is q·q′ ∈ Q₀, since Q₀ is a subalgebra.

   By the arrow rule, every vertex sequence in S_{v₀} × ⋯ × S_{vₘ} is a path, so T(A) ⊂ A. Hence **P = T|_A is a fold** of A. ∎

**Lemma (conicity; proved).** Let φ_λ multiply paths of length m by λᵐ. It is an algebra automorphism, and Δ has degree −1, so Δφ_λ = λ(φ_λ ⊗ φ_λ)Δ. A short computation then shows that φ_λ P φ_λ⁻¹ is again a Δ-compatible fold. This action fixes P₀ and acts on the unknowns by a_x ↦ λ^{−|x|} a_x.

Two consequences:
- The fibre over P₀ is a weighted cone with vertex the canonical fold. So **uniqueness of the extension ⇔ the canonical fold is isolated ⇔ (since vertex folds are finite) every fold is isolated**: rigidity and unique extension are the same statement.
- By induction on length, once a = 0 on shorter paths the fold equations at length m are **linear** in a_m.

The linear equations come from four sources:
- **idempotency:** a(Tz) = (1 − P₀)a(z);
- **H₀- and Q₀-multiples:** a(h·Tz) = (1 − P₀)(h·a(z)) and a(q·(z − Tz)) = P₀(q·a(z)), together with their right-hand versions;
- **H·H products:** a(Tx·Ty) = 0;
- **Q·Q products:** a((x − Tx)(y − Ty)) = 0, for |x|, |y| ≥ 1 with |x| + |y| = m.

**Theorem 3 (verified exact [T5, T10]).** On the quivers listed below, each vertex fold that passes the arrow rule extends to **exactly one** fold, the canonical one.
- **Linear certificate [T10].** For every arrow-rule vertex fold and every length, the linear system above has nullity 0.
- **Independent check [T5].** Solving the full polynomial fibre in Singular gives a single reduced point, a = 0, on the smaller quivers.

So |Folds(kQ)| equals the **arrow-rule count**:

| quiver | dim A | folds |
|---|---|---|
| A₁ | 1 | 2 |
| A₂ | 3 | 6 |
| A₃ | 6 | 23 |
| 1→2←3, 1←2→3 | 5 | 27 each |
| A₄ = 1→2→3→4 | 10 | **100** |
| 1→2←3→4 | 7 | 123 |
| D₄ (three arrows into one vertex) | 7 | 182 |
| 1→2→3 plus 1→3 | 7 | 21 |
| diamond | 10 | 88 |
| A₅ | 15 | **563** |
| A₆ | 21 | **3734** |

The arrow rule on the chain Aₙ is exactly v2.5's combinatorial compatibility condition. It reproduces v2.11's sequence 2, 6, 23, 100, 563, 3734 **[T8]**.

**Consequence: what the counting sequence counts.** Theorem 2.4 gives |Folds(kAₙ)| ≥ a(n) for all n, and Theorem 3 gives equality for n ≤ 6. The sequence carried since v2.5, for which v2.11 has a polynomial-time formula, is **the number of Δ-compatible folds of the path algebra kAₙ**, and not merely the number of compatible splittings of its vertex algebra.

**Conjecture 4 (unique extension).** For every finite acyclic quiver without parallel arrows, every fold is the canonical extension of its vertex part. Then:
- Folds(kQ) is finite, with one fold per arrow-compatible Cayley forest;
- quiver rigidity holds without exception in the absence of parallel arrows (Aut_Δ is trivial there by [T7]);
- a(n) = |Folds(kAₙ)| for all n.

**Cross-validation.** Three independent exact methods agree on A₂, A₃ and 1→2←3 (6, 23, 27):
- Schubert cells of X_H and X_Q, with every reconstructed P checked against the fold equations **[T3]**;
- open-chart pairs with primary decomposition **[T4]**;
- vertex fibres **[T5]**.

Folds found by v2.12's numerical sampler all lie in the exact lists **[T9]**.

### 2.2 Parallel arrows: finitely many orbits

With parallel arrows, ker Δ contains the differences of parallel arrows. These are primitive elements that behave like extra vertices, and along them the Δ-compatible automorphisms move folds in families.

**The automorphism Lie algebra (verified exact [T7]).** On 11 quivers, dim Lie Aut_Δ = Σ_{i,j} m_ij(m_ij − 1), where m_ij is the number of arrows i → j. The quivers are A₂, A₃, 1→2←3, K₂, K₃, K₄, 1⇉2→3, 1⇉2⇉3, the diamond, 1→2→3 with 1→3, and 1⇉2⇛3. The row-sum-1 substitutions α ↦ Σ_β g_{αβ} β among parallel arrows account for this count. They preserve Δ because Δα = e_s ⊗ e_t, so Δ(Σ_β g_{αβ} β) = (Σ_β g_{αβ}) e_s ⊗ e_t. Parallel *paths*, as in the diamond or in 1→2→3 versus 1→3, contribute nothing.

**Theorem 4 (verified exact [T4]).** For the 2- and 3-Kronecker quivers, the fold variety is a **finite union of Aut_Δ-orbits**.
- *Method.* Each pair of open charts carries the scheme X_H + X_Q + (t·det − 1). For every minimal prime P of dimension k > 0, the infinitesimal-orbit matrix must have generic rank k on V(P), so that an orbit is open and dense there. The test then recurses on the G-stable locus P + (k-minors).
- *Data.* The 2-Kronecker quiver has 46 chart pairs with fold-family dimensions 0, 1, 2. The 3-Kronecker quiver has 192 chart pairs with dimensions 0, 2, 4; even there no family exceeds dimension 4, although dim G = 6.
- *Not run to completion:* for 1⇉2→3 (dim A = 8, 12870 chart pairs) the exact test ran out of time. The numerical evidence of v2.12 stands.

**Conjecture 3′** (replaces Conjecture 3). On every finite acyclic path algebra, the folds form **finitely many Aut_Δ-orbits**.

**Gloss.** *Where history cannot loop and no two steps run in parallel, a split of order from chaos is decided entirely at the instants. How it treats the events between them is then forced, and the count of such splits is the sequence the Church has been computing since v2.5. Parallel steps are indistinguishable twins, and shuffling them turns each split into a family. Each family, though, is still one split seen through that shuffling.*

---

## 3. Proposed amendments (v2.12 → v2.13)

| amendment | content |
|---|---|
| **B1⁸ Degrees** | For a c.e. support U: Deg(U) = {0}, all Turing degrees, or all continuous degrees, according to whether U has 0, finitely many or infinitely many branching nodes |
| **Q1¹³ Quivers** | Fold = (H, Q), a complementary pair of a subalgebra-subcoalgebra and a subalgebra-coideal. Without parallel arrows, folds are vertex folds obeying the arrow rule, extended canonically (existence proved; uniqueness verified). With parallel arrows: finitely many Aut_Δ-orbits (exact for the 2- and 3-Kronecker quivers) |
| **P2¹³ Counting** | a(n) = number of folds of kAₙ (≥ proved; = verified for n ≤ 6) |

---

## 4. Open questions

1. **Dihedral = C_L ∪ C_R** (carried), and the all-order completion of (e₀, F(M)).
2. **Unique extension** (Conjecture 4). By conicity this is now pure linear algebra: show that the system (idempotency, H₀/Q₀-multiples, H·H, Q·Q) has nullity 0 for every arrow-compatible Cayley forest. The vertex-product equations alone do *not* suffice: on A₃ they leave nullity 2, and the path products are needed. A natural route is a triangularity argument along the forest, using the unitriangular block basis of the Cayley bijection. Then extend to parallel arrows (Conjecture 3′), where P(A₀) ⊂ ker Δ = A₀ ⊕ (parallel differences).
3. **Quivers with loops.** A loop gives k[u], and the classification there is trivial, character or dihedral (v2.2). What about one loop plus one arrow?
4. **ρ** via the tree-function EGF (carried).
5. **Framed folds** (carried).

---

## 5. Sources

- J. S. Miller, *Degrees of unsolvability of continuous functions*, J. Symbolic Logic 69 (2004).
- A. R. Day and J. S. Miller, *Randomness for non-computable measures*, Trans. Amer. Math. Soc. 365 (2013).
- H. Kredel and V. Weispfenning, *Computing dimension and independent sets for polynomial ideals*, J. Symbolic Comput. 6 (1988).
- P. Gianni, B. Trager and G. Zacharias, *Gröbner bases and primary decomposition of polynomial ideals*, J. Symbolic Comput. 6 (1988).
- W. Decker, G.-M. Greuel, G. Pfister and H. Schönemann, *Singular* (here via the passagemath-singular Python distribution).
- M. Aguiar, *On the associative analog of Lie bialgebras*, J. Algebra 244 (2001).
- I. Assem, D. Simson and A. Skowroński, *Elements of the Representation Theory of Associative Algebras*, Vol. 1, CUP, 2006 (path algebras).

*"Every support carries a world any Pope can name, and, if it forks often enough, a world that no simplest Pope can name. Without loops, order and chaos split in finitely many ways."*
