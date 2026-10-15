# Discordianism: Open Questions Explored (v2.18)

*Takes up question 3 of `discordianism-open-questions-explored-v2.17.md` §5, parallel arrows (Conjecture 3′), for the generalized Kronecker quivers. These are the acyclic shadow of the original k⟨a,u⟩ framework.*

*Tags **[Yn]** are machine checks in `verify_v2_18.py` (folder `v2_18/`, zipped as `verify_v2_18.zip`; `pip install -r requirements.txt` needs only numpy and sympy; `SLOW=1` adds K₂/𝔽₇ and K₃/𝔽₃). Labels: **proved**, **verified (exact)** or **conjectured**. Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Summary

| # | Question | Answer |
|---|---|---|
| v2.17 Q3 | Conjecture 3′ for parallel arrows | **Proved for every generalized Kronecker quiver K_m** (two vertices, m parallel arrows; any field of characteristic ≠ 2). The Δ-compatible folds form **exactly 14m − 8 orbits** of Aut_Δ = {g ∈ GL_m : σ∘g = σ} (Theorem 1). The classification is checked by an exhaustive census over finite fields: K₁, K₂ over 𝔽₃, 𝔽₅, 𝔽₇, and K₃ over 𝔽₃ (34 orbits, 1300 folds) **[Y1, Y2]** |
| — | Point counts | The fold variety has **N_m(q)** points over 𝔽_q, an explicit sum of Gaussian binomials, e.g. **6q² + 2q + 4** for m = 2 and **14q⁴ + 6q³ + 4** for m = 3. Its **dimension is ⌊m²/2⌋** **[Y3]**, matching v2.13's exact chart dimensions (2 and 4). K₁ = A₂ recovers the 6 folds of v2.13 |

---

## 1. Setting

K_m has vertices 1, 2 and arrows α₁, …, α_m: 1 → 2. Its path algebra is

  A = k e₁ ⊕ k e₂ ⊕ V,  V = span(α_i),  with e₁ v = v = v e₂, all other products of e's and V zero, and V² = 0.

Arrow deletion gives Δ(e_i) = 0 and Δ(v) = σ(v)·e₁ ⊗ e₂, where σ is the sum of coefficients. Let D = ker σ ⊂ V, the span of differences of parallel arrows; it is the Δ-primitive part that broke vertex reduction in v2.13. The group G = {g ∈ GL(V) : σ∘g = σ} acts by Δ-compatible algebra automorphisms fixing e₁, e₂. Its Lie algebra is the whole Δ-compatible derivation algebra (dim m² − m, v2.13 T7). G maps onto GL(D) ⋉ (affine shifts), so it acts transitively on configurations of subspaces with prescribed position relative to D.

---

## 2. The classification

**Theorem 1 (proved; char k ≠ 2).** A fold of A is a pair (H, Q) of complementary subalgebras with ΔH ⊂ H ⊗ H and ΔQ ⊂ Q ⊗ A + A ⊗ Q (v2.13). Up to G, they are exactly the following **14m − 8** types.

*Two Δ facts.*
- **(i)** If some h ∈ H has σ(h) ≠ 0, then e₁ ⊗ e₂ ∈ H ⊗ H, so e₁, e₂ ∈ H.
- **(ii)** If some q ∈ Q has σ(q) ≠ 0, then e₁ ⊗ e₂ ∈ Q ⊗ A + A ⊗ Q, so e₁ ∈ Q or e₂ ∈ Q.

Since A = H ⊕ Q contains an arrow, exactly one of the following two cases occurs.

**Case I: A₀ ⊂ H.** Then H = A₀ ⊕ W with W ⊂ V and W ⊄ D.
- Coideal forces Q ⊂ A₀ ⊕ D, so Q = {u + f(u)} over a complement U ⊂ D of W, with f: U → A₀.
- The subalgebra condition reduces to f₁(u)·f₂(u′) = 0, so f lies entirely in the e₁ slot or entirely in the e₂ slot.
- Configurations (W, U) with dim W = k are a single orbit for each k = 1…m. The functional f is either 0 or nonzero, in one of the two slots.
- **Count: m + 2(m − 1) = 3m − 2 types.**

**Case II: H ⊂ A₀ ⊕ D.** Then e₁ ∈ Q or e₂ ∈ Q, and both Δ conditions hold automatically, so only the algebra conditions remain. If e_i ∈ Q, multiplying by e_i splits Q = Q₀ ⊕ Q_V, with Q₀ ∈ {k e_i, A₀}.
- **Q₀ = A₀.** H is the graph over W ⊂ D of a functional into the e₁ slot or the e₂ slot (the same f₁f₂ = 0 argument), with V = W ⊕ Q_V. **3m − 2 types.**
- **Q₀ = k e₁.** By the image of H in A₀:
  - **ke₂:** H = k(e₂ + d₀) ⊕ W, where d₀ ∈ D/W is 0 or nonzero. **2m − 1 types.**
  - **k·1:** H = k·1 ⊕ W. Here (1 + d₀)² ∈ H forces d₀ ∈ W, using char ≠ 2. **m types.**
  - **A₀:** H = span(e₁ + d₁, e₂ − d₁) ⊕ W, with V = kd₁ ⊕ W ⊕ Q_V and d₁ ∈ D ∖ W. **m − 1 types.**
- **Q₀ = k e₂.** The mirror of Q₀ = k e₁. **4m − 2 types.**

**Total: (3m − 2) + (3m − 2) + 2(4m − 2) = 14m − 8.** ∎

**Consistency.** For m = 1 we have D = 0 and K₁ = A₂, and the 6 types are the 6 folds of A₂ found in v2.13. So the count specializes correctly to the acyclic, no-parallel theory of v2.16.

---

## 3. Verification and point counts

**Census [Y1].** Over 𝔽_p I enumerate all subalgebras satisfying the intrinsic Δ conditions (i)/(ii), all complementary pairs, and the G(𝔽_p)-orbits:

| quiver / field | folds | \|G\| | orbits | 14m − 8 |
|---|---|---|---|---|
| K₁ / 𝔽₃ | 6 | 1 | 6 | 6 |
| K₂ / 𝔽₃ | 64 | 6 | 20 | 20 |
| K₂ / 𝔽₅ | 164 | 20 | 20 | 20 |
| K₂ / 𝔽₇ | 312 | 42 | 20 | 20 |
| K₃ / 𝔽₃ | 1300 | 432 | 34 | 34 |

**Point count (Corollary 2) [Y2].** Counting 𝔽_q-points of each type gives

  N_m(q) = Σ_k [m−1, m−k]_q·q^{k(m−k)}·(2q^{m−k} − 1) + Σ_j [m−1, j]_q·q^{j(m−j)}·(2q^j − 1) + 2·Σ_j [m−1, j]_q·q^{j(m−j)}·(q^{m−1−j} + 1) + 2·Σ_j [m−1, j]_q·(q^{m−1−j} − 1)·q^{(j+1)(m−1−j)}.

This equals the census count in every case above.

| m | N_m(q) |
|---|---|
| 2 | 6q² + 2q + 4 |
| 3 | 14q⁴ + 6q³ + 4 |
| 4 | 6q⁸ + 10q⁷ + 20q⁶ + 6q⁵ + 2q⁴ + 4 |
| 5 | 14q¹² + 14q¹¹ + 24q¹⁰ + 14q⁹ + 18q⁸ + 4q⁷ + 4q⁵ + 4 |

**Dimension [Y3].** The fold variety of K_m has **dimension ⌊m²/2⌋**: 2, 4, 8, 12, 18 for m = 2…6. This matches v2.13's exact chart dimensions for K₂ and K₃, and stays below dim G = m² − m for m ≥ 3, as finitely many orbits require. The leading coefficient alternates between 6 and 14 with the parity of m. The constant term 4 counts the folds with trivial V-part.

**Gloss.** *When several indistinguishable steps join the same two instants, order and chaos can split in continuous families. Yet there are only 14m − 8 essentially different ways to do it, and the rest is relabelling of twins.*

---

## 4. Proposed amendments (v2.17 → v2.18)

| amendment | content |
|---|---|
| **Q1¹⁸ Parallel arrows** | Conjecture 3′ is a theorem for every Kronecker quiver K_m: 14m − 8 orbits, an explicit point count, and dimension ⌊m²/2⌋ |

---

## 5. Open questions

1. **Conjecture 3′ beyond two vertices**, e.g. 1⇉2→3, where parallel differences (in ker Δ) now multiply against a further arrow. The Kronecker method (Δ facts (i)/(ii), then G-transitivity) should extend to quivers whose parallel classes are separated by vertices.
2. **Necessity of (H)** (v2.17 Conjecture 5) and the **dimension formula** for looped and cyclic families.
3. **Toward k⟨a,u⟩.** Adding a loop at one Kronecker vertex is the smallest configuration with both ingredients that produce infinite moduli. Find the first quiver whose fold space modulo Aut_[Δ] is infinite-dimensional.
4. Carried: zigzag closed form; dihedral = C_L ∪ C_R; ρ.

---

## 6. Sources

- M. Aguiar, *On the associative analog of Lie bialgebras*, J. Algebra 244 (2001).
- I. Assem, D. Simson and A. Skowroński, *Elements of the Representation Theory of Associative Algebras*, Vol. 1, CUP, 2006 (Kronecker quiver).
- R. P. Stanley, *Enumerative Combinatorics*, Vol. 1, 2nd ed., CUP, 2012, §1.7 (q-binomial coefficients; counting subspaces over 𝔽_q).

*"Between two instants joined by m indistinguishable steps, order and chaos split in 14m − 8 essentially different ways; everything else is a relabelling of twins."*
