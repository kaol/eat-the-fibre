# Discordianism: Open Questions Explored (v2.8)

*Answers the open questions of `discordianism-open-questions-explored-v2.7.md` §5.*

*Tags **[Ln]** are machine checks in `verify_v2_8.py` (folder `v2_8/`, zipped as `verify_v2_8.zip`). Labels: **proved**, **verified** (exact computation to a stated order or window) or **conjectured**. Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Summary

| # | Question (v2.7 §5) | Answer |
|---|---|---|
| 1 | Prove the identification classification for every word length | **Proved.** Identification-type folds are exactly R ∪ L. Consistency on words with two u's alone forces it, through one functional equation, an edge dichotomy and an infinite descent **[L1–L3]** |
| 2 | One-sided structure for the dihedral fold | **A module hierarchy** (verified to order 2) **[L4]**. Every dihedral component is a bimodule over B_core = ⟨Q₀, aQ₀, Q₀a⟩. The *letter* a is absorbed on both sides by the run-wise fold, on one side by the one-sided families, and on **neither** side by C_±. This is the dihedral analogue of "two-sided ideal → left or right ideal" |
| 3 | What cuts C_± down from 6 to 2 dimensions | Partial. At first order each C_± absorbs a **root-tagged letter** on each side, (u² − r²)·a on the left and a·(u² − r′²) on the right, with {r, r′} = {b ± √s}. At second order these letters disappear **[L5]** |
| 4 | Counting sequence; quiver moduli; Borelness | Carried |

---

## 1. Identification-type folds are exactly R ∪ L

**Operator form.** After conjugation the Hodge is k⟨a⟩. By v2.7 Theorem 1, P(x·u·y) for x, y ∈ k[a] is

  **B(x, y) = Σ C_ij · Dⁱx · Dʲy,**

where D is the backward shift (aᵐ ↦ aᵐ⁻¹, D1 = 0). In particular:
- the R family is B(x, y) = x·F(D)y;
- the L family is B(x, y) = F′(D)x·y.

**Two-u consistency.** For generators, (x u z₁ − B(x, z₁))(z₂ u y − B(z₂, y)) ∈ Q. So

  Φ(x, z₁, z₂, y) = B(B(x, z₁)z₂, y) + B(x, z₁B(z₂, y)) − B(x, z₁)B(z₂, y)

must equal P(x u z₁z₂ u y). That means **Φ depends only on the product z₁z₂**.

**Exponential coordinates.** Let e_r = (1 − ra)⁻¹. Then:
- D e_r = r e_r, so B(e_s, e_t) = C(s, t)·e_s e_t;
- products are partial fractions, e_r e_q = (r e_r − q e_q)/(r − q).

**Theorem 1 (proved; every step machine-checked [L1–L3]).** A fold with Hodge k⟨a⟩ has C(s, t) − C(0, 0) ∈ k[[t]] (family R, Podge a left ideal) or ∈ k[[s]] (family L, Podge a right ideal). By v2.7 Theorem 2, both families consist entirely of folds. So **identification-type folds = R ∪ L**.

*Proof.*
1. **One equation [L1].** Compare Φ(e_s, e_σ, e_τ, e_t) with Φ(e_s, 1, e_σe_τ, e_t). The e_σ and e_τ coefficients vanish identically. The e_s and e_t coefficients share a single numerator N(s, σ, τ, t), bilinear in values of C.
2. **Edge dichotomy [L2].** Setting t = τ = 0 leaves

   (C(s,0) − C(s,σ))·(C(s,0) − C(σ,0)) = 0

   in the integral domain k[[s, σ]]. So either C(s, σ) = C(s, 0), which is family L, or C(s, 0) is constant. The mirror normalisation (z₂ = 1, then s = σ = 0) gives (C(0,t) − C(0,τ))·(C(0,t) − C(τ,t)) = 0: either family R, or C(0, t) is constant.
3. **Interior case [L3].** If C is in neither family, both edge restrictions are constant. Translate so the constant is c = 0 (u ↦ u − c is in G), and write C = st·k(s, t). With τ = 0 the equation becomes

   E₂(k) = σ²(t−s)k(σ,t)k(s,σ) + t²(s−σ)k(σ,t)k(s,t) + s²(σ−t)k(s,σ)k(s,t) = 0.

   Let k_d be the lowest nonzero homogeneous part of k; it must satisfy E₂(k_d) = 0. Two identities hold for every n ≥ 2:
   - E_n(q) at σ = 0 equals st·q(s,t)·(t^{n−1}q(0,t) − s^{n−1}q(s,0)). A function of t alone equals a function of s alone, and both vanish at 0, so both are 0, and st divides q.
   - E_n(st·q′) = sσt·E_{n+1}(q′).

   So the degree drops by 2 at every step. That is impossible for a nonzero k_d, hence k = 0 and C = c lies in R ∩ L, a contradiction. ∎

**Consequences.**
- The moduli are cut out entirely by the **quadratic** (two-u) conditions; longer words add nothing.
- v2.7's monomial obstruction ideal, with C₀ⱼC_{i0} and C_ij² in it, is the visible shadow of steps 2 and 3.
- Up to conjugation, an identification fold is a **cyclic one-sided module**: k⟨a⟩ ≅ A/Q, with u acting by c + F(D) on one side.

**Gloss.** *An order that merely identifies can be disturbed from the left or from the right, with any appetite, but a disturbance from both sides at once is not an order at all. That needs no long words to see: two events already decide it.*

---

## 2. The dihedral module hierarchy

**The run-wise Podge is a ⟨a, Q₀⟩-bimodule (proved).** It consists of the words with some Q₀-run.
- Left multiplication by a or by h·a (h ∈ H₀) adds a new first run and leaves that Q₀-run intact.
- Left multiplication by q ∈ Q₀ either merges into the Q₀-run (Q₀ is a subalgebra) or leaves it untouched.
- Right multiplication is the mirror image.
- It is *not* closed under H₀ alone: h·(q a ⋯) breaks the run.

For the identification fold, ⟨a, Q₀⟩ = A, so this is the two-sided ideal of v2.7.

**Hierarchy (verified to order 2 [L4]).** I computed, for points on each component, all elements x of span{h·a·h′} ⊕ k[u] (low levels) with P(x·q) = 0, respectively P(q·x) = 0, for all q in the Podge. The absorbers may carry ε-corrections x₀ + εx₁ + ε²x₂.

| component | letter a absorbed on the left | on the right | B_core = ⟨Q₀, aQ₀, Q₀a⟩ absorbed on both sides |
|---|---|---|---|
| run-wise | ✓ | ✓ | ✓ |
| left-only family (edges h·a·q) | ✗ | ✓ | ✓ |
| right-only family (edges q·a·h) | ✓ | ✗ | ✓ |
| C₊, C₋ | ✗ | ✗ | ✓ |

The identification analogue: R is a left ideal (a absorbed on the left), L is a right ideal, run-wise is two-sided, and *there is no two-sided component*. The dihedral fold has the extra layer C_±: it **keeps only the core**, with the letter lost on both sides.

**First-order root letters (verified [L5]).** At first order, with x held fixed, each C_± point does absorb a *root-tagged letter* on each side:
- (u² − r²)·a + (Q-terms) on the left;
- a·(u² − r′²) on the right.

Here u² − r² is the Hodge function vanishing at r and its Hodge mirror −r, and {r, r′} = {b + √s, b − √s}, swapped between C₊ and C₋. This matches v2.7's root assignment. At second order, even with corrections, no root letter survives; only B_core remains. So the root assignment is a first-order shadow of the two-sided components, not their defining module structure.

**Gloss.** *Identified order is a one-sided ideal of chaos. Mirrored order keeps a core that chaos can never pry off: the pure chaos runs, and events glued to them. It is free to let go of the bare event a on one side, the other side, or both. Letting go on both sides costs a choice of root.*

---

## 3. Proposed amendments (v2.7 → v2.8)

| amendment | content |
|---|---|
| **M3⁸ Identification folds** | Exactly R ∪ L: cyclic one-sided modules with u acting by c + F(D). **Proved**, by two-u consistency alone |
| **M3⁸′ Dihedral hierarchy** | Every component is a B_core-bimodule. The letter a is absorbed on 2 sides (run-wise), 1 side (C_L, C_R) or 0 sides (C_±) |

---

## 4. New open questions

1. **All orders for the hierarchy.** Does each dihedral component keep its absorber profile to all orders? Conversely, is C_± *characterised* as "B_core-bimodule, no letter absorbed", or does it need more?
2. **Cyclic-module form for the dihedral one-sided families.** The right-only family is a left ⟨a, Q₀⟩-module. Is A/Q ≅ ⟨H₀, a⟩ a cyclic left ⟨a, Q₀⟩-module, with a closed form analogous to ρ(u) = F(D)?
3. **The framed branch in module language.** The framed run-wise Podge is a k⟨a⟩-bimodule, since Q₀ = k. Which side does each of C₀, C_L, C_R keep?
4. Carried: the counting sequence (a(10) = 25665390); quiver moduli; Borelness.

---

## 5. Sources

- L. Guo, *An Introduction to Rota–Baxter Algebra*, International Press, 2012 (splittings into subalgebras).
- M. Aguiar, *On the associative analog of Lie bialgebras*, J. Algebra 244 (2001) (ε-bialgebras).
- M. Gerstenhaber, *On the deformation of rings and algebras*, Ann. of Math. 79 (1964) (deformations and obstructions).

*"Two events suffice to forbid chaos on both sides; two roots suffice to permit it."*
