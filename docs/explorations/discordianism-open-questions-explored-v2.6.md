# Discordianism: Open Questions Explored (v2.6)

*Answers the open questions of `discordianism-open-questions-explored-v2.5.md` §5.*

*Tags **[Wn]** are machine checks in `verify_v2_6.py` (folder `v2_6/`, zipped as `verify_v2_6.zip`, with the fold engine and C enumerators). Labels: **proved**, **verified** (exact computation in a stated window) or **conjectured**. Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Summary

| # | Question (v2.5 §5) | Answer |
|---|---|---|
| 1 | The last density class | **Closed. The density lemma is proved for every degree pair.** When gcd(d, e) = 1, Bézout gives a node or cusp, and that yields an annihilator **[W1]**. When gcd(d, e) ≥ 2, a weighted top form plus Hajós's lemma shows fg ∉ k[f] + k[g] **[W2]**. So v2.3's Conjecture 4 is a **theorem**: splittings of k[u] are trivial, character or dihedral **[W3]** |
| 2 | Dihedral rigidity | **Soft, not rigid.** There are non-conjugation deformations at every output level **[W7]**. One-sided families integrate **[W8]**. Two-sided deformations live on two Galois-conjugate components over ℚ(√s) **[W9]**. Byproducts: a general a-degree-1 tangent formula **[W4]**, and the identification fold is not rigid either once non-graded deformations are allowed **[W5]** |
| 3 | The sequence 2, 6, 23, 100, … | **Combinatorial characterization proved** (an adjacency rule on leaf types) **[W10]**. Extended to n = 9: …, 28777, 246164, 2395445. No small recurrence fits. The closed form is open |
| 4 | Quiver moduli; Borelness | Carried |

---

## 1. The density lemma: complete

**Theorem 1 (proved; checks [W1–W3]).** Let char k = 0, deg f = d ≥ 2 and deg g = e ≥ 2. Then **k[f] + k[g] = k[u] ⟹ d = e = 2**, and then the centres differ.

Work over k̄. Spanning commutes with base change, so this loses nothing. Let V = k[f] + k[g].

**Case gcd(d, e) = 1.**
- Let Φ_f(x, y) = (f(x) − f(y))/(x − y). It has degree d − 1. Its top form is ∏_{ζ ∈ μ_d∖1}(x − ζy), so its points at infinity are [ζ : 1 : 0] with ζ ∈ μ_d∖1. Likewise for Φ_g with μ_e∖1.
- Since μ_d ∩ μ_e = {1}, the curves Φ_f = 0 and Φ_g = 0 share no point at infinity. By Bézout, all (d−1)(e−1) ≥ 1 of their intersections are affine; a common component would give infinitely many.
- An affine intersection point (p, q) gives an annihilator of V:
  - if p ≠ q, it is a **node** of the image curve u ↦ (f, g), and λ = ev_p − ev_q kills V;
  - if p = q, then f′(p) = g′(p) = 0, a **cusp**, and λ = d/du|_p kills V.
- In both cases λ(u) ≠ 0.

The resultant degree (d−1)(e−1) is confirmed for five pairs. v2.5's residual example has an exact witness: k[u³] + k[u⁴ + u] is killed by **λ(h) = h(−1) − h(−ω)** **[W1]**.

**Case gcd(d, e) = g ≥ 2.**
- If the map u ↦ (f, g) is not birational onto its image, generic fibres contain two points p ≠ q, a node again. So assume it is birational, and let P be the minimal polynomial with P(f, g) = 0.
- Give x weight d and y weight e. The weighted top form of Res_u(f − x, g − y) is ±(x^{e′} − y^{d′})^g, where d = gd′ and e = ge′ **[W2]**.
- Suppose fg = A(f) + B(g). Then F = xy − A(x) − B(y) lies in (P), so top(F) is divisible by (x^{e′} − y^{d′})^g.
- top(F) has at most three monomials. Dehomogenizing by w = y^{d′}/x^{e′}, divisibility says w = 1 is a root of multiplicity ≥ g of a polynomial with at most three terms.
- **Hajós's lemma:** a k-term polynomial has nonzero roots of multiplicity ≤ k − 1. (Proof: the moments Σ cᵢ aᵢʲ for j < k vanish, a Vandermonde system; the trinomial determinant is jk(k − j).)
- So g = 2 and top(F) is a genuine trinomial xy + αx^{(d+e)/d} + βy^{(d+e)/e}. That forces d | e and e | d, hence d = e = 2. Otherwise **fg ∉ V**.

Checked: fg ∉ V up to degree 48 for (3,6), (4,6), (3,3) and (2,4); fg ∈ V for the dihedral (2,2) pair **[W2]**.

**Corollary (v2.3 Conjecture 4, now proved).** Every splitting of k[u] is trivial, a character splitting, or dihedral, k[(u−a)²] ⊕ ((u−b)² − s)·k[(u−b)²]. This combines the reduction of v2.4 §1 with Theorem 1 **[W3]**.

*Remarks.*
- In the coprime case the annihilators have finite support: point evaluations and derivatives, not the formal series of v2.4.
- Spanning forces (f, g): A¹ → A² to be a closed embedding: no nodes, no cusps. The Abhyankar–Moh–Suzuki theorem then already forces d | e or e | d, which is an independent proof of the coprime case.

**Gloss.** *Two orders whose periods share no factor always either confuse two distinct events (a node) or share a blind spot (a cusp). Orders whose periods share a factor are betrayed by their own product, which neither can express. Only two mirrors escape both traps.*

---

## 2. Deformations of deconstructing folds

### 2.1 Structure (proved)

Fix the splitting p of A₀ = k[u] = H₀ ⊕ Q₀. Since ker Δ = k[u], compatibility with letter deletion forces, by induction on the number of a's:

  P(x₀ a x₁ ⋯ a x_n) = Σ over subsets of the a's that are kept, of c_{k₁}(segment) a c_{k₂}(segment) a ⋯

Here c₀ = p, and c_k: A₀^{⊗(k+1)} → k[u] is the correction attached to a segment from which k a's were deleted. So a fold over p is a sequence (c₁, c₂, …). The run-wise extension is the case where all c_k vanish.

**Proposition 2 (a-degree-1 tangent; proved).** Let K = ker(p⊗p). The conditions on c₁ are *linear*:
- idempotence: c₁(H₀⊗H₀) ⊂ Q₀ and c₁(K) ⊂ H₀;
- RB(−1): c₁ restricted to H₀⊗H₀ is an H₀-bimodule map into A₀/H₀, and c₁ restricted to K is a Q₀-bimodule map into A₀/Q₀.

Hence

  **T¹ = Hom_{H₀}(H₀⊗H₀, A₀/H₀) ⊕ Hom_{Q₀}(K, A₀/Q₀),**

and conjugation θ_h (a ↦ a + h) contributes c_h(x, y) = p(x)·h·p(y) − p(xhy).

If 1 ∈ H₀, the first summand is generated by c₁(1⊗1) ∈ Q₀, and conjugation realizes all of it. So **a fold with unital Hodge can change its Hodge only by conjugation**, and every new modulus sits on the Podge side.

| fold | Hom_{Q₀}(K, A₀/Q₀) beyond conjugation | consequence |
|---|---|---|
| identification (H₀ = k, Q₀ = (u)) | (K/m²)*, 2-dimensional: **P(ua) = α, P(au) = β** (constants) | *non-graded* moduli, invisible to v2.4's graded rigidity test |
| framed (H₀ = (u), Q₀ = k) | all linear maps K → H₀: the **edges** P(uᵐa), P(auᵐ), P(a) | v2.5's edge moduli, plus a Hodge-side **P(uau) = uau + κ** |
| dihedral | infinite, growing by 2 per unit of output shift δ | next section |

The window counts match these predictions: 2 for identification and 10 for framed **[W4]**.

### 2.2 Method, and three artifacts it removed

Integrability is tested order by order in ε, with exact arithmetic modulo 2⁶¹ − 1, on a window of input levels. The basis of k[u] is adapted to the splitting: for the dihedral case, (u−a)^{2k} and ((u−b)² − s)(u−b)^{2k}, which makes p a coordinate projection and truncation exact.

Three artifacts, each of which produced wrong verdicts in this round until removed:
1. **Correction window.** Second-order corrections of conjugations need an output shift of 2δ.
2. **Skipped equations.** A large output shift D pushes products out of the window. Silently dropping those equations manufactured freedom; for a while even the framed fold looked unobstructed.
3. **Basis indices.** Indices depend on the truncation, so vectors must not be mixed across windows.

The final runs are **skip-free**: no equation is ever dropped. So every obstruction reported is a genuine necessary condition, and an "ok" is evidence up to the stated order.

The integrator was calibrated on v2.5's framed classification **[W6]**:
- points on C₀ (uniform), on C_L (arbitrary left) and α₁ alone integrate;
- α₁ + β₂, γ + α₁, and α₁ = β₁ = 1 (others 0) are obstructed.

The last point satisfies v2.5's three degree-4 quadrics. So those quadrics do **not** by themselves cut out C₀ ∪ C_L ∪ C_R, and v2.5's "vanish exactly" overstated it. The classification stands, enforced by higher-level equations.

### 2.3 Answers

**Identification fold (verified [W5]).** P(ua) = α alone and P(au) = β alone integrate to order 4. Together they are obstructed at order 2, with quadric αβ = 0. So v2.4's "rigid" holds **only in the graded category**. Non-graded, the identification fold has a small C_L ∪ C_R: two lines.

**Dihedral fold (verified [W7–W9]).** It is **not rigid**.
- **Tangent.** In the window, there are 10 non-conjugation directions at δ = 1 and 12 at δ = 2 **[W7]**.
- **Hodge mode.** The Hodge-only direction is obstructed by itself (h² is an obstruction), as Proposition 2 predicts.
- **One-sided families.** Left-only and right-only families, with the Hodge mode removed, integrate at a generic point. A generic two-sided point is obstructed **[W8]**.
- **Two-sided components.** Two-sided deformations survive on **two Galois-conjugate components defined over ℚ(√s)**. Points on them integrate, and nudging any coordinate obstructs **[W9]**.
  - The defining quadratic has discriminant in s·ℚ² for (a, b, s) = (0,1,2), (0,1,3), (0,2,2) and (1,0,5); the discriminants are 512, 1728, 2048 and 8000.
  - So the two components correspond to the two points **u = b ± √s** where the Podge's defining quadric (u − b)² − s vanishes.

| fold | one-sided moduli | two-sided moduli |
|---|---|---|
| identification | C_L, C_R: one constant each | none (αβ = 0) |
| framed | C_L, C_R: any edge sequence | C₀: uniform edges plus γ, through conjugation |
| dihedral | C_L, C_R: polynomial and constant edges | **C₊, C₋**, one per root b ± √s |

*Also found:* tangent directions with c₁ = 0 but c₂ ≠ 0, i.e. **interior moduli at a-degree 2**, for the unital-Hodge folds. Examples: P(aua) = κ for identification; 27–200 such directions in dihedral windows. The framed fold has none. Their integrability is untested.

**Gloss.** *With two mirrors, chaos may lean on either edge of an event as it likes. To lean on both edges at once, it must commit to one of the two points where the Podge's mirror-pair vanishes, and the two commitments are Galois twins. Only conjugation moves the Hodge itself: order can be renamed, not bent.*

---

## 3. Compatible splittings: the leaf-type rule

**Theorem 3 (proved; verified by brute force on all forests n ≤ 8 [W10]).** Recall the forest picture of v2.5: rooted forests on the points, with each tree's root colour free and colours alternating. Classify each point as one of:
- **E**: a Podge leaf (S_v = ∅);
- **HL**: a Hodge leaf (S_v = {v});
- **P1**: a Podge vertex whose only child is a leaf c (S_v = {c});
- **B**: anything else.

Then the splitting is arrow-compatible **iff** every consecutive pair (i, i+1) contains an E, or is (HL, HL), or is (P1, P1) with consecutive children c, c + 1.

*Proof.*
- S_v = ∅ iff v's subtree has no Hodge vertex. Colours alternate, so this means v is a Podge leaf.
- S_v is a singleton iff exactly one Hodge vertex lies below v and it is a leaf. That is the HL case or the P1 case.
- The mixed pairs (HL, P1) and (P1, HL) would need a child equal to i or i+1, which is impossible.

Counts for n = 1…9: **2, 6, 23, 100, 563, 3734, 28777, 246164, 2395445**. The n = 9 count took a 147 s exhaustive run over 4.3·10⁸ forests.
- No linear recurrence of order ≤ 3 with polynomial coefficients of degree ≤ 2 (order 1–2 for degree 2) fits a(n) or a(n)/n!.
- The compatible fraction a(n)/2(n+2)^{n−1} decays by a factor of roughly 0.45 per step.

**Gloss.** *An arrow survives the split only between two leaves of order, between two single-child chaos nodes whose children stand in step, or where one end is a bare leaf of chaos.*

---

## 4. Proposed amendments (v2.5 → v2.6)

| amendment | content |
|---|---|
| **M3⁶ Splittings of k[u]** | Trivial, character or dihedral: **proved** (density lemma complete) |
| **M3⁶′ Fold structure** | A fold over p is a sequence of corrections c_k. T¹ = Hom_{H₀}(H₀⊗H₀, A₀/H₀) ⊕ Hom_{Q₀}(K, A₀/Q₀). With unital Hodge, the Hodge moves only by conjugation |
| **M3⁶″ Moduli** | Identification: C_L ∪ C_R (non-graded). Framed: C₀ ∪ C_L ∪ C_R. Dihedral: C_L ∪ C_R ∪ C₊ ∪ C₋, with C_± over ℚ(√s). The "rigid" of v2.4 means graded-rigid |
| **v2.5 correction** | The degree-4 framed quadrics do not cut out C₀ ∪ C_L ∪ C_R alone (α₁ = β₁ = 1 passes them) |
| **P2⁶ Path algebras** | Compatibility is the leaf-type adjacency rule |

---

## 5. New open questions

1. **All-order closed forms for the dihedral components.** Find the analogue of v2.5's generated Podge′ for C_L, C_R and C_±. In particular, does C_± come from evaluating at b ± √s?
2. **Interior moduli.** Do the a-degree-2 directions (P(aua) = κ, …) integrate? Does every unital-Hodge fold have new moduli at every a-degree?
3. **Unifying the two-sided components.** Is framed C₀ the one-root degeneration of dihedral C_±? What plays the role of b ± √s for identification, which has none?
4. **The sequence 2, 6, 23, 100, 563, 3734, 28777, 246164, 2395445.** Find a transfer-matrix or generating function over leaf-type words.
5. Carried: moduli on quivers; language-level Borelness.

---

## 6. Sources

- Bézout's theorem and intersection at infinity: W. Fulton, *Algebraic Curves*, Ch. 5.
- Hajós's lemma (a k-term polynomial has nonzero roots of multiplicity ≤ k − 1): G. Hajós, 1953; see also A. Schinzel, *Polynomials with Special Regard to Reducibility*, Cambridge UP, 2000.
- Embeddings of the line in the plane: S. S. Abhyankar and T.-T. Moh, J. reine angew. Math. 276 (1975); M. Suzuki, J. Math. Soc. Japan 26 (1974).
- Deformation theory and quadratic obstructions: M. Gerstenhaber, *On the deformation of rings and algebras*, Ann. of Math. 79 (1964).

*"Coprime orders confuse a point; commensurate ones are betrayed by their product. Two mirrors split the world, and chaos, leaning on both edges, must pick a root."*
