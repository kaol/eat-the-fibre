# Discordianism: Open Questions Explored (v2.7)

*Answers the open questions of `discordianism-open-questions-explored-v2.6.md` §5.*

*Tags **[Kn]** are machine checks in `verify_v2_7.py` (folder `v2_7/`, zipped as `verify_v2_7.zip`). It reuses the v2.6 fold engine and adds `ident.py`, an exact solver for identification-type folds that works at every number of a's. Labels: **proved**, **verified** (exact computation to a stated size) or **conjectured**. Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Summary

| # | Question (v2.6 §5) | Answer |
|---|---|---|
| 1 | Closed forms for the moduli | **Done for the identification branch, at every number of a's.** Folds with Hodge k⟨a⟩ are classified by one constant series C(s, t) **[K1]**. Their moduli are **R ∪ L**: the Podge is a **left ideal** or a **right ideal** **[K2–K4]**. R is the family ρ(u) = c + F(D), for an arbitrary series F |
| 2 | Interior moduli | **Obstructed.** P(aua) = κ and every other interior constant C_ij (i, j ≥ 1) is tangent but has C_ij² in the obstruction ideal **[K4]** |
| 3 | Meaning of the ℚ(√s) components | **Root assignment.** On each two-sided dihedral component, left-edge corrections vanish at one root of (u − b)² = s and right-edge corrections at the other **[K5]**. This is necessary but not sufficient **[K6]** |
| 4 | The counting sequence | One more term, a(10) = **25665390**. Still no recurrence **[K7]** |
| 5 | Quiver moduli; Borelness | Carried |

---

## 1. Identification-type folds: the complete picture

By v2.6 Proposition 2, a fold whose Hodge contains 1 can change its Hodge only by conjugation. So for the identification branch, fix the **Hodge k⟨a⟩** and ask which Podges complement it.

**Generators (proved).**
- Since P(aⁱuaʲ) ∈ k⟨a⟩, the Podge Q contains g_ij = aⁱuaʲ − κ_ij(a).
- Every word containing u is a product of the leading words aⁱuaʲ, so ⟨g⟩ + k⟨a⟩ = A by induction on the number of u's.
- Hence **Q = ⟨g_ij⟩ exactly when ⟨g⟩ ∩ k⟨a⟩ = 0**.
- Coideal property: the generators handle Δ, and products preserve it, because Δ(qq′) = (q⊗1)Δq′ + Δq(1⊗q′).

**Theorem 1 (the C-series; proved, checked [K1]).** The coideal condition on the generators is the recursion

  Δκ_ij = Σ_{i′} a^{i′} ⊗ κ_{i−1−i′, j} + Σ_{j′} κ_{i, j′} ⊗ a^{j−1−j′}.

Let K(s, t; a) = Σ κ_ij sⁱtʲ, and read Δ on k[a] as the divided difference (κ(x) − κ(y))/(x − y). Then the recursion becomes

  K(x)(1 − sx)(1 − tx) = K(y)(1 − sy)(1 − ty),

so K(s, t; a)(1 − sa)(1 − ta) does not depend on a. That gives

  **Σ g_ij sⁱtʲ = (1 − sa)⁻¹ (u − C(s, t)) (1 − ta)⁻¹,  C ∈ k[[s, t]] constant.**

In general form, solving for κ_ij with i + j ≤ 3 gives a 10-dimensional solution space, one constant per C_ij **[K1]**. In this language:
- C₀₀ is the centre;
- C₀₁ = α, i.e. P(ua) = α;
- C₁₀ = β, i.e. P(au) = β;
- C₁₁ = κ, i.e. P(aua) = κ.

The residual symmetry a ↦ a + λ acts by Möbius substitution, s ↦ s/(1 + λs) and likewise for t.

**Theorem 2 (one-sided components; proved, checked [K2]).** Let C − C₀₀ = F(t) depend on t alone. Then P(w) = ρ(w)·1, where ρ is the representation of A on k[a] with:
- ρ(a) = multiplication by a;
- **ρ(u) = c + F(D)**, where D is the backward shift aᵐ ↦ aᵐ⁻¹ and D·1 = 0.

Proof:
- Q = {w : ρ(w)·1 = 0} is a **left ideal**, so a subalgebra, and A = k⟨a⟩ ⊕ Q.
- It contains every g_ij, since ρ(aⁱuaʲ)·1 = aⁱ·F(D)aʲ = κ_ij.
- So it equals ⟨g⟩, and Δ-compatibility follows from Theorem 1.

Checked for a random F on all words of length ≤ 6: idempotence, RB(−1) and Δ-compatibility hold **[K2]**. The mirror family L, with C depending on s alone, gives **right-ideal Podges**. The run-wise fold, C = c, is the two-sided ideal and lies in R ∩ L.

**Theorem 3 (nothing else; verified [K3, K4]).** Write the obstruction ideal I_N for the condition ⟨g⟩ ∩ k⟨a⟩ = 0, computed exactly over all products of generators up to length N.
- I₆ is **monomial**. Its reduced Gröbner basis is C₀₁C₁₀, C₀₂C₁₀, C₀₃C₁₀, C₀₁C₂₀, C₀₂C₂₀, C₀₁C₃₀, the products of C₁₁ with C₀₁, C₀₂, C₁₀ and C₂₀, **C₁₁²**, and the products of C₁₂ and C₂₁ with C₀₁ and C₁₀.
- Every interior C_ij with i, j ≥ 1 and i + j ≤ 4 is unobstructed at length 2(i+j+1) − 1, but **C_ij² ∈ I** at length 2(i+j+1).
- Mixed restrictions stay monomial.

Hence the radical forces all interior C_ij = 0 and C₀ⱼC_{i0} = 0, so the moduli are **R ∪ L** to these lengths. In particular, v2.6's interior modulus P(aua) = κ **does not integrate**: it is tangent, and κ² is an obstruction.

*Explicit witnesses (proved).* Every obstruction comes from a word that factors into generators in two ways.
- **Interior, with only C_ij = κ ≠ 0 and i ≤ j.** Take w = aⁱuaʲuaⁱ⁺ʲ. Then g_ij·g_{0,i+j} = w − κ·ua^{i+j}, which gives P(w) = 0. Also g_{i0}·g_{j,i+j} = w − aⁱu·κ_{j,i+j}, where κ_{j,i+j} = κ·a^{(j−i)+i} = κ·aʲ, so P(w) = κ·P(aⁱuaʲ) = κ². Hence κ² ∈ Q ∩ k⟨a⟩. The case i > j is the mirror image. This explains why the obstruction appears at length 2(i+j+1).
- **Mixed edges, with only C₀₁ = α and C₁₀ = β.** g₀₁² = (ua − α)² = uaua − 2α·ua + α². But uaua = u·g₁₁ + u·κ₁₁, and κ₁₁ = (α + β)a. So P(g₀₁²) = (α + β)α − 2α² + α² = **αβ**. The other edge's constant leaks through the interior generator. The general C₀ⱼC_{i0} follows the same pattern; the witnesses up to length 6 are listed by the checks.

In the presence of all variables, the radical statement rests on the observed monomiality, which is verified through length 8.

**Corollary.** Up to conjugation, an identification-type fold is the kernel of a cyclic module structure on k⟨a⟩: either chaos acts on the left, ρ(u) = c + F(D), or on the right. This refines v2.6's "C_L ∪ C_R: two lines" to two infinite-dimensional families, each given by a power series in the shift.

**Gloss.** *Chaos may swallow the events that follow it, with any appetite F, or the events that precede it, but never both. An order that is merely identified is a one-sided ideal of chaos: disorder is absorbed from one side only.*

---

## 2. The dihedral two-sided components: root assignment

v2.6 found two Galois-conjugate two-sided components C_± over ℚ(√s).

**Proposition 4 (verified [K5]).** Take points on the two components at Leq = 3 (m₀ = 1, l₃ = 7) and evaluate every low-level edge correction at the roots r_± = b ± √s:

| component | left edges (h·a·q) vanish at | right edges (q·a·h) vanish at |
|---|---|---|
| first | b − √s | b + √s |
| second | b + √s | b − √s |

So **each two-sided component is an assignment of the two roots of the Podge's quadric to the two sides**, and Galois conjugation swaps the assignment. Diagonal assignments do not occur.

This is the dihedral counterpart of Theorem 2. The run-wise dihedral Podge lies in both identification Podges J_{r_±} = ker(u ↦ r_±), because (u−b)² − s vanishes at both roots. A two-sided deformation appears to tilt the left edge toward one of these and the right edge toward the other.

**Not sufficient (verified [K6]).** The linear space of corrections vanishing at an assigned pair of roots is 6-dimensional modulo conjugation, with Hodge corrections removed, and a generic point of it is obstructed. The components themselves are 2-dimensional. The root assignment is a necessary shadow, not the definition.

**Gloss.** *With two mirrors, chaos can press on both edges only by facing the two points where the Podge vanishes, one point per side. Swapping the two points is the only symmetry left.*

---

## 3. The counting sequence

a(n) for n = 1…10: **2, 6, 23, 100, 563, 3734, 28777, 246164, 2395445, 25665390**. The n = 10 count was exhaustive over 1.03·10¹⁰ forests, in 12 slices that sum to 2·12⁹.

Fit tests with all 10 terms: no linear recurrence fits a(n) or a(n)/n! among order 1 with coefficients of degree ≤ 2, order 2 with degree ≤ 1, or order 3 with constant coefficients. The compatible fraction a(n)/2(n+2)^{n−1} keeps falling by a factor of about 0.45 per step (0.00559 → 0.00249). The closed form remains open. The structural description is v2.6's leaf-type rule.

---

## 4. Proposed amendments (v2.6 → v2.7)

| amendment | content |
|---|---|
| **M3⁷ Identification folds** | Hodge k⟨a⟩ up to conjugation. The Podge is a left ideal, ker(ρ(·)·1) with ρ(u) = c + F(D), or a right ideal. No interior moduli |
| **M3⁷′ Dihedral two-sided folds** | Edge corrections are root-assigned: left at one root of (u−b)² = s, right at the other |
| **v2.6 correction** | "C_L ∪ C_R: two lines" for identification becomes two infinite-dimensional families. The interior modulus P(aua) = κ is obstructed |

---

## 5. New open questions

1. **Proof of Theorem 3 at all lengths.** Find an explicit word that produces C_ij², and the mixed C₀ⱼC_{i0}, for every i, j. The candidate is a word with two u's whose middle a-block can be split two ways.
2. **One-sided structure for the dihedral and framed branches.** Is dihedral C_R the kernel of a cyclic action of A on the Hodge ⟨H₀, a⟩, the analogue of ρ(u) = F(D)?
3. **Full description of C_±.** Root assignment cuts 6 dimensions down to 2. What are the remaining equations?
4. **The counting sequence.** Generating function; a transfer matrix over leaf-type words.
5. Carried: quiver moduli; language-level Borelness.

---

## 6. Sources

- Rota–Baxter operators of weight −1 and splittings into two subalgebras: L. Guo, *An Introduction to Rota–Baxter Algebra*, International Press, 2012.
- Infinitesimal (ε-)bialgebras and the coderivation recursion: M. Aguiar, *On the associative analog of Lie bialgebras*, J. Algebra 244 (2001).
- Deformation theory and quadratic obstructions: M. Gerstenhaber, Ann. of Math. 79 (1964).

*"Chaos absorbs from one side, or it faces two roots."*
