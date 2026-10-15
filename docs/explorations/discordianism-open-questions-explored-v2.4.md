# Discordianism: Open Questions Explored (v2.4)

*Answers the open questions of `discordianism-open-questions-explored-v2.3.md` §6.*

*Tags **[Yn]** are machine checks in `verify_v2_4.py` (9/9 pass); [Xn] refer to v2.3. Every claim is labelled **proved**, **verified** or **conjectured**. Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Summary

| # | Question (v2.3 §6) | Answer |
|---|---|---|
| 1 | Density lemma: does k[f] + k[g] = k[u] force two quadratics? | **Proved when the centres differ**, by an exact construction of an annihilating residue functional **[Y1]**. This settles the classification of k[u]-splittings (v2.3 Conjecture 4) except for one sub-case: equal centres with 1/d + 1/e < 1 |
| 2 | Are all deconstructing folds G-conjugates of run-wise extensions (Conjecture 7)? | **No — it fails on the framed branch.** The identification fold is *rigid*: its tangent space is exactly the 1-dimensional conjugation orbit. The framed fold has **2c − 1** deformation directions in degree ≤ c, two new per degree (left and right edges) **[Y7]**. Independent edge twists solve every equation through total degree 5 **[Y8]** |
| 3 | Arrow compatibility for n ≥ 4; quivers with cycles or parallel arrows | **Sufficient as well as necessary for n ≤ 5**, via a *shift-transport* extension: 6, 23, 100, 563 folds **[Y3]**. Splittings of kⁿ number **2(n+2)^{n−1}**, Cayley's count of 2-rooted forests **[Y2]**. The Kronecker quiver has **continuous moduli** of folds **[Y4]** |
| 4 | A non-Borel grid world? | **Per program, yes.** Productive regions can be any G_δ set **[Y5]**, so a single unguarded coin program can name a non-Borel (Σ¹₁-complete) set. Whether the language-wide union is Borel stays **open** |
| 5 | Does social coupling sharpen the clock or lock it? | **Both, depending on relaxation.** If the culture relaxes faster than the clock ticks, coupling sharpens it, until the 3-state Potts transition (K ≈ 4 ln 2) locks it. If it relaxes slower, coupling only adds inertia **[Y6]** |

---

## 1. The density lemma

**Dual formulation (proved).** Define a functional λ on A = k[u] by λ(h) = Σ_m b_m·[u^{m−1}]h. This is a residue pairing with the formal differential φ = Σ b_m u^{−m} du, and every linear functional on k[u] has this form.

- k[f] + k[g] = k[u] ⟺ no λ ≠ 0 kills both k[f] and k[g].
- λ kills k[f] ⟺ λ(f^j) = 0 for all j ≥ 0.

Constraint j involves only b₁,…,b_{dj+1}, with **pivot** b_{dj+1} (coefficient lc(f)^j ≠ 0). The system is *exactly triangular*: no truncation is needed.

*Correction of my own first attempt.* A dimension count of truncated kernels only reproves the trivial count. It cannot see cancellations between high-degree elements, which was the original difficulty. The triangular structure is what avoids it.

**Theorem 1 (proved; verified [Y1]).** Let char k = 0, deg f = d ≥ 2, deg g = e ≥ 2, (d, e) ≠ (2, 2), and suppose the centres differ. (The centre of f is a_f = −(coefficient of u^{d−1})/(d·lc f).) Then **k[f] + k[g] ≠ k[u]**.

*Proof.* Choose b₁, b₂, … in increasing order:
- **Pivot of exactly one family** (m ≡ 1 mod d or mod e, not both): solve that constraint for b_m.
- **Shared pivot** (m ≡ 1 mod lcm(d, e)): both constraints bear on b_m. Eliminating b_m leaves one condition on earlier variables. In it, b_{m−1} has coefficient (m−1)(a_f − a_g) ≠ 0, because [u^{dj−1}]f^j = j·(subleading coefficient). And m−1 ≡ 0 is a pivot of neither family, so b_{m−1} is *reserved* to satisfy the shared constraint.
- **Remaining free indices**: these have density (1 − 1/d − 1/e + 1/lcm) − 1/lcm = 1 − 1/d − 1/e, which is positive. Set one of them to 1.

The result is a nonzero λ that kills both subalgebras.

For (2, 2) every free index is reserved, so the construction forces λ = 0. This is consistent with the dihedral splittings. The exact construction, checked to order 90 for random f and g of degrees (2,3), (3,3), (2,4), (3,4), (3,5) and (4,4), produced a nonzero annihilator every time; for (2,2) none exists **[Y1]**.

**Consequence for v2.3 Conjecture 4 (proved, except one sub-case).** In a nontrivial splitting k[u] = B₁ ⊕ B₂:
1. B₁ and C = k ⊕ B₂ have finite codimension in some k[f] and k[g]. (Lüroth's theorem, plus finiteness of normalization for affine curves.)
2. k[f] + k[g] = k[u] and k[f] ∩ k[g] = k.
3. By Theorem 1, either (d, e) = (2, 2), or f and g share a centre.
4. If (2, 2), the codimensions add up exactly, so B₁ = k[(u−a)²] and C = k[(u−b)²] with a ≠ b. This is the dihedral family.

**Remaining open sub-case.** Equal centres with 1/d + 1/e < 1, for example f = u², g = u³ + u. There the reserved variable must be found deeper in the expansion. The obstruction is the first index at which the series f^{1/d} and g^{1/e} disagree, and they must disagree somewhere, or k[f] ∩ k[g] ≠ k.

Evidence that no splitting arises here:
- spanning fails computationally (v2.3 [X3]);
- exact annihilator solutions to order 72 restrict to a nonzero space at order 24 for three equal-centre pairs, with stable dimensions 8, 10 and 8 **[Y9]**.

**Gloss.** *Two orders can share the world only if each is a mirror. Any richer symmetry (degree ≥ 3) is too sparse to cover what the other misses.*

---

## 2. Conjecture 7: are all deconstructing folds conjugates?

**Graded solve (verified).** All total-degree-graded deconstructing folds up to total degree 3 were solved exactly. There are four families:

| fold | parameters | prediction |
|---|---|---|
| zero | none | ✓ |
| identity | none | ✓ |
| identification | one parameter, a ↦ a + λu | ✓ exactly the G-orbit θ_λ P̂ θ_λ⁻¹ |
| framed | four free parameters at degree 3 (P(a), P(ua), P(auu), P(uua)) | conjugation predicts one |

**Infinitesimal test (verified [Y7]).** Linearize the full equation system (idempotence, RB(−1), Δ-compatibility) at a fold, take the tangent space, and project it onto the coordinates of words of degree ≤ c. Ranks:

| fold | T = 3 | T = 4 | T = 5 |
|---|---|---|---|
| identification | 1, 1, 1 | 1, 1, 1, 1 | 1, 1, 1, 1, 1 |
| framed | 1, 3, 5 | 1, 3, 5, 7 | 1, 3, 5, 7, 9 |

The truncation degree T does not change the low-degree ranks, so these are not top-degree artifacts.
- The identification fold is **rigid**. Its only deformation is the conjugation θ_λ, as Conjecture 7 predicts.
- The framed fold has **2c − 1** directions in degree ≤ c: one at degree 1, then **two new directions per degree**. Conjugation supplies only one in the graded setting.

**Nonlinear test (verified [Y8]).** Impose independent left and right edge twists, P(au) = −λ_L·u² and P(ua) = −λ_R·u². The full graded system then has exact solutions through total degree 5, with residual ~10⁻¹⁶:

| twist | (λ_L, λ_R) | solvable through degree 5 |
|---|---|---|
| conjugation | (1, 1) | ✓ |
| two-sided | (1, 2) | ✓ |
| one-sided | (1, 0) | ✓ |

Conjugation produces only λ_L = λ_R. My explicit guess for the twisted Podge (rewriting boundary letters only) fails, so the true family is subtler than a letter substitution.

**Verdict.** **Conjecture 7 is false as stated** (infinitesimally proved at degree ≤ 5; integrable through degree 5). Deconstructing folds are G-conjugates of run-wise extensions on the identification branch, which is rigid. The framed branch carries **edge moduli**: deformations acting independently on the left and right ends of words, at every run length. The full symmetry of deconstructing folds is larger than the automorphism group G.

*Why the framed fold is soft and the identification fold is rigid.* The identification Podge is a two-sided ideal, which leaves nothing to bend. The framed Podge is a subalgebra but not an ideal: its words are "open" at an a-edge on the left, the right, or both, and each kind of edge can be deformed separately.

**Gloss.** *The framed order is held in place by differences on both sides, and chaos may lean differently against each side. The order that identifies everything is rigid; the order that frames events is soft at its edges.*

---

## 3. Path algebras: sufficiency, Cayley, and parallel arrows

**Theorem 2 (necessity proved for all n; sufficiency verified for n ≤ 5 [Y3]).** A splitting of the vertex algebra kⁿ extends to a deconstructing fold on UT_n **iff** it is arrow-compatible.

The extension is **shift-transport**:

  P(E_ij) = Σ_s (∏_{v=i}^{j} c_{v+s, v}) · E_{i+s, j+s}, where P(e_v) = Σ_w c_{w,v} e_w

Paths are carried along by the same vertex shift as their endpoints. This is needed because compatibility allows P(eᵢ) ∝ e_k with k ≠ i. Two UT₅ splittings (P(e₀) = −e₃, P(e₁) = −e₄) defeat the diagonal extension but not the shift-transport one.

| n | splittings of kⁿ | arrow-compatible = folds |
|---|---|---|
| 1 | 2 | 2 |
| 2 | 8 | 6 |
| 3 | 50 | 23 |
| 4 | 432 | 100 |
| 5 | 4802 | 563 |

**Theorem 3 (verified for n ≤ 5 [Y2]; bijection conjectured).** The number of splittings of kⁿ into two subalgebras is **2(n+2)^{n−1}**. This is Cayley's count of forests on n + 2 labelled vertices with two specified roots, k·N^{N−k−1} with N = n+2 and k = 2.

Conjectured bijection: the points of kⁿ are the n non-root vertices, and the two roots are *Hodge* and *Podge*. Each point hangs, through its blocks, from one root or the other.

**Parallel arrows (verified [Y4]).** For the Kronecker quiver (two arrows α, β: 1 → 2), Δ(α) = Δ(β), so α − β ∈ ker Δ. The kernel is span{e₁, e₂, α−β}, and every fold preserves it.

The folds come in **continuous families**: the exhaustive solve gives components of dimension 0, 1 and 2. One example is the line P(α) = P(β) = (1−t)α + tβ.

**The bridge.** The free grid algebra is the path algebra of one vertex with two loops, and its ker Δ = k[a − b] is generated by the difference of the parallel loops.

| quiver | ker Δ | folds |
|---|---|---|
| UT_n (no parallel arrows) | the vertices | finitely many |
| free algebra (all arrows parallel) | k[a − b] | dihedral and framed families |
| Kronecker | vertices + α − β | continuous moduli |

**Gloss.** *Parallel paths are where deconstruction goes blind (it cannot tell α from β), and that blindness is exactly where order gets freedom.*

---

## 4. The Borel frontier

**Theorem 4 (construction; illustrated [Y5]).** For unguarded coin programs, the productive region can be **any G_δ set**.

*Construction.* Take survival functions S_k(p) = ε_k + Σ_n 2^{−n} max(0, 1 − k·2^n·|p − q_n|), which decrease in k. They tend to D(p) = Σ_n 2^{−n}[p = q_n], which is positive exactly at the rationals.
1. A loop continues with probability S_k/S_{k−1}·(1 − δ_k).
2. Each such probability is continuous and bounded away from 0 and 1 on a compact interior parameter set, so it is simulable by a Keane–O'Brien Bernoulli factory.
3. The loop fails to terminate with positive probability exactly at the rationals.

Hence the productive region is the irrationals (in the illustration, survival stays ≥ 2^{−n} at rationals and is 0 at irrationals) **[Y5]**. With two parameters, R can be any G_δ subset of the square.

**Corollary (proved).** A single program can name a **non-Borel** set:
1. Let the law encode p₁ injectively and continuously, for example i.i.d. Bernoulli(p₁).
2. Let productivity depend on (p₁, p₂) through a G_δ set R whose projection is analytic but not Borel (Suslin).
3. The program's named set is then a homeomorphic copy of that projection.

**Still open.** The *language-wide* set N is the union over all programs, and the language contains the "ungated" version of the program above, which names the whole family. Whether N is Borel for a Turing-complete coin language is therefore not settled by per-program examples.

**Gloss.** *Every Pope can name a world that no Borel description captures, yet the Church's collective naming may wash it out.*

---

## 5. Clocks under social coupling

**Model.** A 3-state culture with logit choice, P(f) ∝ exp(K·m_f + H·[f = fast]), feeds the C₃ "successor of the most-filed" policy of v2.3. The true-tick rate is the fraction of periods in which the fast lane advances by exactly one.

**Results [Y6]**, at weak field H = 0.15:

| sweeps per tick | K = 0 | K = 2 | K = 3.5 |
|---|---|---|---|
| 2 (slow relaxation) | 0.66 | 0.55 | 0.02 |
| 10 (fast relaxation) | 0.64 | 0.92 | 0.02 |

The exploration run gives more detail. With fast relaxation, coupling sharpens the tick rate from 0.78 to 0.95 up to K = 2.5, then the clock locks at K ≈ 3. That matches the mean-field 3-state Potts transition at K_c = 4 ln 2 ≈ 2.77, where the culture's consensus stops following the lane.

**Theorem-shaped summary (verified numerically).** Coupling sharpens a data-driven clock iff the culture's relaxation time is shorter than the tick, and only below the consensus transition. Above it, the clock stops: the culture becomes the order, and the bureaucracy reads only itself.

**Gloss.** *A culture can keep Eris's calendar only if it forgets faster than the seasons turn.*

---

## 6. Proposed amendments (v2.3 → v2.4)

| amendment | content |
|---|---|
| **M3⁗ Splittings of k[u]** | Trivial, character or dihedral, **proved** except for the equal-centre sub-case |
| **M3⁗′ Rigidity** | The identification branch is rigid: its only deformations are G-conjugation. The framed branch has edge moduli, 2 per degree. Conjecture 7 is withdrawn |
| **P1′ Path algebras** | Arrow compatibility ⟺ fold, via shift-transport (n ≤ 5). Splitting count 2(n+2)^{n−1}. Parallel arrows give continuous moduli |
| **B1⁗ Borel** | Per-program named sets can be non-Borel. Language-level Borelness is open |
| **Γ3 Clocks** | Coupling helps a clock only when relaxation is faster than the tick and coupling is below the Potts transition |

---

## 7. New open questions

1. **The equal-centre case** of the density lemma. Find the reserved variable from the first disagreement of f^{1/d} and g^{1/e}.
2. **Edge moduli.** Integrate the 2-per-degree edge deformations of the framed fold to all degrees, and identify the symmetry that generates them. Candidates: one-sided coideal deformations, or a groupoid acting on word ends. Is the dihedral branch rigid like identification, or soft like framed?
3. **The Cayley bijection.** Construct it explicitly, and find the count of *arrow-compatible* splittings (2, 6, 23, 100, 563, …).
4. **Moduli of folds on a quiver.** Is the dimension of the fold variety governed by the parallel-arrow structure, e.g. by dim ker Δ minus the number of vertices?
5. **Language-level Borelness** for Turing-complete coin grids.

---

## 8. Sources

- Keane & O'Brien, *A Bernoulli factory*, ACM TOMACS 4 (1994).
- Kechris, *Classical Descriptive Set Theory* (1995), §14 (Suslin: analytic sets as projections of closed sets).
- Wu, *The Potts model*, Rev. Mod. Phys. 54 (1982): mean-field 3-state transition at K = 4 ln 2.
- Stanley, *Enumerative Combinatorics* Vol. 2, §5.3: rooted forests, k·N^{N−k−1}.
- Lüroth's theorem; normalization of affine curves, e.g. Eisenbud, *Commutative Algebra* §4.6.

*"Two mirrors split the world; a third would leave gaps. Eris's calendar runs only for cultures that forget."*
