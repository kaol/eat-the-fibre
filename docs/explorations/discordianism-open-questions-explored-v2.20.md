# Discordianism: Open Questions Explored (v2.20)

*Takes up question 1 of `discordianism-open-questions-explored-v2.18.md` §5: Conjecture 3′ beyond two vertices. The first case is 1⇉2→3, a parallel pair followed by a further arrow.*

*Tags **[Vn]** are machine checks in `verify_v2_20.py` (folder `v2_20/`, zipped as `verify_v2_20.zip`; `pip install -r requirements.txt`, and `SLOW=1` adds the local analysis). Labels: **proved**, **verified (exact)**, **numerical** or **conjectured**. Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Summary

| # | Question | Answer |
|---|---|---|
| v2.18 Q1 | Finitely many orbits for 1⇉2→3? | **Supported, modulo Aut_Δ ⋊ grading (Conjecture 3″).** Three pieces of evidence: <br>• **Census over 𝔽₃ (exact):** 219 folds in 71 orbits of G(𝔽₃) **[V1]**. <br>• **Local analysis (exact tangents, numerical local dimension):** 161 of the folds lift to ℚ. At 31 of them the local dimension exceeds the orbit dimension, but at a generic nearby fold the two agree in every case. These are fixed-point degeneracies, not moduli **[V3]** |
| — | A structural correction | **With parallel arrows, folds need not preserve path length.** Some folds send β ↦ β + δβ, with δ = α₁ − α₂ **[V2]**. The filtration folds do preserve is the **Δ-filtration** F_k = ker Δ^{(k+1)}, in which δβ has degree 1 although it has length 2. This is the same phenomenon as u in k⟨a,u⟩. A length-staged enumeration of mine silently lost two folds because of it |
| — | Limits | 𝔽₅ was out of reach (about 143 s per liftable kernel fold). The local analysis covers only the 161 folds that lift from 𝔽₃. Generic orbits are measured under Aut_Δ together with grading; whether Aut_Δ alone suffices was not separated |

---

## 1. Method: kernel fibration

For A = k(1⇉2→3), with basis e₁, e₂, e₃, α₁, α₂, β, α₁β, α₂β:
- **The kernel.** K = ker Δ = span(e₁, e₂, e₃, δ) is a subalgebra (Δ is an ε-derivation), and every fold preserves it. So P|_K is an *algebra* fold of K. There are 384 such folds over 𝔽₃.
- **The complement.** On span(α₁, β, α₁β, α₂β), Δ is injective. So each P(x) equals the unique Δ-preimage of (P ⊗ P)Δx, plus an unknown K-component.
- **Enumeration.** These K-components are enumerated in two stages, and every condition is checked only once all the columns it involves are known. The final stage checks idempotency and both subalgebra conditions in full; Δ-compatibility holds by construction.

The enumeration is complete because each step is forced. A first version staged by path length and assumed P(A_{≤1}) ⊂ A_{≤1}. That assumption is false, as the example above shows, and that version found 217 folds. The corrected count is **219** (cross-checked: the 2 extra folds satisfy every fold equation exactly).

**Symmetries.** Aut_Δ has the 2-dimensional Lie algebra of v2.13 T7: α_i ↦ α_i + c_i·δ, extended to α_iβ. Over 𝔽₃ its group is the 6 column-sum-1 substitutions. Length grading (x ↦ λ^{|x|}x) preserves Δ up to a scalar, so it also maps folds to folds.

---

## 2. Local structure

At each of the 161 folds that lift to ℚ (2 ↦ −1):
- the **Zariski tangent dimension** (exact Jacobian of the fold equations);
- the **orbit dimension** under Aut_Δ, and under Aut_Δ plus grading (exact);
- the **local dimension** (v2.12's ε-scaling test, numerical).

| (tangent, Aut orbit, Aut + grading orbit, local) | folds |
|---|---|
| (0, 0, 0, 0) | 12 |
| (1, 1, 1, 1) | 32 |
| (2, 2, 2, 2) | 64 |
| (3, 2, 2, 2) | 20 |
| (2, 1, 1, 1), (5, 1, 1, 1) | 4 |
| (1, 0, 0, 1) | 11 |
| (2, 1, 1, 2) | 8 |
| (3, 1, 1, 2), (3, 1, 1, 3) | 8 |
| (5, 1, 1, 4) | 2 |

Rows where the tangent exceeds the local dimension (for example (3, 2, 2, 2)) are obstructed first-order directions, as at the Kronecker folds of v2.12. The last four rows (29 folds) have local dimension above the orbit dimension. At each such fold I stepped to a generic nearby fold on the local family and measured again (31 flagged in that run, because the sampling is random) **[V3]**:

| at the fold (local, orbit) | at a generic nearby fold (local, orbit) | count |
|---|---|---|
| (1, 0) | (1, 1) | 11 |
| (2, 1) | (2, 2) | 12 |
| (3, 1) | (2, 2) | 6 |
| (4, 1) | (2, 2) | 2 |

In every case the generic local dimension equals the generic orbit dimension. The excess at the fold itself comes from fixed points of the symmetry, such as cone vertices of grading orbits and intersections of orbit closures, not from moduli.

**Gloss.** *When twin steps are followed by a further step, a split may borrow from the twins' difference and attach it to what follows (β ↦ β + δβ). Still, every family of splits found is a relabelling of twins combined with a change of the clock's scale.*

---

## 3. Proposed amendments (v2.19 → v2.20)

| amendment | content |
|---|---|
| **Q1²⁰ Filtration** | With parallel arrows, folds preserve the Δ-filtration, not path length |
| **Q2²⁰ 1⇉2→3** | Census over 𝔽₃ (219 folds, 71 orbits). Locally finitely many Aut_Δ ⋊ grading orbits at all 161 lifted folds (numerical). Conjecture 3″ supported |

---

## 4. Open questions

1. **Prove the 1⇉2→3 case**, ideally by a Kronecker-style hand classification using the Δ-filtration.
2. **Separate Aut_Δ from grading**: does Conjecture 3′ (Aut_Δ alone) hold here, or is grading needed?
3. A **faster census** (𝔽₅, 𝔽₇) to get the point-count polynomial, as for Kronecker.
4. Carried: the combinatorial criterion for unobstructed violations; Kronecker plus a loop; zigzag; dihedral; ρ.

---

## 5. Sources

- M. Aguiar, *On the associative analog of Lie bialgebras*, J. Algebra 244 (2001).
- I. Assem, D. Simson and A. Skowroński, *Elements of the Representation Theory of Associative Algebras*, Vol. 1, CUP, 2006.
- W. Decker, G.-M. Greuel, G. Pfister and H. Schönemann, *Singular* (via passagemath-singular).

*"When twin steps are followed by a further step, a split may borrow the twins' difference and carry it forward; even so, every family of splits found is a relabelling of twins combined with a change of the clock's scale."*
