# Discordianism: Open Questions Explored (v2.61): the λ = 1 split lifts form an affine-linear bundle

*Program 3, item 1 of v2.60 §4: the split lifts with λ = 1 and μ ≠ 0, the only part of the split kind that could carry
moduli, for Q₀ → Q₁ = Q₀ + ℓ: 2→4.
- **The result (char 0, exact).** At every fold P̄ of Q₀, the λ = 1 lifts are the solutions of a linear system in
  μ ∈ B* (dim B = 7) whose coefficients are polynomials in P̄. So the λ = 1 split locus is an affine-linear bundle over
  Folds(Q₀), piece by piece.
- **Its pieces.** It has 72 irreducible components:
  - 60 are the identity lifts, with fibres that are single points;
  - the other 12 have fibres of dimension 0, 1 or 2;
  - generically on each component of Q₀ the fibre is a point or empty.

This is the fibre structure a leaf-extension proof needs. The kQ-only, J-only and mixed kinds remain.*

---

## 0. Results

The predictions are in logs/v2_61/PREDICTIONS.md. 69 and 70 were committed while their runs were in progress. 71 was
added after a test on two strata and before the full run.

| # | Statement | Status |
|---|---|---|
| 69 | Over the generic point of every Q₀ component the μ-fibre is affine, of dim ≤ 2 | **held**, more strongly. All 118 generic fibres are affine of dim 0: 60 are the identity lift, 1 is a unique non-identity lift (100,010,001 comp 4), and 57 are empty (`mu_generic.txt`) |
| 70 | 𝔽₅ enumeration: the same picture as 𝔽₃ | **held**: 21 classes. 17 lines through 0 and 1 plane through 0, all spanned by μ-vectors (14 of the lines and the plane by coordinate functionals p*), plus 3 isolated points without the identity lift, at both q = 3 and q = 5 (`mu_{3,5}.txt`) |
| 71 | On all 21 strata of Q₀ the λ = 1 lift ideal in ℚ[vs, μ] is generated in μ-degree ≤ 1 | **held**: 21 of 21 (`mu_linear.txt`) |

## 1. The objects

From v2.60 §2: over a fold P̄ of Q₀, a split lift with λ = 1 is P(μ) = P̄ ⊕ P_J. Here
- P_J(e₄) = e₄;
- P_J(pℓ) = φ(p)ℓ + μ(p)e₄;
- φ(p) = P̄(p) + Σ μ(suf)P̄(pre), and φ must map B into B.

P(μ) is affine in μ. M(P̄) is the set of μ for which P(μ) is a fold of kQ₁.

## 2. Four computations, from evidence to theorem

1. **Over 𝔽₃ and 𝔽₅, every μ** (`analysis/v2_61_mu.py`). All 3⁷ (respectively 5⁷) μ are tested over each of the 407
   (513) G-orbit representatives of Q₀, with a staged fold test.
   - M(P̄) is always an affine subspace.
   - The μ ≠ 0 lifts sit over 21 census classes: 17 lines through 0 (the identity lift), 1 plane through 0, and 3
     single points where there is no identity lift.
   - This accounts exactly for v2.60's 23 classes with μ ≠ 0: each line is one class, the plane gives 3 (two axes and
     the open part), and each point gives 1.

   The first version of the script dropped the single-point fibres: it discarded |M| ≤ 1 instead of M ⊆ {0}. The
   3 missing classes exposed the bug.
2. **Exactly over ℚ at the integral points** (`analysis/v2_61_mu_exact.py`). 388 (from the q = 3 dump) and 486 (q = 5)
   representatives are folds of ℚQ₀. At each, the full fold system of P(μ) in ℚ[μ] has a reduced Gröbner basis of
   degree ≤ 1, so M(P̄) is an affine space even scheme-theoretically. The dimensions are as over 𝔽_q.
3. **At the generic point of each component** (`analysis/v2_61_mu_generic.py`). Each of Q₀'s 118 components is a
   graph over an independent set t. The fold system of P(μ) is generated over ℚ(t) by Singular, and its reduced basis
   has degree ≤ 1 for all 118. The fibre dimensions are as in the table above.

   The Singular fold-equation generator was first validated against an exact integer fold test on 120 matrices, 63 of
   them folds; the two tests agreed on all 120. An early version used `unitmat` without loading its library, silently
   dropped the kernel-closure conditions, and failed this validation.
4. **At every point of every stratum** (`analysis/v2_61_mu_linear.py`). In ℚ[vs, m₀..m₆], let L be the ideal generated
   by:
   - the stratum's equations;
   - the conditions φ(B) ⊆ B;
   - the fold equations of μ-degree ≤ 1;
   - every remainder of μ-degree ≤ 1 left when a higher-degree equation is reduced modulo L, added until stable.

   Each such remainder lies in the fold ideal I. On all 21 strata every equation of μ-degree ≥ 2 then reduces to 0
   modulo L, so I = L.

   **Theorem (char 0).** For every fold P̄ of Q₀ over a field of characteristic 0, the λ = 1 split lifts of P̄ to Q₁
   are the solutions of a linear system A(P̄)μ = b(P̄), with A and b polynomial in P̄. In particular M(P̄) is an affine
   subspace of B*. Specialising generators of I = L to P̄ gives generators of the fibre's ideal. ∎

   The linear equations are not all visible as fold equations of μ-degree 1. For example, on 100,000,100, m₁ and m₂
   arise only as remainders of quadratic equations. So the saturation step is needed.

## 3. The bundle (`analysis/v2_61_mu_bundle.py`, `logs/v2_61/mu_bundle.txt`)

Decomposing L over ℚ (minAssGTZ, with the overflow-free primdec of v2.55) gives, for each component W of the λ = 1
split locus, its dimension and the dimension of its image in the stratum. The difference is the fibre dimension.

| | components |
|---|---|
| strata with no λ = 1 lift (1 ∈ L: for example 010,010,000) | 4 strata |
| identity-lift components (μ ≡ 0 on W; fibre a point) | 60, over exactly the 60 Q₀ components with P̄(B) ⊆ B generically |
| non-identity components | 12: fibre dimension 0 (2), 1 (9), 2 (1, in 100,010,001) |
| total | 72 |

So the moduli question for item 1 comes down to 12 components:
- **Fibre dimension 0:** a section over its image, with no extra moduli.
- **Fibre dimension 1:** lines. On 𝔽_q these are lines through the identity lift, swept by the ℓ-weighting. One
  nonzero orbit in each, as the census shows.
- **Fibre dimension 2:** a plane, on which the census sees exactly the three nonzero orbits of a 2-torus acting by
  independent characters.

The finiteness of orbits on these 12 is already part of v2.57's theorem for all of Q₁. What is new is the reason: the
fibres are affine of dimension ≤ 2 with torus actions, not an accident of a large computation.

## 4. Where Program 3 stands

| item (v2.60 §4) | status |
|---|---|
| split, λ = 0 | the trivial lift; proved (v2.60) |
| split, λ = 1, μ = 0 | the identity lift iff P̄(B) ⊆ B; proved (v2.60) |
| split, λ = 1, μ ≠ 0 | affine-linear bundle with fibres of dim ≤ 2, everywhere in char 0, for Q₀ → Q₁ (this round). A general proof of linearity for any sink leaf is not yet found |
| kQ-only, J-only | next: the same treatment, with the off-diagonal block as the unknown |
| neither | 16 classes in 6 vertex blocks; to be handled through the vertex block |

**What would generalise.** The theorem of §2.4 is a computation for Q₀. The census pattern behind it (every M(P̄) an
affine space of dim ≤ 2) can be tested at once on the eight smaller leaf extensions of v2.60, using census data
(`v2_61_mu.py` generalised). Proving linearity for every sink leaf would need to explain the remainder mechanism of
§2.4 by hand.

## Files

- analysis/v2_61_mu.py, v2_61_mu_exact.py, v2_61_mu_generic.py, v2_61_mu_linear.py, v2_61_mu_bundle.py
- logs/v2_61/:
  - PREDICTIONS.md;
  - mu_3.txt, mu_5.txt;
  - mu_exact_3.txt, mu_exact_5.txt;
  - mu_generic.txt, mu_linear.txt, mu_bundle.txt.
