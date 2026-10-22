# v2.44 predictions (committed before any computation)

Quiver 1⇉2⇉3 + 1→3, merged orbits (G, local torus, shear U). v2.34's invariants at a fold P:
- b = dim T_P Folds − dim([Der, P] ∩ T): tangent directions not produced by any derivation;
- e = dim([Der, P] ∩ T) − dim(census symmetry directions): derivation directions outside the census group.

Both are invariant under the census group and the shear (Δ-compatible automorphisms), so one representative per
merged orbit suffices.

8. (sanity) At points of the 18 inner-curve families (curve points t = 1, 2 over 𝔽₅), e ≥ 1: the family tangent
   [ad x, P] is a derivation direction outside the census group.
9. In the two dim-7 image groups, the per-(b, e) orbit counts at q = 3, 5, 7 grow only in classes with e = 0.
   Together they grow by 4 per unit of q, 2 per group. That is, the 4 missing families are not orbits of automorphism
   groups to first order: along them the splitting itself, its behaviour, varies.

## Outcomes of 8 and 9: both **held** (`be_orbits_L.txt`)

- 8: points of the inner families have (b, e) ∈ {(0, 1), (0, 2)}.
- 9: in the L target group only the classes (1, 0) (6, 8, 10) and (2, 0) (8, 10, 12) grow, by q each.

Exploration (not predicted): at rep 5107 (class (1, 0), q = 5), T_P has dimension 8 and [Der, P] ∩ T has dimension 7.
The kernel-moving idempotent line P₀ + t·N, with rank-one integer N on rows (e₃, b₁, b₂) and columns
(a₁, a₂, c, a₁b₁ or a₁b₂, a₂b₂), is Δ-compatible and idempotent over ℚ[t]. Its points meet exactly q orbits at
q = 3, 5, 7 (`lines_5107.txt`).

## Predictions 10–12 (committed before the exact checks and the other searches)

10. For that line (the first N in lines_5107.txt):
    (a) ker P(t) is a subalgebra for all t, exactly over ℚ[t], so P(t) is a fold over ℚ for every t;
    (b) N ∉ [Der(kQ), P₀] over ℚ, by an exact rank computation. Hence in characteristic 0 the curve is not contained in
        one Aut(kQ)-orbit: Aut(kQ) is an algebraic group with Lie algebra Der(kQ), so a curve inside the orbit of P₀
        would have tangent in [Der, P₀]. Its members realise infinitely many behaviours.
11. The mirror R group has the same structure: an image- or kernel-moving line from a (1, 0) representative that meets
    q orbits at q = 3, 5, 7.
12. The (2, 0)-class family of each group is also an idempotent line from one of its representatives, found by the
    same search.

## Outcomes of 10–12

- 10 **held** (`exact_5107_N0.txt`):
  - (a) P(t) = P₀ + t·N is a fold over ℚ[t]: idempotent, image and kernel closed, Δ-compatible.
  - (b) dim Der(kQ) = 22 over ℚ, rank [Der, P₀] = 22, and with N the rank is 23. So N ∉ [Der, P₀].
- 11 **held as stated**: the R group has kernel-moving lines meeting q orbits at q = 3, 5, 7 (`driver_R_10.txt`).
- 12 **held as stated**: the L group's (2, 0) representatives give such lines (`driver_L_20.txt`).
- **Caveat on 11–12:** lines are abundant. 32 distinct lines were found from 16 representatives, overlapping in
  orbits, so lines are not in one-to-one correspondence with the 4 families. Coverage (`coverage.txt`): the lines
  meet q + 1 of the L group's (1, 0) orbits and q + 3 of its (2, 0) orbits at q = 3, 5, 7. Coverage of the R group is
  partial (6, 7, 9 and 8, 8, 9), because only 6 representatives were searched. Every line found moves the kernel and
  fixes the image.

## Prediction 13 (committed before the Lean build)

13. The line family from rep 5107 (first N in `lines_5107.txt`) passes the reflection fold check
    (`check c D lineP = true`, one `decide +kernel` per row). So it is a fold over every commutative ring for every t.
