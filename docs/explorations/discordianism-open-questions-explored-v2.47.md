# Discordianism: Open Questions Explored (v2.47): where the second family comes from

*v2.46 proved that the folds with the L group's fixed image A form an irreducible 7-dimensional variety V_A. Its
generic behaviours form a one-parameter family. But the census grows by 2q in that group, and all lines share the
image A. This round accounts for the second family: a divisor of V_A on which the symmetry orbits drop by one
dimension. It also proves the mirror theorem for the R group's image B. Three predictions were committed in advance:
two held (27, 28) and one failed (26).*

---

## 0. Summary

| # | Statement | Status |
|---|---|---|
| 1 | **Mirror theorem.** For image B (R group, rep 6118): V_B is irreducible of dimension 7, and dim T V_B = 7 at P₂. The admissible weights are span{(1, 0, 0), (0, 1, 1)}, the L group's with classes a and b swapped. The slice has rank 6, and 7 with N, at t = 1, 2, 3, and N ∉ [Der, P₂]. So v2.46's dominant-slice theorem holds for B | exact over ℚ (prediction 27 held) |
| 2 | The line through the (2, 0) representative 5103 is generic in V_A (rank 6 + 1 at t = 2, 3; 5 + 1 at t = 1). The (2, 0) points are special points on generic lines | exact over ℚ (prediction 26 failed) |
| 3 | **The rank-drop divisor.** D = {P ∈ V_A : rank [𝔥′, P] ≤ 5} has dimension 6, a divisor of V_A, and at its generic points the stabiliser grows by exactly 1. So H′ has 5-dimensional orbits on D, and D carries a one-parameter family of its own | Singular over ℚ (prediction 28 held) |

## 1. The count

In the L group, V_A splits into two strata that both have one-dimensional quotients:
- **the generic stratum:** orbits of the A-stabiliser H′ of dimension 6, quotient of dimension 7 − 6 = 1. This gives
  the class (b, e) = (1, 0), q orbits;
- **the divisor D:** orbits of dimension 5, quotient of dimension 6 − 5 = 1. This gives the class (2, 0), q orbits.

Together they give the 2q growth of the group. By the mirror theorem the same holds for the R group, giving the four
non-inner families of v2.42.

**Caveats.**
- The match between "one-dimensional quotient" and "q census orbits" is a heuristic. Over 𝔽_q, orbits can split or
  merge.
- Whether D is irreducible was not computed: each chart gives dimension 6, but the components were not counted. The
  count (one extra family per group) predicts one component, or several identified by symmetry.

## 2. What is proved, what is computed

| | status |
|---|---|
| line families are folds | Lean (`line_isFold`) |
| V_A, V_B irreducible of dimension 7, geometrically | Singular over ℚ + smooth rational point (v2.46) |
| dominant slice; one-parameter behaviour family on V_A and V_B | proved in characteristic 0 (v2.46; B this round) |
| divisor D of dimension 6 with 5-dimensional orbits | Singular over ℚ (this round) |
| 873 + 22q = 18 inner families + 2 × (generic stratum + divisor) | count computed at q = 3, 5, 7; matching heuristic |

## 3. Next

1. The components of D (`minAssGTZ` on the eliminated ideal), and the behaviour statement on D: does D also carry
   one-parameter behaviours, or position moduli only?
2. Lean closed forms for the 18 inner families of 1⇉2⇉3 + 1→3.
3. Job 14 (GPU).

Files: `analysis/v2_47_divisor.py`, `analysis/v2_46_{slice,exact_checks}.py` (parsers generalised), `logs/v2_47/`.
