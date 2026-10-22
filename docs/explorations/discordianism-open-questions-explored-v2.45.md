# Discordianism: Open Questions Explored (v2.45): the variety of folds with a fixed image

*v2.44 found lines of folds on 1⇉2⇉3 + 1→3 along which the behaviour changes, all with a fixed 7-dimensional image
and a tilting kernel. This round studies the variety V_A of all folds with a given image A. For the image of those
lines it is irreducible of dimension 7 over ℚ and smooth at the base point. Automorphisms preserving A sweep out 6 of
those 7 directions at generic points, which leaves (at least) a one-dimensional family of behaviours. Seven
predictions were committed before the computations. Three held as stated, two held in part, and two failed. The
failures corrected a wrong reading of the tangent invariants.*

---

## 0. Summary

| # | Statement | Status |
|---|---|---|
| 1 | V_A = {P₀ + N : N ∈ Hom(ker P₀, A), P₀ + N a fold} for A = im P₀, rep 5107 (L group). 35 unknowns, 995 quadratic equations. Singular: one minimal associated prime, dimension 7, containing P₀. The tangent dimension at P₀ is 7, so V_A is smooth there | computed exactly over ℚ (prediction 16 held) |
| 2 | Rep 5088 (a constant orbit with (b, e) = (3, 0)) has the same image A. It is a special point of the same V_A | computed exactly (prediction 19 held) |
| 3 | At generic points of V_A (the v2.44 line at random t), the directions [Der_A, P] ∩ T_P V_A have dimension 6 of 7, so m = 1. At the special point 5088, m = 3 | computed mod 9999991 (prediction 20 held; predictions 17, 18 failed) |
| 4 | Behaviours with image A form a family of dimension ≥ 1 (characteristic 0). The first-order count says exactly 1 at generic points, assuming the automorphism orbits reach their tangent dimension, which is not checked | ≥ 1 from v2.44's exact rank argument; = 1 is a first-order statement |
| 5 | dim T V_A = 7 at every representative of both dim-7 groups. m = b at most L-group representatives but not in the R group | computed (prediction 14 failed; 15 held in L) |

## 1. The picture

Fix the 7-dimensional subalgebra A, the common image. The folds with image A are the complements K = (1 − N)·ker P₀
that are subalgebras and make P Δ-compatible. Over ℚ they form one irreducible 7-dimensional variety. Automorphisms
of kQ preserving A act on it with orbits of dimension at most 6 (first order). So, generically, V_A is a
one-parameter family of behaviours, each spread over a 6-dimensional orbit.

The v2.44 lines are transversals to those orbits. The census sees this as q orbits of a growing class.

The earlier reading of m was wrong. I had treated m = 3 at a constant orbit as a tangent overcount (prediction 17),
but V_A is smooth of dimension 7 at that point too. The point is special: its automorphism orbit is smaller. m is the
codimension of the orbit tangent at one point, not a moduli dimension. Only its generic value, 1, is a moduli count.
- The L group's (1, 0)/(2, 0) split and the R group's spread of m values (1–3) are statements about orbit
  dimensions at special points.
- They are not counts of families.

## 2. Status of 1⇉2⇉3 + 1→3

| | |
|---|---|
| orbit count under the full Δ-compatible group | 873 + 22q at q = 3, 5, 7 |
| 18 families, inner curves, one behaviour each | closed forms; computed at three primes |
| 4 families in the dim-7 groups: kernel tilts with fixed image | one line exact over ℚ and a fold in Lean (v2.44) |
| behaviours with the image A of that line | V_A irreducible, dimension 7; generic behaviour dimension ≥ 1 (= 1 to first order) |

Still open:
- how the 4 families are distributed over images. Is it one V_A per family, or two families inside one V_A (the
  L group grows by 2q)?
- the corresponding statement for the R group.

## 3. Next

1. The images: list the distinct images A among the growing orbits of each group, up to the census group, and run
   V_A for each.
2. The orbit dimensions, not just tangents: integrate the 6 stabiliser directions at a generic point, or compute the
   stabiliser group of A, to turn "= 1 to first order" into a statement.
3. Job 14 on the GPU machine.

Files:
- `analysis/v2_45_{fixed_image,singular_VA,generic_m}.py`;
- `logs/v2_45/` (predictions and outcomes, Singular outputs).
