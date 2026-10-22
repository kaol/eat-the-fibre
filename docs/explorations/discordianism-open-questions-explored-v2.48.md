# Discordianism: Open Questions Explored (v2.48): the divisor is irreducible and also carries behaviours

*v2.47 found the rank-drop divisor D of the fixed-image variety V_A on 1⇉2⇉3 + 1→3 (L group). This round shows that D
has a unique divisorial component and that its one-parameter family is a family of behaviours, like V_A's generic
one. Both pre-registered predictions held.*

---

## 0. Results

| # | Statement | Status |
|---|---|---|
| 1 | The rank-drop locus {P ∈ V_A : rank [𝔥′, P] ≤ 5} has exactly one component of dimension 6, the same prime in the three relevant charts, plus components of dimension 5 | Singular over ℚ (`minAssGTZ`), equality of primes by mutual reduction (prediction 29) |
| 2 | P₁ = P₀ + N (the v2.44 line at t = 1) is a smooth point of that component, with dim T_{P₁}D = 6 | exact Jacobian rank over ℚ (prediction 30) |
| 3 | dim(T_{P₁}D ∩ [Der(kQ), P₁]) = 5, so exactly one tangent direction of D at P₁ is not produced by any derivation | exact rank over ℚ (prediction 30) |

## 1. Consequence

Rank conditions are open, so the counts of item 3 hold on a dense open subset of the component D (irreducible). On
that subset, each Aut(kQ)-orbit meets D in a set of dimension at most 5, because its tangent lies in [Der, P] ∩ T D.
So D is not swept by finitely many orbits of dimension ≤ 5 and any lower-dimensional pieces: **the behaviours on D
form a family of dimension at least 1.** The census-group orbits on D have dimension 5 (v2.47), so the count is
exactly one census parameter, at least one behaviour parameter.

Exactness ("= 1") would follow as in v2.46 from a dominant slice inside D. That is not attempted here.

## 2. The picture for 1⇉2⇉3 + 1→3

| family type | number | behaviour |
|---|---|---|
| inner curves Inn(1 + t·x)·P₀ | 18 | one behaviour each (inner automorphisms) |
| generic stratum of V_A, V_B (kernel tilts) | 2 (L, R) | one-parameter families of behaviours (proved, char 0, v2.46–47) |
| divisor D of V_A (and, by mirror, of V_B) | 2 | families of behaviours, ≥ 1 parameter (this round; B by mirror symmetry, not computed separately) |

The count is 873 + 22q at q = 3, 5, 7. The correspondence of each piece with q census orbits remains the heuristic
of v2.47.

## 3. Next

1. The divisor of V_B, computed directly rather than by symmetry.
2. Lean closed forms for the 18 inner families on this quiver.
3. Job 14 (GPU).

Files: `analysis/v2_48_divisor_components.py`, `logs/v2_48/`.
