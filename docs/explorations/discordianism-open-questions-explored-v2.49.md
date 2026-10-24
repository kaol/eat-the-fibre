# Discordianism: Open Questions Explored (v2.49): the inner families in Lean, and the R divisor computed directly

*v2.48 computed the rank-drop divisor of the fixed-image variety V_A (L group) of 1⇉2⇉3 + 1→3 and inferred the R
group's by mirror symmetry. This round computes D_B directly. It also generates Lean files for the 18 inner families,
so every family piece of that quiver now has a machine-checked or exact-ℚ certificate. Both pre-registered predictions
held (logs/v2_49/PREDICTIONS.md).*

---

## 0. Results

| # | Statement | Status |
|---|---|---|
| 1 | The 18 inner families P_t = Inn(1 + t·x)·P₀·Inn(1 + t·x)⁻¹ of 1⇉2⇉3 + 1→3 (x = αδ_a + βδ_b + γδ_aδ_b) are folds over ℚ[t], with Φ_t = Inn(1 + t·x) invertible | exact sympy check in the generator, all 18 (prediction 31) |
| 2 | Families 1 and 2: `famk_isFold` (a fold over every commutative ring, every t) and `famk_one_orbit` (Φ_t P₀ = P_t Φ_t, Φ_t a multiplicative bijection) | Lean, `decide +kernel` reflection (prediction 31). Family 1: 754 s; family 2: FAM2 |
| 3 | D_B = {P ∈ V_B : rank [𝔥′, P] ≤ 5} (𝔥′ for the 6118 line) has a unique 6-dimensional component: the same prime in charts 0, 2, 3; the other charts have incidence dimension 5 | Singular over ℚ, `minAssGTZ` and mutual reduction (prediction 32) |
| 4 | The 6118 line P₀ + t·N meets D_B at t = −1 (rank 5 there; rank 6 at t ∈ {±2, ±3, 1, ±1/2, 1/3}); P₋₁ is a smooth point of D_B with dim T = 6 and dim(T ∩ [Der, P]) = 5 | exact rank over ℚ (prediction 32) |

## 1. Consequences

Item 4 gives for D_B the argument v2.48 gave for D_A: rank conditions are open and the component is irreducible. So on
a dense open subset of D_B, each Aut(kQ)-orbit meets D_B in dimension ≤ 5. **The behaviours on D_B form a family of
dimension ≥ 1.** The mirror inference of v2.48 is now a direct computation.

One difference from the L group: on the L line the rank-drop point was t = 1. On the R line it is t = −1. This is
an observation about the two chosen lines (their N are normalised independently); it has not been checked whether the
mirror map takes one normalisation to the other.

Items 1–2 turn the inner families from census observations (they meet q orbits at q = 3, 5, 7) into closed forms over
ℤ[t]. Every member is a fold over every commutative ring, and each family is one inner orbit. Two of the 18 files were
kernel-checked here. The other 16 are generated and registered (build.sh, EatTheFibre.lean, Audit.lean) and need about
5 GB each.

## 2. The picture for 1⇉2⇉3 + 1→3

| family type | number | certificate |
|---|---|---|
| inner curves Inn(1 + t·x)·P₀ | 18 | Lean files (2 built, 16 pending); exact ℚ[t] check for all 18 |
| generic stratum of V_A, V_B (kernel tilts) | 2 | dominant slice, char 0 (v2.46; mirror v2.47) |
| divisors D_A, D_B | 2 | unique component, behaviour family of dimension ≥ 1 (v2.48, this round) |

## 3. For the Strix machine

To build the remaining 16 family files:

```
cd lean && for k in $(seq 3 18); do lake build EatTheFibre.Gen.two_pairs_short_fam$k; done
```

(Lean, not ROCm, so no `HSA_OVERRIDE_GFX_VERSION` is needed.) Afterwards the full `./build.sh` runs the audit with the new
`#print axioms` lines.

## 4. Next

1. Exactness on D_A and D_B: is there a dominant slice inside D, as in v2.46, giving exactly one behaviour parameter?
2. Job 14 (pairs_in_short_tail over 𝔽₅), when its results arrive.

## Files

- lean/gen/gen_short_families.py → lean/EatTheFibre/Gen/two_pairs_short_fam{1..18}.lean
- analysis/v2_48_divisor_components.py (new env `T`: test point P₀ + T·N)
- logs/v2_49/: PREDICTIONS.md, divisor_B.txt, components_B.txt, lean_fam1.log, lean_fam2.log

---

## Addendum (2026-10-24): all 18 family files built

The remaining 16 files (fam3–fam18) were built on the Strix machine, followed by a full `./build.sh`. The run had no
errors and ended "Build completed successfully (9118 jobs)". The audit now has 184 results, up from 148: 183 depend
only on standard axioms ([propext, Classical.choice, Quot.sound] or a subset) and one depends on none. All 36 new
results (`famk_isFold`, `famk_one_orbit`, k = 1..18) depend on [propext, Classical.choice, Quot.sound].

So the inner-curve row of §2 is complete: **all 18 inner families of 1⇉2⇉3 + 1→3 are machine-checked.** Each is a fold
over every commutative ring for every t, and each is one inner orbit. Excerpt: logs/v2_49/strix_build_audit.txt.
