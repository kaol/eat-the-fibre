# v2.59 predictions (committed before the runs)

Program 2 for Q₁ = 1⇉2⇇3 + 1→3 + 2→4, after v2.58 did Q₀. The integral structure of Q₁ is already checked
(logs/v2_58/integral_q1.txt: Smith form 1, saturated K-basis, det B = −1, group identities over ℤ).

63. The certified rerun of Q₁ (now recording presolve pivots, and the constants of the case-split steps) reproduces
    v2.57: 90 blocks, the same component counts, all passing; every presolve pivot is ±1, and the case-split constants
    of block 0000,0000,0000,0001 involve no prime other than possibly 2.
64. Every specialisation claim of Q₁ holds at every prime (no exceptional prime), as for Q₀; candidate primes, if any,
    are 2 or 3. The one risk is computational: the coverage chain of 1000,0100,0010,0001 (38 components) did not finish
    over ℚ in v2.57, and may stay unresolved.
65. The vertex blocks of Q₁: the analogue of v2.52's system eliminates to exactly the 90 census blocks over ℚ, and the
    same holds over every 𝔽_p (elimination over ℤ with unit leading coefficients and a radical certificate).

## Outcome

- **63 held.** The rerun gives 90/90 passes with v2.57's per-block component counts. Pivots are ±1 except in the
  case-split block (±1/2 and a constant 2, from monic normalisation); the integral mode (FOLD_INTEGRAL=1) makes every
  pivot and recorded constant ±1 there.
- **64 held.** 5413 claims over 139 strata/leaves; candidates 2 (five blocks) and 3 (one block), all holding by direct
  check; the coverage chain of 1000,0100,0010,0001 finished (pairwise intersections).
- **65 held.** Elimination: the 90 census blocks exactly, over ℚ; over ℤ all leading coefficients 1 with a radical
  certificate; 90 points over 𝔽₂ and 𝔽₃.
