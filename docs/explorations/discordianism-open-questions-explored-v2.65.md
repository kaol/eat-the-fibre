# Discordianism: Open Questions Explored (v2.65): the uniform bound for Q₀ + 3→4

*v2.64 proved that Q₂ = 1⇉2⇇3 + 1→3 + 3→4 has finitely many census orbits in characteristic 0. The GPU then found
1773 classes at q = 3, 5 and 7. This round applies the uniformity machinery of v2.58–v2.59 to Q₂ unchanged, with no new
code. The certificates specialise to every prime, and there is no exceptional characteristic. So Conjecture 3ᶜ holds for
Q₂ as well, under the same trusted Singular computations as for Q₀ and Q₁ (v2.58 §4). That makes three unobstructed
quivers with a uniform bound: Q₀, Q₁ and Q₂.*

---

## 0. Results

The predictions are in logs/v2_65/PREDICTIONS.md, committed before the runs.

| # | Statement | Status |
|---|---|---|
| 81 | Unimodular integral structure; pivots and case-split constants ±1; the 99 vertex blocks are the only ones over every field | **held** |
| 82 | Every specialisation claim holds at every prime; candidates at most 2 and 3, each passing a direct check | **held**: 6423 claims in 128 strata and leaves; candidates only 2, in 6 blocks, and all pass there; full cross-check at 2 on every stratum, none unresolved |

**Theorem (Conjecture 3ᶜ for Q₂).** For every prime power q, the census count of folds of 𝔽_q(1⇉2⇇3 + 1→3 + 3→4) is
at most B ≈ 10^65. The proof is that of v2.58 §4 with Q₂'s certificates. The observed count is 1773 at q = 3, 5, 7.

## 1. The checks

| step | result | file |
|---|---|---|
| Smith form of Δ; Kmat; det [Kmat \| e_C] | elementary divisors 1; saturated; det = −1 | `integral.txt` |
| census group and weightings over ℤ | all identities hold (Φ multiplicative and Δ- and Δ_c-compatible over ℤ[g]; φ_w; Φφ_w = φ_wΦ) | `integral.txt` |
| presolve pivots, case-split constants (99 blocks, including the integral case split of 1000,010m,010m,0000) | 64 321 pivots, all ±1; constants ±1 | v2.64 certificates |
| vertex blocks over ℤ (688 equations, 88 unknowns) | strong basis of the elimination ideal: 143 elements, all leading coefficient 1; 133 point-ideal generators with ρ^k ∈ E_ℤ, k ≤ 2; the point ideal's ℤ-basis has leading coefficients 1 and vdim 99 | `vertex_integral.txt` |
| direct elimination over 𝔽₂, 𝔽₃ | 99 points each | `vertex_char{2,3}.txt` |
| specialisation (DIM, DEF, COV claims) | 6423 claims; candidate 2 in 1000,0000,0000,0000, 1000,0100,0010,0000, 1000,0000,0000,0001, 1000,0100,0010,0001, 1000,0001,0001,0001 and 100m,0000,0000,0000; all hold over 𝔽₂ | `goodprimes_{aa,ab,rerun}.txt` |
| Bézout bound | Σ D_V^{n_V}·max o_Z ≤ 10^56.8; stabilisers in A¹⁴ of degree ≤ 4 (paths of length ≤ 2), so 4¹⁴ ≈ 10^8.4; **B ≈ 10^65.2** | `bound.txt` |

## 2. Status of the three programs

| program | Q₀ | Q₁ = Q₀ + 2→4 | Q₂ = Q₀ + 3→4 |
|---|---|---|---|
| 1. exact char 0 finiteness | v2.54 | v2.57 | v2.64 |
| 2. uniform in q, all characteristics | v2.58 | v2.59 | v2.65 |
| census (constant in q) | 366 (q = 3, 5, 7, 11) | 1220 (q = 3, 5, 7) | 1773 (q = 3, 5, 7) |
| 3. leaf extension from Q₀ | — | all four kinds explained (v2.60–v2.63) | the relative criterion passes; the mixed kind has rank > 1 and the kQ-preserving λ = 0 part has families (v2.64) |

The per-quiver procedure is now routine:
1. vertex elimination;
2. the char 0 pipeline with certificates, case splits where the presolve runs out of memory;
3. good primes;
4. Lang and Bézout.

For the "if" direction of Conjecture 3ᶜ in general, Program 3's structural question remains: why should step 2 always
pass on an unobstructed quiver?

## Files

- logs/v2_65/:
  - PREDICTIONS.md;
  - integral.txt;
  - vertex_integral.txt, vertex_char{2,3}.txt;
  - goodprimes_{aa,ab,rerun}.txt;
  - bound.txt.
