# Discordianism: Open Questions Explored (v2.59): the uniform bound for 1⇉2⇇3 + 1→3 + 2→4

*v2.58 carried out Program 2 for Q₀. This round does the same for Q₁ = Q₀ + 2→4, the 90-block quiver of v2.57. The
census count over 𝔽_q is bounded by one constant for every q, with no exceptional characteristic. Along the way:
- the vertex-block theorem of v2.52 extends to Q₁ over every field: exactly the 90 census blocks;
- the case-split presolve of v2.57 gets an integral mode, so that it is valid over ℤ;
- the good-prime check now tests each candidate prime only on the claims it is a candidate for.

With v2.57 and v2.58, both Q₀ and Q₁ now satisfy Conjecture 3ᶜ as theorems, modulo the trusted Singular computations
listed in v2.58 §4.*

---

## 0. Results

The predictions are in logs/v2_59/PREDICTIONS.md, committed before the runs.

| # | Statement | Status |
|---|---|---|
| 63 | The certified rerun of Q₁, now with pivots and case-split constants, reproduces v2.57; pivots ±1; the case-split constants involve no prime but possibly 2 | **held**. 90 / 90 pass, with the same component count in every block (`tally_q1.txt`). Pivots are ±1 in 89 blocks. In the case-split block, pivots ±1/2 and a constant 2 came from monic normalisation. The new integral mode removes them: all pivots and constants are ±1 |
| 64 | Every specialisation claim of Q₁ holds at every prime; candidates at most 2 and 3; the coverage of 1000,0100,0010,0001 may stay unresolved | **held**, risk avoided. 5413 claims in 139 strata and leaves. Candidates: 2 (5 blocks) and 3 (1 block), and every candidate claim holds there by direct computation. The ℚ-coverage chain of 1000,0100,0010,0001 finished this time, built pairwise, where v2.57's single `intersect` had timed out |
| 65 | The vertex blocks of Q₁ are exactly the 90 census blocks over ℚ and over every 𝔽_p | **held**. Elimination of 704 equations in 108 unknowns gives exactly the 90 census blocks. Over ℤ: every leading coefficient is 1, the 133 point-ideal generators all have ρ^k ∈ E_ℤ with k ≤ 2, and the point ideal has a ℤ-basis with leading coefficients 1 and vdim 90. Direct elimination over 𝔽₂ and 𝔽₃ gives 90 points each |

**Theorem (Conjecture 3ᶜ for Q₁).** For every prime power q, the census count of folds of 𝔽_qQ₁ is at most
B ≈ 10^96.8. The census count is the number of local-torus classes. The proof is v2.58 §4 verbatim, with Q₁'s
certificates. The observed count is 1220 at q = 3, 5, 7.

## 1. The vertex blocks of Q₁ in every characteristic

`analysis/v2_58_vertex_integral.py pairs_in_short_tail … --arrow 1,3,4,5` builds the v2.52 system for Q₁:
- P_K (6 × 6) idempotent;
- products closed on K, written with the integral rows of B⁻¹;
- the Δ-conditions of one arrow from each class: a, b, c and the new d = 2→4.

That is 704 equations in 108 unknowns, eliminated to the 16 vertex-block entries.

| step | result |
|---|---|
| over ℚ (`vertex_char0_q1.txt`) | zero-dimensional, vdim 105, **90 points, exactly the 90 census blocks** (orientation checked entrywise) |
| over ℤ (`vertex_integral_q1.txt`) | strong basis of the elimination ideal: 142 elements, all leading coefficient 1 |
| radical certificate | the ideal R of the 90 points over ℚ, plus v³ − v for each entry (entries are in {−1, 0, 1}): every one of its 133 generators has ρ^k ∈ E_ℤ with k ≤ 2, by exact reduction over ℤ |
| R over ℤ | a strong basis with every leading coefficient 1 and vdim 90 over ℚ. So for every p the quotient by R mod p has dimension ≤ 90, and V(R mod p) is exactly the reductions of the 90 points |
| direct, over 𝔽₂ and 𝔽₃ | dimension 0, vdim 105, 90 points each (the blocks stay distinct even mod 2) |

v2.58's Q₀ certificate used v² − v, which suffices there because Q₀'s blocks are 0/1 matrices. Q₁ has blocks with
entries −1 (written m), so the script now adds v³ − v.

## 2. The case split over ℤ

v2.57's case-split presolve normalised factors to be monic and accepted any constant pivot. On block
0000,0000,0000,0001 this introduced ±1/2 pivots and a constant 2. At p = 2 that would leave the leaf cover unproved:
replacing 2x − y by x − y/2 is not a valid step mod 2.

**The integral mode (`FOLD_INTEGRAL=1`).**
- Presolve pivots must be ±1.
- Factors and new equations are normalised to primitive integral polynomials, not monic ones.
- Every constant whose primes could break a step mod p is recorded (`BRANCH_CONSTS`, with `BRANCH_TAGS` saying
  where). These are:
  - the content c₀ in e = c₀ ∏ f_i^{k_i};
  - the denominators of factors;
  - the contents removed when making a polynomial primitive;
  - constants that made a branch empty.

**Result.** The block again has 54 leaves (53 pass, 1 point), and every pivot and recorded constant is ±1. Each step
is then an identity over ℤ:
- e = ±∏ f_i^{k_i};
- splitting a product;
- dividing out a known-nonzero factor;
- substituting x = ±(rest).

So the leaves cover the stratum over every 𝔽_p.

In the earlier, non-integral run, the one remaining 2 was a factor's leading coefficient (for example in 2x − y). That
does not matter once nothing is divided by it, and the record was dropped.

**Coverage for a leaf.** A leaf variety is V(rem_L) minus V(∏ nz_L). Its closure lies in ∪ Z_p once
(k·∏nz)^e ∈ (rem_L) for each generator k of the component chain. This is checked by exact reduction modulo a strong
Gröbner basis of rem_L over ℤ (e ≤ 4), so it holds at every p with no candidate prime. Rabinowitsch over ℤ is the
fallback.

## 3. The specialisation run (`logs/v2_59/goodprimes_q1_*.txt`)

| part | strata / leaves | claims | candidate primes | direct checks |
|---|---|---|---|---|
| 86 blocks with positive-dimensional strata (3 more blocks are single points) | 86 | 5 099 | 2 in 1000,0000,0000,0000, 1000,0000,0000,0001, 1000,0100,0010,0000, 1000,0100,0010,0001 and 0000,0100,0000,0001; 3 in 0000,0100,0000,0001 | all hold. Full cross-check at p = 2 on every stratum; it timed out on 1000,0100,0010,0001, where the targeted check at 2 passes (`goodprimes_q1_1000010000100001_targeted.txt`) |
| block 0000,0000,0000,0001, integral case split | 53 leaves (+1 point) | 314 | none | full cross-check at 2: all hold |
| **total** | 139 | **5413** | {2, 3} | **no exceptional prime** |

**Why the candidates appear.** In 1000,0100,0010,0001, the prime 2 comes only from the leading coefficients of two
dimension cuts (components 26 and 38, r = 1). Over 𝔽₂ these dimensions are still ≤ 1.

**What changed in the checker.**
- Each candidate prime is now checked only on the claims it is a candidate for.
- The full check of every claim at the `--always` primes is a separate cross-check. A timeout there is reported as an
  unresolved cross-check, not as an exceptional prime.
- The negative test (every minor × 5) still fails at 5, at the same 9 cuts.
- Presolve pivots and case-split constants not equal to ±1 are reported as candidates that a direct check cannot
  settle, since they invalidate the coordinates themselves. None remain.

**Speed.** v2.59 passes integral strings to Singular without sympy parsing. Some Q₁ strata have 3833 presolved
equations totalling 17 MB of text, which sympy did not parse in 20 minutes.

## 4. The bound (`analysis/v2_58_bound.py`, `logs/v2_59/bound_q1.txt`)

- Σ over strata and leaves of D_V^{n_V} · max_Z o_Z ≤ 10^85.9.
- Stabilisers sit in H̃ ⊂ A¹⁴: 2 classes × 3 coordinates for G, plus 4 classes × 2 for T. Their equations have
  degree ≤ 6, since paths have length ≤ 3. So |π₀| ≤ 6¹⁴ ≈ 10^10.9.

Together, c(q) ≤ 10^96.8 for every q. As for Q₀, the bound says nothing about the true count (1220). It only shows
that a uniform bound exists, in every characteristic.

## 5. Where this leaves the programs

| program | status |
|---|---|
| 1. exact char 0 finiteness for Q₁ | done (v2.57) |
| 2. uniformity in q, including exceptional primes | done for Q₀ (v2.58) and Q₁ (v2.59). No exceptional prime in either. Trusted: Singular's `std` and `eliminate` over ℤ, `minAssGTZ` and `intersect` over ℚ, each cross-checked by direct computation at small primes |
| 3. leaf extension and structural induction | next. The v2.57 observation stands: 6 of Q₁'s blocks restrict to no Q₀ block |

The method is now a complete per-quiver decision procedure for Conjecture 3ᶜ, in all characteristics:
1. the vertex-block elimination over ℤ;
2. the char 0 pipeline with certificates;
3. good primes and direct checks;
4. Lang and Bézout.

What it lacks is a reason why it should succeed on every unobstructed quiver. That is Program 3.

## Files

- analysis/fold_ideal.py: integral mode (`FOLD_INTEGRAL`), and `BRANCH_CONSTS` / `BRANCH_TAGS`.
- analysis/v2_55_pipeline.py: certificates carry the case-split constants.
- analysis/v2_58_goodprimes.py: leaves (inequations), the power-membership coverage, targeted checks, `--skip`,
  fast integral strings.
- analysis/v2_58_vertex_integral.py: v³ − v, and a ℤ-basis of the point ideal.
- logs/v2_59/:
  - PREDICTIONS.md;
  - rerun_q1_s{0,1}.txt, certs_q1_s{0,1}.jsonl.xz, tally_q1.txt;
  - pipeline_q1_0000000000000001_integral.txt, certs_q1_branch_integral.jsonl;
  - vertex_char{0,2,3}_q1.txt, vertex_integral_q1.txt;
  - goodprimes_q1_s0_first.txt, goodprimes_q1_s0.txt, goodprimes_q1_s1.txt, goodprimes_q1_branch.txt,
    goodprimes_q1_1000010000100001_targeted.txt;
  - bound_q1.txt.
