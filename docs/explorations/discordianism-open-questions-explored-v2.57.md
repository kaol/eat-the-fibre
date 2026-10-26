# Discordianism: Open Questions Explored (v2.57): finitely many census orbits on 1⇉2⇇3 + 1→3 + 2→4 (char 0)

*ChatGPT's Phase 2 brief (at commit a781d2b) set three programs. This round carries out Program 1: the v2.54
orbit-dimension criterion on every vertex block of Q₁ = 1⇉2⇇3 + 1→3 + 2→4 (job 14's quiver, 1220 local-torus orbits
at q = 3, 5, 7). All 90 blocks pass. There are no geometric failures, and nothing is left computationally unresolved.
So the folds of k̄Q₁ form finitely many census orbits in characteristic 0. To get there, the presolve had to learn to
split cases: the last block has no linear pivots left after the quadratic equations. A certificate format and an
independent checker were added, and a full certified rerun is in §3. Programs 2 and 3 are planned in §5, with a first
observation for Program 3.*

---

## 0. Results

The predictions are in logs/v2_57/PREDICTIONS.md, committed before the full run.

| # | Statement | Status |
|---|---|---|
| 59 | Every component of every one of the 90 strata of Q₁ is a certified pass, so the connected census groups have finitely many orbits on Folds(k̄Q₁) in char 0 | **held**: 90 of 90 blocks (`logs/v2_57/tally_q1.txt`) |
| 60 | Regression: the control's 000,001,001 still has exactly one geometric failure; pairs_in_short 100,010,001 still has 16 components, all passing | **held**, also in the new case-split mode (`regression.txt`, `regression_branch.txt`) |

The three outcome classes, reported separately as the brief asks:

| outcome | blocks | components |
|---|---|---|
| certified pass | **90** | 372 components in 89 blocks (353 pass for every r, 19 are points) + 54 leaf primes in the case-split block (53 pass, 1 point) |
| geometric failure | 0 | 0 |
| computationally unresolved | 0 | 0 |

Components by dimension, for the 89 blocks run without case splits: 19 of dim 0, 59 of dim 1, 66 of dim 2, 115 of dim 3, 62
of dim 4, 40 of dim 5, 9 of dim 6, 2 of dim 7.

**A consistency check against the census.** Each component carries a dense orbit (it passes at r = dim Z − 1), and
dense orbits of distinct components are distinct. So the number of components of a block should not exceed its
local-torus count. It does not in any of the 89 blocks: the components total 372, and the local orbits total 1220 at
q = 7.

## 1. The theorem

Let k be algebraically closed of characteristic 0 and Q₁ = 1⇉2⇇3 + 1→3 + 2→4.

**Theorem (char 0).** Folds(kQ₁) is the union of 90 vertex-block strata. On every irreducible component Z of every
stratum, the connected group H_Z = Aut_Δ°·T_{W_Z} has finitely many orbits. So Folds(kQ₁) is a finite union of
census orbits.

*Proof.* This is the v2.54 argument, applied block by block:
- **The groups.** Lie H_Z = 𝔥_Z = Der_Δ ⊕ E_{W_Z}. In characteristic 0 group schemes are smooth, so Lie(Aut_Δ°) is the
  algebra of Δ-compatible derivations. Here dim Der_Δ = 4 and there are 4 classes.
- **H_Z preserves Z.** H_Z preserves the W_Z-compatible folds of the stratum (v2.23), and being connected it preserves
  each irreducible component.
- **The criterion.** H_Z has finitely many orbits if dim Z_r ≤ r for every r < dim Z (v2.54 §1).
- **The certificate for each r.** The certificate is a set F of (r+1)-minors of M_Z(P) = ([h_i, P])_i. Each minor is
  chosen at a rational point where it is nonzero, so it cuts the component on which that point lies. Singular then
  shows over ℚ that Z ∩ V(F) has no component of dimension > r. Every minor vanishes on Z_r, so dim Z_r ≤ r.

The points are used only to choose minors. The final dimension statement is exact. ∎

What the theorem does not say: it gives no orbit count, and it says nothing about finite fields. Over 𝔽_q, Q₁'s count
is 1220 at q = 3, 5, 7, but turning char 0 finiteness into a bound uniform in q is Program 2 (§5).

## 2. How the last blocks were finished

After v2.55, the pipeline handled every stratum of Q₀ (d = 13). Q₁ has d = 18 and up to 236 unknowns per stratum,
and the first full run left 9 blocks unfinished. Each fix is kept as an option of `analysis/v2_55_pipeline.py`:

| step | what | effect |
|---|---|---|
| indexed presolve (`presolve_fast`) | substitutions touch only the equations that contain the variable | 856 s → 118 s on 1000,1000,1000,1000 |
| staged build (`staged_presolve`, default) | the quadratic equations (P² = P, Δ) first, then the cubic closure equations in batches of column pairs, presolving after each batch | 81 of 90 blocks done; 9 killed by memory |
| ring arithmetic in the component stage | sympy expressions of P ran out of memory at d = 18 | (needed for all blocks) |
| Markowitz pivoting (`PRESOLVE_RULE=markowitz`) | substitute where (terms of the value) × (equations containing the variable) is smallest | 8 of the 9 (e.g. 0010,0000,0010,0010: 305 s, 1.1 GB) |
| closure equations built in Singular (`CLOSURE_IN_SINGULAR=1`) | — | **failed**: over 40 min and 3 GB on one heavy block; stopped |
| **case-split presolve** (`BRANCH=1`, `fold_ideal.branch_presolve`) | below | the last block, 0000,0000,0000,0001: 54 leaves, 527 s |

**Why the last block resisted.**
- After the quadratic stage it has 60 unknowns and 339 equations with no linear pivot. 60 of the equations are
  monomials and 100 factor.
- The closure batches then add cubic equations that never become linear. They pile up: 2000 equations and 280,000
  terms after a third of the batches, still growing.
- Singular's minAssGTZ on the quadratic ideal alone did not finish in 500 s.

**The case split.**
- After each presolve, every equation that factors is replaced by its radical: the product of its distinct irreducible
  factors. Monomials are always split. Other equations are split when they have at most 12 terms; they are factored in
  the ring of their own variables, because sympy's dense factorisation over all 236 generators was the first
  bottleneck.
- An equation f₁⋯f_k = 0 with k ≥ 2 splits the branch. Branch i has f_i = 0 and f₁, …, f_{i−1} ≠ 0, so the
  branches are disjoint.
- Known-nonzero factors are divided out of later equations. A branch dies when a known-nonzero polynomial becomes 0.
- Each leaf is decomposed in its own forked child. Its ideal is first saturated by the leaf's inequations.

**Why this is sound.**
- Every rewrite strengthens the equations: substitution, radical, dividing out a nonzero factor, adding a factor. So
  each leaf variety lies in the stratum.
- Every point of the stratum lies in the leaf picked by the first vanishing factor at each split. So the stratum is the
  union of the saturated leaf varieties.
- Its components are therefore among the leaf primes. Every leaf prime passes, so every component passes.

**Two failed versions of the case split.** Both are kept in the log of this round.
- Splitting into overlapping branches, with no inequations, kept the factored equation and looped. With the equation
  dropped, it still produced more than 2960 leaves before it was stopped.
- Factoring with sympy's dense representation spent 278 of 300 s in `factor_list`.

**Regression of the new mode.**
- The control's growing blocks 000,001,001 and 001,001,001 still fail exactly once each, at components of dimension 4
  and 6, as in v2.55.
- pairs_in_short 100,010,001 passes on all 122 leaves.

## 3. Certificates and an independent check

The brief asks for reproducible certificates for every component. The pipeline now writes them when `CERT=file.jsonl`
is set:
- **One stratum record per block** (or per leaf): the coordinates, the substitutions, the presolved equations, the
  inequations, the Der_Δ basis and the class scalings.
- **One record per component:** its prime's generators, its dimension, an independent set, W_Z, and for each r the
  minors, each given by its rows (elements of 𝔥_Z), its columns (entries of P) and the determinant.

`analysis/v2_57_verify_certs.py` checks a certificate file without the pipeline's decomposition steps:
- **(a) The minors.** It rebuilds P from `fold_ideal.stratum` and the recorded substitutions, and 𝔥_Z from the record.
  It recomputes every minor from its recorded rows and columns and compares it with the recorded polynomial.
- **(b) Z lies in the stratum.** In Singular, the presolved equations reduce to 0 modulo the prime, and the prime's
  dimension is as recorded.
- **(c) The dimension bounds.** dim std(prime + minors) ≤ r for every r, computed with `std` and `dim` only, so it does
  not depend on `minAssGTZ`, the routine that once overflowed silently (v2.55).
- **(d) Coverage**, optional (`--coverage`). Every generator of the intersection of the primes lies in the radical of
  the stratum's ideal (Rabinowitsch), so the listed primes cover it.

**The checker has been tested.** On pairs_in_short 100,010,001 and 100,100,100 (case-split mode), all 61 cuts pass,
the minors agree and coverage holds. With one minor deleted from a certificate, it reports `CUTFAIL 1 3`.

**The full rerun of Q₁ with certificates** is in §3a. It uses one presolve mode for 89 blocks (Markowitz, batches of 6
column pairs) and the case split for 0000,0000,0000,0001.

### 3a. The certified rerun

All 90 blocks were rerun from scratch with certificates. They ran in two streams, about 3 h on the 2-core container
(`rerun_q1_s0.txt`, `rerun_q1_s1.txt`).

| check | result |
|---|---|
| tally of the rerun (`tally_q1_rerun.txt`) | 90 of 90 pass, 426 components or leaf primes; **every block has the same number of components as in the first run** (the first run used mixed presolve modes) |
| certificates (`certs_q1_s{0,1}.jsonl.xz`) | 90 blocks: 143 stratum records (89 blocks + 54 leaves), 426 component records (406 pass, 20 points) |
| (a) minors recomputed | all agree |
| (b) components lie in their strata, with the recorded dimensions | all, 0 errors |
| (c) dim(prime + minors) ≤ r by `std` | **1190 of 1190 cuts verified**, 0 failures (`verify_certs_q1_s{0,1}.txt`) |
| (d) coverage (the primes cover the stratum) | 138 of 139 non-point strata or leaves verified within 900 s; 4 more are single points; **not verified** for 1000,0100,0010,0001 (38 components): it timed out at 900 s, and a rerun was stopped after 2400 s of Singular (`verify_certs_coverage_q1_1000010000100001.txt`) |

So the dimension bounds of the theorem rest on two separate computations: `minAssGTZ` in the pipeline and `std` in the
checker. The decomposition itself is checked independently, by coverage, everywhere except block 1000,0100,0010,0001. There,
the claim that its 38 primes are all of its components rests on `minAssGTZ` alone, with the overflow guard of v2.55.
The dimension cuts on those 38 primes are verified.

## 4. What changed in the code

- `analysis/fold_ideal.py`:
  - `_equations` (the fold equations by part);
  - `stage1_presolve`, `staged_presolve`, `presolve_markowitz`;
  - `branch_presolve`, with `_radical_factors`, which factors in per-equation subrings and caches results.
- `analysis/v2_55_pipeline.py`:
  - the environment options `FOLD_STAGED`, `FOLD_BATCH`, `PRESOLVE_RULE`, `BRANCH`, `CLOSURE_IN_SINGULAR` and `CERT`;
  - the ring-arithmetic component stage;
  - single-point strata; leaf saturation.
- `analysis/v2_55_invariant.py` rebuilds P from the ring form.
- New:
  - `analysis/v2_57_tally.py`: classifies every block from the logs into pass / fail / unresolved. It reproduces
    v2.55's control table (8 / 6 / 7) and v2.54's (21 / 0 / 0).
  - `analysis/v2_57_verify_certs.py`.

## 5. Programs 2 and 3: plan and a first observation

### Program 2: uniformity in q, written out

The char 0 theorems for Q₀ and Q₁ do not bound the 𝔽_q counts by themselves. The argument to be written next round
has four steps. Each step names what must be computed.

1. **The group over ℤ.** Write each element of a basis of Der_Δ as the derivative of an explicit one-parameter group
   of substitutions defined over ℤ (for example 1 + tD with D² = 0, if that holds). Then the subgroup H generated
   over 𝔽̄_p is smooth and connected, and Lie H contains the reduction of 𝔥_Z.

   This settles the char p issue: there dim H·P ≥ rank [𝔥_Z, P] still holds, because the orbit dimension is
   dim H − dim Stab ≥ dim H − dim Lie Stab. So a rank certificate bounds Z_r in characteristic p as well. It only
   needs a basis of 𝔥 realised by subgroups. Without that, Lie(Aut_Δ) in char p could be larger than the Lie algebra
   of its reduced group.

   *First computation (exact, `analysis/v2_57_derd.py`, `logs/v2_57/derd.txt`).* For both Q₀ and Q₁, Der_Δ is spanned class by
   class by the derivations a ↦ λ(a)·(a′ − a) on each two-arrow class {a, a′}, with λ any linear form on the class,
   extended to paths by the Leibniz rule. The arrows 1→3 and 2→4 are not moved. These are not nilpotent (D² ≠ 0),
   but D(x)D(y) = 0 for all basis elements x, y.

   They are the tangent vectors of the substitution groups x ↦ x + λ(x)(a′ − a), one copy of Aff₁ = 𝔾_a ⋊ 𝔾_m per
   class:
   - the unipotent part, with λ(a′ − a) = 0;
   - the torus part, with λ(a′ − a) = 1.

   Both are defined over ℤ and smooth in every characteristic. It remains to check that they commute with Δ, as
   identities over ℤ[t, s].
2. **Good primes.** The certificates of §3 are polynomial identities over ℚ:
   - presolve and substitutions;
   - I ⊆ P_Z, coverage, and the vanishing of the defect normal forms that define W_Z;
   - the minors;
   - dim(P_Z + F_r) ≤ r.

   Let N be the product of all denominators, together with the primes at which a Gröbner basis of P_Z + F_r over ℤ
   (Singular `ring r = integer`) changes leading terms. For p ∤ N every certificate reduces mod p.

   One more point is needed. A component Z can become reducible mod p. Each geometric component Y of the stratum mod p
   lies in some Z mod p, and its own group H_Y ⊇ H_Z acts on it: H_Y is connected and acts on the W_Y-compatible locus,
   of which Y is a component. Its rank loci are contained in those of Z mod p. So the criterion passes on Y.
3. **From k̄-orbits to 𝔽_q-orbits.**
   - Lang's theorem: a Frobenius-stable orbit H·x contains at most |π₀(Stab x)| orbits of H(𝔽_q).
   - The bound must be uniform in p and x. Stab x is cut out of H by equations of bounded degree in a fixed affine
     space, so a Bézout inequality (Heintz) bounds its number of components independently of p and x.
   - The number of geometric orbits is bounded the same way, by the degrees of the varieties Z ∩ V(F_r).

   This gives an explicit bound on the census count for p ∤ N.
4. **Bad primes.** For each p | N, rerun the pipeline over 𝔽_p:
   - Singular in characteristic p, with the decomposition over 𝔽_p;
   - each 𝔽_p-component is a Galois orbit of geometric components of one dimension;
   - the char p rank criterion of step 1.

   Two outcomes are possible. A pass gives the bound at that p as well. A failure there, if confirmed by an invariant,
   would be a genuine exceptional characteristic. The brief asks for exactly this distinction, and it is not assumed
   away.

Q₀ comes first, then Q₁.

### Program 3: the leaf extension Q₀ → Q₁ = Q₀ + 2→4. A first observation

The brief warns that a fold of Q⁺ need not restrict to a fold of Q. The vertex blocks already show this. The natural
candidate map is the quotient: kQ₀ = kQ₁/(e₄), with the restriction π P ι. Its vertex block is the top-left 3 × 3
corner of Q₁'s vertex block.

Comparing the 90 blocks of Q₁ with the 21 of Q₀ (from the census tables; every Q₁ block is nonempty in char 0 by §1):
- **84 Q₁ blocks have a Q₀ block as their corner.** Every Q₀ block occurs.
- **6 Q₁ blocks have a corner that is not a Q₀ block.** These 6 carry 16 of the 1220 local orbits. In each of them
  P(e₄) has a component on e₁, e₂ or e₃ (column 4 is not supported on e₄ alone), so P does not preserve the ideal (e₄),
  and π P ι is not a fold.

So the extension lemma cannot be proved through a restriction map on all of Folds(Q₁). The candidate structure is:
- a part that preserves (e₄) and fibres over Folds(Q₀);
- a small part with new vertex blocks that has to be handled on its own.

Over the 21 Q₀ blocks the fibres vary a lot. One example: the Q₀ block with 42 orbits (vertex block the identity)
carries Q₁ blocks with 144 local orbits. The leaf's coordinates on a fibre are the entries of P in rows and columns
that involve e₄ and the arrow 2→4. Next: find which of the 84 blocks preserve (e₄) as an ideal, not only on vertices.
Then decide whether the census group acts with finitely many orbits on each fibre over a fixed Q₀ orbit. That would
be a relative version of the v2.54 criterion, with the base orbit fixed and only the leaf coordinates moving.

## Files

- analysis/fold_ideal.py, v2_55_pipeline.py, v2_55_invariant.py (changed); v2_57_tally.py, v2_57_verify_certs.py (new)
- logs/v2_57/:
  - PREDICTIONS.md, with its outcome section;
  - the run logs: pipeline_q1_part0.txt, part1.txt, retry*.txt, *_markowitz*.txt and
    pipeline_q1_0000000000000001_branch.txt (the last block); the earlier attempts are in first_run/ and
    pipeline_q1_singular_closure.txt;
  - regression.txt, regression_branch.txt, regression_singular_closure.txt;
  - tally_q1.txt;
  - the certified rerun: rerun_q1_s0.txt and rerun_q1_s1.txt, with certs_q1_s0.jsonl and certs_q1_s1.jsonl.
