# v2.42 predictions (committed before any computation)

Quiver 1⇉2⇉3 + 1→3 (`two_pairs_short`). Orbits are taken under the full Δ-compatible group found in v2.28: G, the
local torus, and the shear U = {c ↦ c + s·δ_aδ_b}. Counts are 939 (𝔽₃) and 983 (𝔽₅) (logs/unipotent_merge_two_pairs_short.log).

1. The inner-curve search of v2.36 (`analysis/family_curves.py`), run against the merged orbit labels at q = 5, finds
   exactly 22 families: pairwise disjoint curves P_t = Inn(1 + t·x)·P₀·Inn⁻¹, each meeting ≥ 4 merged orbits.
   This is the slope (983 − 939)/(5 − 3) = 22.
2. Their integer closed forms, reduced mod 3, meet exactly 3 merged orbits each, pairwise disjoint. That would make
   939 − 3·22 = 983 − 5·22 = 873 orbits constant.
3. All families lie along x ∈ {δ_a, δ_b, δ_aδ_b}, as on 1⇉2⇉3. No direction outside the span of δ_a, δ_b and their
   product is needed.

Caveat recorded in advance: censuses over prime fields cannot see a family defined only over 𝔽_{p²} and paired with
its Frobenius conjugate. Such a pair contributes finitely many 𝔽_p-points. So "no other families" claims are about
prime-field counts only.

## Outcomes of 1–3 (before the targeted search)

1. **Failed.** 35 curves meet ≥ 3 orbits at q = 5. Reducing all 35 mod 3 separates them:
   - 18 are genuine families: each meets q orbits at both primes, pairwise disjoint.
   - 17 are bounded curves: they meet 3 orbits at both primes.
   - The 18 explain 36 of the growth of 44, so 4 families are missing.
2. **Held for the 18 found** (they meet 3 orbits mod 3, pairwise disjoint). The constant count does not follow, because
   of the 4 missing families: the orbits they miss go 834 → 842.
3. **Failed as written.** Two of the 18 need x = δ_a − 2δ_aδ_b (reps 350 and 351, classes R and L). Whether that is
   equivalent to a pure direction is not checked.

The missing growth sits in exactly two signature groups, an L/R pair with image dim 7 and kernel dim 5. Each goes
23 → 27, i.e. two families each (`growth_3_5.txt`).

## Prediction 4 (committed before the targeted search)

4. Inner curves Inn(1 + t·x)·P₀ with x ∈ rad, *including the arrow c: 1→3 and the paths a_ib_j*, found by a tangent
   computation at liftable representatives of those two groups, give exactly the 4 missing families: 2 L and 2 R,
   each meeting q orbits at q = 3 and 5.

### Outcome of 4 (tangent search) and prediction 5

- 4 **failed**: at every liftable representative of the target groups, no nonzero inner direction x ∈ rad is
  tangent to Folds to first order (`tangent_search.txt`).
- (b, e) on these groups (`be_targets.txt`): b > 0 almost everywhere, also at constant orbits, so it does not
  discriminate. e = 1 at 4 orbits at both primes.

5. (committed before the run) The missing families are orbits of the shears σ_w : c ↦ c + t·w, with
   w ∈ ⟨a_ib_j⟩ \ k·δ_aδ_b. These are automorphisms of kQ but not Δ-compatible. They may be composed with an
   inner δ-direction. Each such curve, from a liftable representative, meets q orbits at q = 3 and 5, and together
   they account for the growth 23 → 27 in both target groups.

### Outcome of 5 and prediction 6

- 5 **failed**: no shear curve σ_w, alone or composed with an inner δ-direction, from any liftable representative of
  the two target groups is a fold curve at all (`shear_search_{0,1}.txt`). The numeric filter was sanity-checked on
  a known family (rep 75 along δ_a passes).
- The vertex block does not carry the parameter: each target group has one block, and it grows by 4
  (`vblock_targets.txt`).
- Two primes cannot distinguish "2 families of q orbits per group" from other q-dependence, e.g. terms in (q − 1)/2.

6. (committed before the q = 7 census of 1⇉2⇉3 + 1→3) The merged orbit count at q = 7 is
   873 + 22·7 = **1027**. Within that:
   - the 18 integer closed forms meet 7 orbits each, pairwise disjoint;
   - the remaining growth is 4 per target group per unit of q.

### Status of 6

Not decided in the container. The q = 7 census was OOM-killed three times:
1. `expand` sized its slices by the int8 images while `pack` makes an int64 copy, 8× larger. Fixed in `fibre_torch.py`.
2. The engine stack held every child chunk of a parent batch. Made lazy in `fibre_torch.py`; this is an improvement,
   but it was not the cause.
3. Storing every fold of every fibre: cgroup OOM at 6.1 GB RSS.

`finish()` materialises all folds and keys them by bytes, so q = 7 needs a machine with more memory. It was added as
run_batch.py job 16 with this prediction (unmerged count 2012, merged 1027). Logs:
`logs/run_reps_two_pairs_short_7{_oom,_oom2,}.log`.
