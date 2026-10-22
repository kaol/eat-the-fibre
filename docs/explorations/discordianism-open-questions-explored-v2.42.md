# Discordianism: Open Questions Explored (v2.42): the families of 1⇉2⇉3 + 1→3: 18 found, 4 missing

*1⇉2⇉3 + 1→3 is the quiver where Theorem 2 fails but Theorem 4 survives. v2.28 found that its Δ-compatible
symmetries include the shear c ↦ c + s·δ_aδ_b. Under the full symmetry group (the census group G, the class weightings,
and the shear) its orbit count goes 939 → 983 from 𝔽₃ to 𝔽₅, a slope of 22. This round looks for the 22 families.*

*The inner-curve search finds 18, cross-checked at both primes. The other 4 are not inner curves, not shear curves and
not composites of the two, and their nature is open. Five predictions were committed before the computations:*
- *1, 3, 4 and 5 failed;*
- *2 held only for the families found;*
- *6, the q = 7 count, is handed to the GPU machine after three out-of-memory failures in the container.*

---

## 0. Summary

| # | Statement | Status |
|---|---|---|
| 1 | The inner-curve search (v2.36), run against orbit labels merged under the full group, finds 35 curves meeting ≥ 3 orbits at q = 5 | computed (`curves_two_pairs_short_5_merged.txt`) |
| 2 | Reduced mod 3, 18 of them meet q orbits at both primes, pairwise disjoint: **18 genuine families**. The other 17 meet 3 orbits at both primes: bounded curves | computed (`crossprime_q{3,5}.txt`) |
| 3 | The 18 families lie along δ_a, δ_b, δ_aδ_b, and in two cases (classes L and R) along δ_a − 2δ_aδ_b | computed |
| 4 | The remaining growth, 834 → 842 orbits, sits in exactly two signature groups: an L/R pair with image dimension 7 and kernel dimension 5, each going 23 → 27 | computed (`growth_3_5.txt`) |
| 5 | In those groups, no nonzero inner direction x ∈ rad is tangent to Folds at any liftable representative. No shear curve c ↦ c + t·w (w ∈ ⟨a_ib_j⟩ ∖ k·δ_aδ_b), alone or composed with an inner δ-direction, is a fold curve | computed (`tangent_search.txt`, `shear_search_{0,1}.txt`) |
| 6 | The vertex block is not the modulus: each target group has one vertex block, which grows by 4 | computed (`vblock_targets.txt`) |

## 1. What the 18 are

They are the analogues of the 1⇉2⇉3 families: one-parameter inner orbits Inn(1 + t·x)·P₀ with image dimension 2 or 3.
- 13 are unclassed (none of the module conditions L, R, K, E holds).
- 3 are class L and 2 class R.

The search found no LRKE family. On 1⇉2⇉3 that was Theorem 2's family 13, and the joining arrow kills Theorem 2
(v2.28).

## 2. What the missing 4 are not

They live in two mirror-image groups of large-image folds, dim im P = 7 out of 12. On 1⇉2⇉3 every family had
dim im P ≤ 3, so nothing in 1⇉2⇉3's picture predicts them. Three ansätze failed:
- inner curves along any direction in rad (tangent computation);
- the non-Δ-compatible shears, which exist only because the paths a_ib_j are parallel to c;
- composites of those shears with inner δ-directions.

v2.34's invariants (b, e) do not discriminate here: b > 0 at most representatives, including those of constant
orbits.

Two primes also cannot rule out another form of growth. For example, terms in (q − 1)/2, such as square classes,
would also add 4 per group between q = 3 and 5. The q = 7 count decides between "two families of q orbits per group"
and such alternatives.

## 3. Engine fixes (from the failed q = 7 runs)

- **`expand` slice size.** The completed-fold expansion by |N| sized its slices by the int8 images. `pack` then makes
  an int64 copy, 8× larger, and asked for 9.2 GB. Slices are now a tenth of that size.
- **Lazy work stack.** The stack now holds a generator per parent batch instead of all its child chunks at once.
  This bounds pending memory to one chunk per column. It was not the cause of the third failure.
- **Remaining limit.** The fibres and `finish()` store every fold: about 1e7 folds at q = 7 (estimated). That needs
  more than this container's 7 GB, so the run is run_batch.py job 16, with the prediction committed.

Engine tests pass (`tests/test_engines.py`, 9 tests).

## 4. A caveat that applies to every census conclusion

A one-parameter family defined over 𝔽_{p²} but not over 𝔽_p comes with its Frobenius conjugate. Frobenius swaps the
two curves, so only their finitely many intersection points are 𝔽_p-rational. Prime-field censuses therefore cannot
see such a pair. Every "no other families" and "bounded" claim in this project means: over prime fields, at the primes
computed. The census engine works over prime fields only; testing 𝔽_{p²} would need extension-field arithmetic in
`fibre_torch.py`.

## 5. Next

1. Job 16 (q = 7) on the GPU machine. Then, if prediction 6 holds, merge with `analysis/v2_42_merge_dump.py` and
   repeat the cross-prime analysis at three primes.
2. Ansätze for the large-image families:
   - torus composites, i.e. individual arrow scalings, as Theorem 3's F1 needed in v2.33;
   - unipotent class substitutions outside G;
   - fitting a polynomial curve through the growing orbits of the two groups directly.
3. Lean closed forms for the 18 found families, as for 1⇉2⇉3 (v2.36–v2.38), once the count is settled.

Files:
- `analysis/v2_42_{merge_dump,crossprime,growth,tangent_search,be,shear_search,vblock}.py`;
- `fibre_torch.py` (two fixes), `run_batch.py` (job 16);
- `logs/v2_42/`, `logs/run_reps_two_pairs_short_7*.log`.
