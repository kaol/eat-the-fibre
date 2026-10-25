# Discordianism: Open Questions Explored (v2.55): the finiteness pipeline on an obstructed control quiver

*v2.54 proved finitely many census orbits on 1⇉2⇇3 + 1→3 by a computational criterion. A criterion is only worth
something if it can fail where it should. This round runs it on the obstructed control 1⇉2⇉3 + 1→3, whose
local-torus count grows with q. In each of the 6 finished blocks whose count grows, it fails at exactly one component,
and in each of those the failure comes with an explicit invariant. All 8 blocks whose count is constant pass. Along the
way a silent Singular overflow was found and worked around; re-running v2.54 with the fix gives identical results.*

---

## 0. Results

The predictions are in logs/v2_55/PREDICTIONS.md.

| # | Statement | Status |
|---|---|---|
| obs. | The control's vertex blocks: 21 at each of q = 3, 5, 7. In 13 of them the local-torus count grows, in 8 it is constant. | census reps (logs/v2_55/block_census_two_pairs_short.txt) |
| 55 | The pipeline fails somewhere on the control | **held**: one failing component in each of the 6 finished growing blocks |
| 56 | Failures occur only in growing blocks; all components of constant blocks pass | **held**: 8 of 8 constant blocks, 59 components, all pass |
| infra | v2.53–v2.54 do not depend on the Singular overflow | **confirmed**: the rerun of all 21 pairs_in_short blocks gives the same 118 components, 111 passes and 7 points |

The 7 largest growing blocks (strata of dimension 6–8, with up to ~60 coordinates) were still in their first
decomposition when the run was stopped. They are not needed for 55 or 56.

## 1. The failures, with invariants

Each failing component Z has no open orbit of the connected census group H_Z: the orbit-dimension criterion fails at
r = dim Z − 1. `analysis/v2_55_invariant.py` makes the failure explicit:
- it parametrises Z, which is a graph over an independent set of coordinates;
- it computes the fields of 𝔥_Z on Z;
- it finds the semi-invariants (invariants of [𝔥_Z, 𝔥_Z] that are eigenvectors of 𝔥_Z);
- an integer relation between their characters gives an invariant rational function κ.

A nonconstant κ proves that the generic orbits have codimension ≥ 1, so the component carries a modulus.

| block | component | dim | invariant κ (in the presolved coordinates of the stratum) |
|---|---|---|---|
| 000,001,001 | 2 | 4 | k₃₄ / (k₁₄² k₅₂) |
| 100,100,000 | 2 | 4 | k₁₃² k₅₀ / k₄₃ (the mirror) |
| 001,001,001 | 5 | 6 | k₃₅ / k₁₄ |
| 100,100,100 | 5 | 6 | k₄₅ / k₁₃ (the mirror) |
| 000,010,010 | 3 | 4 | k₅₄ / k₃₁ |
| 010,010,000 | 2 | 4 | k₅₃ / k₄₁ (the mirror) |

Here k_ij are entries of P_K (P restricted to K = ker Δ) in the integer K basis. The blocks come in mirror pairs,
under the source–sink symmetry of the quiver.

**The growth matches, roughly.** Block 000,001,001 has 25, 27, 29 local-torus orbits at q = 3, 5, 7: slope 1 in q,
which is one family of orbits parametrised by κ. Block 001,001,001 has 37, 41, 45: slope 2 with a single failing
component. Either κ's fibres split into two census orbits, or a second family sits in that component's boundary. The
pipeline stops at the first failure in a component, so it cannot tell which.

## 2. What this says about the method

The criterion has both properties:
- On the unobstructed 1⇉2⇇3 + 1→3 it passes everywhere (v2.54), giving finitely many orbits over k̄.
- On the obstructed control it fails exactly at components that carry an explicit continuous invariant. These lie
  only in vertex blocks whose census count grows.

So the pipeline is a working decision procedure (in characteristic 0, per block) for "finitely many census orbits",
with certificates either way:
- a dimension bound for every r when it passes;
- an invariant rational function when it fails.

That is the candidate reduction ChatGPT's Program 2 asked for. Proving Conjecture 3ᶜ in general would need it to pass
on every unobstructed quiver, which is not a finite computation. Each quiver checked is a theorem for that quiver.

## 3. Infrastructure

- **`analysis/fold_ideal.py`** builds the stratum ideals with sparse polynomials (sympy PolyRing), with a faster
  presolve. Building the equations for one control stratum takes 7 s. The expression-based version had not finished building and presolving in 10 minutes.
- **`analysis/singular_nohilb.py`.** On the control, Singular 4.3.2's `minAssGTZ` overflowed ("overflow at t^24" in
  `hilb`) and silently returned no components. The first attempt reported "0 components" for a nonempty stratum, which
  is how the bug showed up.
  - The fix is a copy of `primdec.lib` in which Hilbert-driven `std(I, hilb, w)` becomes `std(I)` and `groebner(`
    becomes `std(`. It is loaded through `SINGULARPATH`. The ideals are the same; only the Gröbner strategy changes.
  - The pipeline now also treats any `? error` or `overflow` in Singular's output as a failure, instead of trusting an
    empty list.
- **`analysis/v2_55_pipeline.py`** is the general pipeline for any quiver; `v2_55_invariant.py` makes the invariants.

## 4. Next

1. Run the pipeline on job 14's quiver (1⇉2⇇3 + 1→3 + 2→4, d = 18, 1220 orbits) and on 1⇛2⇇3 + 1→3. Each pass is a
   finiteness theorem for that quiver. d = 18 will need the same speedups taken further, or Singular over a finite
   field first as a filter.
2. Finish the 7 largest control blocks, and count the families per growing block against the census slopes.

## Files

- analysis/fold_ideal.py, singular_nohilb.py, v2_55_pipeline.py, v2_55_invariant.py
- logs/v2_55/: PREDICTIONS.md, block_census_two_pairs_short.txt, pipeline_two_pairs_short.txt (14 of 21 blocks),
  pipeline_pairs_in_short.txt (the full rerun), validation_pairs_in_short.txt, invariants_two_pairs_short.txt
