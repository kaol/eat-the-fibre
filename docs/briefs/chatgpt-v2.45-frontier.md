# eat-the-fibre: proof strategies for the v2.45 frontier

Work on the current repository branch and verify its actual state before proceeding. Start with:

- `docs/explorations/discordianism-open-questions-explored-v2.45.md`
- `analysis/v2_45_fixed_image.py`
- `analysis/v2_45_singular_VA.py`
- `analysis/v2_45_generic_m.py`
- the v2.44 line construction and corresponding Lean files

Do not treat finite-field rank observations as characteristic-zero proofs. Do not equate tangent-space dimension with actual orbit-intersection dimension.

## Task 1: Resolve the fixed-image distribution

Identify the four kernel-tilting families found in v2.44. Compute their image subspaces exactly over Q and determine which are equal and which are equivalent under the relevant automorphism group.

For families with a common image, construct the pairwise algebraic orbit-incidence relation. Determine whether each pair represents the same generic behaviour curve or only has finitely many cross-identifications.

Use the exact image comparison as the first experiment. Report the image bases, exact rank computations, and any explicit automorphisms or separating invariants.

## Task 2: Prove the generic tangent statement exactly

For the seven-dimensional fixed-image variety V_A associated to representative 5107, replace the modular rank evidence by an exact generic-rank calculation over the function field of its irreducible component.

Target:

- dim T_P(V_A) = 7 at the generic point.
- dim(T_P(V_A) intersect [Der_A, P]) = 6 at the generic point.
- The resulting first-order value m = 1.

Give explicit rank certificates using polynomial minors, ideal membership, or a reproducible exact computation. Investigate the rank-drop locus and use representative 5088 as a special-point check.

If these generic rank claims fail, identify the exact failed assertion and explain the corrected dimension picture.

## Task 3: Attempt the dominant-slice theorem

Try to construct an explicit six-parameter family of automorphisms g(s) preserving A such that

Phi(s,t) = g(s) P_t g(s)^(-1)

is a fold with image A for every parameter in a nonempty open set.

Here P_t = P_0 + tN is the already verified kernel-tilting line.

Target the following lemma:

There exists a point at which Phi has differential rank 7. Such a result would make Phi dominant onto V_A. Combined with the existing proof that the line intersects each algebra-automorphism orbit in finitely many points, this would establish that the generic behaviour space has dimension exactly one.

The critical obligation is to prove that the conjugated operators remain folds. Six independent tangent vectors do not suffice.

If integration fails, compute the first nonzero higher-order obstruction for the chosen directions. Do not discard the one-dimensional moduli conjecture merely because one proposed slice fails.

## Required output

For each task, report:

1. The exact theorem or intermediate lemma targeted.
2. Existing repository results used, with actual file and theorem names.
3. A complete proof or a clearly identified unresolved step.
4. An exact computational test that could falsify the claim.
5. The smallest next deliverable.

Clearly distinguish proved mathematical statements, exact computer-algebra results, finite-field evidence, and conjectures. Prefer one substantive theorem or counterexample over a large collection of speculative approaches.