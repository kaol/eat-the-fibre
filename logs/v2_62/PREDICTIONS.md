# v2.62 predictions (committed before the sampled run; 3 folds tried first)

Program 3, items 2 (one-sided kinds) for Q₀ → Q₁, at integral folds P̄ of Q₀ (census representatives that are folds
over ℚ), exactly over ℚ: the whole fibre over P̄ of the J-preserving folds (P = [[P̄, 0], [X, P_J]]) and of the
kQ-preserving folds (P = [[P̄, Y], [0, P_J]]), for λ = 0 and λ = 1 separately (analysis/v2_62_fibre.py).
First 3 folds: every fibre a point or an affine space (one of dimension 3, kQ-preserving with λ = 1), except one
J-preserving λ = 0 fibre with two points (the trivial lift and one J-only lift).

72. On a random sample of 80 integral folds of Q₀ (seed 0): every fibre is empty, a finite set of points, or an
    affine space; no fibre has dimension above 4.
73. λ = 0 fibres are finite (points) in both kinds: positive-dimensional one-sided families occur only with λ = 1.
74. (Added before the full run; one block tried: 000,001,001 gives 11 = its census count.) Explicit characteristic-0
    orbit representatives of Q₀ (analysis/v2_62_orbits.py: one rational point of maximal rank on each orbit closure
    found among the certified cut varieties) number exactly the census local-torus count in every one of the 21 blocks,
    366 in total.
75. (Added before the full run; 4 classes tried.) Relative criterion over the 366 class representatives
    (analysis/v2_62_relative.py): every positive-dimensional component of every one-sided fibre (kinds J and Q, λ = 0, 1)
    passes, so Stab(P̄)° has finitely many orbits on each; no fibre component has dimension above 3.
