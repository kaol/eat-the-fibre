# v2.53 predictions (committed before any computation)

Context: v2.52 proved that on Q₀ = 1⇉2⇇3 + 1→3 (char 0) the vertex block takes exactly 21 values. This round stratifies
the folds by vertex block. G-orbit representatives with local-torus labels were saved for q = 3, 5, 7
(logs/reps_pairs_in_short_{3,5,7}.pkl).

49. For each of the 21 vertex blocks, the number of local-torus orbits with that block is the same at q = 3, 5, 7. So the
    366 splits the same way at every prime.
50. Over ℚ, the stratum S_V (folds with vertex block V; unknowns P_K with V fixed, plus P on a complement of K) is a
    variety whose dimension, for every V, equals the largest degree in q of the census counts |S_V(𝔽_q)|. The degree is
    estimated from q = 3, 5, 7 as round(log(N₇/N₅)/log(7/5)).

## Outcomes of 49–50, recorded before 51

49 held: every block has the same local-torus count at q = 3, 5, 7 (total 366).
50 failed for 1 of 21 blocks. Stratum 100,010,001 has dimension 5 over ℚ, while the crude degree estimate from q = 5, 7
gives 4 (log ratio 4.48). Its unique 5-dimensional component is absolutely irreducible (absPrimdecGTZ, extension
degree 1). Its counts 748, 7036, 31732 fit q⁵ + about 6q⁴, so the estimator, not the stratum, was off. The other
20 blocks agree.

51. Every component of maximal dimension of every one of the 21 strata is absolutely irreducible (absPrimdecGTZ
    extension degree 1, one absolute prime per ℚ-prime).
