# v2.57 predictions (committed before the full run)

ChatGPT's Phase 2, Program 1: the v2.54 criterion on Q₁ = 1⇉2⇇3 + 1→3 + 2→4 (job 14's quiver), all 90 vertex blocks.
The pipeline now uses an indexed presolve (fold_ideal.presolve_fast; 118 s instead of 856 s on the 236-unknown
stratum 1000,1000,1000,1000). Test block 1000,1000,1000,1000 (old presolve): 5 components, all passing.

Outcome classes, reported separately: certified pass; geometric failure (rank ≤ r at every sampled point, then an
explicit invariant from v2_55_invariant.py); computationally unresolved (timeout, no rational point, Singular error).

59. Every component of every one of the 90 strata is a certified pass. So the connected census groups have finitely
    many orbits on Folds(k̄Q₁) in characteristic 0.
60. Regression: with the new presolve, the control's growing block 000,001,001 still has exactly one geometric failure,
    and pairs_in_short block 100,010,001 still has 16 components, all passing.
