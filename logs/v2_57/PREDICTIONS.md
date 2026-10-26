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

## Outcome (recorded after the runs)

- **59 held.** All 90 blocks are certified passes (`logs/v2_57/tally_q1.txt`, from `analysis/v2_57_tally.py`):
  89 blocks with 372 components (353 pass every r, 19 points), and block 0000,0000,0000,0001 by the case-split presolve
  (54 leaf primes covering the stratum, 53 pass, 1 point). No geometric failure, nothing unresolved.
- **60 held.** `logs/v2_57/regression.txt`: 000,001,001 fails at exactly one component (dim 4), 100,010,001 has 16
  components, all passing. The case-split mode was regressed too (`regression_branch.txt`): the control's 000,001,001
  and 001,001,001 still fail exactly once each (dim 4 and dim 6), and pairs_in_short 100,010,001 passes on all leaves.
