# v2.55 predictions (committed before any stratum or orbit computation on the control)

Control: 1⇉2⇉3 + 1→3 (two_pairs_short). It is obstructed (consecutive classes), and its local-torus counts grow:
1668, 1840, 2028 at q = 3, 5, 7. Observation, not a prediction: its vertex blocks are 21 at each of q = 3, 5, 7
(logs/v2_55/block_census_two_pairs_short.txt). In 13 of the 21 blocks the local-torus count grows with q; in 8 it is
constant.

55. The v2.53–v2.54 pipeline (QUIVER=two_pairs_short) fails on the control: in at least one stratum, some component has
    no open H_Z-orbit, or some r has dim Z_r > r (reported as rank ≤ r at every sampled point).
56. Failures occur only in strata whose local-torus count grows with q. Every component of the 8 constant-count
    strata passes.
