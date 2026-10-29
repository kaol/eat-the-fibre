# v2.69 predictions (committed before the runs)

Goal: explain the q-independence of the census counts (Conjecture G, v2.67) by Lang's theorem. Survey of the v2.67
trees (7920 nodes over 17 quivers): every 0-dimensional node is a single rational point, and every positive-
dimensional node has its dense orbit through a rational point (6663 found at rational points directly, 1 by
absPrimdecGTZ). The remaining ingredient is the stabiliser. For a node Y with group H_Y = (column-sum-1 substitution
groups of the classes) × (torus of the weightings W(Y)), take a rational point x of the dense orbit and compute the
stabiliser scheme Stab_{H_Y}(x) (equations σ P_x = P_x σ in the group coordinates, determinants inverted) and its
minimal primes over ℚ (analysis/v2_69_stab.py). The identity component is defined over ℚ and contains a smooth
rational point, so the stabiliser is connected iff there is exactly one minimal prime.

100. Every positive-dimensional node of every tree has a connected stabiliser, of dimension dim H_Y − dim Y.
     (A disconnected stabiliser with component group A would split its orbit over 𝔽_q into |H¹(𝔽_q, A)| > 1
     classes, contradicting census = number of nodes; so the census agreement predicts connectedness.)
101. Local weightings: at every node, the weightings compatible on Y (defect condition alone) are all tangent to Y,
     i.e. W_defect(Y) = W(Y). (Otherwise census moves would join orbits of different nodes.)
