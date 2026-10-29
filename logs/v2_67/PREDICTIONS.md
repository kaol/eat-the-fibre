# v2.67 predictions (committed before the runs)

v2.66 left 1⇛2⇇3 + 1→3 open on one component: block 100,001,001, component 1 (dimension 10) of the default run; the same
component reappears in leaves 114 and 116 of the case split. The v2.54 cut (one minor at a time, then a primary
decomposition) ran out of memory: 6 GB here, 51 GB in 16 hours on the user's machine.

New method (analysis/v2_67_orbits.py): orbit stratification by invariant boundaries. For an H-stable prime X
(H the connected group with Lie algebra 𝔥 = Der_Δ + class weightings of the component), with generic orbit dimension
e = dim X, the complement of the dense orbit is the H-core of V(f) ∩ X for one e × e minor f, and its ideal is the
closure of (I_X, f) under the vector fields of 𝔥. Its minimal primes are H-stable; recurse. Finitely many orbits on X
iff every prime met has e = dim X. Every step is checked: tangency of the vector fields to I_X, closure, dims.

91. On control components that passed the v2.54 criterion, the new method also reports finitely many orbits.
92. On component 1 of block 100,001,001 of 1⇛2⇇3 + 1→3 it terminates within the container's memory and reports finitely
    many orbits. This would complete the char 0 theorem for 1⇛2⇇3 + 1→3 (the other 5 components of that block passed
    in v2.66).

Added after the first runs (prediction 92 held: 326 primes, every one with a dense orbit, 264 s; the obstructed control
two_pairs_short 000,001,001 component 2 is reported as having no dense orbit, the passing components of that block
as finite), and after one block of 1⇉2→3 (100,010,001) gave 7 geometric H-orbits against 7 census classes:

93. Census = geometric orbit count. For every vertex block of K₂, K₃, 1⇉2→3, 1⇉2 + 1→3, 1⇛2→3, 1⇛2 + 1→3 and
    1⇉2→3→4 (no case-split blocks, saved q = 3 representatives), the number of geometric H-orbits on the stratum
    (distinct H-stable primes met, each with a dense orbit at a rational point, plus the points of 0-dimensional
    primes) equals the number of local-torus classes in that block at q = 3. Summed: 20, 34, 71, 101, 120, 183, 334.

Added after prediction 93 held on all seven quivers (20, 34, 71, 101, 120, 183, 334: equal block by block), before
running these:

94. The same block-by-block equality for 1→2⇉3 (71), 1⇉2 + 3→2 (101), 1→2⇉3→4 (257), 1⇉2 + 2→3 + 2→4 (414),
    1⇉2 + 1→3 + 4→2 (485), and Q₀ = 1⇉2⇇3 + 1→3 (366, v2.58 certificates, v2.53 block census). For Q₁, Q₂ and
    1⇉2→3→4→5 on the blocks without a case split.

Results so far (before this addition): 93 held on all seven quivers; 94 held on all six quivers named, after two
script fixes (a Singular `size` vs `ncols` bug that dropped defect equations, and resampling when a random point is
special); the triple quiver's q = 3 census (torch engine, 46 s) gives 72 classes in block 100,001,001, equal to the 72
geometric orbits found there before the census was computed. Before running the other 20 blocks:

95. 1⇛2⇇3 + 1→3: the geometric orbit count equals the census class count in every one of its 21 blocks; total 752.
