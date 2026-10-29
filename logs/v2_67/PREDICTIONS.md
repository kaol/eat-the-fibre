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
