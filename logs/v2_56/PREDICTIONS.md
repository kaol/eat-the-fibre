# v2.56 predictions (committed before the q = 5 validation finished and before any q = 7 census)

Engineering round: the order portfolio runs batched (many fibres per breadth-first pass, per-fibre row caps), and
finish()'s local-torus test is linear in w (one product instead of a W·d⁵ einsum per representative).

57. The batched portfolio plus the new finish() reproduce job 14 (q = 5) in all fields except t. Already checked at
    q = 3: equal to job 13 in all 21 fields, and all 17748 fibres equal as sets to the single-order run.
58. Job 17 (1⇉2⇇3 + 1→3 + 2→4 over 𝔽₇, to run on the GPU machine) gives orbits_GxTloc = 1220.
