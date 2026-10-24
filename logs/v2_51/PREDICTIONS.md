# v2.51 predictions (committed before looking at any profile output)

Context: job 14 (1⇉2⇇3 + 1→3 + 2→4 over 𝔽₅) ran out of memory inside one fibre. analysis/v2_51_profile.py records,
per source fibre, the engine's frontier width after each column (after the N-quotient).

37. At q = 3, the frontier is concentrated in few fibres: fewer than 1% of the 12210 source fibres ever reach a width
    of 1000 rows.
38. In the heaviest q = 3 fibres the peak width occurs at a middle column (strictly between kdim + 1 = 7 and d = 18),
    and the peak is at least 10× the fibre's final leaf count (N-orbit representatives).
39. At q = 5, at least one source fibre exceeds 2·10⁶ rows at some column (the cap), consistent with the 44 M children
    recorded at one column in fibre_torch's docstring.
