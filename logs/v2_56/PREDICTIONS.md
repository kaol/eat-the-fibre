# v2.56 predictions (committed before the q = 5 validation finished and before any q = 7 census)

Engineering round: the order portfolio runs batched (many fibres per breadth-first pass, per-fibre row caps), and
finish()'s local-torus test is linear in w (one product instead of a W·d⁵ einsum per representative).

57. The batched portfolio plus the new finish() reproduce job 14 (q = 5) in all fields except t. Already checked at
    q = 3: equal to job 13 in all 21 fields, and all 17748 fibres equal as sets to the single-order run.
58. Job 17 (1⇉2⇇3 + 1→3 + 2→4 over 𝔽₇, to run on the GPU machine) gives orbits_GxTloc = 1220.

## Outcome of 57
57 held. The q = 5 rerun (batched portfolio, linear local-torus test) equals job 14 in every field except t. On the
2-core CPU the portfolio took 3010 s, against about 2400 s one fibre at a time: batching does not pay on 2 CPU cores,
and the 5 fibres over the first-pass cap repeat their work. finish() took about 110 s, against about 900 s.

## Outcome of 58 (GPU run on the Strix machine, 2026-10-25)
58 held: job 17 gives orbits_GxTloc = 1220 (N = 3,390,442 folds, 3215 G-orbits, t = 18,472 s).
- The portfolio took 17,936 s. The first order handled 23,597 fibres in 9,288 s. The second order handled the 5
  deferred fibres (932, 34319, 34361, 38367, 81564) in 8,648 s. No fibre was left for the depth-first engine.
- Block census (logs/v2_56/block_census_tail.txt): 90 vertex blocks at each of q = 3, 5, 7, and the local-torus count
  of every block is the same at all three primes.
