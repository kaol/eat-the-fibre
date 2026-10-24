# v2.51 predictions (committed before looking at any profile output)

Context: job 14 (1⇉2⇇3 + 1→3 + 2→4 over 𝔽₅) ran out of memory inside one fibre. analysis/v2_51_profile.py records,
per source fibre, the engine's frontier width after each column (after the N-quotient).

37. At q = 3, the frontier is concentrated in few fibres: fewer than 1% of the 12210 source fibres ever reach a width
    of 1000 rows.
38. In the heaviest q = 3 fibres the peak width occurs at a middle column (strictly between kdim + 1 = 7 and d = 18),
    and the peak is at least 10× the fibre's final leaf count (N-orbit representatives).
39. At q = 5, at least one source fibre exceeds 2·10⁶ rows at some column (the cap), consistent with the 44 M children
    recorded at one column in fibre_torch's docstring.

## Added after the q = 3 profile and the column-order experiments on three fibres (before the runs below)

The order experiments (fibres 5973, 14916, 5965 at q = 3) suggest filling the paths through the joining arrow first
(c·b₁, c·b₂, right after a₁, b₁, c and before the tail arrow t = 2→4 and the paths a_i·t, b_j·t):
ETF_CORDER = 4,6,8,14,15,9,10,12,11,13 (basis indices).

40. With that order, the full q = 3 census of pairs_in_short_tail in the container reproduces job 13 exactly:
    N = 25622, orbits_G = 1443, orbits_GxGr = 1281, orbits_GxT = 1281, orbits_GxTloc = 1220.
41. With that order, the q = 5 profile of every source fibre stays below the 2·10⁶-row cap.

## Added after the q = 5 profile with the first order (stopped at 13,215 of 17,434 fibres), before the runs below

Fibre 16019 at q = 5 still exploded with ETF_CORDER = 4,6,8,14,15,9,10,12,11,13 (503,274 rows when the time cap hit) and
with 4,6,8,14,15,9,10,11,12,13. The order 4,6,8,14,15,9,12,16,13,17,10,11 (each path followed at once by its extensions:
a₁, b₁, c, cb₁, cb₂, t, b₁t, cb₁t, b₂t, cb₂t, a₁t, a₂t) finished it in 19 s with peak 627.

42. With the order 4,6,8,14,15,9,12,16,13,17,10,11: the full q = 3 census reproduces job 13 exactly, and every q = 5
    source fibre finishes within 120 s with peak frontier below 10⁶ rows.
43. Job 14 itself (q = 5) then completes in the container, and orbits_GxTloc = 1220, the q = 3 value (the run_batch
    note's prediction: the local-torus count is constant in q, as Conjecture 3ᶜ requires for this unobstructed quiver).

## Outcome of 42 (recorded before 44), and the next prediction

42 failed at q = 3: with the path-tree order alone, four q = 3 fibres (17663, 14964, 15009, 17727) exploded at columns
13–16, fibres that the first order handles easily. No single order suits every fibre. fibre_torch.portfolio now tries
the orders per fibre (breadth-first with a row cap, first success kept, converted to the main column basis). Check at
q = 3 with ETF_CORDERS = path-tree;first and cap 5·10⁴: 12206 fibres by the first order and 4 by the second. The census is
identical to job 13 in all 21 fields, and all 17748 fibres are identical as sets of folds to the single-order run.

44. Job 14 with ETF_CORDERS = 4,6,8,14,15,9,10,12,11,13;4,6,8,14,15,9,12,16,13,17,10,11 and cap 3·10⁶: every q = 5 source
    fibre is completed by one of the two orders (none left for run()), and the census completes in the container.
    (Prediction 43, orbits_GxTloc = 1220, stands for this run.)
