# Discordianism: Open Questions Explored (v2.43): GPU results: joined ends stay at 366; 1⇉2⇉3 + 1→3 is not linear

*Four census jobs ran on the Strix Point GPU machine (logs/results_strix.jsonl). Two pre-registered constancy
predictions for unobstructed quivers held, including 366 orbits over 𝔽₁₁. The q = 7 census of 1⇉2⇉3 + 1→3 broke
v2.42's linear prediction. It also exposed a bug: run_batch.py had never passed `--save-reps` on, so the
shear-merged count is still pending.*

---

## 0. Results

| quiver | q | local-torus orbits | prediction | |
|---|---|---|---|---|
| 1⇉2⇇3 + 1→3 (joined converging ends) | 3, 5, 7, **11** | 366, 366, 366, **366** | constant (Conjecture 3ᶜ) | held |
| 1⇛2⇇3 + 1→3 | 3, **5** | 752, **752** | constant (v2.28 P2) | held |
| 1⇉2⇇3 + 1→3 + 2→4 | **3** | **1220** | — (𝔽₅ still to run) | — |
| 1⇉2⇉3 + 1→3 (consecutive classes, ends joined) | 3, 5, **7** | 1668, 1840, **2028** | 2012 (linear, v2.42 P6) | **failed** |

## 1. Joined converging ends: four primes

1⇉2⇇3 + 1→3 now has the same count, 366, at q = 3, 5, 7 and 11. Its unjoined version 1⇉2⇇3 grows (Theorem 3), and
the joining arrow removes every family-carrying inner curve (v2.36). Two caveats apply:
- constancy over prime fields at four primes is evidence, not proof;
- v2.42 §4 applies: prime-field censuses cannot see a family that exists only over 𝔽_{p²}.

1⇛2⇇3 + 1→3 is also constant at 752 over 𝔽₃ and 𝔽₅.

## 2. 1⇉2⇉3 + 1→3 is not linear in q

The orbit counts at q = 3, 5, 7 are:

| count | q = 3 | q = 5 | q = 7 | interpolating quadratic |
|---|---|---|---|---|
| under the census group and class weightings | 1668 | 1840 | 2028 | 2q² + 70q + 1440 |
| G-orbits | 2824 | 8148 | 19808 | 792q² − 3674q + 6718 |
| folds N | 77 576 | 2 656 570 | 29 947 356 | — |

Three points always fit a quadratic, so the quadratic column records the interpolant, not a law. What is established
is that the first count is not linear. This revisits v2.42's picture of 22 one-parameter families.
- **A possible source of a q² term.** The shear U = {c ↦ c + s·δ_aδ_b} is in the full Δ-compatible group but not in
  the census group. If U acts with orbits of size q along a one-parameter family, the census group sees q·q orbits
  where the full group sees q. Whether that is the explanation, and whether the merged count is linear (v2.42
  predicted 873 + 22·7 = 1027), needs the representatives.
- **The four missing families of v2.42**, in the two image-dimension-7 groups, are the natural suspects for the
  non-linear part. This is not checked.

## 3. Tooling fix

`run_batch.py` accepted `--save-reps` and silently dropped it, which I should have caught when I recommended the flag
for jobs 12–16. Now:
- it passes `--save-reps logs/reps_QUIVER_p.pkl` to pcensus.py;
- with the flag, a job counts as done only once that file exists, so a finished job is redone from its checkpoint
  without recomputing fibres.

This was tested with a stub pcensus.py: the flag is forwarded for job 16, and the job is skipped without it.

## 4. Next

1. `python3 run_batch.py 16 --save-reps` (with `HSA_OVERRIDE_GFX_VERSION=11.0.0` on the Strix Point). This resumes
   from `checkpoints/two_pairs_short_7.pkl`, so only the orbit stage reruns. Send back
   `logs/reps_two_pairs_short_7.pkl`.
2. Then: the shear merge (`analysis/v2_42_merge_dump.py`), the merged count against 1027, and the cross-prime family
   analysis at three primes.
3. Job 14 (1⇉2⇇3 + 1→3 + 2→4 over 𝔽₅), to compare with job 13's 1220.
