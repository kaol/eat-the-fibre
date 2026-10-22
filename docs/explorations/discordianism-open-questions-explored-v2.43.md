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

## 2a. Update (same day): the representatives arrived; the merged count is linear

With `logs/reps_two_pairs_short_7.pkl` from the rerun:

| | q = 3 | q = 5 | q = 7 | |
|---|---|---|---|---|
| merged with the shear U (full Δ-compatible group) | 939 | 983 | **1027** | = 873 + 22q; **v2.42 prediction 6 held** |
| the 18 inner-curve families of v2.42 | 3 orbits each | 5 each | **7 each** | pairwise disjoint at every prime |
| the 17 bounded curves | 3 each | 3 each | **3 each** | |
| remaining after the 18 | 885 | 893 | **901** | = 873 + 4q |
| each of the two dim-7 image groups | 23 | 27 | **31** | = 17 + 2q |

So the non-linear part of the unmerged count is exactly what the shear absorbs. Under the full Δ-compatible group,
1⇉2⇉3 + 1→3 has 22 one-parameter families at three primes:
- 18 are inner curves Inn(1 + t·x)·P₀ with integer closed forms;
- 4 lie in the two mirror-image groups with dim im P = 7, two per group. They are not inner curves or shear curves
  (v2.42), and their form is still open.

Prediction 7 (the q = 7 checks above) held, but it was committed after the run: a gitignore error stopped the commit
and the chained command ran on (logs/v2_42/PREDICTIONS.md, process note).

## 3. Tooling fix

`run_batch.py` accepted `--save-reps` and silently dropped it, which I should have caught when I recommended the flag
for jobs 12–16. Now:
- it passes `--save-reps logs/reps_QUIVER_p.pkl` to pcensus.py;
- with the flag, a job counts as done only once that file exists, so a finished job is redone from its checkpoint
  without recomputing fibres.

This was tested with a stub pcensus.py: the flag is forwarded for job 16, and the job is skipped without it.

## 4. Next

1. Find the 4 large-image families. Candidates:
   - torus composites with rational coefficients, as Theorem 3's F1 needed (v2.33);
   - class substitutions outside G;
   - or families that are not single automorphism orbits at all, i.e. moduli of the splitting itself. Here
     Aut(kQ) is larger than Inn ⋊ Sub (shears), so v2.33's behaviour census would need extending to decide this.
2. Lean closed forms for the 18 inner-curve families, with rigidity and disjointness, as for 1⇉2⇉3 (v2.36–v2.38).
3. Job 14 (1⇉2⇇3 + 1→3 + 2→4 over 𝔽₅).
