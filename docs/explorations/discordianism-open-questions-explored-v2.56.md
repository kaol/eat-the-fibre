# Discordianism: Open Questions Explored (v2.56): the census engine, ready for job 14's quiver over 𝔽₇

*Engineering round, at the user's request: make the GPU program fit for another attempt on 1⇉2⇇3 + 1→3 + 2→4.
Job 14 (q = 5) was finished in the container in v2.51, so the new run is the next prime, q = 7: job 17 in
run_batch.py. Two changes, both checked to reproduce jobs 13 and 14 exactly. Prediction 57 held; 58 waits for the
GPU run.*

---

## 1. Changes

**The order portfolio is batched** (`fibre_torch.portfolio`, `bfs_group`). v2.51 tried the column orders one fibre at
a time. Now each order runs breadth-first over batches of fibres:
- 256 fibres by default (`ETF_PGROUP`), with each frontier row tagged by its fibre.
- A fibre whose own frontier exceeds the cap is dropped from the batch and moves on to the next order. The cap is
  `ETF_PCAP0` for every order but the last and `ETF_PCAP` for the last.
- If the whole batch exceeds `ETF_PGCAP` (default twice the cap), the fibre with the largest frontier is dropped
  until it fits.

So memory stays bounded per batch, and the GPU gets large launches instead of one small fibre at a time. Fibres that
no order manages fall back to the depth-first engine.

**finish()'s local-torus test is linear in the weighting.**
- *Before:* for each G-orbit representative, P was tested against all (q−1)^classes weightings w with two einsums
  over (W, d, d, d). The second einsum ran naively, at cost W·d⁵ per representative.
- *Now:* the defect Σ_m P[m,l] Δ_w[m] − P Δ_w[l] Pᵀ is linear in w. The pipeline computes the per-class defects E_c
  once, then tests every w with one matrix product (W × classes)·(classes × d³).
- *Effect at q = 5:* 900 s → 110 s. At q = 7 there are 1296 weightings instead of 256, so the old version would have
  taken hours.

## 2. Checks

| run | result |
|---|---|
| q = 3, orders reversed, first-pass cap 5000 (so both orders are used: 12205 + 5 fibres) | equal to job 13 in all 21 fields; all 17748 fibres equal as sets to the single-order run |
| q = 5, the job 14 settings | equal to job 14 in every field except t. Portfolio 3010 s on the 2-core container (not faster there than one at a time); finish() 110 s |
| q = 7 setup | kernel stage: 61.9 M subspaces, 533,832 kernel folds, 92,572 torus orbits, 23,602 source fibres, \|N\| = 294; 83 s on the container CPU |
| q = 7 sample | 512 random source fibres, all done by the first order: 219 s on the container CPU (about 2.8 h for all sources at that rate, not counting the heavy fibres) |

## 3. Running job 17

```
export HSA_OVERRIDE_GFX_VERSION=11.0.0
git pull
python3 run_batch.py 17 --save-reps
```

Job 17 carries:
- the two column orders, `ETF_PCAP0 = 10⁶` (so the fibres whose frontier spikes only at the c column stay in the
  first order);
- `ETF_PCAP = 2·10⁷` for the last order (about 6.5 GB of int8 frontier at most per batch);
- `ETF_MEM = 512` MiB, to keep the depth-first fallback well inside 64 GB if it is ever needed.

It writes results.jsonl, batch.log and logs/reps_pairs_in_short_tail_7.pkl; please send all three back. The reps file
lets the v2.55 pipeline and the block census read the q = 7 orbits.

**Prediction 58:** orbits_GxTloc = 1220, as at q = 3 and q = 5.

## Files

- fibre_torch.py: `bfs_group` and the batched `portfolio`; env `ETF_PCAP0`, `ETF_PCAP`, `ETF_PGROUP`, `ETF_PGCAP`.
- pcensus.py: the linear local-torus test in `finish()`.
- run_batch.py: job 17.
- logs/v2_56/: PREDICTIONS.md, results and run logs of the q = 3 and q = 5 checks.
