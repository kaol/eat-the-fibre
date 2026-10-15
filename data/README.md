# Raw data: the v2.23 batch

| file | content |
|---|---|
| `checkpoints.tar.gz` | the 10 checkpoints of the v2.23 batch (`run_batch.py` jobs 1–10), as produced on the user's machine (torch engine, ROCm on an AMD Radeon 880M/890M iGPU, 2026-10-14/15). Extracts to `checkpoints/*.pkl`, 44 MB |
| `results.jsonl` | the JSON lines `pcensus.py` wrote for those runs (one per job, v2.24 fields included) |

A checkpoint is the raw census: `{name, p, nreps, fibres: {torus-orbit index: [int16 (d, d) arrays]}}`. Each fibre over a
representative kernel fold is stored in the new basis. Every fold and every statistic is derived from these by
`pcensus.finish` (materialization via the grading, orbit labelling, torus merging). Nothing in `results.jsonl` needs to
be trusted on its own.

## Use

```bash
tar xzf data/checkpoints.tar.gz -C data/                                   # -> data/checkpoints/*.pkl (git-ignored)
python3 refinish.py --dir data/checkpoints --out /tmp/refinished.jsonl     # recompute every statistic, no fibre work
python3 pcensus.py two_pairs 5 --checkpoint COPY.pkl --no-orbits --save-folds logs/folds_two_pairs_5.pkl   # all folds
```

Pass a *copy* of a checkpoint to `pcensus.py`, because it rewrites the checkpoint file on exit. `refinish.py` passes the
file itself, but the content written back is the same.

## Checks done on receipt (v2.24)

- **Loading safety.** `pickletools` scan: each pickle references only `numpy._core.numeric._frombuffer` and
  `numpy.dtype`, and no other globals.
- **Completeness.** Every checkpoint has all its fibres (`len(fibres) == nreps`).
- **Reproduction.** `refinish.py` on all ten, with the current code, gives JSON identical to `results.jsonl` in every
  field except the wall time `t`.

| job | N | G-orbits | with grading | G × T_Q | local torus |
|---|---|---|---|---|---|
| 1⇉2, 2→3, 2→4 / 𝔽₅ | 3912 | 456 | 456 | 456 | 414 |
| 1⇉2→3→4 / 𝔽₇ | 5257 | 334 | 334 | 334 | 334 |
| 1⇉2, 1→3, 4→2 / 𝔽₇ | 14928 | 615 | 485 | 485 | 485 |
| 1→2⇉3→4 / 𝔽₇ | 3276 | 257 | 257 | 257 | 257 |
| 1⇉2, 1→3 / 𝔽₁₁ | 8474 | 137 | 101 | 101 | 101 |
| 1⇛2, 1→3 / 𝔽₅ | 128640 | 207 | 183 | 183 | 183 |
| 1⇛2→3 / 𝔽₅ | 54687 | 120 | 120 | 120 | 120 |
| 1⇉2→3⇉4 / 𝔽₃ | 11008 | 1158 | 1158 | 1158 | 1158 |
| 1⇉2→3→4→5 / 𝔽₃ | 5405 | 1770 | 1770 | 1770 | 1770 |
| 1⇉2⇉3 / 𝔽₅ | 251550 | 1112 | 790 | 754 | 754 |
