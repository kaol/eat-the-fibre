# Discordianism: Open Questions Explored (v2.51): job 14 diagnosed and finished in the container

*Job 14 (1⇉2⇇3 + 1→3 + 2→4 over 𝔽₅) ran for days on the GPU machine and then ran out of memory inside one fibre. This
round finds the cause: the column order. The fibre search fills one column of P at a time. Under the fixed order, some
fibres reach the column that would prune them only after a wide stretch of unconstrained columns, so the frontier
explodes. Which columns cause this depends on the fibre, so no single order suits every fibre. The new engine driver
tries a short list of orders per fibre and keeps the first one that stays small. With it, job 14 finished in the
2-core container in under an hour, and its local-torus count is 1220, the same as over 𝔽₃.*

---

## 0. Results

The predictions are in logs/v2_51/PREDICTIONS.md; some were added in stages, each before its runs.

| # | Statement | Status |
|---|---|---|
| 37 | Fewer than 1% of the 12210 q = 3 source fibres reach a frontier of 1000 rows | held: 11 fibres (0.09%) |
| 38 | In the heaviest q = 3 fibres, the peak is at a middle column and at least 10× the final leaf count | held: fibre 5973 had more than 1.3 M rows at column 12 (time cap) against 5 leaves; fibre 14916 had 31,214 rows at column 10 against 1 leaf |
| 39 | At q = 5, under the default order, some fibre exceeds 2·10⁶ rows | not tested: the default order was never profiled at q = 5 |
| 40 | With order O₁ (paths through c before the tail arrow), the q = 3 census equals job 13 | held: all 21 fields equal; 100 s on CPU, against 666 s for the original run on the Strix GPU |
| 41 | With O₁, every q = 5 fibre stays below 2·10⁶ rows | not decided: fibre 16019 reached 503,274 rows when the time cap stopped it |
| 42 | With order O₂ (each path followed at once by its extensions), q = 3 equals job 13 and every q = 5 fibre finishes | **failed**: four q = 3 fibres explode under O₂ that O₁ handles easily |
| 43 | Job 14's orbits_GxTloc = 1220, the q = 3 value | **held** |
| 44 | With the portfolio (O₁, then O₂), every q = 5 fibre is finished by one of the two orders and the census completes in the container | **held**: 11,273 fibres by O₁ and 2 by O₂ (16019, 17763); none left over |

The orders, as basis indices, are O₁ = 4,6,8,14,15,9,10,12,11,13 and O₂ = 4,6,8,14,15,9,12,16,13,17,10,11. In paths:

| order | sequence |
|---|---|
| default | a₁, b₁, c, t, a₁t, a₂t, b₁t, b₂t, cb₁, cb₂, cb₁t, cb₂t |
| O₁ | a₁, b₁, c, cb₁, cb₂, t, a₁t, b₁t, a₂t, b₂t, cb₁t, cb₂t |
| O₂ | a₁, b₁, c, cb₁, cb₂, t, b₁t, cb₁t, b₂t, cb₂t, a₁t, a₂t |

## 1. Diagnosis

The engine fills the columns X[:, c] = P·b_c in a fixed order. Each column has p^{dim K} = p⁶ candidates. A condition
is checked only once every column it involves is filled. When the columns that close a constraint (a product such as
P(a₂)·P(t), or a longer path's Δ-lift) come late, every candidate in between survives. The frontier then grows by
factors of about p³ per column and collapses later.

Examples from the profile:
- **Fibre 5973 at q = 3**: 4,696 rows → more than 1.3 M rows at b₁t under the default order. Under O₁ it peaks at
  3,314 rows and finishes in 9 s.
- **Fibre 16019 at q = 5**: 1,871 rows → 503,274 rows at the second t-path under O₁. Under O₂ it peaks at 627 rows and
  finishes in 19 s.
- **Fibres 14964, 15009, 17663, 17727 at q = 3**: under O₂ they explode at columns 13–16, where O₁ keeps them in the
  hundreds.

The depth-first engine bounds memory per stacked batch, not per fibre. With 18 columns, a fibre that explodes in the
middle columns can keep a gigabyte-scale batch alive at many levels at once. That is how job 14 ran for days, slowly
filling 64 GB.

## 2. The fix

All three are in the repository:
- `pcensus.setup(name, p, order=…)` and the environment variable `ETF_CORDER` choose the column order.
  - A checkpoint now records its column indices, and resuming under another order is refused, because fibres are
    stored in the column basis.
- `fibre_torch.portfolio`, enabled by `ETF_CORDERS = o₁;o₂;…`, works one fibre at a time. It runs a breadth-first
  search with the engine's own `level` and `quotient` under o₁, then o₂, and so on.
  - Each attempt is capped at `ETF_PCAP` rows (default 3·10⁶). Every order but the last is also capped at `ETF_PTIME`
    seconds (default 300).
  - The first success is kept and converted to the main column basis: X_main = X_o · B_o⁻¹ · B_main.
  - Fibres that no order manages are left for the old depth-first `run()`.
- `run_batch.py` jobs can carry environment variables; job 14 now carries the two-order portfolio.

Correctness checks:
- At q = 3, with the orders reversed and the cap lowered to 5·10⁴ so that both orders are used (12,206 fibres and 4
  fibres), the census equals job 13 in all 21 fields.
- All 17,748 fibres are identical as sets of folds to those of the single-order run.

## 3. Job 14

| | q = 3 (job 13) | q = 5 (job 14, this round) |
|---|---|---|
| folds N | 25,622 | 442,960 |
| G-orbits | 1,443 | 2,129 |
| G × grading-torus orbits | 1,281 | 1,387 |
| **local-torus orbits** | **1,220** | **1,220** |
| run time | 666 s (GPU) | about 47 min (2-core CPU, with an earlier interrupted part) |

1⇉2⇇3 + 1→3 + 2→4 is unobstructed: the converging pair has joined outer ends, plus a tail. Its local-torus count is
constant across q = 3 and q = 5, as Conjecture 3ᶜ predicts. The coarser counts grow, and the local class torus absorbs
the growth, as on the other joined-ends quivers. This is finite-field evidence, not a proof. G-orbit representatives
are saved in logs/reps_pairs_in_short_tail_5.pkl.

**For the Strix machine: job 14 is done and there is no need to rerun it.** If it is ever rerun, `run_batch.py 14`
uses the portfolio automatically.

## 4. Next

1. ChatGPT's Program 2: the twist-fibre theorem for the joined-ends configuration 1⇉2⇇3 + 1→3. v2.29 already proved in
   Lean that Theorem 3's families fail on it, in the Δ condition of the arrow 1→3. The open part is the generic lemma
   (L): over k̄, every twist fibre is a finite union of census-symmetry orbits.
2. The 5-dimensional components of the divisors of V_A and V_B (v2.50).
3. Profile the other heavy censuses (for example two_pairs_short at q = 7) to see whether column order is their bottleneck too.

## Files

- analysis/v2_51_profile.py: per-fibre column profile; options `--only`, `--cap`, `--tcap`; env `ETF_CORDER`, `KCACHE`.
- pcensus.py: `corders`, `column_order`, `setup(…, order)`, the checkpoint guard, and the portfolio hook.
- fibre_torch.py: `bfs_fibre`, `portfolio`, `ETF_PCAP`, `ETF_PTIME`.
- run_batch.py: per-job environment variables.
- logs/v2_51/:
  - PREDICTIONS.md;
  - profile_3_all.jsonl, profile_3_O7_all.jsonl, and profile_5_O5_all.jsonl (partial, 13,215 fibres);
  - results_tail3_{O5,portfolio}.jsonl and results_tail5.jsonl, with their run logs.
- logs/reps_pairs_in_short_tail_5.pkl
