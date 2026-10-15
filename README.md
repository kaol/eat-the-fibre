# eat-the-fibre

Census of Δ-compatible folds of path algebras over finite fields, plus the mathematical project it serves: the
Discordianism Mu–Nu axioms and their open-question explorations (`docs/`, start at `docs/README.md`), with the machine
checks for every round (`verify/`) and the analysis scripts used from v2.23 on (`analysis/`).


The v2.23 fold census (`pcensus.py`, `run_batch.py`) with a second fibre engine, `fibre_torch.py`. It runs on a GPU
through PyTorch (ROCm/HIP on AMD, e.g. the Radeon 890M in Strix Point, or CUDA), or on the CPU through torch.
Same mathematics and same output JSON as the process-pool engine, with a shared checkpoint format.
Design and correctness argument: `DESIGN.md`.

Besides batching, the engine solves the conditions that are linear in the K-component exactly (batched Gaussian
elimination mod p). It then enumerates only that solution space instead of all p^dimK candidates. Survivors are
identical (checked per parent against the reference). The number of checked candidates drops by 2–3 orders of
magnitude at p = 3, and by more at larger p.

## Setup on Strix Point (Linux, ROCm)

```bash
python3 -m venv ~/etf && . ~/etf/bin/activate
pip install numpy scipy sympy pytest hypothesis
# PyTorch with ROCm: pick the newest ROCm index listed at https://pytorch.org/get-started/locally/
pip install torch --index-url https://download.pytorch.org/whl/rocm6.4
python3 -c "import torch; print(torch.cuda.is_available(), torch.cuda.get_device_name(0), torch.cuda.get_arch_list())"
```

The iGPU is gfx1150 (RDNA 3.5). If `get_arch_list()` does not include `gfx1150` and kernels fail with "invalid device
function", either install a wheel built for gfx1150 (AMD's TheRock nightlies) or try
`export HSA_OVERRIDE_GFX_VERSION=11.0.0`. Whatever you choose, run the validation below before the batch. The engine
refuses to start if device float32 matmul is not exact (`Plan.selftest`), and `ETF_XCHECK=1` compares every step with
the numpy reference.

The iGPU shares system RAM. `ETF_MEM` (MiB, default 1024) bounds one working set. Raise it if the UMA frame
buffer/GTT allows, and lower it if you hit HIP out-of-memory errors.

## Validate

```bash
python3 -m pytest -q tests                                     # solver property test + pool ≡ torch, whole JSON
ETF_XCHECK=1 python3 pcensus.py K3_srcout 3 --engine torch     # per-parent cross-check on the GPU (~1 min)
```

## Run

```bash
python3 run_batch.py --list        # jobs and predictions
python3 run_batch.py               # all jobs, torch engine on the GPU; resumable (Ctrl-C safe)
python3 run_batch.py 3 5           # only jobs 3 and 5
ENGINE=pool JOBS=16 python3 run_batch.py      # the original process-pool engine
DEVICE=cpu python3 run_batch.py               # torch engine on the CPU
```

Single job: `python3 pcensus.py two_pairs 5 --engine torch --checkpoint ck.pkl --out results.jsonl`.
Send back `results.jsonl` and `batch.log`. At the end of every run the engine logs, per column,
`[parents, lifted, candidates, children]`. These numbers show where the time goes.

| env | default | meaning |
|---|---|---|
| `ENGINE` | `torch` | `run_batch.py`: `torch` or `pool` |
| `DEVICE` | auto | torch device: `cuda` (also ROCm), `cpu`, `cuda:1`, … |
| `ETF_MEM` | 1024 | MiB per working set (parent batch / candidate chunk) |
| `ETF_GROUP` | 4096 | fibres per checkpoint group |
| `ETF_XCHECK` | 0 | 1 = compare lift/children with `pcensus.lift`/`check` for every parent (slow) |
| `JOBS`, `PCENSUS_*` | | pool engine, as before |

## After a batch: re-finish from checkpoints (v2.24)

v2.24 adds `orbits_GxTloc`: G-orbits merged along local class-torus moves (see
`docs/explorations/discordianism-open-questions-explored-v2.24.md` §2). Finished checkpoints already hold every fibre, so the
new statistics need only the orbit stage:

```bash
python3 refinish.py                # every complete checkpoint in checkpoints/ -> results_v224.jsonl, with predictions
python3 refinish.py pair_tgtfork   # just one
```

## Files

- `gcensus.py`: quiver data, path algebra, `Field` (unchanged reference).
- `pcensus.py`: setup, kernel folds, torus reps, pool engine, `build_G`, orbit statistics (G, G × grading, G × T_Q, local torus), CLI. Changes from the v2.23 version (`verify/v2_23/pcensus.py`):
  `--engine/--device`, main split into `pool_fibres`/`finish`, and fibres sorted canonically before materialization.
- `fibre_torch.py`: the torch engine (`Plan`, `lift`, `Parents`, `prefilter`/`solve_affine`/`rref`, `check`,
  `level`, `run`).
- `run_batch.py`: the v2.23 batch, using the torch engine by default. It now also prints `orbits_GxT` and `orbits_GxTloc`.
- `refinish.py`: recompute orbit statistics from finished checkpoints (no fibres recomputed).
- `analysis/`: exact tangent spaces, local dimension, class-torus and orbit checks (`orbitcheck.py`); run from anywhere, they
  read `logs/folds_QUIVER_p.pkl` (`pcensus.py ... --save-folds logs/folds_QUIVER_p.pkl`).
- `docs/`, `verify/`: the project documents and the per-round verifiers (see `docs/README.md`).
- `tests/`: `test_solve_affine.py` (hypothesis vs brute force) and `test_engines.py` (whole-pipeline equivalence and
  cross-engine checkpoint resume).
- `logs/results_here.jsonl`: reference results from the CPU runs.

Note: `transversal.lifted`/`lifted_transversal` depend on which orbit representative is chosen (see `DESIGN.md` §4).
Older CPU results vary between runs for that reason. They are now deterministic, but they are still not orbit
invariants.
