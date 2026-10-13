# eat-the-fibre: GPU port of the v2.23 fold census

Dense reference for AI readers. Human quickstart is in `README.md`.

## 0. Object

`pcensus.py QUIVER p` enumerates the Δ-compatible folds P of kQ over F_p (idempotent P, image and kernel subalgebras,
ΔP = (P⊗P)Δ). Kernel fibration: K = ker Δ. P|_K is an algebra fold of K (`kernel_folds`), one rep per grading-torus
orbit. Over each rep the fibre is filled column by column in the new basis (K first, then the complement in length
order): `X[:, c] = L(X_{<c}) + k`, with L the unique Δ-preimage in the complement and k ∈ K. A masked check prunes
after every column. The hot loop is `extend = check ∘ enumerate ∘ lift` (the Le Chat version runs it in a process pool).
`run_batch.py` sequences the long jobs. This repo adds a second fibre engine, `fibre_torch.py`, behind
`pcensus.py --engine torch`. `run_batch.py` uses it by default.

## 1. Mu side: the plan is static

| Mu element | content | data-dependent? |
|---|---|---|
| `setup(name,p)` | basis, M, DL, K, B, Binv, Cidx, DC, Linv | no |
| column schedule | the word c_kdim … c_{d−1} | no |
| per-column Δ-slice | nonzeros of DL[Cidx[c−kdim]], and the static assert that the lift needs only columns < c (`Plan.col`) | no |
| condition set at c | C1 = P b_i idempotent, C2 = P b_i·P b_j ∈ H, C3 = (1−P)b_i·(1−P)b_j ∈ ker, V1/V2 = U_i·W, W·U_j, V3 = W·W, V4 = idempotency of W; all masked by "Binv h has no support ≥ c+1" | layout no; masks per parent |

**Arrow condition holds.** Neither the schedule nor the condition layout depends on the data, so
`census = level_{d−1} ∘ … ∘ level_kdim` is an Arrow composition. Each `level_c : Frontier_c → Frontier_{c+1}` is pure. The
frontier (leaf-data measure) is the only thing that flows. Consequence: one batch can mix parents from **all** fibres at
the same column (level-synchronous BFS across fibres). The pool engine's per-fibre DFS tasks (STRIDE, MAXB) are a
Nu-side scheduling artifact, not structure, and they do not carry over.

**level_c is a hylomorphism**: ana (candidates from a seed) followed by cata (filter). The filters are predicates that
combine by ∧, a meet-semilattice (commutative, associative, idempotent, unit ⊤). Two laws carry the design:
- *reordering/fusion* (commutativity, associativity): conditions may run in any order or block, which justifies batching all conditions into a few
  tensor ops and evaluating cheap ones first.
- *absorption* `pre ∧ chk = chk` whenever chk ⇒ pre. The prefilter's equations are a subset of the check's conditions,
  encoded exactly. So `filter chk ∘ enumerate(sol(pre))` = `filter chk ∘ enumerate(L + K)`: identical survivors, less
  work. This is the whole soundness argument for §2.

## 2. The prefilter: the leaf-data measure, changed

Write k = Kmat·κ, κ ∈ F_p^kdim. The brute force enumerates all p^kdim values of κ (that is `allK`). Classify each
condition h ↦ [P h = h] (or [P h = 0]) by how h depends on κ:

- **constant** (C1–C3; H = 0): the residual is `cc_c·(L + Kκ) − (T − R0)`. If cc_c ≠ 0 the condition *pins*
  Kκ = v_q (Kmat has full column rank). If cc_c = 0 it is a consistency bit.
- **affine with κ-independent mask** (V1, V2, V4 with h = h0 + Hκ; eligible iff (Binv H)[c:] = 0 and
  (Binv h0)[c+1:] = 0). Then cb[c] = c0[c] is constant, the mask is constant, and the residual
  `X_{<c}(c0_{<c} + CH_{<c}κ) + (L + Kκ)c0_c − target` is affine in κ.
  V4 is always eligible, because Binv·Kmat = [I; 0] and c ≥ kdim.
- **nonlinear** (V3 = W·W, and ineligible V1/V2): left to the check.

Per parent this is a linear system over F_p with ≤ d + (4c+1)·d rows and kdim unknowns. Batched RREF over all parents
gives either ∅ or κ0 + span(N), p^f points with f = kdim − rank. The check then runs on those points only.
Measured on samples (probe, F_3): two_pairs needs 1 510 candidates against 675 054 for brute force and 1 264 true
survivors; K3_srcout needs 950 against 125 388. At p = 5 the ratio grows as 5^(kdim−f). The workload shifts from
"p^kdim tiny checks per parent" to "one small batched elimination per parent": GPU-shaped either way, and far less work.

## 3. Nu side: the driver

`fibre_torch.run` is a Mealy machine. State: (stack of (column, batch), fibres dict). Step: pop a batch. If it is
larger than `Plan.parent_batch(c)`, split it and push the pieces. Otherwise push `level_c(batch)`, or record it if c = d.
Output: progress lines, plus a checkpoint after each group of `ETF_GROUP` fibres. Productivity and termination: the
measure (d − c, |batch|) decreases lexicographically on every step, because a split happens only when |batch| > nb ≥ 1.
Memory is bounded by the batch budgets (`ETF_MEM`), not by fibre sizes. Inside a level, candidates are processed
in flat chunks of `Plan.cand_batch(c)`, so a single parent with large p^f still fits.

**Peer equivalence (bisimulation).** Engines `pool` and `torch` agree at the strongest available observation,
`children(X, c)` as a set, per parent and per column. `ETF_XCHECK=1` asserts this against `pcensus.lift`/`pcensus.check`
over all of K. Both engines write the same checkpoint (`{name, p, nreps, fibres: {r: [int16 (d,d)]}}`), so a run may
switch engine mid-way (`tests/test_engines.py::test_checkpoint_crosses_engines`). The observable fold is `finish()`.
It now sorts every fibre canonically, so the JSON is a function of the fibre *sets* only. Engine- and
schedule-independence of the output is therefore exact, and `tests/test_engines.py` compares whole JSON lines.

## 4. Constraint levels

- **L1 exact arithmetic (substrate + tool).** Residues are float32 in [0, p). Every contraction goes through `Fp.mm`,
  which asserts `len·(p−1)² < 2^24`. `Plan` asserts `d²(p−1)² < 2^24` (the longest contraction is the d² Δ-tensor in
  the lift). TF32 is disabled, and `Plan.selftest` compares a max-length device matmul with int64 numpy. This turns the
  one silent failure (an inexact matmul, e.g. TF32/xf32, giving wrong residues) into a loud one at startup.
  Integer matmul is unavailable on CUDA/ROCm; fp64 runs at 1/16–1/32 rate on RDNA, so float32 is the right carrier.
- **L2 engine contract.** XCHECK (above), and `tests/test_solve_affine.py`: hypothesis checks `solve_affine` against
  brute force on random, rank-deficient and inconsistent systems over F_2, F_3, F_5, F_7.
- **L3, interpretive (inherited from the reference, flagged).** `transversal.lifted` and `lifted_transversal` count
  properties of the *chosen representative* of each orbit (whether its residues lift to an exact Q-fold with p−1 ↦ −1).
  They are not orbit invariants. In the pool engine they depended on task-completion order:
  `logs/results_here.jsonl` has two_pairs/F_3 with lifted 343 vs 336, and K3_after/F_3 with 62 vs 59. Canonical
  ordering makes them deterministic, but they are still not invariants: read them as "for the lexicographically
  first representative". `transversal_modp`, N, the fibre histogram and the orbit counts are invariant.

## 5. Structural axioms: working set (Mode 1)

Kept: CP, CL, RSing, LC, CR, RS, TSA. Dropped as not load-bearing for an engine port: SD, CDB, PA, CI. Each design
choice is traced to an axiom (dominant-pattern audit):

| choice | axiom trace | rejected dominant pattern |
|---|---|---|
| one CLI and one post-processing path; the engine is a flag | RSing: one canonical implementation of the math and of `finish` | a forked `gpu_census.py` with copied orbit code |
| the prefilter is a strict subset of the check, never a replacement | RS (absorption law, §1), CP (soundness = subset, testable) | a GPU-only "faster check" with its own semantics |
| `Plan` (static) / `Parents` (per-batch constants) / `level` (pure) / `run` (driver) | CR (Arrow composition), CL (static data in one record) | scheduler interleaved with math, as in the pool loop |
| float32 residues with an asserted bound and a self-test | LC (removes a class-2 trap: silent TF32 inexactness), CP (error at construction) | "use fp32 and hope", or int64 einsum (unsupported on GPU) |
| shapes in comments next to every tensor op | CL, TSA (the shape is the type) | untyped tensor plumbing |
| canonical sort of fibres before materialization | RSing (one output per fibre set), and it exposes the L3 issue | order-dependent output |

Remaining weak spots: TSA ~5 (the reference's terse names `cc`, `sel`, `ver`, `Tc` are kept for 1:1 traceability
with `pcensus.check`), CP ~5 (shapes are asserted only by use, as in Python generally).

## 6. Further work (not needed for correctness)

- Incremental elimination: fold row blocks into a (kdim+1)-row reduced state (meet of affine subspaces) to cut the
  `Parents`/prefilter working set by ~|rows|/kdim and allow larger parent batches.
- Orbit stage on GPU: canonical form = lexicographic min of the G-conjugates (a quotient map) for all folds at once.
  Not needed at current sizes; orbit labelling is per orbit, not per fold.
- V3 (W·W) is quadratic in κ. When the nullspace is 1-dimensional it becomes a univariate quadratic over F_p, solvable
  in closed form. Only worth it if measured candidates/children stays well above 1.
