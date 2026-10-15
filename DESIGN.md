# eat-the-fibre: GPU port of the v2.23 fold census

Dense reference for AI readers. Human quickstart is in `README.md`.

## 0. Object

`pcensus.py QUIVER p` enumerates the Δ-compatible folds P of kQ over F_p (idempotent P, image and kernel subalgebras,
ΔP = (P⊗P)Δ). Kernel fibration: K = ker Δ. P|_K is an algebra fold of K (`kernel_folds`), one rep per grading-torus
orbit. Over each rep the fibre is filled column by column in the new basis (K first, then the complement in length
order): `X[:, c] = L(X_{<c}) + k`, with L the unique Δ-preimage in the complement and k ∈ K. A masked check prunes
after every column. The hot loop is `extend = check ∘ enumerate ∘ lift` (the v2.23 version, `verify/v2_23/pcensus.py`, runs it in a process pool).
`run_batch.py` sequences the long jobs. This repo adds a second fibre engine, `fibre_torch.py`, behind
`pcensus.py --engine torch`. `run_batch.py` uses it by default. v2.25 adds the kernel stage on the GPU
(`kernel_torch.py`, §6), the H = G × T symmetry of the fibration (transport and stabilizer quotient, §7), a vectorized
finish (§8) and a kernel cache (§6).

Pipeline (the stages are an Arrow composition; only the frontier flows):

    setup ─▶ kernel stage ─▶ symmetry plan ─▶ fibres (sources only) ─▶ transport ─▶ finish
             §6, memoized    §7             §1–3 + §7 quotient       §7          §8

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

### Measured (torch engine on 4 CPU threads, this container; the pool engine with 4 workers)

| job | torch | pool | per column: candidates checked / children (worst column) |
|---|---|---|---|
| pair_tgtfork / F_3 | 52 s | 564 s | c=10: 179 650 / 64 906 (brute force: 91 286·3^5 = 22.2 M) |
| K3_after / F_3 | 69 s | 706 s | c=8: 784 280 / 248 704 (brute force: 74 130·3^5 = 18.0 M) |
| two_pairs / F_3 | 55 s | 1260 s (logged) | c=9: 364 130 / 95 268 (brute force: 35 630·3^6 = 26.0 M) |

The JSON output is identical between engines; for two_pairs it matches the logged CPU run except `transversal.lifted`
and `lifted_transversal` (§4). The deep columns carry the frontier mass, where the remaining nonlinear conditions (V3,
ineligible V1/V2) leave ≈ 2–4 candidates per child. That is the next lever if needed (§9). On a GPU the same tensors
run with larger batches; raise `ETF_MEM` there.

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
  `logs/results_here.jsonl` has two_pairs/F_3 with lifted 343 vs 336, and K3_after/F_3 with 62 vs 59. With
  canonical (lexicographically first) representatives, two_pairs/F_3 gives lifted 597 and lifted_transversal **21**,
  where the logged runs give 23. So lifted_transversal is representative-dependent too. Canonical ordering makes these
  fields deterministic, but they are still not invariants: read them as "for the lexicographically first
  representative". An invariant version would count orbits with *some* member that lifts; that is a change in the
  mathematics, not in the port, so it is not made here. `transversal_modp`, N, the fibre histogram and the orbit counts are invariant.

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

## 6. Kernel stage (v2.25): `kernel_torch.py`

The v2.23 kernel stage (`kernel_folds`, kept as the reference) loops in Python over every subspace of K = F_p^kdim and
tests each with ranks. For two_pairs/F_7 that is Σ_k [6 k]_7 = 61 946 920 subspaces, single-threaded: hours, which is
what stalled job 11. Mu–Nu reading: the subspace set is a static index set, so the Arrow condition holds and the stage
is a hylomorphism, `quotient ∘ join ∘ filter ∘ unrank`, with nothing data-dependent in its shape.

- **Mu (index set).** Gr(k, K) = ⊔_π Schubert cells: the RREFs with pivot set π form F_p^f, f = #free entries (row i,
  column j ∉ π, j > π_i). Index t ∈ [0, p^f) unranks to free values in `itertools.product` order (first free entry =
  most significant digit). So (k, π, t) order is the reference order, and a cell is one batched tensor.
- **filter (cata, ∧-semilattice).** S is a subalgebra ⇔ r_a·r_b ∈ span S for all k² row pairs. Because S is in RREF,
  membership is linear: v ∈ span S ⇔ v − v[π]·S = 0. No ranks are needed. The k² predicates commute under ∧, so they
  run one after another on the survivors (absorption/reordering, as in §1): the first sees the batch, the rest see
  ~nothing. MK (K's structure constants) is sparse (10 of 216 nonzero for 1⇉2⇉3), so a product is nnz gathered terms.
- **join.** Folds of K ↔ complementary pairs (H, Q): P = Bt[:, :k]·(Bt⁻¹)[:k], Bt = [H; Q]ᵀ. Invertibility and inverse
  both come from one batched Gauss–Jordan (`fibre_torch.rref(…, ncols=n)` on [Bt | I]). Order: k, then H, then Q.
- **quotient (torus).** All conjugates A_l PK A_l⁻¹ come from one GEMM, vec(PK)·(A_l ⊗ A_l⁻ᵀ). One `np.unique` over the
  byte rows of Kf and of the conjugates maps each conjugate to its index. Orbit label = min index, so the reps and their
  l-lists are the reference loop's.
- **Peer equivalence, strongest observation.** `kernel_stage(…, 'torch')` and `(…, 'ref')` return *equal lists* (same
  folds, same order, same reps and l-lists). Order matters: checkpoints key fibres by rep index, and the
  representative-dependent fields (§4) depend on materialization order. `tests/test_kernel.py` checks 7 quivers, and
  `tests/test_engines.py` runs pool+`--kernel ref` against torch+`--kernel torch` through the whole pipeline.
- **Memoization (L2: seed purity).** (Kf, reps) is a pure function of (p, M, DL, Kmat, Binv). The cache
  `CHECKPOINT.kernel.npz` (or `--kernel-cache`) stores it under `kernel_digest`, a sha256 of exactly that seed plus a
  version number. The version manifest is in the seed, so a stale cache is detected, never served. A resumed run skips the
  stage.

| quiver / p | subspaces | subalgebras | kernel folds | torus reps | reference | `kernel_torch` (4 CPU threads) |
|---|---|---|---|---|---|---|
| two_pairs / 3 | 56 632 | 573 | 20 864 | 10 984 | 12.6 s | 0.3 s |
| two_pairs / 5 | 3 583 232 | 1 617 | 196 788 | 49 984 | ≈ 10 min (est.) | 5 s |
| two_pairs / 7 | 61 946 920 | 3 429 | 915 048 | 153 480 | hours | 29 s (13 s subalgebras, 6 s pairs, 10 s torus) |

## 7. Symmetry: H = G × T acts on the fibration

G (column-sum-1 substitutions inside parallel classes) and the grading torus T commute. Every h ∈ H is an algebra
automorphism with Δh = (h⊗h)Δ (the column sums are 1), so h maps folds to folds. K is H-stable (h|_K = hK). Hence the
kernel fibration is **H-equivariant**: fibre(h·PK) = h·fibre(PK)·h⁻¹. The fibres were computed as if the leaf-data space
had no quotient structure. It has one (Mu–Nu: a leaf-data quotient), and v2.25 uses it twice. Both uses are exact, never
heuristic. In the new basis, X = P B, so h acts by X ↦ Lh X Rh with Lh = h and Rh = Binv h⁻¹ B.

1. **Transport (orbits of torus reps).** `symmetry_plan` computes Mx[r, q] = rep index of hK_q PK_r hK_q⁻¹. The H-orbit
   of PK_r is H·PK_r, so source(r) = min of row r. Only sources are computed. Every other rep gets
   fibre(r) = h·fibre(source) (`transport`, which asserts that X'[:, :kdim] = Kmat PK_r). Mu–Nu: lossiness *with
   retraction*. The section is the transport map.
2. **Kernel quotient of the frontier (torch engine).** N = ker(H → H|_K) fixes every PK, so it acts on every fibre
   (|N| = 294 for two_pairs/F_7: |H| = 10 584, |H|_K| = 36). At column c, the subgroup
   N_c = {n : Rn[c:, :c] = 0} acts on column prefixes: (n·X)[:, :c] = Ln X[:, :c] Rn[:c, :c]. For n ∈ N_c,
   completions(n·X) = n·completions(X). `quotient` keeps one frontier member per (fibre, N_c-orbit). The orbit key is
   the lexicographically smallest image of columns kdim..c−1 (the K columns are fixed by N), packed into int64 words of
   base-p digits. `expand` applies N to the leaves, and `run` deduplicates each fibre (orbits that were split across
   batches expand to the same folds). Invariant: fibre = N·completions(frontier).
   - **Global-L1 trap, avoided.** The engine skips conditions "already verified at the parent" (`ver`, inherited from
     `pcensus.check`). That is a global invariant of the frontier: every node passed every condition determined so far.
     An *image* n·X never passed those checks, and the masks are basis-dependent, so feeding images into the BFS could
     leak non-folds. So `quotient` keeps genuine frontier members (the first of each orbit) and uses images only as keys.
     Merging only on *equal* images is sound for any subset of N_c. Exactness of the orbit key is a performance
     property, not a correctness one.
   - **Measured and rejected: the full stabilizer.** Stab(PK) ⊋ N on symmetric kernel folds (up to all of H). Quotienting
     by Stab(PK)_c shrank two_pairs/F_5 further (col 7: 7 508 instead of 10 928 kept). But the cost is |Stab_c| images
     per node, up to 10⁴ on the symmetric fibres that carry most of the mass, and on F_7 that was slower overall.
     The kernel N is one group for all fibres and costs ≤ |N| images per node.
   - `ETF_QUOTIENT=0` disables it (`--no-symmetry` disables transport and quotient).

| quiver / p | torus reps | H-orbits (sources) | frontier with the N-quotient (torch) |
|---|---|---|---|
| two_pairs / 5 | 49 984 | 14 206 | col 7 kept 10 928 of 87 258; col 9 candidates 4 889 (327 948 before); col 10 parents 2 689 (184 968) |
| two_pairs / 7 | 153 480 | 27 994 | col 7 kept 31 892 of 471 538; col 9 candidates 13 471; col 10 parents 6 591 |
| K3_after, K3_srcout / 5 | 13 720 | 598 | |
| single-pair quivers | — | = torus reps (H\|_K = T) | |

Verification:
- On all 10 committed checkpoints (`data/checkpoints`), transport from the sources reproduces *every* stored fibre as a
  set (156k fibres; 4 quivers have nontrivial transport).
- Whole-JSON equality with transport + quotient: two_pairs/F_5 and two_pairs_sep/F_3 equal `data/results.jsonl`, and
  two_pairs/F_3 equals its `--no-symmetry` run. `tests/test_engines.py::test_symmetry_transport` compares K3_srcout/F_3
  with `--no-symmetry`.
- two_pairs/F_7: the G-kernel quotient (N ∩ G, 49 elements, an earlier variant) and the H-kernel quotient (294) give
  identical fibres for all 153 480 reps.
- The pool engine has transport but no quotient, so it stays the plain reference.

## 8. Finish, vectorized

`finish` was per-fold Python (≈ 90 s of 96 s for two_pairs/F_5). Now:
- materialization: one int8 (N, d, d) array; φ_l is diagonal, so φ_l P φ_l⁻¹ = P ∘ (diag φ ⊗ diag φ⁻¹), batched per
  fibre; uniqueness by `np.unique` over byte rows (`fold_keys`).
- `torus_preserving`: Folds = set(lab) is a union of G-orbits, and φ_w commutes with G (one scalar per parallel class).
  So φ_w preserves Folds iff it maps every G-orbit representative into it: 1 506 checks per w instead of N = 1 682 084
  for two_pairs/F_7. (`np.isin` on byte-row keys was tried first. It falls back to a quadratic path: 20 min at F_7.)
- local class torus: compatibility with Δ_w for all w at once (two einsums per G-orbit rep).
- `orbit_label`/`orbit_stats`: g P g⁻¹ as batched matmul instead of a 5-index einsum (14× per orbit), with orbit keys
  from `np.unique` over byte rows.
- `transversal`: `q_fold` as batched integer contractions. The Q-rank is `rank_Q`, fraction-free elimination in Python
  ints (exact; hypothesis-tested against sympy), replacing `sympy.Matrix.rank` (28 ms per call).
- `build_G` is memoized in ctx (the symmetry plan and finish share it).

Measured (4 CPU threads): two_pairs/F_5 finish 89 s → 14 s; two_pairs/F_7 finish (N = 1 682 084, |G| = 1 764)
≈ 40 s. `refinish.py` on all 10 committed checkpoints reproduces `data/results.jsonl` in every field but t.

## 9. Further work (not needed for correctness)

- Full-stabilizer quotient (§7, rejected as implemented). It would pay if orbit keys came from a few generators of Stab_c
  (union-find over generator moves inside the frontier: sound by the equal-image argument) and `expand` closed the
  leaves under generators of Stab(PK). The cost would then be ~|gens| per node and per fold, instead of |Stab_c| and |Stab|.

- Incremental elimination: fold row blocks into a (kdim+1)-row reduced state (meet of affine subspaces) to cut the
  `Parents`/prefilter working set by ~|rows|/kdim and allow larger parent batches.
- Orbit stage on GPU: canonical form = lexicographic min of the G-conjugates (a quotient map) for all folds at once.
  Not needed at current sizes (§8: ≈ 40 s at F_7); orbit labelling is per orbit, not per fold.
- V3 (W·W) is quadratic in κ. When the nullspace is 1-dimensional it becomes a univariate quadratic over F_p, solvable
  in closed form. Only worth it if measured candidates/children stays well above 1.
