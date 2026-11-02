# Discordianism: Open Questions Explored (v2.78): the heavy runs; 1⇛2→3⇉4 grows; one chain quiver proved finite

*The user ran v2.78's heavy jobs (`logs/v2_78/heavy.sh`) on their own machines: an AMD Strix Point laptop and a
Radeon Pro W6600. Two engine changes made that possible. The quotient step now budgets its 64-bit copies. And the
engine recovers from allocation failures by halving its batch budgets; fault injection leaves every count unchanged.*
- *Three predictions held, with q = 5 for the first time on a counterexample quiver.*
- *1⇛2→3⇉4 is not flat (2236 → 2240). That refutes the parallel-class exception in Conjecture 3ᶜ″. Theorem 5's
  family does not explain the growth: all four vertex-part families of the path version have no modulus there.*
- *1⇛2 + 1→3 + 4→2 and 1→2⇛3→4 are finite, uniformly in q and in every characteristic. Those are the first "chain"
  quivers of size 3 proved finite beyond 1⇛2→3→4. The other two driver runs each have one or two blocks left over.*

---

## 0. Results

The predictions are in logs/v2_78/PREDICTIONS.md, committed before the runs.

| # | Statement | Status |
|---|---|---|
| 138 | 1⇛2, 2→3, 2→4 at q = 5: 714, only the identity block moves | **held**: 708, 710, 714. Identity block 111, 113, 117, all 147 other blocks flat. This settles v2.74's predictions 119 (it grows) and 120 (706 + 2(q − 1)), at the three q computed |
| 139 | 1⇛2→3, 3→4, 3→5 at q = 3: 4903, only the identity block moves | **held**: 4899 → 4903; identity block 249 → 253; the other 766 blocks flat |
| 140 | 1⇛2→3⇉4 flat: 2240 at q = 3 | **failed**: 2236 → 2240, all in the identity block (338 → 342) |
| 141 | 1⇛2, 2→3→5, 2→4 at q = 5: 3435 | **pending** (still running: fibre group 6 of 13 after 19 hours) |
| 142 | the v2.66 driver proves the four chain quivers finite | **partial**: complete for 1⇛2 + 1→3 + 4→2 and 1→2⇛3→4; three blocks open in the other two (§3) |

## 1. Census (user's machines; `logs/v2_78/census_*`, block splits `*.blockcensus.txt`, `*.blockdiff.txt`)

| quiver | q = 2 | q = 3 | q = 5 | blocks that move |
|---|---|---|---|---|
| 1⇛2, 2→3, 2→4 | 708 | 710 | 714 | identity only: 111, 113, 117 = 109 + 2(q − 1) |
| 1⇛2→3, 3→4, 3→5 | 4899 | 4903 | — | identity only: +4 (v2.76's four families) |
| 1⇛2→3⇉4 | 2236 | 2240 | — | identity only: +4 |
| 1⇛2, 2→3→5, 2→4 | 3417 | 3423 | running | identity +4, e₅-killed block +2 (v2.77) |

**1⇛2, 2→3, 2→4.** The two fork families give at least 2(q − 1) classes in the identity block (v2.75). The identity
block is exactly 109 + 2(q − 1) at q = 2, 3, 5, and every other block is constant. So, at these three q, the census is
the finite part plus exactly the two families:

  census(q) = 706 + 2(q − 1).

## 2. 1⇛2→3⇉4 grows

v2.77 found no modulus in the path family on this quiver: its 2-class absorbs the cross-ratio. That was not a proof of
finiteness, and the census now moves: +4, all in the identity block.

**What the families do and do not give.** `analysis/v2_76_vertexpart.py` with the paths b·c₁ and b·c₂ finds the same
four vertex-part families as on 1⇛2→3, 3→4, 3→5 (X = {1}, {1, 3} and their complements). Run on each of the four,
`analysis/v2_79_tangent.py` gives:
- the member is a fold identically;
- dim 𝔥 = 11 and rank T = 7;
- **0 moduli**, at two random members each.

So the growth comes from somewhere else. There are two possibilities:
- **another family**, outside the path version: Q would be infinite, and the right conjecture counts paths
  (Conjecture 3ᶜ‴, logs/v2_79/PREDICTIONS.md 144);
- **a q = 2 effect**: four geometric orbits whose 𝔽₂-points are missing or merged, for example orbits with disconnected
  stabilisers. Q would then be finite, and the census at q = 2 would undercount (Conjecture G would fail at q = 2 for
  this quiver).

**How to decide.** Run the v2.55 pipeline on the identity block. It decides finiteness in characteristic 0, component
by component. It is running here (its primary decomposition had not finished after an hour), and it is also job
`pairid` in `logs/v2_79/heavy2.sh`. The q = 5 census (job `pair5`) would show which, if any, blocks keep moving.

## 3. The v2.66 driver on the four chain quivers (`logs/v2_78/drv/`)

| quiver | vertex blocks | pipeline | specialisation | status |
|---|---|---|---|---|
| 1⇛2 + 1→3 + 4→2 | 123 | 123 PASS, 513 components | 6473 claims; candidate primes 2, 3, all checked directly; no exceptional prime; bound 10^36.5 | **finite, uniformly in q, every characteristic** |
| 1→2⇛3→4 | 100 | 100 PASS, 248 components | 3779 claims; candidates 2, 3, 5, each checked directly on its claims; no exceptional prime; bound 10^102.5. On the identity block the redundant full cross-check at p = 2 (v2.58's `--always`) timed out | **finite, uniformly in q, every characteristic** |
| 1⇛2→3 + 4→2 | 148 | 147 PASS, 1 UNRESOLVED (0001,0100,0010,0001) | 6056 claims, no exceptional prime | one block open |
| 1⇛2→4 + 1→3 | 109 | 107 PASS, 2 UNRESOLVED (0000,0100,0010,0001; 0000,0000,0000,00m1) | 4747 claims (candidates include 17 and 193, checked directly), no exceptional prime | two blocks open |

**The unresolved blocks.** All three are "FAILED in both modes": the pipeline did not finish within the driver's
6000-second limit, either directly or with the case split. They are not failures of the finiteness criterion. A rerun
here of 0000,0000,0000,00m1 (1⇛2→4 + 1→3) also did not finish its presolve within 15 minutes.
`logs/v2_79/heavy2.sh drivers2` reruns both drivers with a long limit; RESUME keeps the finished blocks.

The census counts are all flat at q = 2, 3: 893, 419, 783 and 640. So these four fit Conjecture 3ᶜ″, and also the
path-count version 3ᶜ‴, which agrees with 3ᶜ″ on them.

## 4. Engine changes (`fibre_torch.py`)

- **Quotient batches** now count pack's int64 copies. They overran `ETF_MEM` about threefold: a 5 GB request with
  `ETF_MEM` = 8192 on the laptop.
- **Out-of-memory recovery.**
  - On an allocation failure (CUDA/HIP, or the CPU allocator) the engine halves every batch budget, down to 16 MiB,
    frees the device cache and redoes the parent batch.
  - Children produced twice disappear in the existing per-fibre deduplication.
  - The prefilter's row reduction had also overrun: 6 GB with `ETF_MEM` = 4096 on the 8 GB W6600.
  - **Checks.** `ETF_FAULT_EVERY` and `ETF_FAULT_MAX` simulate failures; expand failures were injected separately.
    The counts are unchanged: 708 and 710 (1⇛2, 2→3, 2→4), and 414 (1⇉2, 2→3, 2→4).

## 5. Next

1. **1⇛2→3⇉4: finite or not.** Finish the identity-block pipeline (here, or `heavy2.sh pairid`). If a component fails,
   find its family. Run the q = 5 census (`pair5`).
2. **The open driver blocks** (`heavy2.sh drivers2`).
3. **141**, when the q = 5 run of 1⇛2, 2→3→5, 2→4 finishes.

## Files

- **fibre_torch.py:** the batch budgets and out-of-memory recovery (§4).
- **analysis/:** v2_79_tangent.py (infinitesimal moduli with a general vertex part).
- **logs/v2_78/:** the user's census, block and driver outputs; PREDICTIONS.md; heavy.sh.
- **logs/v2_79/:** PREDICTIONS.md (143, 144); vertexpart and tangent outputs for 1⇛2→3⇉4; heavy2.sh.
