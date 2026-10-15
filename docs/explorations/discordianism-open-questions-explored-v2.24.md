# Discordianism: Open Questions Explored (v2.24)

*Reads the v2.23 batch, which the user ran on a GPU port of the census (`fibre_torch.py`, written by another instance; repository `kaol/eat-the-fibre`). Four predictions made before the runs are confirmed at new primes. One prediction fails, and the failure is explained by a **local** form of the v2.23 class torus. Conjecture 3\* is restated accordingly.*

*Scripts: `pcensus.py` (root; new fields `orbits_GxTloc`, `local_torus_sizes`), `refinish.py` (recomputes orbit statistics from finished checkpoints), `analysis/orbitcheck.py` (exact tangent spaces at one lifted member of every G-orbit, plus a generic-point local-dimension test). Logs are in `logs/`. Raw data: `data/` (the batch checkpoints and results). Labels: **proved**, **verified (exact)**, **numerical**, **conjectured**, **pending**. Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Summary

| # | Question | Answer |
|---|---|---|
| batch | Predictions from smaller primes | **Confirmed** (exact census): 1⇉2, 1→3, 4→2 over 𝔽₇ (14928 folds, 615 G-orbits, 485 with grading); 1⇉2→3→4 over 𝔽₇ (5257, 334); 1→2⇉3→4 over 𝔽₇ (3276, 257); 1⇉2, 1→3 over 𝔽₁₁ (8474, 137, 101) |
| batch | 1⇉2, 2→3, 2→4 over 𝔽₅ | **Prediction failed**: 456 G-orbits (and 456 with grading), against 428 over 𝔽₃. G-orbits = 400 + 14(q − 1): 14 one-parameter families that the grading does not see |
| — | Why | **A local class torus.** For a fold P, the weightings w with P compatible with Δ_w = Σ_c w_c Δ_c form a linear subspace W(P) containing (1, …, 1). Every w ∈ W(P) with nonzero entries moves P to a fold (proved). W is constant on G-orbits (proved). The census's T_Q (weightings preserving *all* folds) is only the part common to every P, and here it is too small. Exactly 28 G-orbits over 𝔽₃ have r_T = r_G + 1 = tangent dimension, all split compatible. Merging along local torus moves gives **414 over 𝔽₃ and 414 over 𝔽₅** **[exact; the 𝔽₅ value was predicted before the re-finish]**, while the global T_Q count stays at 456 over 𝔽₅ |
| v2.23 Conj. 3\* | Restated | **Conjecture 3\*ˡ: Folds(kQ) has finitely many classes under Aut_Δ and local class-torus moves.** It holds in every count except one: **1⇉2⇉3 gives 716 over 𝔽₃ but 754 over 𝔽₅** (every fold there admits the whole torus, so local and global agree). The extra 38 orbits are lower-dimensional orbits with connected-looking stabilizers (§2a). Two readings remain: orbits defined over ℚ(i) (𝔽_q-points only when 4 ∣ q − 1), or one-parameter families of lower-dimensional orbits, which would refute 3\*ˡ. **𝔽₇ decides: 716 or 792** (batch job 11) |
| — | K₃ + arrow over 𝔽₅ | 1⇛2→3: 54687 folds, **120** orbits as over 𝔽₃ (G = G × grading). 1⇛2, 1→3: 128640 folds; G-orbits 191 → 207 = 175 + 8(q − 1); **183** with grading at both primes. So 8 families, absorbed by the grading |
| — | New data over 𝔽₃ | 1⇉2→3⇉4: 11008 folds, 1158 G-orbits. 1⇉2→3→4→5: 5405 folds, 1770 G-orbits. Equal counts with grading prove nothing about families (see 1⇉2, 2→3, 2→4), so the r_T test is still to run |
| — | Engine | The GPU port passes its tests here (7/7, CPU torch) and reproduces the CPU censuses. Its exact affine prefilter cuts candidates by 2–3 orders of magnitude. 1⇉2, 2→3, 2→4 over 𝔽₅, which stalled at 85 CPU-minutes in v2.22, took 2.6 GPU-hours |

---

## 1. The batch

| job | folds N | G-orbits | with grading | prediction | status |
|---|---|---|---|---|---|
| 1⇉2, 2→3, 2→4 / 𝔽₅ | 3912 | 456 | 456 | 428, 428 | **failed** (explained in §2) |
| 1⇉2, 1→3, 4→2 / 𝔽₇ | 14928 | 615 | 485 | 14928, 615, 485 | confirmed |
| 1⇉2→3→4 / 𝔽₇ | 5257 | 334 | 334 | 5257, 334 | confirmed |
| 1→2⇉3→4 / 𝔽₇ | 3276 | 257 | 257 | 3276, 257 | confirmed |
| 1⇉2, 1→3 / 𝔽₁₁ | 8474 | 137 | 101 | 8474, 137, 101 | confirmed |
| 1⇛2, 1→3 / 𝔽₅ | 128640 | 207 | 183 | — | 8 families, absorbed by grading |
| 1⇛2→3 / 𝔽₅ | 54687 | 120 | 120 | — | constant |
| 1⇉2→3⇉4 / 𝔽₃ | 11008 | 1158 | 1158 | — | new |
| 1⇉2→3→4→5 / 𝔽₃ | 5405 | 1770 | 1770 | — | new |
| 1⇉2⇉3 / 𝔽₅ | 251550 | 1112 | 790 | — | 754 under G × T (716 over 𝔽₃); see §2a |

Families show up as G-orbit counts linear in q. Each one-parameter family of G-orbits contributes q − 1 orbits over 𝔽_q:

| quiver | G-orbits | families | orbits of the right group |
|---|---|---|---|
| 1⇉2, 1→3 | 97 + 4(q − 1) | 4 | 101 (grading) |
| 1⇉2, 1→3, 4→2 | 459 + 26(q − 1) | 26 | 485 (grading) |
| 1⇛2, 1→3 | 175 + 8(q − 1) | 8 | 183 (grading) |
| 1⇉2, 2→3, 2→4 | 400 + 14(q − 1) | 14 | 414 (local class torus) |

---

## 2. The local class torus

Let Δ_c delete only the arrows of class c (parallel classes, single arrows included), and Δ_w = Σ_c w_c Δ_c. The class weighting φ_w satisfies Δφ_w = (φ_w ⊗ φ_w)Δ_w (v2.23 §1).

**Lemma 2 (proved).** If a fold P satisfies Δ_w P = (P ⊗ P)Δ_w for a weighting w with all w_c ≠ 0, then φ_w P φ_w⁻¹ is a fold.
*Proof.* As in v2.23's Lemma, which only used compatibility with Δ_w itself: Δ(φPφ⁻¹) = (φ ⊗ φ)Δ_w P φ⁻¹ = (φ ⊗ φ)(P ⊗ P)Δ_w φ⁻¹ = (φPφ⁻¹ ⊗ φPφ⁻¹)Δ. ∎

The condition is linear in w, so W(P) = {w : Δ_w P = (P ⊗ P)Δ_w} is a linear subspace, and it contains (1, …, 1) because P is a fold.

**Lemma 3 (proved).** W(gPg⁻¹) = W(P) for g ∈ G.
*Proof.* g substitutes arrows within each class with column sums 1, so Δ_c g = (g ⊗ g)Δ_c for every c, and hence Δ_w g = (g ⊗ g)Δ_w. ∎

So "P ~ φ_w P φ_w⁻¹ for w ∈ W(P) ∩ (k^×)^classes" together with G-conjugation generates an equivalence on folds that can be computed on G-orbit representatives. `pcensus.py` now does this (`orbits_GxTloc`), asserting that each move lands in the fold set.

**Why the global T_Q was too small.** T_Q = ⋂_P W(P) ∩ torus. In 1⇉2, 2→3, 2→4 (three classes) some folds are not split compatible, so T_Q over 𝔽₃ is only 4 of 8 weightings. The 28 family orbits are split compatible and admit all 8:

| quiver / 𝔽₃ | G-orbits | with grading | with T_Q | with the local torus | weightings admitted per G-orbit |
|---|---|---|---|---|---|
| 1⇉2, 1→3 | 105 | 101 | 101 | 101 | 2 (10 orbits), 4 (95) |
| 1⇉2→3 | 71 | 71 | 71 | 71 | 4 (all) |
| 1⇉2⇉3 | 814 | 728 | 716 | 716 | 4 (all) |
| 1⇉2, 2→3, 2→4 | 428 | 428 | 428 | **414** | 4 (64), 8 (364) |

**Local test [exact; numerical part pending] (`analysis/orbitcheck.py`, `logs/orbitcheck_tgtfork_3.log`).** Every one of the 428 G-orbits of 1⇉2, 2→3, 2→4 has a member that lifts to an exact ℚ-fold. Tabulated by (r_G, r_T, tangent dimension), the torus enlarges the orbit at exactly **28 orbits, (2, 3, 3), all split compatible**. There tangent = r_T, so the fold variety is locally one orbit of Aut_Δ ⋊ T. Grading alone gives r_H = r_G there (transversal count 0), which is why the grading-based census could not see them. At many other orbits the tangent space is larger than r_T (e.g. 82 of type (0, 0, 1) and 14 of type (1, 1, 2)). The generic-point local-dimension test there was stopped by its 2-hour limit before reporting. So this check shows that the torus accounts for the 28 family orbits, but not that no other moduli exist. The 𝔽₅ count tests that globally, and it came out as predicted: **414**.

**A calibration lesson.** On 1⇉2, 1→3, the first version of the generic-point test flagged one class (local dimension 3, orbit rank 2). The flagged points lay about 0.01 from a degenerate fold, where orbit-rank singular values are O(distance) (e.g. 2.0, 0.013, 0, 0). Replacing the single 1e-2 step with a projected random walk (12 tangent steps of 0.05) gives well-conditioned ranks (smallest retained singular value ≥ 0.017) and no flags. That agrees with the 101-orbit count at four primes.

**Gloss.** *When a step forks after the twins, the two branches may keep separate clocks, but only for the splits that do not tie the branches together. A split that binds them must keep one clock for both.*

## 2a. 1⇉2⇉3 over 𝔽₅: 754, not 716

Every fold of 1⇉2⇉3 admits all class weightings (over 𝔽₅, local torus sizes are 16 for all 1112 G-orbits), so H = G × T²,
of dimension 6, acts on the whole fold set. Its orbit count is 716 over 𝔽₃ and 754 over 𝔽₅. `analysis/stabtype.py`
tabulates, for every H-orbit, its size against its orbit dimension r (rank of the Lie-algebra action mod p).
Logs: `logs/stabtype_two_pairs_{3,5}.log`.

| r = orbit dim | dim Stab | \|Stab(𝔽₃)\| | orbits / 𝔽₃ | \|Stab(𝔽₅)\| | orbits / 𝔽₅ | change |
|---|---|---|---|---|---|---|
| 0 | 6 | 144 | 23 | 6400 | 23 | |
| 1 | 5 | 72, 48 | 35, 32 | 1600, 1280 | 35, 32 | |
| 2 | 4 | 36, 24, 16 | 20, 82, 43 | 400, 320, 256 | 22, 82, 43 | **+2** among orbits of size (q − 1)² (\|Stab\| = q²(q − 1)²) |
| 3 | 3 | 12, 8 | 54, 196 | 80, 64 | 62, 196 | **+8** among orbits of size q(q − 1)² (\|Stab\| = q(q − 1)²) |
| 4 | 2 | 4 | 179 | 16 | 207 | **+28** among orbits of size q²(q − 1)² (\|Stab\| = (q − 1)²) |
| 5 | 1 | 2 | 40 | 4 | 40 | |
| 6 | 0 | 1 | 12 | 1 | 12 | |

The orbit types match one for one between the primes. Each stabilizer's point count fits a connected group of the
stated dimension: for example (q − 1)² for the 2-dimensional stabilizers, 4 and then 16. So a finite stabilizer whose
𝔽_q-points depend on q (such as μ₄ × μ₄) does **not** explain the change. Every extra orbit is a *non-generic* orbit
(r = 2, 3, 4 < 6), with 2, 8 and 28 more of them over 𝔽₅.

Two readings fit both primes:
- **Galois.** Some H-orbits are defined over ℚ(i) and not over ℚ, so they have 𝔽_q-points only when 4 ∣ q − 1 (true for
  5, false for 3). Then there are finitely many orbits over the algebraic closure, and **𝔽₇ gives 716**.
- **Families.** One-parameter families of r-dimensional orbits inside lower strata of the fold variety. Each family
  contributes about q − 1 orbits, and **𝔽₇ gives 792** (per type: 24, 70, 235). Conjecture 3\*ˡ would then be false.
  v2.23's generic-point test could not have seen this: it walks into the top stratum, where everything checked was a
  single orbit.

(Any other continuous family also grows with q. Any other Galois effect would show up as a dependence on q mod n.)

---

## 3. Point counts

| quiver | N(3) | N(5) | N(7) | N(11) | formula |
|---|---|---|---|---|---|
| 1⇉2, 1→3 | 378 | 1190 | 2682 | 8474 | 4q³ + 25q² + 10q + 15 (four primes) |
| 1⇉2, 1→3, 4→2 | 1932 | 6410 | 14928 | — | 26q³ + 115q² + 45q + 60 (three primes) |
| 1⇉2→3→4 | 1057 | 2749 | 5257 | — | 102q² + 30q + 49 (three primes) |
| 1→2⇉3→4 | 696 | 1742 | 3276 | — | 61q² + 35q + 42 (three primes) |
| 1⇉2, 2→3, 2→4 | 1280 | 3912 | — | — | open: v2.22's quadratic guess 126q² + 28q + 62 assumed no families; the 14 families add a q³ term |
| 1⇛2→3 | 5439 | 54687 | — | — | open |
| 1⇛2, 1→3 | 10410 | 128640 | — | — | open |
| 1⇉2⇉3 | 16376 | 251550 | — | — | open |

---

## 4. Proposed amendments (v2.23 → v2.24)

| amendment | content |
|---|---|
| **Q1²⁴ Local class torus (proved)** | Lemma 2 and Lemma 3: W(P) is a linear space of weightings, constant on G-orbits, and its torus part moves P to folds |
| **Q2²⁴ Conjecture 3\*ˡ** | Finitely many classes of folds under Aut_Δ and local class-torus moves. This replaces 3\* (global T_Q). 1⇉2, 2→3, 2→4 has 456 orbits under G × T_Q over 𝔽₅ against 428 over 𝔽₃, but 414 under the local torus at both, which refutes 3\* in its global form. Open for 1⇉2⇉3 (𝔽₇ decides) |
| **Q3²⁴ Families** | Families of G-orbits are recognized exactly by r_T > r_G at a lifted member, with r_T computed from W(P); counts over 𝔽₃ alone cannot see them |

---

## 5. Pending (user's machine)

1. **Done.** `refinish.py` on the batch checkpoints (now committed as `data/checkpoints.tar.gz`): 1⇉2, 2→3, 2→4 / 𝔽₅ gives
   orbits_GxTloc = **414** as predicted, and all others keep their grading counts. 1⇉2⇉3 / 𝔽₅ gives **754, not 716**
   (§2a). All ten re-finish here to identical JSON (`data/README.md`).
2. **1⇉2⇉3 over 𝔽₇** (`python3 run_batch.py 11`): 716 (Galois) or 792 (families; 3\*ˡ fails).
3. `analysis/orbitcheck.py` on 1⇉2→3⇉4 and 1⇉2→3→4→5 over 𝔽₃, and `analysis/orbitcheck.py pair_tgtfork 3 1` to finish
   the generic-point test (stopped here at 2 hours).

## 6. Open questions

1. A combinatorial description of W(P): which weightings a fold admits, read off from which classes it mixes beyond ker Δ (v2.23 §1).
2. Prove Conjecture 3\*ˡ for one parallel pair.
3. Carried: unobstructed violations; Kronecker plus a loop; zigzag; dihedral; ρ.

## 7. Sources

- M. Aguiar, *On the associative analog of Lie bialgebras*, J. Algebra 244 (2001).
- I. Assem, D. Simson and A. Skowroński, *Elements of the Representation Theory of Associative Algebras*, Vol. 1, CUP, 2006.

*"When a step forks after the twins, the branches may keep separate clocks, except where a split binds them; then one clock serves both."*
