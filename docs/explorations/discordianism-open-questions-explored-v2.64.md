# Discordianism: Open Questions Explored (v2.64): a second leaf on Q₀, and what it falsifies

*v2.60–v2.63 built a picture of the leaf extension Q₀ → Q₁ = Q₀ + 2→4. This round tests that picture on a second
extension of the same base: Q₂ = Q₀ + 3→4 (1⇉2⇇3, 1→3, 3→4).
- **New theorem (char 0).** Q₂ has finitely many census orbits. The direct pipeline certifies all 99 vertex blocks.
- **The relative criterion carries over.** Over Q₀'s 366 classes, every positive-dimensional fibre of the J- and
  kQ₀-preserving kinds passes.
- **Two parts of the Q₁ picture do not carry over:**
  1. The mixed kind is no longer confined to cornerless blocks, and it is mostly not rank one: 43 mixed components,
     37 of rank > 1, 6 of them in blocks whose corner is a Q₀ block.
  2. The kQ₀-preserving λ = 0 lifts form positive-dimensional families: 83 components of dimension 1–3, the Lemma B
     idempotent families of v2.63. Stab(P̄) still has finitely many orbits on them.

So a general leaf lemma must handle the mixed kind and Lemma B through the stabiliser, not by fibre finiteness.*

---

## 0. Results

The predictions are in logs/v2_64/PREDICTIONS.md, committed before any computation on Q₂.

| # | Statement | Status |
|---|---|---|
| 77 | The census count of Q₂ is the same at q = 3, 5, 7 | **open**: 1773 at q = 3 (8346 s on one core). q = 5 and 7 are queued as GPU jobs 18 and 19 in run_batch.py |
| 78 | The direct char 0 pipeline certifies every block of Q₂ | **held**: 99 of 99 blocks, 594 components or leaf primes, no failure. One block (1000,010m,010m,0000) ran out of memory with both the default and the Markowitz presolve, and needed the integral case split (34 leaves) |
| 79 | As for Q₁: mixed components only in cornerless blocks, all of rank ≤ 1 | **failed**: 43 mixed components in 12 blocks; 37 have rank > 1; 6 of the 12 blocks have a Q₀ corner (for example 0000,0100,0010,0001) |
| 80 | The relative criterion passes over the 366 classes, and every λ = 0 fibre is finite | **half held**: every positive-dimensional fibre component passes (205 components). J-preserving λ = 0 fibres are finite (Lemma A). But kQ₀-preserving λ = 0 fibres include 77 components of dimension 1, 5 of dimension 2 and 1 of dimension 3 |

## 1. Q₂ in characteristic 0

- **Vertex blocks** (`analysis/v2_58_vertex_integral.py pairs_in_short_c4 … --arrow 1,3,4,5`, `vertex_char0.txt`).
  Eliminating 688 equations in 88 unknowns gives a zero-dimensional ideal with 99 rational points, the vertex blocks
  (vdim 119).
- **Pipeline** (`pipeline_aa.txt`, `pipeline_ab.txt`, `pipeline_rerun.txt`, `tally.txt`). Every component of every
  block passes the v2.54 criterion, with certificates in `certs_*.jsonl.xz`.

**Theorem.** Folds(k̄Q₂) is a finite union of census orbits, for k̄ of characteristic 0.

**Census.** 1773 local-torus classes at q = 3. Q₂ has the same kernel algebra as Q₁: K⁺ = K₀ × k e₄, with 32904
kernel folds at q = 3 in both, but more classes (1773 against 1220).

## 2. The four kinds on Q₂ (`analysis/v2_63_kinds.py` with LEAF_QUIVER=pairs_in_short_c4, `kinds.txt`)

| kind on the whole component | Q₁ = Q₀ + 2→4 | Q₂ = Q₀ + 3→4 |
|---|---|---|
| split | 234 | 195 |
| J-preserving only | 98 | 135 |
| kQ₀-preserving only | 88 | 221 |
| mixed, rank ≤ 1 | 6 (6 blocks) | 6 |
| mixed, rank > 1 | 0 | 37 |
| blocks with mixed components | 6, all cornerless | 12, of which 6 have a Q₀ corner |

The mixed components are concentrated in blocks where row 4 of the vertex block is (0, −1, 0, 1), 0m01 in the
block notation: P(e₂) and P(e₄) share e₄-components. This is vertex-level mixing between the leaf and vertex 2, which
is not where the leaf is attached; the leaf hangs at vertex 3. One mixed component sits in a block with no vertex-level
mixing at all (0000,0100,0010,0001), so the mixing there happens on paths.

## 3. The relative criterion for Q₀ → Q₂ (`analysis/v2_62_relative.py` with LEAF_QUIVER=pairs_in_short_c4, `relative.txt`)

| kind, λ | one point | finite | positive-dimensional components (all pass) |
|---|---|---|---|
| J, 0 | 160 | 477 points | none (Lemma A) |
| J, 1 | 129 | 238 points (114 empty) | none |
| Q, 0 | 220 | 241 points | dim 1: 77, dim 2: 5, dim 3: 1 |
| Q, 1 | 77 | 290 points (114 empty) | dim 1: 120, dim 2: 2 |

The kQ₀-preserving λ = 0 families are the idempotent families of Lemma B made real. In Q₂, P(e₄) = −ε with ε an
idempotent of im P̄ ∩ K₀, and ε can move along e_i + t(a − a′)-type lines. On each such family the census group's
stabiliser still has finitely many orbits: substitutions within a class rescale the differences a − a′.

## 4. What this changes for Program 3

**Survives on both extensions:**
- the four kinds and the compression lemma;
- the trivial lift;
- Lemma A (J-preserving λ = 0 is finite: characters);
- Lemma C (rank-one structure);
- the relative criterion over the base classes, with every positive-dimensional fibre component passing.

**Does not survive:**
- the confinement of the mixed kind to cornerless blocks;
- "mixed means rank one";
- finiteness of the kQ₀-preserving λ = 0 fibres.

**Smallest counterexample.** The brief asks for the smallest counterexample to each candidate lemma. For "the mixed
kind of a sink-leaf extension of an unobstructed quiver is rank one", it is K₂ → K₂ + 1→3 (1⇉2, 1→3), from v2.63's
census table: 5 mixed classes of rank 2 at q = 3. For "kQ-preserving λ = 0 lifts are finite over each base fold", Q₀ → Q₂ is
the first example found, through the Lemma B families; the small extensions were not searched for this.

**Where the leaf hangs matters.** At vertex 2, the common target of both parallel classes, the leaf mixes only
through rank-one idempotent folds. At vertex 3, the source of one class and the target of the joining arrow, it mixes in rank > 1, partly through paths. A general leaf lemma will have to reason about the attaching vertex's position in the
parallel-class incidence pattern, which is what the brief's "generalization target" anticipated.

**The method still works.** For both extensions, finiteness over the base classes holds part by part, and the mixed
kind is handled directly. As a proof procedure it is now:
1. the base classes, explicit (v2.62);
2. the relative criterion on the compressible kinds;
3. the direct pipeline only on the blocks that carry mixed components.

Step 3 is still per quiver. For Q₂ the mixed components sit in 12 of the 99 blocks.

## Files

- pcensus.py: the quiver pairs_in_short_c4.
- run_batch.py: jobs 18 and 19 (q = 5, 7), predicting 1773.
- analysis/v2_62_relative.py and v2_63_kinds.py: parametrised by LEAF_QUIVER.
- logs/v2_64/:
  - PREDICTIONS.md;
  - results.jsonl, census.txt, reps/;
  - vertex_char0.txt, blocks.txt, blocks_census_format.txt;
  - pipeline_{aa,ab,rerun}.txt, certs_{aa,ab,rerun}.jsonl.xz, tally.txt;
  - kinds.txt, kinds_rerun.txt;
  - relative.txt.
