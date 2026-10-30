# Discordianism: Open Questions Explored (v2.73): the case-split blocks; Conjecture G complete on 17 quivers

*v2.70–72 proved Conjecture G (census count = geometric orbit count) at every q, except on three vertex blocks whose
strata needed the v2.57 case split: Q₁ 0000,0000,0000,0001, Q₂ 1000,010m,010m,0000 and 1⇉2→3→4→5
00000,00000,0010m,00000,00000. This round closes them. The leaf decomposition gets its own ℤ-certificates (A–E
below). Where those fail (only at p = 2) the round falls back on characteristic-2 trees inside the leaf charts, and,
for Q₁'s block, on the integral case split v2.59 had already computed.

**Theorem.** For each of the 17 quivers (and their opposites) and every prime power q, the 𝔽_q local-torus census
count equals the number of geometric orbits:
- K₂, K₃: 20, 34;
- the pair quivers: 71, 71, 101, 101, 334, 257, 414, 485;
- the 1⇛2 quivers: 120, 183;
- Q₀: 366;
- 1⇛2⇇3 + 1→3: 752;
- Q₁: 1220;
- Q₂: 1773;
- 1⇉2→3→4→5: 1770.*

---

## 0. Results

The predictions are in logs/v2_73/PREDICTIONS.md, committed before the runs.

| # | Statement | Status |
|---|---|---|
| 112 | The case split is valid at every p for Q₂'s and 1⇉2→3→4→5's blocks (constants ±1); for Q₁'s block (v2.57 run) at every p ≠ 2 | **held**. v2.59's integral case split of Q₁'s block (constants ±1) then covers p = 2 |
| 113 | B, C, D, E have candidates in {2, 3} after direct checks | **held, sharper**: only 2 |
| 114 | Conjecture G on these blocks for all p ≥ 5, and p = 3 where the direct checks pass; counts 17, 15, 3 | **held, stronger**: at every p including 3 (no candidate) and 2 (§3) |

| block | leaves | H-stable roots | count | candidate primes (ℤ-certificates) | p = 2 |
|---|---|---|---|---|---|
| 1⇉2→3→4→5 00000,00000,0010m,00000,00000 | 5 | 2 (+1 fixed point) | 3 | **none** | covered by ℤ |
| Q₂ 1000,010m,010m,0000 | 34 | 11 roots (+1 fixed point) | 15 | 2 (C1–C4 on three trees) | char-2 trees: 15 |
| Q₁ 0000,0000,0000,0001 (v2.57 split) | 54 | 14 (+1 fixed point) | 17 | 2 (branch constants ±2, C1–C4 on eight trees) | integral split + char-2 trees: 17 |

## 1. The certificates for a case-split block

The v2.57 case split is a decision tree. At each node it either:
- substitutes at a linear pivot,
- replaces an equation by its radical, or
- splits f₁⋯f_k = 0 into disjoint branches (f_i = 0 with f_1, …, f_{i−1} ≠ 0).

A leaf is a chart, with:
- coordinates vs;
- every other stratum coordinate a polynomial in vs (the composite substitutions);
- remaining equations rem and inequations nz.

The stratum is the union of the leaves' points. v2.67 showed, in characteristic 0, that the stratum's components are
the H-stable leaf primes, and counted orbits on the union of their trees, identified by full-coordinate key.

At a prime p, Lemma L for the block needs (`analysis/v2_73_leafcheck.py`, `v2_73_split.py`, `v2_73_combine.py`):

| | claim | certificate |
|---|---|---|
| A | the case split covers S_p | every recorded branch constant (contents, factor denominators, pivots, emptying constants) prime to p |
| B | a leaf's points lie in its leaf prime mod p | per generator ρ of the leaf prime, an integer in (rem) + (1 − T·∏nz·ρ) over ℤ |
| C | a leaf prime that is not H-stable lies in an H-stable root mod p | the root's full-coordinate key, pulled back to the leaf's chart, reduces to 0 modulo a strong ℤ-basis of the leaf prime |
| D | the roots' trees specialise | C1–C4 (v2_68_uniform) and the node certificates (v2_71_gcert, no exempt prime), in the root's leaf chart |
| E | nodes of different trees are distinct mod p | same key = same node; different keys of equal dimension separated at a point, in full coordinates, combined per point with the node's other candidates |

**Why this suffices.**
1. **The roots cover.** By A, B and C, S_p is covered by the reductions R_p of the H-stable roots, each H_p-stable
   (the coaction argument over ℤ).
2. **Every orbit sits in a tree.** An orbit lies in some R_p. By Lemma C, applied in R's chart (the chart's fields
   are the true fields' chart components), it is the top orbit of a node of R's tree.
3. **The count.** D and E then give one class per distinct node.
4. **Fixed point leaves are roots too.** They are isolated 0-dimensional components. Their points must be fixed by
   G and by every compatible weighting move, which is checked as an identity over ℤ.

**Results** (`logs/v2_73/<block>/combined.json`).
- **Checks A–C:** every check in B and C is exact over ℤ, with no candidate prime, on all 93 leaves of the three
  blocks (5 + 34 + 54).
- **Candidates:** 2 is the only prime, and it enters only in two ways:
  - from C1–C4 (Rabinowitsch integers);
  - from v2.57's branch constants on Q₁'s block (claim A, on all 54 leaves).
- **Counts:** 3, 15 and 17, the v2.67 counts.

## 2. Two bugs in the first pass, both fixed

- **A root lost.** A root whose leaf prime is the zero ideal of its chart has the empty key ''. The driver's
  `if rk:` dropped it, so two 1-dimensional roots of Q₂'s block were missing. That showed up as leaves that no root
  contained, and a count of 13.
- **A point counted twice.** A fixed point leaf and a tree's point node can be the same point under keys of
  different formats (count 16). Points are now merged by their coordinates.

## 3. p = 2

**Q₂'s block.** The three trees whose C1–C4 list 2 also fail them directly over 𝔽₂ (radical membership), as most
blocks did in v2.72. So, as there, characteristic-2 trees, this time inside leaf charts (`analysis/v2_73_split_p.py`):
- **Roots:** the minimal primes over 𝔽₂ of each H-stable root's leaf prime, in its chart. The chart is a closed
  graph and the root's key contains the graph equations, so R₂ lies in it.
- **Nodes:** Lemma L at each node, with v2_71_charp_cert and Lemma N.
- **Deduplication** (`analysis/v2_73_dedupe_p.py`), across the charts in full coordinates mod 2. Two certified nodes
  are one when each one's point lies on the other's closure: on the other chart's graph, with the node's
  generators vanishing.

The result: 34 certified nodes in 11 charts collapse to 15 orbit closures, including the fixed point. That is the
count.

**Q₁'s block.** v2.57's split divides by ±2, which is not a valid step mod 2. v2.59 had already redone this block
in integral mode, with pivots ±1, primitive normalisation and every recorded constant ±1
(`logs/v2_59/certs_q1_branch_integral.jsonl`). On that certificate:
- the leaf union gives 54 leaves and 17 orbit closures (`logs/v2_73/q1i_leaf`);
- A, B and C hold exactly, with no candidate;
- the characteristic-2 trees in its 14 root charts have 42 certified nodes, deduplicating to 16 orbit closures plus
  the fixed point, 17 in all.

## 4. The theorem

**Theorem G (complete).** Let Q be one of the 17 quivers of v2.66–67, or its opposite, and q any prime power. Then
the local-torus census count of Folds(𝔽_q Q) equals the number of geometric orbits of the local census groups.
These are the counts listed above.

*Proof.* Lemma L on every vertex block at every prime.
- **Blocks without a case split:** v2.70 (p ≥ 5, with the v2.71 rerun), v2.71 (p = 3) and v2.72 (p = 2).
- **The three case-split blocks:** §1 at every p with no candidate prime; §3 at p = 2.

Vertex blocks are the same list in every characteristic (v2.66) and pairwise distinct for every p (v2.72). ∎

**Where this leaves Conjecture G.** It is proved for every unobstructed quiver of the census list. The proof is
automatic: v2.66 (vertex blocks, components), then v2.67 (orbit trees), then v2.70–73 (certificates, char-p trees,
leaf charts). So the next step is to run the whole chain on quivers beyond the list.

## 5. Next

1. **New quivers.** The census list's natural continuations: K₄; 1⇛2→3→4; the K₃ variants in pcensus.py; Q₀ with
   two tails. The chain is now automatic end to end.
2. **Conjecture S** (vertex blocks = folds of Q^s). It is verified on 24 quivers and still unproved.
3. **Program 3** (the leaf extension Q₀ → Q₁ as a relative orbit problem). The case-split structure of this round
   is its first worked example.

## Files

- **analysis/:**
  - v2_73_leafcheck.py (A, B, C, point leaves, per leaf);
  - v2_73_split.py (one block: roots, trees, all ℤ-certificates);
  - v2_73_fullpts.py (full coordinates of certified points);
  - v2_73_combine.py (E, the block's primes and count);
  - v2_73_split_p.py (char-p trees in leaf charts);
  - v2_73_dedupe_p.py (char-p deduplication across charts, and the count).
- **Changed:**
  - v2_67_orbits.py (point leaves keep their setup and borrow the group data, under NO_POINT_EXIT);
  - v2_68_uniform.py and v2_71_charp_tree.py (leaf charts; a given root ideal).
- **logs/v2_73/:**
  - PREDICTIONS.md;
  - pair_after3/, q2/, q1/ (v2.57 split), q1i/ and q1i_leaf/ (integral split): the leaf checks, trees, certificates,
    combined.json, char-2 trees and dedupe_p2.json.
