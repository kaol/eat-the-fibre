# Discordianism: Open Questions Explored (v2.72): Conjecture G for every q

*v2.71 left p = 2, expecting it to need a new formulation: vertex blocks have entries in {0, 1, −1}, and 1 ≡ −1
mod 2. It does not.
- **No block collides mod 2.** On none of the 17 quivers do two vertex blocks become equal mod 2.
- **The census at q = 2 matches.** It was never run before on these quivers. It equals the char 0 geometric orbit
  count on every one of them.
- **Lemma L at p = 2,** by the p = 3 machinery of v2.71. Characteristic 2 stresses it much harder than 3:
  - the char 0 trees fail to specialise at 2 on 213 blocks (of the 1307 with a tree);
  - on those blocks the char-2 trees have up to 500 nodes;
  - most of these nodes carry no orbit (Lemma N).
  The certified nodes reproduce the char 0 count on every block.

**Theorem G.** For every prime power q and each of the 17 quivers (and their opposites), the 𝔽_q local-torus census
count equals the number of geometric orbits, over all blocks except the three case-split blocks. Conjecture G is
proved for these quivers, up to those three blocks.*

---

## 0. Results

The predictions are in logs/v2_72/PREDICTIONS.md, committed before the runs.

| # | Statement | Status |
|---|---|---|
| 109 | On K₂, K₃, Q₀, 1⇛2⇇3 + 1→3 (no −1 entries), census(q = 2) = char 0 count | **held**: 20, 34, 366, 752 |
| 110 | On the other 13 quivers, census(q = 2) differs from the char 0 count and is smaller | **failed**: equal on all 13. No two vertex blocks of any quiver agree mod 2 (the −1 entries never sit where a block with +1 also exists), so nothing merges |
| 111 | census(q = 2) = number of geometric orbits in characteristic 2, by Lemma L at p = 2 wherever it finishes | **held** on every block. 1211 blocks go through the char 0 trees, with 3094 direct node checks at 2. The other 213 go through char-2 trees: 6905 nodes, 3940 certified and 2965 with no orbit, matching the char 0 count on every block |

**The census at q = 2** (`pcensus.py Q 2`, local-torus classes), against the char 0 count (v2.67):

| quiver | K₂ | K₃ | 1⇉2→3 | 1→2⇉3 | 1⇉2, 1→3 | 1⇉2, 3→2 | 1⇛2→3 | 1⇛2, 1→3 | 1⇉2→3→4 | 1→2⇉3→4 |
|---|---|---|---|---|---|---|---|---|---|---|
| q = 2 | 20 | 34 | 71 | 71 | 101 | 101 | 120 | 183 | 334 | 257 |

| quiver | 1⇉2, 2→3, 2→4 | 1⇉2, 1→3, 4→2 | Q₀ | 1⇛2⇇3 + 1→3 | Q₁ | Q₂ | 1⇉2→3→4→5 |
|---|---|---|---|---|---|---|---|
| q = 2 | 414 | 485 | 366 | 752 | 1220 | 1773 | 1770 |

Every entry equals the char 0 count, including the case-split blocks' 17, 15 and 3.

## 1. Lemma L at p = 2

The same three steps as v2.71, with no prime exempt (`analysis/v2_72_all.py`, `v2_72_tally.py`):
1. **Node certificates** (v2_71_gcert.py with GCERT_EXEMPT empty): every residual prime, 2 included, is checked
   directly at a point where it comes only from the G1 or S leading coefficients.
2. **C1–C4 directly over 𝔽_2** on the trees whose ℤ claims list 2.
3. **The char-2 route** (v2_71_charp_tree.py, v2_71_charp_cert.py) for the blocks where C1–C4 fail at 2.

**Two changes to the char-p route, both found on characteristic-2 failures.**
- **Separation is not required.** v2.71 asked a certified node's point to avoid every other node of the same
  dimension. That is redundant on a char-p tree:
  - its nodes are distinct 𝔽_p-primes;
  - G1 makes each certified node the closure of its orbit;
  - so distinct certified nodes have distinct orbits.
  (On the char 0 route separation is needed, because distinct ℚ-primes can reduce to the same ideal mod p.) It
  mattered on K₃: a point in the dense orbit of one node also lay on a non-H-stable node of the same dimension.
- **Points in positive-dimensional fibres.** Over 𝔽_2, every assignment of the independent coordinates can land in
  an empty or positive-dimensional fibre. This happened at a 4-dimensional, degree-10 node of 1⇛2→3: its 16 fibres
  over 𝔽_2 have dimension −1, 1 or 2, never 0. The point search now keeps fixing coordinates until the fibre is
  finite.

Neither change affects v2.71. All its separation checks passed, and all its point searches found points.

**The char-2 route** (`logs/v2_72/charp/<quiver>/summary_p2.json`). It covers the blocks whose char 0 tree fails
C1–C4 directly at 2:

| quiver | blocks | char-2 nodes | certified | without orbit (Lemma N) | blocks = char 0 count |
|---|---|---|---|---|---|
| K₃ | 6 | 66 | 34 | 32 | 6 of 6 |
| 1⇛2→3 | 21 | 238 | 118 | 120 | 21 of 21 |
| 1⇛2, 1→3 | 27 | 385 | 183 | 202 | 27 of 27 |
| pair_before, pair_after | 1, 1 | 7, 7 | 7, 7 | 0, 0 | both |
| 1⇉2→3→4, 1→2⇉3→4 | 2, 3 | 19, 27 | 19, 27 | 0, 0 | all |
| 1⇉2, 2→3, 2→4 | 1 | 53 | 32 | 21 | 1 of 1 |
| Q₀ | 12 | 313 | 254 | 59 | 12 of 12 |
| 1⇛2⇇3 + 1→3 | 21 | 2808 | 752 | 2056 | 21 of 21 |
| Q₁ | 40 | 1177 | 978 | 199 | 40 of 40 |
| Q₂ | 70 | 1726 | 1450 | 276 | 70 of 70 |
| 1⇉2→3→4→5 | 8 | 79 | 79 | 0 | 8 of 8 |
| **total** | **213** | **6905** | **3940** | **2965** | **213 of 213** |

No node is inconclusive and no hard failure occurs anywhere (the tally: `v2_72_tally.py 2`).

**1⇛2⇇3 + 1→3 is the extreme case.**
- Every one of its 21 blocks needs a char-2 tree, 2808 nodes in all.
- The certified nodes number 752, the whole census of the quiver.
- The largest tree, block 000,000,000, has 500 nodes: 106 certified and 394 without orbit.

**What makes characteristic 2 hard.** The closures are only ξ-stable, so the char-2 trees fill up with components
that are not H-stable. Lemma N disposes of every one of them: each has a point where the orbit rank exceeds its
dimension. What remains is exactly the char 0 stratification: the certified nodes match the char 0 tree count block
by block, and so does the q = 2 census.

## 2. The theorem

**Theorem G.** Let Q be one of the 17 quivers (or its opposite) and q any prime power. Then the local-torus census
count of Folds(𝔽_q Q), over the vertex blocks other than the case-split block, equals the number of geometric orbits
of the local census groups. In particular:
- K₂, K₃, the eight pair quivers, the two 1⇛2 quivers, Q₀ and 1⇛2⇇3 + 1→3 have census counts 20, 34, 71, 71, 101,
  101, 120, 183, 334, 257, 414, 485, 366 and 752 for **every** q;
- Q₁, Q₂ and 1⇉2→3→4→5 have 1203, 1758 and 1767 on their other blocks for every q.

*Proof.*
- **p ≥ 5:** v2.70, with the v2.71 rerun.
- **p = 3:** v2.71.
- **p = 2:** §1, by the same three routes.

Vertex blocks are the same list in every characteristic (v2.66), and pairwise distinct mod 2 (checked on all 17
quivers). ∎

**What is still open for these quivers.** The three case-split blocks (v2.57). The census agrees on them at q = 2
and 3 (17, 15, 3), and v2.67 proved the char 0 count, but no certificate covers them at any p yet.

## 3. Next

1. **The case-split blocks.** Specialise the v2.57 pivot case split, whose pivots are ±1, and then the leaf trees,
   by the same chain.
2. **New quivers** beyond the census list: the chain v2.66 → v2.67 → v2.70–72 is now automatic and gives Conjecture G
   at every q.
3. **Conjecture S, and Program 3.**

## Files

- **analysis/:**
  - v2_72_all.py (one quiver at a prime P: node certificates with no exempt prime, C1–C4 directly at P);
  - v2_72_tally.py.
- **Changed:**
  - v2_71_gcert.py (an empty exempt set is allowed);
  - v2_71_charp_cert.py (separation recorded, not required; points in positive-dimensional fibres);
  - v2_71_charp_all.py (tally and output paths configurable).
- **logs/v2_72/:**
  - PREDICTIONS.md;
  - census_<quiver>_2.{json,txt};
  - <quiver>/ (node certificates at 2, *.u2.json);
  - charp/<quiver>/ (char-2 trees, certificates, summary_p2.json);
  - tally_p2.json.
