# Discordianism: Open Questions Explored (v2.71): Conjecture G for every odd q

*v2.70 proved Conjecture G (census count = geometric orbit count) for every q = p^k with p ≥ 5. This round does
two things.
- **An erratum to v2.70, and the rerun that repairs it.** The separation primes had been combined unsoundly. All
  7920 node certificates are recomputed with the sound combination; Theorem G₅ stands.
- **p = 3.** Lemma L is now proved at p = 3, by two routes.
  - On 1390 of the 1424 blocks outside the case splits, the char 0 trees specialise to 3.
  - On the other 34, they do not. There, new characteristic 3 trees carry the argument. They need one new idea: a
    node with no orbit in it.

**Theorem G_odd.** For every odd prime power q and each of the 17 quivers (and their opposites), the 𝔽_q
local-torus census count equals the number of geometric orbits, over all blocks except the three case-split blocks.
Only p = 2 remains open.*

---

## 0. Results

The predictions are in logs/v2_71/PREDICTIONS.md, committed (with the erratum) before the runs.

| # | Statement | Status |
|---|---|---|
| 105 | The sound combination of separation primes changes no conclusion of v2.70 | **held**: every node is certified for every p ≥ 5; no hard failure |
| 106 | At p = 3, every dense and point node of the char 0 trees is certified (directly where needed) | **held** on the blocks whose trees specialise. 1925 node checks were done directly in characteristic 3. On the 34 blocks whose trees do not specialise, the char 0 nodes are replaced by char-3 trees (§3) |
| 107 | The direct C1–C4 check at p = 3 passes on all but a few of the 34 trees that list 3; it fails on 1⇛2⇇3 + 1→3 block 100,001,001 | **failed as stated**: it fails on *all* 34. In each, some radical membership C3 does not hold mod 3: the boundary has extra components there |
| 108 | Conjecture G at every q = 3^k for every quiver whose trees pass 107 | **held, stronger**: for all 17 quivers. All 34 blocks are certified in characteristic 3, each with exactly its char 0 count |

## 1. Erratum to v2.70, and the rerun

**The flaw.** v2_70_gcert.py computed:
- the certificate candidates per point, intersected over points;
- separately, the separation primes, also intersected over points;
- and then took the union of the two.

**Why it is unsound.** A prime that fails one point's certificate and another point's separation survives neither
intersection, yet no single point certifies it. The sound set is ∩_x (cand(x) ∪ sep(x)).

**`analysis/v2_71_gcert.py` fixes it**, with three further changes:
- **Separation per point.** At each point, separation is computed against every other node of the same dimension,
  by evaluating the node ideals in the polynomial ring (fast). It is added to that point's candidates.
- **Only 2 is exempt.** Sampling continues while some residual prime q ≥ 3 has no point at which q comes from the
  G1 or S leading coefficients alone. Points whose coordinates have denominators divisible by such a q are skipped
  before the expensive work. The cap is 40 points.
- **Direct checks at every residual prime q ≥ 3**, as in v2.70 §2 (G1 with a slice chosen mod q; the stabiliser's
  minimal primes over 𝔽_q).

All 7920 nodes were recomputed. **No prime ≥ 3 remains at any node.** 3084 direct checks were made, at 15 different
primes: 3 (1925), 5 (670), 7 (295), 11 (89), 13 (42), and a few each up to 59. Each check is made at a prime prime
to the block's C1 leading coefficients, as the tally verifies.

## 2. p = 3 on the trees that specialise

Lemma L at p = 3 needs the node certificates at 3 (§1) and the tree claims C1–C5 at 3.
- **C5** (coverage) has no candidate at all.
- **C1–C4** list 3 on 34 trees, through Rabinowitsch integers divisible by 3. On the others the ℤ claims hold at 3.

So on the other 1390 blocks (the 1424 blocks outside the case splits, 117 of them point strata, minus these 34),
Lemma L holds at p = 3 with the char 0 tree.

On the 34 trees, `v2_68_uniform.py --p 3` checks C1–C4 directly over 𝔽_3. It fails on every one: some generator of
the children's intersection is not in the radical of J_Y mod 3. The boundaries acquire extra components in
characteristic 3, as v2.68 found for one of them.

## 3. p = 3 on the other 34 blocks: char-3 trees, and nodes without orbits

**The char-3 trees** (`analysis/v2_71_charp_tree.py`, generalising v2_68_charp.py). Roots are the minimal primes
of the stratum ideal over 𝔽_3. Each root Z gets its own group H_Z = G × T_{W_Z}:
- W_Z is the space of weightings compatible on Z, computed over 𝔽_3 and checked to be a partition span.
- v2.68 assumed one W for all roots, which failed on a Q₀ block whose roots differ.
- Stability under H_Z holds by the v2.68 argument. Z is a component of the w-compatible locus of the stratum, and
  that locus is stable under the connected torus T_w.

Children are the minimal primes of the closure of I_Z + (minor) under the fields, over 𝔽_3. Lemma C (v2.68) then
says that every H_Z-orbit is a top orbit *contained in* some node, with no coverage or specialisation claims needed.

**Lemma L at a char-p node** (`analysis/v2_71_charp_cert.py`). For a node Z of dimension e, search 𝔽_3-points x
of Z. A point x **certifies** Z if all of the following hold over 𝔽_3:
- W_x (the weightings compatible at x) is a partition span; set H_x = G × T_{W_x} ⊇ H_Z;
- the orbit matrix of H_x at x has rank e;
- G1: the graph of a slice of H_x (pivots mod 3) eliminates to I_Z. So Z is the closure of H_x·x, and is
  geometrically irreducible;
- S: the stabiliser has one minimal prime, of dimension dim H_x − e, and is smooth at 1;
- x lies on no other node of dimension e.

*Then O_Z := H_x·x is a single census class over every 𝔽_{3^k}.* H_Z is normal in H_x. Every H_Z-orbit inside
H_x·x is a translate of the open one, so H_x·x = O_Z. The census moves at points of O_Z lie in H_x, and Lang
applies.

**A new phenomenon: nodes without orbits.** In characteristic p the closures are only ξ-stable, so a node need not
be H-stable. Some char-3 nodes contain no orbit of their own dimension. They are detected by:

*Lemma N.* If the orbit matrix has rank > e = dim Z at some point of Z, then no orbit of dimension e lies in Z.

*Proof.*
1. The rank is lower semicontinuous, so the generic rank on Z is > e.
2. Suppose a dense orbit O ⊆ Z existed. It would meet the open set where the rank exceeds e.
3. But at every point of O, the rank is at most dim O = e, since the image of Lie(H) lies in the tangent space of
   the orbit. Contradiction. ∎

Such a node contributes no class. Its orbits of lower dimension lie in its children, by Lemma C. A 0-dimensional node
is handled the same way:
- it contributes a class when its point is rational, fixed by G (an identity in g over the residue field) and fixed
  by every weighting move compatible at it;
- it contributes none when G moves it.

**Results** (`logs/v2_71/charp/<quiver>/summary_p3.json`):

| quiver | blocks | char-3 nodes | certified | without orbit | blocks = char 0 count |
|---|---|---|---|---|---|
| 1⇛2→3 (K3_after) | 2 | 30 | 26 | 4 | 2 of 2 |
| 1⇉2, 2→3, 2→4 (pair_tgtfork) | 1 | 32 | 32 | 0 | 1 of 1 |
| Q₀ | 2 | 93 | 87 | 6 | 2 of 2 |
| 1⇛2⇇3 + 1→3 | 18 | 836 | 679 | 157 | 18 of 18 |
| Q₁ | 8 | 377 | 357 | 20 | 8 of 8 |
| Q₂ | 3 | 154 | 142 | 12 | 3 of 3 |
| **total** | **34** | **1522** | **1323** | **199** | **34 of 34** |

Every node is either certified or without orbit; none is
inconclusive.

The extra char-3 nodes all turn out to be orbit-free, and the certified ones reproduce the char 0 count exactly.
Examples:
- triple block 000,000,000: 156 char-3 nodes = 106 certified + 50 without orbit, and its char 0 count is 106;
- 100,010,001: 125 = 85 + 40;
- Q₁ block 1000,0100,0010,0001: 122 = 102 + 20.

An example of how a node can contain no orbit: K3_after 000,010,001 has a char-3 node Z2.1.2 that is a linear 4-space.
Its 81 𝔽_3-points have orbit rank 5 (54 points), 4 (18) and 3 (9). The rank-5 points prove there is no 4-dimensional
orbit inside it. The rank-4 points lie in the orbits of other nodes (they fail separation).

**Two bugs met on the way.** Both are in the new code and both are fixed:
- an 𝔽_3-point search that was too shallow (one node needed a second batch of points);
- a vdim computed in a ring that still contained the group variables, which made a rational point look irrational.

## 4. The theorem

**Theorem G_odd.** Let Q be one of the 17 quivers of v2.70 (or its opposite) and q a power of an odd prime p. Then
the local-torus census count of Folds(𝔽_q Q), over the vertex blocks other than the case-split block, equals the
number of geometric orbits of the local census groups (the v2.67 count).
- For K₂, K₃, the eight pair quivers, the two 1⇛2 quivers, Q₀ and 1⇛2⇇3 + 1→3, this is the whole census:
  20, 34, 71, 71, 101, 101, 120, 183, 334, 257, 414, 485, 366, 752, for every odd q.

*Proof.* There are three cases.
- **p ≥ 5:** v2.70 with the §1 rerun.
- **p = 3, on the 1390 blocks of §2:** Lemma L, with the char 0 tree, the ℤ claims (3 not a candidate), and the node
  certificates of §1 at 3.
- **p = 3, on the 34 blocks of §3:** the char-3 tree and Lemma L at its nodes, with Lemma N for the nodes without
  orbits.

Vertex blocks are the same list in every characteristic (v2.66), and distinct for odd p. ∎

**p = 2 is genuinely different.** The vertex block entries 1 and −1 coincide mod 2, so blocks merge, and the
statement itself has to be reformulated before it can be checked. The census was never run at q = 2 or 4 on these
quivers.

## 5. Next

1. **p = 2:** find the right statement first (which blocks merge, and what the census counts there), then the
   char-2 trees (v2.68 has one, with 326 nodes for the triple's block 100,001,001).
2. **The three case-split blocks:** specialise the v2.57 pivot case split.
3. **New quivers**, through the now automatic chain: v2.66 → v2.67 → v2.70/71.
4. **Conjecture S, and Program 3.**

## Files

- **analysis/:**
  - v2_71_gcert.py (node certificates, separation per point, residual primes ≥ 3 checked directly);
  - v2_71_all.py (one quiver: gcert on every block, C1–C4 directly at 3 where listed);
  - v2_71_charp_tree.py (char-p trees, one group per root);
  - v2_71_charp_cert.py (Lemma L at char-p nodes, Lemma N);
  - v2_71_charp_all.py (the char-3 route for the blocks failing at 3);
  - v2_71_tally.py.
- **logs/v2_71/:**
  - PREDICTIONS.md;
  - <quiver>/<block>.gc.json, *.u3.json;
  - charp/<quiver>/ (trees, certificates, summary_p3.json);
  - tally.json.
