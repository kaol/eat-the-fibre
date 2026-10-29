# Discordianism: Open Questions Explored (v2.67): orbit stratification, and the census counts geometric orbits

*v2.66 left 1⇛2⇇3 + 1→3 open on one component. The v2.54 cut ran out of memory on it: 6 GB here, 51 GB over
16 hours on the user's machine. This round replaces the cut with a different computation, and that computation does
more than finish the case.
- **A new finiteness method (proved correct; `analysis/v2_67_orbits.py`).** Each step cuts with the group-invariant
  hull of one minor, and returns the actual orbit closures. On the stuck component it finishes in 9 seconds for the
  first step and under two minutes in all.
- **1⇛2⇇3 + 1→3 is finite in characteristic 0 (theorem).** Every one of its 21 vertex blocks has finitely many
  orbits.
- **The census counts geometric orbits.** On 17 quivers, 1427 vertex blocks, the number of geometric orbits in each
  block equals the census's 𝔽₃ local-torus class count in that block, block by block, 8072 classes in all. The
  census's q-independence is now explained at the level of the count itself: it counts the orbits of an algebraic
  group action in characteristic 0.*

---

## 0. Results

The predictions are in logs/v2_67/PREDICTIONS.md. 91–92 were committed before any run. 93–97 were added in four
later commits, each before the runs it predicts.

| # | Statement | Status |
|---|---|---|
| 91 | The new method agrees with the v2.54 criterion on controls | **held**: finite on every passing control component. On the obstructed 1⇉2⇉3 + 1→3, block 000,001,001, it finds no dense orbit on exactly the component v2.55 failed |
| 92 | It finishes component 1 of block 100,001,001 of 1⇛2⇇3 + 1→3, and finds finitely many orbits | **held**: 72 orbit closures for the whole block, 102 s (the component alone: 264 s before deduplication) |
| 93 | Geometric orbits = census classes, block by block, on 7 quivers | **held** on all 7 |
| 94 | The same on 6 more, including Q₀ | **held** on all 6, after two fixes (§3) |
| 95 | The same for 1⇛2⇇3 + 1→3, 21 blocks, total 752 | **held** |
| 96 | The same for Q₁ and Q₂ on their blocks without a case split | **held**: 1203 and 1758 |
| 97 | Case-split blocks: Q₁ 17, Q₂ 15; 1⇉2→3→4→5 equal everywhere, total 1770, the case-split block giving 7 | **held as stated for Q₁, Q₂ and the total; failed for the "7"**. The "7" came from a version of the leaf method that was wrong (§3). The corrected method gives 3, and the census, computed afterwards, has 3 |

## 1. Orbit stratification by invariant boundaries

**Setting.** As in v2.54, take one vertex-block stratum with coordinates vs (after presolve). On it acts H, the
connected group with Lie algebra 𝔥 = Der_Δ ⊕ (class weightings). Its vector fields in the coordinates are:
- ξ_D(y_{r,j}) = [D, P]_{r, c_j};
- ξ_D(k_{ij}) = (Kl [D, P] Kmat)_{ij}.

These are checked, not assumed. The script verifies ξ_D(I_X) ⊆ I_X on every input component; this passes on all of
them.

**Local weightings.** At an H-stable prime Y, 𝔥_Y = Der_Δ ⊕ {weightings w with Σ_c w_c · (Δ_c-defect) ≡ 0 on Y and
ξ_w tangent to Y}. This is v2.23's local torus at the generic point of Y. It grows down the tree.

**Lemma O.** Let X be an H_X-stable prime of dimension e. Let f be an e × e minor of the orbit matrix (rows ξ_D for
D ∈ 𝔥_X, columns the entries of [D, P]) with f ∉ I_X. Let J be the closure of I_X + (f) under all ξ_D. Then:
1. X has a dense orbit, on each geometric component;
2. X minus that orbit is contained in V(J);
3. the minimal primes of J are H-stable, and each has dimension below e.

*Proof.*
- **Dense orbit.** f ∉ I_X means orbit dimension e = dim X at a generic point, so the orbit is open in each geometric
  component of X.
- **Only one such orbit.** Two orbits of dimension e in one geometric component would both be dense, so they would
  meet, so they are equal. Every other point of X therefore has orbit dimension below e, and lies in V(f).
- **The rest is in V(J).** X minus the dense orbit is H-stable. So for such x, f(exp(tD)·x) ≡ 0. In characteristic 0,
  differentiating this at t = 0 shows every iterated ξ-derivative of f vanishes at x.
- **The children are H-stable.** J is ξ-stable by construction. In characteristic 0 its minimal primes are ξ-stable
  too (Seidenberg), so they are H-stable, H being connected.
- **They are smaller.** Each minimal prime has dimension below e, because f ∉ I_X. ∎

**The algorithm.** Recurse on the minimal primes of J, and identify equal primes by their reduced Gröbner basis.
- **Each node is one orbit closure.** By Lemma O, every node is the closure of exactly one orbit per geometric
  component. Every orbit is reached, by induction on dimension.
- **The count.** The number of geometric orbits is the sum, over the nodes, of their numbers of geometric components:
  - 1 when the minor was found at a rational point (a smooth rational point means geometric irreducibility);
  - otherwise counted with absPrimdecGTZ (needed once, and it gave 1);
  - vdim for 0-dimensional primes.
- **What this replaces.** It replaces "dim Z_r ≤ r for all r" (v2.54): finiteness is the statement that every node
  has a dense orbit. The difference is in what gets decomposed. v2.54 decomposed V(I_X + f), whose extra components
  pile up from step to step. Here J cuts out exactly the union of lower orbits.

**Two cases needed care, and both mattered (§3).**
- **Case-split leaves.** A leaf's coordinates come with substitutions that hold only on the leaf. So a leaf prime is
  H-stable only if the vector field also matches on the eliminated coordinates. The script now checks this, and
  skips leaf primes that fail it, since those are not components of the stratum.
- **Point leaves.** A single-fold leaf counts as its own orbit only if H fixes it. Otherwise it lies in an orbit
  counted in another leaf.

`analysis/v2_67_leafunion.py` combines the leaves of a case-split block. It identifies orbit closures across leaves
by their reduced Gröbner basis in all coordinates.

## 2. 1⇛2⇇3 + 1→3 in characteristic 0

Block 100,001,001 has component 1 of dimension 10 in 30 coordinates, with 𝔥 of dimension 10.
- **First boundary.** The first minor has 11,669 characters. Its invariant hull closes in 5 rounds and splits into
  five primes of dimension 9; this takes 9 s.
- **The whole tree.** It has 72 distinct orbit closures over the 6 components of the block, every one with a dense
  orbit, in 102 s.
- **The other 20 blocks.** The other 20 blocks run in at most a minute each.

**Theorem.** Folds(k̄(1⇛2⇇3 + 1→3)) is a finite union of census orbits for k̄ of characteristic 0. It has 752 orbit
classes of the local census groups, the census count at q = 3 and 5.

**Uniformity in q is not yet proved for this quiver.** The specialisation argument of v2.58 is written for v2.54 cut
certificates. It covers the other 20 blocks (v2.66) but not the new certificates of block 100,001,001. §4 sketches
what it needs.

## 3. The census counts geometric orbits (`analysis/v2_67_count.py`, `count_*.txt`)

For every quiver with saved census representatives, the procedure is:
1. run the orbit stratification on every block, all components together;
2. count the geometric orbits;
3. compare with the q = 3 local-torus classes of that block.

| quiver | blocks | geometric orbits = census classes |
|---|---|---|
| K₂, K₃ | 6, 6 | 20, 34 |
| 1⇉2→3, 1→2⇉3 | 23, 23 | 71, 71 |
| 1⇉2, 1→3 and 1⇉2, 3→2 | 27, 27 | 101, 101 |
| 1⇛2→3, 1⇛2, 1→3 | 23, 27 | 120, 183 |
| 1⇉2→3→4, 1→2⇉3→4 | 100, 100 | 334, 257 |
| 1⇉2, 2→3, 2→4 and 1⇉2, 1→3, 4→2 | 148, 123 | 414, 485 |
| Q₀ = 1⇉2⇇3 + 1→3 | 21 | 366 |
| 1⇛2⇇3 + 1→3 | 21 | 752 |
| Q₁ = Q₀ + 2→4 | 90 (1 case split) | 1220 (= 1203 + 17) |
| Q₂ = Q₀ + 3→4 | 99 (1 case split) | 1773 (= 1758 + 15) |
| 1⇉2→3→4→5 | 563 (1 case split) | 1770 (= 1767 + 3; census computed here, 56 s with the torch engine) |

**No mismatch in any of the 1427 blocks.**

**What it means.** The 𝔽_q census counts, flat in q for these quivers, are the number of geometric orbits of the
local census groups in characteristic 0. So the q-independence the project has observed since v2.27 is now
explained at the level of the count itself.

**What is still missing for a proof of "census = orbits for every q":**
- Each geometric orbit should be defined over the prime field.
- It should give exactly one class over 𝔽_q. This follows from Lang's theorem when stabilisers are connected and the
  local tori are handled as the census handles them.

The computation shows the equality but does not yet prove why. Stated as a conjecture:

**Conjecture G.** For an unobstructed Q and every q, the census local-torus class count equals the number of
geometric orbits of the local census groups on Folds(k̄Q), in characteristic 0.

**Bugs found on the way, all caught by mismatches against the census.** The census was the check, so each fix below
is argued for its own sake, not tuned to the answer:
1. **A Singular `size` bug.** `size` (non-zero generators) was used where `ncols` was meant, which dropped some
   defect equations. The local tori came out too small, and K₂ gave 4 instead of 20.
2. **Unlucky sample points.** A dense orbit was missed at special random points, in two blocks of the double fork.
   The script now resamples, with a wider range.
3. **Leaf primes that are not H-stable.** Leaf primes were treated as H-stable when the field was tangent only in the
   leaf's own coordinates (Q₁'s case-split block gave 89). The fix checks the eliminated coordinates too (§1). It
   then gave 17.
4. **Point leaves counted twice.** H does not fix them; they lie in orbits counted elsewhere. Q₂'s block gave 16; the
   fix (§1) gives 15.

The "7" of prediction 97 was produced between fixes 2 and 3, so it was wrong and is recorded as failed.

## 4. What is next

1. **Uniformity for the new certificates.** For each node, the claims to specialise are:
   - f ∉ I_X;
   - the closure J;
   - its minimal primes;
   - the absolute irreducibility of the nodes.

   The candidate bad primes would come from the leading coefficients of strong Gröbner bases over ℤ of J, as in
   v2.58. The ideals involved are the invariant hulls, which are small, so this should be cheaper than v2.58 was.
   This would close 1⇛2⇇3 + 1→3 uniformly.
2. **Conjecture G.** Prove it from Lang's theorem by computing the stabilisers. The orbit tree gives a rational point
   in every orbit (or says when it doesn't), so the component groups of the stabilisers are computable.
3. **Conjecture S** (v2.66) and Program 3, as before. The orbit tree gives Program 3 what it lacked: explicit orbit
   closures to compare across a leaf extension.

## Files

- **analysis/:**
  - v2_67_orbits.py (Lemma O; `--leafmode`);
  - v2_67_count.py (per-block comparison with the census);
  - v2_67_leafunion.py (case-split blocks).
- **Changed:** v2_53_block_census.py reads pickles with a numpy-only unpickler (some were uploaded).

- **logs/v2_67/:**
  - PREDICTIONS.md;
  - count_*.txt, count/<quiver>/<block>.{txt,json} (the orbit trees);
  - leafunion/;
  - orbits/ (the first runs on the stuck component);
  - census_triple_3.txt, census_pair_after3_3.txt;
  - *.blockcensus.txt.
