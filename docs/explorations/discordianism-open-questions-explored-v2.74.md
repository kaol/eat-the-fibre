# Discordianism: Open Questions Explored (v2.74): beyond the census list; a counterexample to Conjecture 3ᶜ at multiplicity 3

*With Conjecture G complete on the 17 quivers of the census list (v2.73), this round ran the chain on five quivers
beyond it: K₄ and four quivers with a class of three parallel arrows. The first step, the census itself, found
something new. 1⇛2, 2→3, 2→4 has 708 classes at q = 2 and 710 at q = 3. The difference sits in one vertex block.
That block has two 7-dimensional components on which every orbit has dimension at most 6.
- **The invariants.** Each component carries an explicit H-invariant rational function, a ratio of two products of
  coordinates, certified over ℤ.
- **The consequence.** The quiver has infinitely many orbits, although it has a single parallel class.
- **Conjecture 3ᶜ fails at multiplicity 3.** Its multiplicity-2 version 1⇉2, 2→3, 2→4 is finite (414 classes for
  every q, v2.73).*

---

## 0. Results

The predictions are in logs/v2_74/PREDICTIONS.md. 115–117 were committed before any run; 118–120 were added, as an
addendum, after the census and before the checks that decide them.

| # | Statement | Status |
|---|---|---|
| 115 | K₄: census 48 = 14·4 − 8 at q = 2, 3; Conjecture G with count 48 at every q | **census held** (48, 48). For odd q, Conjecture G for K₄ is v2.18's Theorem 1 (14m − 8 orbit types over any field of characteristic ≠ 2). p = 2 not yet certified |
| 116 | the four K3-class quivers are flat at q = 2, 3 and match their trees block by block | **failed** for 1⇛2, 2→3, 2→4 (708 vs 710). Held for 1⇛2→3→4 (570 at q = 2, 3, 5), 1→2⇛3→4 (419, 419) and 1⇛2, 1→3, 4→2 (893, 893) |
| 117 | the chain completes on all five without new methods; Conjecture G at every q | **failed**: 1⇛2, 2→3, 2→4 is infinite (§2). 1⇛2→3→4 completed the finiteness stage (finite, uniformly in q, every characteristic); the rest is pending (§3) |
| 118 | on components 12, 13 of block 1000,0100,0010,0001 the exact generic orbit rank is 6 < 7: infinitely many orbits | **held** |
| 119 | the census grows with q: the q = 5 count exceeds 710 | **pending**: the q = 5 census (5160 H-orbit fibres) had not finished after 66 minutes; it is checkpointed. Growth itself follows from the invariants (§2) |
| 120 | the census is 706 + 2(q − 1) (714 at q = 5) | **pending** (same run) |

## 1. The census beyond the list (`pcensus.py`, torch engine)

| quiver | q = 2 | q = 3 | q = 5 |
|---|---|---|---|
| K₄ | 48 | 48 | (2.5·10⁶ kernel folds; not run) |
| 1⇛2→3→4 (K3_after2) | 570 | 570 | 570 |
| 1→2⇛3→4 (K3_mid) | 419 | 419 | — |
| 1⇛2, 2→3, 2→4 (K3_tgtfork) | **708** | **710** | (running) |
| 1⇛2, 1→3, 4→2 (K3_srcboth) | 893 | 893 | — |

The torch engine does each q = 2, 3 census in seconds to minutes. The pool engine had not finished 1⇛2→3→4 at q = 3
after 17 minutes; torch took 14 seconds.

**Where the 2 classes come from.** I ran the census split by vertex block at q = 2 and 3, normalising the blocks mod
2: the printed symmetric lifts write 1 as −1 at p = 2. All 148 blocks agree except the identity block
1000,0100,0010,0001, which has 111 classes at q = 2 and 113 at q = 3. At q = 3 that block has 19,630 folds in 158
G-orbits, merged to 113 by the local torus.

## 2. 1⇛2, 2→3, 2→4 has infinitely many orbits

**The block** 1000,0100,0010,0001 (`analysis/v2_55_pipeline.py`, `logs/v2_74/tgt/`):
- presolved to 1338 equations in 45 coordinates;
- 29 components; Der_Δ has dimension 6, with 3 classes.

27 components pass the v2.54 criterion. Components 12 and 13, both of dimension 7, do not: the orbit rank is ≤ 6
at every sampled point.

**Exact generic rank** (`analysis/v2_74_genrank.py`). Fraction-free Gaussian elimination of the orbit matrix runs
over R/I_Z, which is a domain, so the result is the rank at the generic point with nothing sampled. It uses the
largest possible group, 𝔥 = Der_Δ + all three class weightings (dimension 9).
- components 12 and 13: generic orbit rank **6**, dimension 7;
- the control component 29: rank 4 = dimension 4.

In characteristic 0 the orbit dimension at x is the rank at x, and the generic rank is the maximum. So every orbit in
components 12 and 13 has dimension ≤ 6, and each component is a union of infinitely many orbits.

**The invariants** (`analysis/v2_74_orbitclosure.py`, `v2_74_invariant.py`). The closure of the orbit of a random
rational point is I_Z plus one condition. Eliminating the orbit map of a 6-dimensional slice at three points gives
I_Z + (a·y0_8·y6_7 + b·y0_7·y6_8, …), with (a : b) varying with the point. So

  φ₁₂ = (y0_7 · y6_8) / (y0_8 · y6_7) on component 12,  φ₁₃ = (y3_7 · y6_8) / (y3_8 · y6_7) on component 13.

The y-coordinates are entries of P in the complement columns.

Certified for both components:
- **Field-invariance:** for every field ξ of 𝔥 (integral rows), ξ(N)·D − N·ξ(D) ∈ I_Z, by exact reduction over ℚ
  and over ℤ (a strong basis with leading coefficients ±1);
- **φ is defined and non-zero:** neither N nor D lies in I_Z;
- **φ is non-constant:** three points give three different values.

H is connected, so in characteristic 0 field-invariance is H-invariance. The identity
N(h·x)·D(x) − N(x)·D(h·x) ∈ I_Z ⊗ O(H) holds over ℚ, hence over ℤ (both sides are saturated) and mod every p.
So for every q the census classes on Z(𝔽_q) lie in level sets of φ, and their number grows with q.

**The count.** Each component's generic level sets are single orbit closures (the eliminations above). That
suggests one class per value of φ, with a fixed number of exceptional values: 706 + 2(q − 1), which gives 708 and 710.
This is checked against q = 5 once that census finishes.

**Conjecture 3ᶜ.** It predicts finiteness exactly when there are no consecutive parallel classes and no
converging/diverging pair of classes with unjoined ends. 1⇛2, 2→3, 2→4 has a single parallel class, so 3ᶜ predicts
finite, and it is infinite. Its multiplicity-2 version 1⇉2, 2→3, 2→4 is finite, with census 414 at every q
(v2.73). **So Conjecture 3ᶜ is false as stated: finiteness depends on the multiplicities, not just on the shape of
the class graph.**
- **The invariant.** It is a cross-ratio between two arrow paths through the class's target vertex: the coordinates
  y0_·, y3_· are paths through the class, and y6_· involves the tail arrows.
- **The likely mechanism.** At multiplicity 3 the class's substitution group has dimension 6 (AGL₂-like). With the
  fork, the compatible folds have one modulus too many for it.

## 3. The chain on the other quivers

- **1⇛2→3→4: finite, uniformly in q, in every characteristic** (the v2.66 procedure, `logs/v2_74/K3_after2/`):
  - **integral structure:** Smith form 1, det B = 1;
  - **vertex blocks:** 100 points. One leading coefficient is 2; the direct decomposition over 𝔽₂ gives the same
    100 points, linear and injective (`vertex_badprimes.json`);
  - **pipeline:** 100 of 100 blocks pass, 378 components or leaf primes;
  - **specialisation:** 4986 claims, candidate primes 2 (15 claims), 3 (4) and 5 (1), all checked directly; no
    exceptional prime;
  - **Bézout bound:** 10^67.8.

  The census is 570 at q = 2, 3, 5. Conjecture G for it (trees, certificates) is next.
- **1→2⇛3→4, 1⇛2, 1→3, 4→2.** The driver started on 1→2⇛3→4; it resumes with RESUME=1. Their census is flat at q = 2, 3, but 1⇛2, 2→3, 2→4 shows that a family
  can hide behind small q only when its count is constant there. Its own census did move between q = 2 and 3, so
  flatness is still informative.
- **K₄ at p = 2.** Its zero block is expensive (over 15 minutes in Singular without finishing); deferred.

## 4. Next

1. **Explain the family.** Find a characteristic-free statement for the family on 1⇛2, 2→3, 2→4 (a Theorem 2/3/4-type
   family), and find which multiplicity and shape combinations produce it. Test candidates: 1⇛2, 2→3, 2→4, 2→5;
   1⇛2→3, 2→4; and the K₄ versions. This replaces Conjecture 3ᶜ with a multiplicity-aware version.
2. **Finish the chain** on 1⇛2→3→4 (specialisation, trees, certificates), then 1→2⇛3→4 and 1⇛2, 1→3, 4→2.
3. **K₄ at p = 2.**

## Files

- **analysis/:**
  - v2_74_genrank.py (exact generic orbit rank);
  - v2_74_orbitclosure.py (orbit closures of points);
  - v2_74_invariant.py (certifying an invariant over ℚ and ℤ).
- **logs/v2_74/:**
  - PREDICTIONS.md;
  - census_<quiver>_<q>.{json,txt};
  - K3_tgtfork.blockcensus23.txt;
  - tgt/ (the identity block: certificates, tree, generic ranks, orbit closures, invariants);
  - K3_after2/ (v2.66 driver).
