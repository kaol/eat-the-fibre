# Discordianism: Open Questions Explored (v2.70): Conjecture G for every p ≥ 5

*v2.69 proved Conjecture G (census count = geometric orbit count) for p ≫ 0 by spreading out, with N not computed.
This round makes it effective.
- **Certificates at every node.** For every node of every orbit tree (17 quivers, 1307 block trees, 7920 nodes), each step of Lemma L has a
  ℤ-certificate. Each gives the finite set of primes at which that step could fail:
  - G1: the orbit closure is irreducible mod p;
  - S: the stabiliser is smooth and connected mod p;
  - W: the census weighting moves are the expected ones;
  - separation of distinct nodes;
  - Lemma C: every orbit lies in some node.
- **Two cheap arguments.** Each makes a certificate affordable:
  - a reduced Gröbner basis over ℚ suffices for the orbit closure, because one inclusion is automatic;
  - the stabiliser equations are affine-linear in the substitution coordinates (at all 6696 nodes), so connectedness
    mod p reduces to a torus.
- **Theorem.** For every prime p ≥ 5 and each of the 17 quivers, the 𝔽_p local-torus census count of every vertex block
  equals the number of geometric orbits (v2.67). Three blocks are excluded: one case-split block each in Q₁, Q₂ and
  1⇉2→3→4→5. So the census counts (20, 34, 71, …, 752) hold for all p ≥ 5 on 14 of the quivers, and up to those
  three blocks on the other 3.
- **Two gaps in v2.69, now closed.** Lemma L silently used two facts:
  - every W(Y) is the span of a partition of the classes, which makes the census's *linear* weighting moves a torus;
  - at the 0-dimensional nodes, the weighting moves fix the point.
  Both are now checked.*

---

## 0. Results

The predictions are in logs/v2_70/PREDICTIONS.md, committed before the runs.

| # | Statement | Status |
|---|---|---|
| 102 | G1, S, W hold at every node over ℚ; after intersecting over points, the candidates lie in {2, 3, 5} | **failed as stated**: every check holds at every node (6664 dense, 1256 points), but intersecting over up to 10 points still left 5 at 81 nodes and 7 at 13 (all 13 in 1⇛2⇇3 + 1→3). Direct checks in characteristic 5 and 7 (§2) cleared all 94; final candidates {2, 3} |
| 103 | The Lemma C claims C1–C5 hold on every tree, with candidates in {2, 3} | **held**: C1–C4 give 2 and 3 (Rabinowitsch integers, empty boundaries); C5 gives no candidate; no hard failure |
| 104 | N_Q has prime factors in {2, 3, 5} for every quiver | **held, sharper**: N_Q has prime factors in {2, 3} for every quiver |

| quiver (v2.66/67 name) | trees + point strata | dense nodes | point nodes | census = orbits for p ≥ 5 | direct checks |
|---|---|---|---|---|---|
| K₂, K₃ | 6 + 0, 6 + 0 | 14, 28 | 6, 6 | 20, 34 | — |
| 1⇉2→3, 1→2⇉3 (pair_before, pair_after) | 21 + 2 each | 48 each | 21 each | 71, 71 | — |
| 1⇉2, 1→3 and 1⇉2, 3→2 (pair_srcout, pair_tgtin) | 27 + 0 each | 76 each | 25 each | 101, 101 | — |
| 1⇛2→3, 1⇛2, 1→3 (K3_after, K3_srcout) | 21 + 2, 27 + 0 | 97, 158 | 21, 25 | 120, 183 | 0, 3 at p = 5 |
| 1⇉2→3→4, 1→2⇉3→4 (pair_after2, pair_mid) | 93 + 7, 80 + 20 | 234, 157 | 93, 80 | 334, 257 | — |
| 1⇉2, 2→3, 2→4 and 1⇉2, 1→3, 4→2 (pair_tgtfork, pair_srcboth) | 128 + 20, 123 + 0 | 266, 379 | 128, 106 | 414, 485 | — |
| Q₀ (pairs_in_short) | 21 + 0 | 348 | 18 | 366 | — |
| 1⇛2⇇3 + 1→3 (triple_pair_in_short) | 21 + 0 | 734 | 18 | 752 | 74 at 5, 13 at 7 |
| Q₁ (pairs_in_short_tail) | 86 + 3 | 1123 | 77 | 1203 + (case-split block, 17 in char 0) | 2 at 5 |
| Q₂ (pairs_in_short_c4) | 98 + 0 | 1672 | 86 | 1758 + (case-split block, 15) | 2 at 5 |
| 1⇉2→3→4→5 (pair_after3) | 501 + 61 | 1206 | 500 | 1767 + (case-split block, 3) | — |
| **total** | 1307 + 117 | **6664** | **1256** | | 94 |

- **Point strata.** Each is a single fold, the reduction of one rational point in every characteristic (v2.66), so it
  contributes 1.
- **Agreement with v2.67.** The census column (dense + points + point strata) reproduces the v2.67 counts exactly.
- **Where the primes come from.**
  - 2 and 3 come from all four kinds of certificate.
  - 5 and 7 came only from the leading coefficients of the ℚ-basis of the slice graph (G1) and of the stabiliser's
    strong ℤ-basis (S). They come back at every point and every slice tried, so they are structural to these
    presentations, not to the geometry: the direct checks pass.

## 1. The census torus (a gap in v2.69)

The census (pcensus.py, v2.24 onwards) defines the local torus at a fold P as follows:
- the moves are the class scalings by s ∈ (𝔽_p^×)^{classes} that satisfy Σ_c s_c E_c(P) = 0;
- here E_c(P) is the Δ_c-defect.

So the move set is the 𝔽_p-points with non-zero coordinates of the *linear* space W_P = ker E(P), not the points of a
torus.

The two agree exactly when W_P is the span of the indicator vectors of a partition of the classes. In that case
s ∈ W_P with all s_c ≠ 0 means s is constant and non-zero on each block, which is T_{W_P}(𝔽_p).
- **Positive-dimensional nodes.** Checked: W(Y) is a partition span at all 6664 dense nodes (its reduced row echelon
  form is a set of block indicators).
- **Points of the dense orbit.** At a point x of the dense orbit, W_x = W(Y) (v2.69, prediction 101). At x_p the same
  holds when rank_p E(x) = rank_ℚ E(x).
- **Constancy along the orbit.** W_P is constant along G × T_W-orbits: substitutions commute with class scalings, and
  a block-constant scaling rescales E_c uniformly on each block.

So at every point of O_{Y,p} the census moves are exactly H_Y(𝔽_p). v2.69 used this without saying so.

**0-dimensional nodes** (single rational points x) need a separate check, since W_x need not be a partition span and
the census might move x to another point. Two identities, checked over ℤ:
- σ_g P_x = P_x σ_g, identically in g (G fixes x);
- φ_s P_x = P_x φ_s for s = Σ_k λ_k w_k, identically in λ, with (w_k) a ℤ-basis of W_x.

Then, whenever rank_p E(x) = rank_ℚ E(x), every census move at x_p fixes x_p.

## 2. The certificates (`analysis/v2_70_gcert.py`)

Fix a node Y with group H_Y. Let x be a rational point of its dense orbit, so the node's minor does not vanish at x.

**den.** The primes in the denominators of x, of P_x and of the stratum parametrisation. Away from them, x_p and the
equations below reduce correctly.

**W.** The defect matrix E(x) has rows indexed by the defect entries and columns by the classes.
- Checked: its kernel is W(Y), and W(Y) is a partition span.
- Candidates: the primes of its Smith invariants.

**S (stabiliser).** The equations σ_{g,t} P_x = P_x σ_{g,t} come with z·∏det(g_c) = 1 and t_k u_k = 1.
- **Checked:** at all 6696 nodes the equations are affine-linear in g. So K = Stab ∩ G is an affine subspace minus a
  determinant hypersurface: irreducible and connected, *in every characteristic*.
- **Candidates:** the leading coefficients of a strong ℤ-Gröbner basis for the block order (g, z) > (t, u). For p
  prime to them:
  - the basis reduces to a Gröbner basis mod p;
  - the (t, u)-part is p-saturated, hence equal to the ℤ-lattice ideal of the image torus T′;
  - T′ is connected over ℚ (v2.69), so its character lattice is saturated and T′_p is a torus.
- **Hence Stab_p is connected.** Its kernel to T′_p is K_p, and its image is T′_p.
- **Smoothness.** The Hilbert function of the basis gives dim Stab_p = dim H − dim Y. The Smith invariants of the
  Jacobian at 1 give the same tangent dimension. So Stab_p is smooth, and Lang applies.

**G1 (orbit closure), over ℚ only.** Take a slice Σ ∋ 1 of e = dim Y coordinates of (g, t): the pivots of the
differential of h ↦ h·x at 1. Compute the reduced Gröbner basis G of the graph ideal Γ of φ: Σ → Y over ℚ (block
order Σ-side > y), and check that its y-part is I_Y.

*Claim.* If p does not divide the leading coefficients of the primitive integral forms of the elements of G (nor the
den primes), then Y_p = closure(φ_p(Σ_p)). In particular Y_p is geometrically irreducible, with x_p in its dense
orbit.

*Proof.*
1. **G is p-integral.** Divide each element by its leading coefficient. The elements are then monic and
   p-integral.
2. **(Ḡ) contains Γ_p.** Every input generator reduces to 0 modulo G by p-integral division steps, so Γ_p ⊆ (Ḡ).
3. **Ḡ is a Gröbner basis.** S-polynomials of monic p-integral elements reduce to 0 by p-integral steps.
4. **One inclusion.** By elimination, ker φ_p* = Γ_p ∩ 𝔽_p[y] ⊆ (Ḡ) ∩ 𝔽_p[y] = (Ē). The p-integral reduced basis
   of I_Y,ℚ generates I_Y ⊗ ℤ_(p), so (Ē) = I_Y ⊗ 𝔽_p.
5. **The other inclusion holds for free.** The generators of I_Y vanish on the orbit: φ*(I_Y) = 0 is a p-integral
   identity modulo the monic relations z·det = 1 and tu = 1.
6. **Conclusion.** So ker φ_p* = I_Y ⊗ 𝔽_p. This kernel is prime, being the kernel of a map into the domain O(Σ_p).
   ∎

The strong ℤ-basis of Γ is never computed, and that matters: for the 10-dimensional root of 1⇛2⇇3 + 1→3 it ran out
of memory (4 GB in 4 minutes). The ℚ computation on a slice takes seconds.

**Separation.** Two nodes Y ≠ Y′ of the same dimension: x ∈ O_Y does not lie on Y′ over ℚ.
- Candidates are the primes dividing q(x) for every integral generator q of I_{Y′}.
- Away from them, x_p ∉ Y′_p, so O_{Y,p} ≠ O_{Y′,p}.
- Two point nodes: the primes dividing all their coordinate differences.

**Points.** A prime is kept as a candidate for Y only if it is a candidate at every point tried. Points are added
until no candidate ≥ 5 survives, up to 10 points.

**Direct checks at a residual prime.** Some primes q ≥ 5 survive the intersection over points. Take a point where q
comes only from the leading coefficients of G1 or S.
- **G1 directly over 𝔽_q.** Choose a slice by pivots *mod q*: a slice chosen over ℚ can lose transversality mod q,
  which happened at one Q₁ node, where its image mod 5 had dimension e − 1. Then the elimination ideal of the graph
  ideal over 𝔽_q must equal I_Y ⊗ 𝔽_q. That is a valid comparison because C1 gives I_Y a strong ℤ-basis with unit
  leading coefficients; the tally checks this for every direct check.
- **S directly over 𝔽_q.** The stabiliser ideal over 𝔽_q must have exactly one minimal prime, of dimension
  dim H − dim Y. The identity component is defined over 𝔽_q, so one minimal prime means connected. Smoothness comes
  from the Jacobian at 1, which is not divisible by q at that point.

94 nodes needed this: 81 at q = 5, 13 at q = 7. Every one passed.

**Lemma C** (v2.68, `v2_68_uniform.py`, `v2_68_cover.py`, now on every tree): claims C1–C5 over ℤ. They give that
every orbit of S_p is the top orbit of some node Y_p.

## 3. The theorem

**Lemma L (v2.69), restated with its hypotheses.** Let p be a prime and B a vertex block. Suppose the following
hold at p:
- for every node Y of B's tree: G1, S, W and separation (dense nodes), or the point checks (0-dimensional nodes);
- for the tree: C1–C5.

Then the local-torus census count of B over 𝔽_q, for any power q of p, equals the number of nodes.

*Proof.*
1. **Every census class lies in some O_{Y,p}.** By C1–C5 and Lemma C (v2.68), every H-orbit of S_p is a top orbit
   of some Y_p. By G1, Y_p is irreducible, so it has exactly one top orbit, O_{Y,p} = H_p·x_p.
2. **Distinct nodes give distinct sets.** By separation, the O_{Y,p} are pairwise distinct for distinct nodes. For
   points this is the coordinate check; across dimensions the orbit dimensions differ.
3. **Census classes stay inside O_{Y,p}.** By W and §1, the census moves at every point of O_{Y,p} are H_Y(𝔽_q),
   and at a point node they fix the point.
4. **Exactly one class per node.** By S, Stab_{x_p} is smooth and connected. Lang gives an 𝔽_q-point of O_{Y,p} and a
   single H_Y(𝔽_q)-orbit on O_{Y,p}(𝔽_q). ∎

The vertex blocks are the same list in every characteristic, and the stratum parametrisations specialise (v2.58,
v2.66, unit leading coefficients). Blocks with distinct entries in {0, 1, −1} stay distinct for p ≠ 2.

**Theorem G₅.** Let Q be one of the 17 quivers of §0 (or its opposite), and q a power of a prime p ≥ 5. Then the
local-torus census count of Folds(𝔽_q Q) over the vertex blocks other than the case-split block equals the number of
geometric orbits of the local census groups (the v2.67 count).
- For K₂, K₃, the eight pair quivers, the two 1⇛2 quivers, Q₀ and 1⇛2⇇3 + 1→3, this is the whole census: 20, 34, 71,
  71, 101, 101, 120, 183, 334, 257, 414, 485, 366, 752 for every such q.
- For Q₁, Q₂ and 1⇉2→3→4→5 it is the whole census except one block each.

*Proof.* Lemma L at every block. The candidate primes of all 7920 node certificates and of the 1307 tree claims lie
in {2, 3}. The 94 nodes where a ℤ-certificate left 5 or 7 were checked directly in that characteristic. ∎

**What is still conjectural.** Conjecture G at p = 2 and 3, and on the three case-split blocks.
- **p = 2 and 3.** The census agrees there numerically (q = 3 for all 1427 blocks), but the certificates fail at
  these primes, as they did in v2.68.
- **The case-split blocks.** Their leaves come from v2.57's pivot case split, which would itself have to be
  revalidated mod p.

**How this sits with the census data.** The census was run at q = 3 for every block, and at q = 5 (and 7) for some
quivers (Q₀: 366 at q = 3, 5, 7; 1⇛2⇇3 + 1→3: 752 at q = 3, 5). Theorem G₅ covers q = 5 and 7 and says nothing new at q = 3, the one prime where every block was compared.
The q = 3 agreement is therefore the evidence that p = 3 is not actually bad.

## 4. Next

1. **p = 2 and 3.** Two routes:
   - redo Lemma L in characteristic p from scratch, as v2.68 §3 did for Lemma C: char-p trees exist so far only for
     block 100,001,001 of 1⇛2⇇3 + 1→3;
   - or explain the census agreement at q = 3 structurally.
   At p = 2 the vertex blocks themselves merge (1 ≡ −1), so the statement must be reformulated first.
2. **The three case-split blocks**: specialise the v2.57 case split (its pivots are ±1, which suggests a clean
   ℤ-statement), then the leaf trees.
3. **New quivers** beyond the census list, now with the full chain v2.66 → v2.70 automatic.
4. **Conjecture S, and Program 3.**

## Files

- **analysis/:**
  - v2_70_gcert.py (G1, S, W, separation, point nodes, per block);
  - v2_70_all.py (one quiver: gcert, C1–C4, C5 for every block).
  - v2_70_tally.py (hard failures, candidate primes, direct checks).
- **Changed:** v2_68_cover.py (a root equal to the whole stratum).
- **logs/v2_70/:**
  - PREDICTIONS.md;
  - all_<quiver>.txt;
  - <quiver>/<block>.{gc,u,cv}.json, with the points used stored in the newer gc files;
  - tally.json.
