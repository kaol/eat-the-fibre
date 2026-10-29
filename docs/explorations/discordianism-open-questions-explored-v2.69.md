# Discordianism: Open Questions Explored (v2.69): connected stabilisers, and Conjecture G for p ≫ 0

*v2.67 observed Conjecture G on every block: the 𝔽_q census count equals the number of geometric orbits, for every q
the census reached. v2.68 proved this number is finite in every characteristic, but said nothing about why the census
does not see the orbits split over 𝔽_q. This round adds the missing piece and turns the observation into a theorem
for all large characteristics.
- **Stabilisers.** At a rational point of every dense orbit, in every orbit tree (17 quivers, every block, including
  the leaves of the three case-split blocks), the stabiliser in the node's group H_Y was computed exactly over ℚ. It
  is connected every time, of dimension dim H_Y − dim Y.
- **Local weightings.** At every node, the weightings that satisfy the compatibility condition on Y are already
  tangent to Y. So the census move group at a point of a dense orbit is exactly H_Y.
- **Theorem G (p ≫ 0).** For each of the 17 quivers there is an N such that, for every prime power q prime to N, the
  census count equals the geometric orbit count. The proof is Lang's theorem twice, plus spreading out.
- **A bug and its fix.** The first pass reported 12 dimension mismatches in one pair_mid block. The cause was that
  the torus relations t·u = 1 had been dropped. They are restored, and the rerun is clean.*

---

## 0. Results

The predictions are in logs/v2_69/PREDICTIONS.md, committed before the runs.

| # | Statement | Status |
|---|---|---|
| 100 | Every positive-dimensional node has a connected stabiliser, of dimension dim H_Y − dim Y | **held**: 6696 of 6696 nodes connected, no dimension mismatch, no node without a rational point in its dense orbit |
| 101 | At every node, the defect-compatible weightings equal W(Y) (all tangent) | **held**: W_defect > W(Y) at 0 nodes |

| quiver (v2.66/67 name) | blocks with a tree | dense nodes | connected |
|---|---|---|---|
| K2, K3 | 6, 6 | 14, 28 | all |
| pair_before, pair_after, pair_mid | 21, 21, 80 | 48, 48, 157 | all |
| pair_srcout, pair_tgtin | 27, 27 | 76, 76 | all |
| K3_after, K3_srcout | 21, 27 | 97, 158 | all |
| pair_after2, pair_tgtfork, pair_srcboth | 93, 128, 123 | 234, 266, 379 | all |
| Q₀ (pairs_in_short) | 21 | 348 | all |
| 1⇛2⇇3 + 1→3 (triple_pair_in_short) | 21 | 734 | all |
| Q₁ (pairs_in_short_tail) | 86 + leaves of 0000,0000,0000,0001 | 1123 + 16 | all |
| Q₂ (pairs_in_short_c4) | 98 + leaves of 1000,010m,010m,0000 | 1672 + 14 | all |
| 1⇉2→3→4→5 (pair_after3) | 501 + leaves of 00000,00000,0010m,00000,00000 | 1206 + 2 | all |
| **total** | 1307 + 3 case-split | **6664 + 32 = 6696** | **all** |

- **The main trees.** The 6664 dense nodes are exactly the positive-dimensional nodes of the v2.67 survey (6663
  found at rational points directly, 1 through absPrimdecGTZ).
- **The leaves.** Each case-split block's leaf trees, taken once per orbit closure, give 16, 14 and 2 dense orbits.
  With the one H-fixed point leaf in each, that is the 17, 15 and 3 closures of v2.67.
- **Zero-dimensional nodes** are single rational points fixed by H_Y. Their stabiliser is H_Y, which is connected.
- **The 117 blocks without a tree** are point strata: a single fold, H-stable and hence fixed. Together,
  1307 + 117 + 3 = 1427 blocks.

## 1. The computation (`analysis/v2_69_stab.py`)

Take a node Y of an orbit tree with a dense orbit O_Y.
- **The group.** H_Y = ∏_c G_c × T_{W(Y)}:
  - G_c is the group of invertible m_c × m_c matrices with column sums 1, for each parallel class c of size m_c. It
    acts by a_j ↦ Σ_i g_ij a_i. Its Lie algebra is Der_Δ, which the script asserts: dim G = dim Der_Δ.
  - T_{W(Y)} is the torus of the weightings: t ↦ (s_c) with s_c = ∏_k t_k^{w_{k,c}}, where the w_k are the primitive
    integral vectors of a basis of W(Y).
  - Both act multiplicatively on paths, which gives a d × d matrix σ_{g,t} polynomial in g, t, u (with u = t⁻¹).
- **The point.** A rational point x of Y at which the node's e × e minor is non-zero, so x lies in O_Y.
- **The stabiliser.** The scheme defined by:
  - σ_{g,t} P_x − P_x σ_{g,t} = 0, with denominators cleared;
  - z · ∏_c det(g_c) − 1 = 0;
  - t_k u_k − 1 = 0.
- **Its minimal primes over ℚ** (minAssGTZ) and their dimensions.

**Why one minimal prime means connected.**
- The identity component Stab° is geometrically irreducible and contains the rational point 1, so it is defined
  over ℚ.
- Every irreducible component over ℚ is a union of Galois orbits of geometric components, and Stab° is
  Galois-stable. So Stab° is by itself one ℚ-component.
- Hence one minimal prime means Stab = Stab°.

The dimension check, dim Stab = dim H_Y − dim Y, confirms that x lies in the dense orbit and that the group was set up
correctly. This check is what caught the bug.

**The bug.** In the first pass, the line that cleared denominators also swallowed the relations t_k u_k = 1 into a
comment. The u_k were then free, which adds exactly len(W) to every dimension: 6 instead of 3 with three weightings,
4 instead of 2 with two.
- It showed in one block, pair_mid 1000,0100,0010,0001; whatever else ran with the faulty line is covered by the
  rerun below.
- One of its nodes even appeared disconnected (two components of dimension 7). This was an artefact of the free u.
- After the fix the whole batch was rerun from scratch: all 17 quivers, plus the case-split leaves.

## 2. From the computation to the census: Lang twice

Fix a quiver and a block, and let k = 𝔽̄_p. The v2.67 tree gives, in characteristic 0:
- a finite set of nodes Y, each with a group H_Y and an open dense orbit O_Y ⊆ Y;
- the O_Y partition the stratum S, since every orbit is the dense orbit of exactly one node.

Suppose the following facts hold over 𝔽_p:

| | fact over 𝔽_p | char 0 source |
|---|---|---|
| (a) | the tree specialises: the nodes Y_p, with O_{Y,p} open and dense and the O_{Y,p} partitioning S_p | v2.67 tree |
| (b) | each Y_p is defined over 𝔽_p and geometrically irreducible | a smooth rational point on a ℚ-irreducible Y (the point x above: orbits are smooth) |
| (c) | Stab_{H_Y}(x_p) is connected and smooth | §1 (smooth in char 0 by Cartier) |
| (d) | at every point of O_{Y,p}, the compatible weightings are W(Y) | prediction 101 |

**Lemma L.** Under (a)–(d), the census count over 𝔽_q (q a power of p) equals the number of nodes.

*Proof.*
- **Census classes stay inside one O_Y.** Census moves are the substitutions G(𝔽_q) and the weightings compatible
  at the current point. G and the class torus commute, so the compatible set is constant on G·T-orbits. Hence it is
  constant on O_Y, and by (d) it is W(Y). So the moves available at any point of O_Y(𝔽_q) lie in H_Y(𝔽_q), and they
  keep the point in O_Y. Distinct nodes give distinct classes.
- **Every node has a census point.** O_Y ≅ H_Y / Stab_x is a homogeneous space over 𝔽_q under a connected group (by
  (b), O_Y is defined over 𝔽_q). By Lang's theorem it has an 𝔽_q-point.
- **One class per node.** Let y ∈ O_Y(𝔽_q). The H_Y(𝔽_q)-orbits on O_Y(𝔽_q) are in bijection with
  ker(H¹(𝔽_q, Stab_y) → H¹(𝔽_q, H_Y)). Stab_y is a k-conjugate of Stab_x, so it is connected by (c). It is also
  smooth for p ∤ N (below). Hence H¹(𝔽_q, Stab_y) = 1 by Lang.
- **The torus as the census sees it.** The census scales class c by s_c ∈ 𝔽_q^* along the integral weightings. The
  image of T_{W(Y)}(𝔽_q) is contained in these moves, so the census class of y contains H_Y(𝔽_q)·y = O_Y(𝔽_q). By
  the first point it is no larger. ∎

**Theorem G (p ≫ 0).** For each of the 17 quivers of §0, there is an integer N > 0 such that for every prime power q
prime to N, the census count of every block equals its geometric orbit count. So the census count of the quiver
equals Σ_blocks (geometric orbits) for all such q.

*Proof.* Every char 0 fact used in (a)–(d) is a finite collection of statements about finitely presented ℚ-schemes:
- the node ideals and their dimensions;
- the closure ideals and minors;
- the orbit maps H_Y → Y, x;
- the stabiliser and its single minimal prime;
- the weighting defects.

Spread these out over ℤ[1/N] for a suitable N. For p ∤ N:
- the fibres keep their dimensions, irreducibility, and the containments of the tree (EGA IV 9.7.7 for geometric
  irreducibility and connectedness of fibres);
- they keep their rank conditions;
- the smooth rational point x_p of Y_p gives (b);
- the stabiliser group scheme over ℤ[1/N] has geometrically connected fibres, giving (c).

Lemma L applies. ∎

**What this explains.** The observed q-independence, census = geometric orbits wherever the census was computed (q = 3 on all 1427 blocks of v2.67,
flat in q since v2.27), is now a theorem for all large p. The theorem is not effective: N
is not computed. The v2.68 machinery (strong ℤ-bases, Rabinowitsch integers) would compute it. The new ingredient is
the stabiliser over ℤ: its strong ℤ-basis, and a check that the single component stays single mod p.

**Why the census agreement had predicted it.** If a stabiliser had component group A ≠ 1, then for q with Frobenius
acting trivially on A, O_Y(𝔽_q) would split into as many classes as A has conjugacy classes. The census, which
counts one class per node, would then be larger. Prediction 100 was the census agreement read backwards.

## 3. Next

1. **Effective N.** Compute the stabiliser ideals over ℤ, as v2.68 did for the trees, and check (a)–(d) at the
   candidate primes directly. That would give Conjecture G for every q, and with it a closed formula for each census
   count.
2. **New quivers** beyond the census list (K₄, 1⇛2⇛…), through the now automatic pipeline: the v2.66 driver, then
   the v2.67 trees, the v2.68 specialisation, and v2.69 stabilisers.
3. **Conjecture S, and Program 3.**

## Files

- **analysis/:**
  - v2_69_stab.py (one block's stabilisers and local weightings);
  - v2_69_stab_all.py (one quiver);
  - v2_69_stab_leaves.py (case-split blocks, one entry per orbit closure).
- **logs/v2_69/:**
  - PREDICTIONS.md;
  - stab_<quiver>.txt with stab/<quiver>/ (per-block JSON; the Singular inputs are regenerated, not committed);
  - leaves_*.txt with leaves/<block>/.
