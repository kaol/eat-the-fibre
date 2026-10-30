# Discordianism: Open Questions Explored (v2.77): Theorem 5 in Lean on every quiver with a full fork; the six classes accounted for

*Three threads from v2.76.*
1. *Lean. Theorem 5 is now kernel-checked on every acyclic quiver that contains 1⇛2, 2→3, 2→4 as a full subquiver:
   - both fork families are folds over every commutative ring;
   - they are rigid under every symmetry that acts on the subquiver as the census does;
   - the two families are disjoint.*
2. *1⇛2, 2→3→5, 2→4 grows by 6 from q = 2 to q = 3, and all six are now accounted for. Four sit in the identity block:
   the four vertex-part families of v2.76's path version. Two sit in the block that kills the far sink 5: there the
   quiver behaves like 1⇛2, 2→3, 2→4, and the arrow version supplies the two families.*
3. *The infinitesimal test separates the cases. On the four quivers covered by Theorem 5′, the path family has exactly
   one modulus. On 1⇛2→3⇉4, where the paths pass through a 2-class, it has none.*

---

## 0. Results

The predictions are in logs/v2_77/PREDICTIONS.md, committed before the runs.

| # | Statement | Status |
|---|---|---|
| 135 | infinitesimal moduli of the path family: 0 on 1⇛2→3⇉4, exactly 1 on the four Theorem 5′ quivers | **held** (two random members each, exact over ℚ) |
| 136 | 1⇛2, 2→3→5, 2→4: the identity block grows by exactly 4; the other 2 extra classes sit in a single other block | **held**: identity block 278 → 282, block 10000,01000,00100,00010,00000 111 → 113, every other block of 758 unchanged |
| 137 | 1⇛2→3, 3→4, 3→5 grows by at least 4 from q = 2 to q = 3 | **pending**: the q = 3 census was at fibre group 10 of 12 after 105 minutes and dies when the turn ends. Theorem 5′ already gives at least 4(q − 1) classes in the identity block |

## 1. Theorem 5 on every quiver with a full fork (Lean)

`lean/EatTheFibre/ForkGeneral.lean` follows v2.30's `General.lean` (Theorems 3 and 4 on every acyclic quiver).

**The structure `Fork Q`.** Q is a finite quiver with an acyclicity potential, containing:
- vertices i, j, k, k′ with k ≠ k′;
- three distinct arrows a₁, a₂, a₃ : i→j;
- arrows x : j→k and x′ : j→k′;
- fullness: these five are the only arrows of Q with both ends in {i, j, k, k′}.

Nothing else is assumed. Q may have more arrows into or out of these vertices, longer paths, and other classes.

**The results.** All four depend only on propext, Classical.choice and Quot.sound
(`logs/v2_77/lean_audit_forkgeneral.txt`).
- **`fork_general`.** Push P_t and P⁻_t (v2.75, `Gen/K3_tgtfork.lean`) into kQ, zero on every other path. Both are folds
  of kQ, for every commutative ring and every t.
- **`fork_general_rigid`, `forkM_general_rigid`.** Let h be any endomorphism of kQ that acts on the span of the
  subquiver's paths as the census symmetry substK g wa wb wc. If h intertwines P_t with P_t′, then t = t′. No
  invertibility is assumed. By fullness, the three classes of the subquiver are classes of Q, so every class
  substitution and weighting of Q acts this way.
- **`fork_forkM_general_disjoint`.** No such h joins a member of one family to a member of the other, over a
  nontrivial ring.

**The proof.**
1. The census table of 1⇛2, 2→3, 2→4 is its path algebra (`decide +kernel` on 15³ structure constants, both for
   products and for Δ).
2. The subquiver gives an injective quiver morphism, hence an `Embeds`.
3. `isFold_push` leaves one condition: P·Δ(e_l)·Pᵀ = 0 for every path l off the image.
4. **Off the image.** Write a deletion as l = p·c·q with p, q image paths. Then c joins two vertices of the subquiver,
   so by fullness it is one of the five arrows, and l is itself an image path, a contradiction (`Fork.delta_out`).
5. **Rigidity.** It restricts to v2.75's `fork_rigid_GT` through the push (`Fork.restrict`).

**What this covers.** Of the quivers shown infinite in v2.75–76, the Lean theorem covers:
- 1⇛2, 2→3, 2→4;
- 1⇛2 with three tails (take any two);
- 5→1⇛2, 2→3, 2→4;
- 1⇛2, 2→3, 2→4, 5→2;
- 1⇛2, 2→3→5, 2→4.

It does not cover four of them:
- 1⇛2, 2→3, 2→4, 1→3 and 1⇛2, 2→3, 2→4, 3→4: the extra arrow breaks fullness;
- the K₄ variant: the fourth parallel arrow breaks fullness;
- 1⇛2→3, 3→4, 3→5: there is no fork of arrows.

For those, Theorem 5 or 5′ stands on the paper proof and the computer checks of v2.75–76.

## 2. The six classes of 1⇛2, 2→3→5, 2→4

The census was rerun with saved representatives at q = 2 and q = 3 (3417 and 3423, as in v2.75). It was then split by
vertex block, identifying blocks mod 2 (`analysis/v2_77_blockdiff.py`). Two of the 758 blocks change:

| block | q = 2 | q = 3 | |
|---|---|---|---|
| identity | 278 | 282 | +4: the four vertex-part families of the path version (π = a3·a5, π′ = a4; X = {1}, {1, 3} and complements) |
| 10000,01000,00100,00010,00000 (e₅ killed) | 111 | 113 | +2 |

The second block has exactly the counts of 1⇛2, 2→3, 2→4's identity block (111, 113). That is no coincidence:

**Killing a successor-closed set.** Let K be a set of vertices closed under successors. Every path touching K ends in
K. So a fold of the full subquiver on the other vertices, extended by zero, is a fold of Q: deleting an arrow from a
path ending in K leaves a suffix ending in K, which P kills. This is the same `isFold_push` mechanism as §1.

With K = {5}, the arrow a3 ends at a sink of Q ∖ K. So Theorem 5 itself (arrow tails a3, a4, both families) gives the
+2 there. `analysis/v2_77_killsink.py` checks this directly: the arrow family with every path touching vertex 5 sent to
0 is a fold of 1⇛2, 2→3→5, 2→4, identically.

So the census growth of this quiver from q = 2 to q = 3 is 6 = 4 + 2, all explained. The count is therefore at least
3417 + 6(q − 2) for q ≥ 2, if nothing else grows at larger q.

## 3. Infinitesimal moduli (`analysis/v2_77_tangent.py`)

At a random rational member P of the path family, the test compares two spaces, exactly over ℚ:
- the tangent of the H-orbit, T = ⟨[D, P] : D ∈ 𝔥⟩, where 𝔥 is spanned by the derivatives at 1 of the class
  substitutions and weightings;
- the four family directions ∂P/∂λ, ∂P/∂m, ∂P/∂λ′, ∂P/∂m′.

| quiver | dim 𝔥 | rank T | moduli |
|---|---|---|---|
| 1⇛2, 2→3, 2→4 | 9 | 6 | 1 |
| 1⇛2, 2→3→5, 2→4 | 10 | 6 | 1 |
| 1⇛2→3, 3→4, 3→5 | 10 | 6 | 1 |
| 1⇛2, 2→3, 2→4, 3→4 | 10 | 6 | 1 |
| 1⇛2→3⇉4 | 11 | 7 | **0** |

On 1⇛2→3⇉4 the 2-class's substitution group adds one orbit direction, and that absorbs the cross-ratio. This is the
mechanism behind Conjecture 3ᶜ″'s exception for branchings through a parallel class. It is not a proof of finiteness:
other families could still carry moduli. The census of this quiver at q = 3 is still out (v2.76).

## 4. Census

| quiver | q = 2 | q = 3 | |
|---|---|---|---|
| 1⇛2, 2→3→5, 2→4 | 3417 | 3423 | rerun with saved representatives; block split in §2 |
| 1⇛2→3, 3→4, 3→5 | 4899 | unfinished (group 10 of 12 after 105 min) | infinite by Theorem 5′ |

## 5. Next

1. **Theorem 5′ in Lean,** for paths. The push argument needs a full subquiver, so the path version needs its own
   construction: Π₀ + N is the identity off the subquiver, not zero.
2. **1⇛2→3⇉4.** Decide it. Its q = 3 census did not finish in 76 minutes. Alternatively, run the identity-block
   pipeline with exact generic ranks.
3. **The finite half of 3ᶜ″.** Run the v2.66 driver on 1→2⇛3→4 (it resumes, but its pipeline needs 3.5 GB, so it
   cannot run beside a census here).
4. **Exact counts** at q = 5, which need a bigger machine.

## Files

- **analysis/:**
  - v2_77_tangent.py (infinitesimal moduli);
  - v2_77_blockdiff.py (census difference by vertex block);
  - v2_77_killsink.py (families with a successor-closed vertex set killed).
- **lean/:**
  - EatTheFibre/ForkGeneral.lean;
  - EatTheFibre.lean, Audit.lean, build.sh and README updated.
- **logs/v2_77/:**
  - PREDICTIONS.md;
  - tangent_*.json;
  - census_* (the representatives are not committed: *.pkl is ignored);
  - K3_tgtfork_long.blockcensus23.txt, K3_tgtfork_long.blockdiff.txt;
  - killsink_K3_tgtfork_long.json;
  - lean_forkgeneral.log, lean_audit_forkgeneral.txt.
