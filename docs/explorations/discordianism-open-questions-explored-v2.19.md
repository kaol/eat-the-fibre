# Discordianism: Open Questions Explored (v2.19)

*Takes up question 2 of `discordianism-open-questions-explored-v2.18.md` §5: the necessity of (H), v2.17 Conjecture 5, and the dimension of the looped and cyclic fold families.*

*Tags **[Zn]** are machine checks in `verify_v2_19.py` (folder `v2_19/`, zipped as `verify_v2_19.zip`; `pip install -r requirements.txt`, and `SLOW=1` adds L = 4 and the 4-vertex sample). Labels: **proved**, **verified (exact)** or **conjectured**. Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Summary

| # | Question | Answer |
|---|---|---|
| v2.17 Conj. 5 | Is (H) necessary for uniqueness? | **No: refuted.** On the 4-vertex quiver 0→2, 0→3, 2→0, 3→1, a vertex fold violates (H) twice, a first-order direction exists, and yet the exact fibre is the single canonical fold at every truncation tested **[Z1]**. (H) remains **sufficient** (v2.17 Theorem 1), and a violation remains **necessary** for non-uniqueness, but violations can obstruct one another |
| — | How common is the obstruction? | **Rare.** On a random sample of 12 cyclic 4-vertex quivers, **27 of 28** (H)-violating vertex folds carry genuine families (exact fibres) **[Z3]** |
| — | Shape of a family | In the basic example (0⇄1, 0→2) the family is **a single vertex-level constant on the violating path**, a(1→0→2) = c·e₁. Every other deviation from the canonical fold is forced through Δ **[Z2]** |
| — | Dimension formula | **Open.** The naive count (|S_c| per violating head–block pair) fits 19 of 28 cases. It overcounts when violations share a closed walk and undercounts in the larger 5-dimensional fibre |

---

## 1. The counterexample

**Quiver:** 0 → 2, 0 → 3, 2 → 0, 3 → 1, a 2-cycle 0 ⇄ 2 with a pendant arrow at each end.

**Vertex fold:** e₀ ↦ e₀, e₁ ↦ −e₀, e₂ ↦ e₂, e₃ ↦ −e₂. Vertex 1 collapses onto 0, and vertex 3 collapses onto 2.

**(H) is violated twice.** In Podge-block terms, 1 is a head whose block contains 0, and 3 is a head whose block contains 2. The paths 0→3→1 (child 0 to head 1) and 2→0→3 (child 2 to head 3) each close a walk through the 2-cycle. The canonical images are T(0→3→1) = 0→2→0 and T(2→0→3) = −(2→0→2).

**Exact result [Z1]:**

| truncation L | linear nullities by length | exact fibre |
|---|---|---|
| 2 | 0, 1 | 1 point (canonical) |
| 3 | 0, 1, 0 | 1 point |
| 4 | 0, 1, 0, 1 | 1 point |

The length-2 equations alone leave one free direction, as both violations predict. Δ-compatibility for length-3 paths running through both violations, such as 2→0→3→1, then kills it. The constraint is linear, since the a's of shorter paths are already zero. **The two violations cancel each other's freedom.** In cyc A (v2.17) only one violation is present and the family survives. In cyc C two violations share a head and both survive (dimension 2). Here the two violations sit at different heads, interlocked by the same 2-cycle.

**Corrected statement.** For loopless quivers without parallel arrows:

  (H) ⇒ unique extension ⇒ ¬(an unobstructed violation),

and the gap between "violation" and "unobstructed violation" is real. Both ends are now exact:
- the first implication is v2.17 Theorem 1;
- non-uniqueness is decided by the full linear system at each length, *including* Δ-compatibility one length up, because the fibre is a cone.

So **uniqueness is decidable by finite linear algebra length by length.** A clean combinatorial characterization is still open.

---

## 2. Families and their dimensions

**Shape (cyc A) [Z2].** Fix the free coordinate to 1 and solve the exact fibre. The only vertex-level correction is on the violating path itself, P(βγ) = c·e₁ − βα. All longer paths deviate from the canonical fold only through the Δ-recursion. In this example, then, a violation contributes **one constant on the path that closes the walk**, and nothing else is free.

**The 4-vertex sample [Z3].** 12 random cyclic 4-vertex quivers contain 28 (H)-violating vertex folds; full fibres at L = 3:

| fibre dimension | 0 | 1 | 2 | 3 | 5 |
|---|---|---|---|---|---|
| violating folds (count) | 1 (the counterexample) | 21 | 3 | 2 | 1 |

The pair count fits 19 of 28. It fails in both directions:
- **overcount:** for example, pairs (1,3), (3,1) and the closed pair (3,3) at one head give 3, but the fibre is 2;
- **undercount:** a fold with two pairs at one head (|S_c| = 2 each) gives 4, but the fibre is 5.

The correct count must see how violating walks overlap, which is exactly what produces the obstruction in §1.

**Gloss.** *When a chaotic instant's descendants can circle back to it, the split usually gains a free constant, but not always: two such circles that pass through each other can each pin the other's constant, and the split stays rigid. Whether the clock can drift is decided step by step, by finite bookkeeping, not by the mere presence of circles.*

---

## 3. Proposed amendments (v2.18 → v2.19)

| amendment | content |
|---|---|
| **Q2¹⁹ Necessity** | v2.17 Conjecture 5 is **withdrawn**. (H) is sufficient; violations are necessary for non-uniqueness but can obstruct one another. Uniqueness is decided by the length-wise linear system including Δ one length up |
| **Q3¹⁹ Families** | A single violation contributes one vertex constant on the violating path (cyc A). The dimension count is open |

---

## 4. Open questions

1. **Combinatorial criterion.** Characterize the unobstructed violations, e.g. by a condition on how the closed walks of different heads overlap.
2. **Dimension formula** for the fibres, consistent with the 28-case table.
3. Carried: parallel classes separated by arrows (1⇉2→3); Kronecker plus a loop; zigzag closed form; dihedral = C_L ∪ C_R; ρ.

---

## 5. Sources

- M. Aguiar, *On the associative analog of Lie bialgebras*, J. Algebra 244 (2001).
- W. Decker, G.-M. Greuel, G. Pfister and H. Schönemann, *Singular* (via passagemath-singular).

*"A chaotic instant whose descendants circle back usually lets the split drift by a constant; two such circles running through each other can pin each other, and the split stays rigid."*
