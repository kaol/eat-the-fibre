# Discordianism: Open Questions Explored (v2.17)

*Answers question 2 of `discordianism-open-questions-explored-v2.16.md` §4: does the uniqueness theorem extend from acyclic quivers to loopless quivers with cycles?*

*Tags **[Xn]** are machine checks in `verify_v2_17.py` (folder `v2_17/`, zipped as `verify_v2_17.zip`; `pip install -r requirements.txt`, and `SLOW=1` adds the looped and 4-vertex comparisons). Labels: **proved**, **verified (exact)** or **conjectured**. Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Summary

| # | Question | Answer |
|---|---|---|
| v2.16 Q2 | Is v2.16's theorem true for loopless quivers with cycles? | **No; the acyclic hypothesis is sharp.** The smallest counterexample is 0⇄1 with 0→2: three vertices, three arrows. It carries a **genuine one-parameter fold family**, and the quiver has no Δ-compatible automorphisms at all **[X1–X3]** |
| — | What exactly does the proof need? | **A sharper theorem (proved).** Uniqueness holds over every vertex fold satisfying **(H)**: no path through S-nonempty vertices joins a Podge head to a vertex of its own block. Acyclic quivers satisfy (H) automatically **[X4]** |
| — | Is (H) also necessary? | **Yes on all data (conjectured in general).** On all 10 loopless cyclic 3-vertex quivers, uniqueness fails exactly at the 4 (H)-violations **[X1]**; also on a random 4-vertex sample, where the certificate fails on exactly the 28 predicted vertex folds **[X6]**. With loops, every (H)-violation is again non-unique. The only other non-unique fibres are translation orbits through looped Hodge blocks **[X5]**. So the loop-collapse families of v2.14–v2.15 are (H)-violations: **loops and cycles break uniqueness by the same mechanism** |

---

## 1. The counterexample

Take Q: α: 0→1, β: 1→0, γ: 0→2 ("cyc A"). Take the vertex fold e₀ ↦ e₀, e₁ ↦ e₁, e₂ ↦ −e₁: vertex 2 collapses onto 1. It passes the arrow rule, and its canonical extension sends γ ↦ −α and βγ ↦ −βα.

**The family (verified exact, stable for L = 3, 4, 5 [X2]).** The fibre over this vertex fold is a line. Its only free coordinate is

  **P(βγ) = c·e₁ − βα,**

with every other correction determined. So the path 1→0→2 is identified with "constant minus the closed path βα at vertex 1": a 2-cycle playing the role of the loop in v2.14's loop-collapse. The other examples:
- cyc B (0⇄1, 2→0) gives the mirror family, dimension 1;
- cyc C (0⇄1, 0⇄2) gives two such families, dimension 2 (L = 3, 4).

**Not an orbit [X3].** The filtered Δ-compatible derivations of cyc A form the zero space (there are no loops to translate and no parallel arrows to shuffle). So the fibre dimension 1 against orbit dimension 0 is a **genuine modulus**. Grading rescaling acts by c ↦ λ⁻¹c, which leaves two orbits, c = 0 and c ≠ 0. So Conjecture 3″ survives.

*Method note.* The linear certificate reports positive nullity on exactly these 4 vertex folds among all arrow-rule vertex folds of the 10 cyclic quivers. That nullity is only a symptom, since it can come from omitting Δ-compatibility one length up (v2.14). The full fibres confirm that all 4 are genuine.

---

## 2. Where the proof breaks, and a sharper theorem

v2.16's proof used acyclicity in one step only: the endpoints of a diagonal Podge element are Hodge. Write C = {c} ∪ ch(c) for a Podge block, and N for the set of paths through vertices with S_v ≠ ∅. The step excludes endpoint pairs (c, d), (d, c) and (c, c) with d ∈ ch(c), because S_c ∋ d ∈ S_d would close a walk. In cyc A the vertex fold makes 2 a Podge head with child 1, and the path 1→0→2 is exactly such a pair. The walk it closes, 1→0→1, exists.

**Hypothesis (H).** No path of positive length in N joins a Podge head c to a vertex of its own block C.

**Theorem 1 (proved).** Let Q be loopless without parallel arrows, cycles allowed. Every vertex fold satisfying the arrow rule and (H) has exactly one Δ-compatible extension, the canonical one.

*Proof.* v2.16's argument goes through verbatim, with (H) in place of acyclicity in the endpoint step. The other steps never used acyclicity: the Hodge half (only loops matter), the bimodule reduction, Podge-leaf factorization, the Hodge expansion and the run factorization. The induction on length works for infinite-dimensional kQ because each A_m is finite-dimensional. ∎

**Corollary.** Acyclic quivers satisfy (H): an (H)-violating path would close a walk through S-sets. This recovers v2.16 **[X4]**.

**Conjecture 5 (necessity).** For loopless quivers without parallel arrows, a vertex fold has a non-canonical extension if and only if it violates (H).
- **Evidence:** all 10 cyclic 3-vertex quivers, where the (H)-violations, the closed-walk predictor and the non-unique fibres are the same 4 vertex folds **[X1, X2]**; and 12 random cyclic 4-vertex quivers **[X6]**, where the linear certificate fails on exactly the 28 vertex folds flagged by the closed-walk predictor and certifies uniqueness on all others. Those 28 were not individually confirmed with full fibres.
- **Mechanism:** a violating path p from a child d to its head c, or the reverse, lets P(p) acquire a vertex-level constant along the closed walk. That is the free coordinate c·e₁ above.

---

## 3. Loops, revisited

In the looped families of v2.14–v2.15, the collapsed vertex s′ becomes a Podge head whose block contains the looped vertex s. The arrows between them are (H)-violating paths. This covers loop + arrow, the 2-cycle with a loop, two loops + arrow, the complete digraph's collapse fibres, and the block collapse of the loop source fork. **On all looped 2-vertex quivers, every (H)-violation has a positive fibre [X5].**

With loops there is exactly one further source of positive fibres. A *looped Hodge block* (B × B ⊂ Q₁) breaks the length-1 Hodge step of the proof, and the result is the H₀ = k·1 "total loop" fibres of the complete digraph. v2.15 identified these as orbits of loop translations, so they are not moduli.

**Unified picture.** Without parallel arrows, a fold is the canonical extension of its vertex part unless a walk closes up inside a Podge block (an (H)-violation, possible through loops or cycles), or a Hodge block is looped (giving translation orbits).

**Gloss.** *A split of order from chaos is pinned at the instants unless some chaotic instant's own descendants can circle back to it. A loop or a cycle provides that circle, and along it the chaos may shift by a constant that no symmetry absorbs. On histories without circles, nothing can drift.*

---

## 4. Proposed amendments (v2.16 → v2.17)

| amendment | content |
|---|---|
| **Q1¹⁷ Uniqueness** | Proved under (H) for loopless quivers, cycles allowed. Acyclic implies (H) |
| **Q2¹⁷ Sharpness** | The smallest cyclic counterexample (0⇄1, 0→2) has a genuine modulus. Necessity of (H) is conjectured, with all data agreeing |
| **Q3¹⁷ Loops** | Loop-collapse families are (H)-violations. Looped Hodge blocks give translation orbits |

---

## 5. Open questions

1. **Prove Conjecture 5**: construct the family from an (H)-violating path in general.
2. **Dimension formula.** Pair counting (|S_c| per violating head–block pair) matches most cases but fails for the collapse through a spectator. What is the right count?
3. **Parallel arrows** (Conjecture 3′): kernel reduction with ker Δ ⊋ A₀.
4. Carried: zigzag closed form; dihedral = C_L ∪ C_R; ρ.

---

## 6. Sources

- M. Aguiar, *On the associative analog of Lie bialgebras*, J. Algebra 244 (2001).
- I. Assem, D. Simson and A. Skowroński, *Elements of the Representation Theory of Associative Algebras*, Vol. 1, CUP, 2006.
- W. Decker, G.-M. Greuel, G. Pfister and H. Schönemann, *Singular* (via passagemath-singular).

*"Order and chaos split rigidly unless a chaotic instant's descendants can circle back to it; a loop or a cycle supplies the circle, and along it the chaos may drift by a constant."*
