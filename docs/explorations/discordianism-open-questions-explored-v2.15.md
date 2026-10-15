# Discordianism: Open Questions Explored (v2.15)

*Answers questions from `discordianism-open-questions-explored-v2.14.md` §4: the Podge half of unique extension (Q1) and the classification of folds on looped quivers (Q2).*

*Tags **[Vn]** are machine checks in `verify_v2_15.py` (folder `v2_15/`, zipped as `verify_v2_15.zip`; `pip install -r requirements.txt`, and `SLOW=1` adds the exhaustive 4-vertex sweep). Labels: **proved**, **verified (exact)** or **conjectured**. Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Summary

| # | Question | Answer |
|---|---|---|
| v2.14 Q1 | Podge half of Conjecture 4 | **Reduced to one lemma, and proved for chains.** The Podge correction ψ is a Q₀-bimodule map (proved), so mixed block components die. At length 1 the diagonal components vanish for every loopless quiver (proved), so **uniqueness holds at length 1 in general**. At length ≥ 2 everything reduces to a *spanning lemma*, verified exhaustively on all 4-vertex DAGs **[V1, V3]**. For the linearly oriented chain the lemma is **proved**, so **\|Folds(kAₙ)\| = a(n) for every n** **[V2]** |
| v2.14 Q2 | Classification on looped quivers | **Richer than loop-collapse.** On 2 vertices every fold family is a loop-collapse family or an automorphism orbit **[V4]**. On 3 vertices three new shapes appear: block collapse, unit spectators, and collapse through a spectator. One of them carries a **genuine modulus with no 2-cycle** **[V5]**. This supersedes v2.14's claim that the 2-cycle-with-loop offset is the only modulus. Conjecture 3″ survives |
| — | Framed folds | **Located.** In the 2-loop quiver (= k⟨a,u⟩) the framed folds of v2.4–v2.5 are the fibre over P₀ = 0, and their moduli are a pure parallel-arrow effect |

---

## 1. Unique extension: the Podge half

**Setting** (v2.14). Q is loopless with no parallel arrows, and a|_H = 0 is already proved. The unknown is ψ = a|_{Q_m}: Q_m → H₀. It satisfies:
- ψ(q·u) = q⋆ψ(u) and ψ(u·q) = q⋆ψ(u) for q ∈ Q₀;
- ψ(Q_i·Q_j) = 0 for i, j ≥ 1;

where **q⋆h := π_H(qh)**.

**Explicit forest data (proved).** Write P₀(e_v) = ±1_{S_v}. Then **S_v = subtree(v) for a Hodge vertex and subtree(v) ∖ {v} for a Podge vertex**, so S_v = ∅ exactly at Podge leaves. This follows from the alternating expansion e_v = Σ_{x ∈ subtree(v)} (−1)^{d(x)−d(v)} 1_{block(x)}.

**Lemma A (bimodule reduction; proved).**
- ⋆ is a Q₀-module structure on H₀: π_H(q·π_H(y)) = π_H(qy), because q·(y − π_H y) ∈ Q₀.
- So ψ is a Q₀-bimodule map, and Q₀ ≅ k^{F2} is semisimple.
- Let ι be the unit of Q₀ and 1_∅ = 1 − ι. Then Q_m = ⊕_{X,Y ∈ F2 ∪ {∅}} 1_X·Q_m·1_Y.
- **ψ vanishes on every mixed component (X ≠ Y)**: 1_C·1_{C′} = 0, and 1_C·ι = 1_C.

**Lemma B (length 1; proved for loopless quivers).** Both diagonal components of Q₁ are zero:
- **1_C·Q₁·1_C = 0.** Write C = {c} ∪ ch(c). An arrow c → d with d a child would need the loop (d, d), since S_c ∋ d and S_d ∋ d. The kernel of π_H on span{e_x : x ∈ C} is k·1_C. Hence an element of Q₁ supported on C × C is of the form 1_C ⊗ g + g′ ⊗ 1_C, and the missing arrows force it to be 0.
- **(1−ι)·Q₁·(1−ι) = 0.** Its vertices are the Hodge children of the Hodge root, and π_H is injective on them.

So ψ = 0 on Q₁. With the Hodge half, **every fold of a loopless quiver without parallel arrows agrees with the canonical extension through length 1**.

**Lemma C (endpoints; proved).** In a diagonal component at length ≥ 1 of an acyclic quiver:
- **the endpoints are never Podge leaves**, because a Podge leaf's Podge block is a singleton and the two endpoints are distinct;
- **the first slot never carries the whole block 1_C**: a vertex of C that starts a path to C cannot also end one, by acyclicity, and π_H is injective on any proper subset of C.

**Spanning lemma (verified).** For m ≥ 2, 1_X·Q_m·1_X ⊂ Σ_{0<i<m} Q_i·Q_{m−i}.
- It holds on every vertex fold of A₃, 1→2←3, A₄, 0→1→2←3, D₄ and the diamond **[V1]**.
- It holds on all 30 nonempty DAGs on 4 vertices **[V3]**.
- It holds on the 5-vertex configuration built to break the "swap across a block" step (0→2, 1→3, 2→4, 3→4; 581 vertex folds).

The spanning lemma, the Q·Q equations and Lemmas A–B together give Conjecture 4.

**Theorem 1 (proved): the chain.** For the linearly oriented Aₙ, every fold is the canonical extension of its vertex part. Hence

  **|Folds(kAₙ)| = a(n) for every n**,

the sequence 2, 6, 23, 100, 563, … whose polynomial-time formula is in v2.11.

*Proof of the spanning lemma on the chain.*
1. Let f ∈ 1_X·Q_m·1_X and split it as f = z + n, where z collects the paths with an interior Podge leaf.
2. A path through an interior vertex v with S_v = ∅ factors at v into two pieces that T kills. So z ∈ Q·Q, and Tn = −Tz = 0.
3. Every vertex of a path in n has S ≠ ∅: interior vertices by construction, endpoints by Lemma C.
4. On the chain, the arrow rule with both S nonempty forces adjacent singletons S_j = {k}, S_{j+1} = {k+1}. So T sends each path of n to ± a single path.
5. In a block C, endpoints (c, d) with d a child would force S_c = S_d = {d}, i.e. a closed path, which the chain does not have. Endpoints (d, d′) are Hodge leaves, and then T fixes the path. In the (1−ι) part the endpoints are Hodge leaves under the root, and again T fixes the path.
6. So T is injective on n's span, hence n = 0. ∎

The intermediate claim in step 2, "every diagonal element is carried by paths through an interior Podge leaf", is machine-checked on A₃–A₅ **[V2]**.

**Consequence.** The polynomial-time formula of v2.11, and its terms up to a(26) = 9706864926713557355624856470, **count Δ-compatible folds of the path algebra kAₙ**. On kAₙ these folds form a finite set, one fold per arrow-compatible Cayley forest.

**Gloss.** *On a single line of history, every way of splitting order from chaos is decided at the instants, and the Church's count of them, carried since v2.5, is now a theorem. On branching histories the same holds through the first step and is checked everywhere it has been examined.*

---

## 2. Looped quivers: the classification is richer

**Method.** For every vertex fold that passes the arrow rule, I computed the exact truncated fibre (all equations involving a_{≤L}, with L = 3). Its dimension is compared with a predictor for *loop-collapse with spectators*:
- the collapse of v2.14 Theorem 3 acts on {s, s′};
- the other vertices keep their own vertex fold, with no component on s;
- the predicted dimension is the number of arrows among {s, s′} other than the loop at s.

**Two vertices [V4].** On all 7 looped quivers every positive fibre matches the predictor, with one exception: the complete digraph with both loops. There, the vertex fold with H₀ = k·1 carries a 1-dimensional fibre. Its canonical Hodge part is k[w], where w is the sum of all four arrows. Since Δw = 1 ⊗ 1, the whole quiver behaves like a single looped vertex. The fibre is the orbit of the simultaneous translation of both loops (fibre 1 = orbit 1). The collapse fibres of the complete digraph have dimension 3 against orbit 2, one genuine modulus, as the 2-cycle rule of v2.14 predicts.

**Three vertices.** Six quivers were computed, and three of them are fully explained by the predictor: loop at 0 with the chain 0→1→2, the chain with a loop at the middle vertex, and loops at both 0 and 1 with arrows 0→1→2. Three new shapes remain:

| quiver | vertex fold (P₀ on e₀, e₁, e₂) | fibre | shape |
|---|---|---|---|
| loop at 0, arrows 0→1, 0→2 (and its dual) | e₀, −e₀ − e₂, e₂ | **2** | **block collapse**: 1 folds onto the block {0, 2} |
| same | e₀, e₁, −e₀ − e₁ | **2** | block collapse onto {0, 1} |
| same | e₀, −e₀, 1 (and its mirror) | 1 | collapse with a **unit spectator** |
| 2-cycle 0⇄1 with loop at 0, tail 1→2 | e₀, e₁, −e₀ | 2 | collapse of 2 **through the spectator** 1 |

**A modulus without a 2-cycle [V5].** For the loop source fork, the filtered Δ-derivations are 1-dimensional: the translation of the single loop. Each block-collapse fibre has **dimension 2 against orbit dimension 1**, so it carries a genuine modulus although the quiver has no 2-cycle. Concretely, α ↦ −(u + β) + c₀e₀ + c₂e₂, and the corrections are geometric in c₀ − c₂: the coefficient of u^k·α at vertex 2 is c₂(c₀ − c₂)^k. Grading rescaling adds one more dimension of symmetry, so **Conjecture 3″ (finitely many orbits of Aut_Δ ⋊ grading) survives**. v2.14's sentence "the 2-cycle offset is the only modulus found" is superseded.

**Framed folds, located.** The quiver with one vertex and two loops is k⟨a,u⟩ (v2.14 Proposition 2), and its vertex algebra is k. So P₀ ∈ {0, 1}:
- P₀ = 1 gives the unital-Hodge folds of v2.6–v2.10 (identification, dihedral);
- **P₀ = 0 gives the framed folds of v2.4–v2.5.**

Without parallel arrows the fibre over P₀ = 0 is the single fold P = 0, so the framed edge moduli are a pure parallel-arrow effect.

**Gloss.** *A loop lets a neighbouring instant collapse onto the repeating one. With several branches leaving the loop, a neighbour can instead collapse onto the loop and one of its branches together, and that leaves an offset no symmetry absorbs, without any cycle to explain it. The Church's catalogue of such collapses is open, but every case examined is finite once the clock's scale is allowed to change.*

---

## 3. Proposed amendments (v2.14 → v2.15)

| amendment | content |
|---|---|
| **P2¹⁵ Counting** | **a(n) = \|Folds(kAₙ)\|** (theorem). v2.11's formula counts folds |
| **Q1¹⁵ Uniqueness** | Bimodule reduction and length-1 uniqueness proved for all loopless quivers. General length reduced to the spanning lemma (verified; proved on chains) |
| **Q2¹⁵ Loops** | Loop-collapse with spectators, total loops (H₀ = k·1 on a complete looped digraph), block collapse, unit spectators. Moduli also arise without 2-cycles. Conjecture 3″ stands |

---

## 4. Open questions

1. **The spanning lemma for all acyclic quivers.** The chain proof used the fact that on a line the arrow rule forces singleton S-sets. In general the S-sets are subtrees, and the lemma needs a swap argument across the blocks 1_C of non-leaf Podge nodes. The configuration built to break that swap still satisfies the lemma.
2. **A catalogue of looped folds.** Find the general shape behind block collapse and spectators, plausibly collapse onto a "looped Hodge block" with sign data, and count moduli modulo Aut_Δ.
3. **Parallel arrows plus loops**, where the framed and identification moduli live.
4. Carried: dihedral = C_L ∪ C_R; ρ.

---

## 5. Sources

- M. Aguiar, *On the associative analog of Lie bialgebras*, J. Algebra 244 (2001).
- I. Assem, D. Simson and A. Skowroński, *Elements of the Representation Theory of Associative Algebras*, Vol. 1, CUP, 2006.
- W. Decker, G.-M. Greuel, G. Pfister and H. Schönemann, *Singular* (via passagemath-singular).
- R. P. Stanley, *Enumerative Combinatorics*, Vol. 2, CUP, 1999 (Cayley forests, used through v2.5's bijection).

*"On a single line of history the count of splits is a theorem. With loops and branches, chaos finds offsets that no symmetry absorbs, yet the catalogue stays finite once the clock's scale is allowed to change."*
