# Discordianism: Open Questions Explored (v2.16)

*Answers question 1 of `discordianism-open-questions-explored-v2.15.md` §4, the spanning lemma, and with it **Conjecture 4**. It then draws the consequences.*

*Tags **[Wn]** are machine checks in `verify_v2_16.py` (folder `v2_16/`, zipped as `verify_v2_16.zip`; `pip install -r requirements.txt`, and `SLOW=1` for the 5-vertex sample and n = 7). Labels: **proved**, **verified (exact)** or **conjectured**. Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Summary

| # | Question | Answer |
|---|---|---|
| v2.15 Q1 | Spanning lemma, hence Conjecture 4 | **Proved.** For every finite acyclic quiver without parallel arrows, **every Δ-compatible fold of kQ is the canonical extension P₀^{⊗•} of its vertex part.** The folds therefore correspond bijectively to the Cayley forests that satisfy the arrow rule: they are finite and rigid. Every step of the proof is machine-checked on named quivers, on all 4-vertex DAGs and on 25 random 5-vertex DAGs **[W1–W3]** |
| — | Consequences | v2.12's rigidity conjecture is a **theorem** without parallel arrows. Fold counts are now theorems: chain a(n), zigzag, and **stars: (n+1)^{n−1} + 2(n+1)^{n−2} + 2^{n−1} − 1** (proved; checked **[W5]**). The compatible fraction tends to **1/(2e)** for stars but decays geometrically for chains |
| — | Loops; parallel arrows; dihedral; ρ | Carried |

---

## 1. The theorem

**Theorem 1 (proved).** Let Q be a finite acyclic quiver without parallel arrows, and A = kQ with arrow deletion Δ. Then **P ↦ P|_{A₀} is a bijection** from the Δ-compatible folds of A to the vertex folds P₀ of k^{Q₀} satisfying the arrow rule: supp P₀(e_s) × supp P₀(e_t) ⊂ Q₁ for every arrow s → t whose two sides are both nonzero. The inverse is P₀ ↦ P₀^{⊗•}|_A.

The proof assembles the results of v2.13–v2.15 with one new step.

**Notation.** The vertex fold is a Cayley forest (v2.5), with Hodge blocks F1 and Podge blocks F2. P₀(e_v) = ±1_{S_v}, where S_v = subtree(v) if v is Hodge and S_v = subtree(v) ∖ {v} if v is Podge (v2.15). For a Podge vertex c with children ch(c), its block is C = {c} ∪ ch(c), and

  e_c = 1_C + h(c),  h(c) := −Σ_{d ∈ ch(c)} e_d,  with π_H(h(c)) = π_H(e_c).

For a Hodge vertex, h(v) := e_v.

**Induction** (v2.13). Suppose P = T on lengths < m. Then P = T + a on A_m, with a: A_m → A₀. Write H_m = T(A_m) and Q_m = (1 − T)(A_m).

1. **Hodge half** (v2.14): a|_{H_m} = 0.
2. **Bimodule reduction** (v2.15): ψ = a|_{Q_m}: Q_m → H₀ satisfies ψ(q·u) = π_H(q·ψ(u)) on both sides. So ψ kills the mixed block components. It also kills every product of two or more elements of positive length in ker T, because such a product lies in ker P.
3. **Diagonal components (new).** Let f ∈ 1_X·Q_m·1_X.
   - The paths of f through an interior Podge leaf (S_v = ∅) factor there into two pieces killed by T. Call that part z, so z is a product.
   - The rest, n = f − z, lies on the set N of paths whose vertices all have S_v ≠ ∅.
   - **Endpoints are Hodge** (P3). If one endpoint were c and the other a child d of c, then T's support would contain a closed walk d → ⋯ → d, since S_c ∋ d, S_d ∋ d and the arrow rule links consecutive S-sets. That contradicts acyclicity. The (1−ι) part consists of Hodge vertices.
   - **Hodge expansion.** Write H(p) = ⊗_j h(v_j). Because π_H is injective on Hodge-vertex tensors (triangular: π_H(e_h) = 1_{B_h} + deeper terms), Tn = 0 ⇔ H(n) = 0. Hence **n = Σ λ_p (p − H(p))**.
   - **Runs** (P4). Cut p into maximal runs of consecutive Podge vertices. Expanding ⊗e_c = Rdiff + ⊗h(c) on each run writes p − H(p) as a sum of terms, each containing at least one Rdiff = ⊗e_{c_j} − ⊗h(c_j). Factor each term at its first such run:

     **Lpart ⊗ Rdiff ⊗ Rpart = (Lpart ⊗ 1_{C₁}) · Rdiff · (1_{C_r} ⊗ Rpart).**

     All three factors are path-supported. The neighbours of a run are Hodge vertices u, w with u ∈ S_u, so the arrow rule gives u → d and d → w for every child d. Children-to-children arrows inside the run also come from the arrow rule. All three factors are killed by T: the outer two carry a block slot 1_C, and Rdiff because π_H(h(c)) = π_H(e_c).

   So n is a sum of products of positive-length elements of ker T, and ψ(n) = 0. ∎

The step that was missing in v2.15 is the factorization at *block* slots: put the whole 1_C at the run boundaries, not split at Hodge vertices. A factorization at Hodge vertices leaves H·Q products, which the fold equations do not constrain.

**Machine checks.** The checker verifies four steps exactly for every arrow-rule vertex fold:
- **(P1)** the forest data;
- **(P2)** that H(p) is path-supported;
- **(P3)** that diagonal endpoints are Hodge;
- **(P4)** that the run terms are products of admissible factors and sum exactly to p − H(p).

It covers six named quivers and a hand-built run of two consecutive Podge vertices (x→c→c′→y with children d→d′) **[W1]**, all 30 nonempty DAGs on 4 vertices **[W2]**, and 25 random 5-vertex DAGs, 15840 vertex folds **[W3]**. An independent linear certificate agrees **[W4]**. A matrix-free counter built only from the explicit S_v description reproduces the matrix-based counts **[W6]**.

**Gloss.** *On any branching history without loops or parallel steps, a split of order from chaos is decided at the instants. A run of chaotic instants is replaced by the chaos of their direct descendants, and that substitution is itself a product of chaotic pieces, so it cannot move the split.*

---

## 2. Consequences

**Corollary 2 (rigidity).** Without parallel arrows, Folds(kQ) is a finite set of points, whatever Aut_Δ is (v2.13 found it trivial on every example). So v2.12's rigidity conjecture is a theorem for these quivers. The parallel-arrow case (Kronecker; Conjecture 3′) remains open in general, though it is exact for the 2- and 3-Kronecker quivers.

**Corollary 3 (counts).** |Folds(kQ)| is the number of arrow-rule Cayley forests:

| family | n = 1…6 |
|---|---|
| linear Aₙ | 2, 6, 23, 100, 563, 3734 (= v2.11's a(n), with the formula and terms up to a(26)) |
| zigzag Aₙ | 2, 6, 27, 123, 669, 4371 |
| source star = sink star | 2, 6, 27, 182, 1743, 21640 |

**Proposition 4 (stars; proved; checked by enumeration for n ≤ 7, where n = 7 gives 327743 [W5]).**

  **|Folds(kStarₙ)| = (n+1)^{n−1} + 2(n+1)^{n−2} + 2^{n−1} − 1  (n ≥ 2).**

*Proof.* Let the centre be 0. Every arrow joins 0 to a leaf, and the arrow rule S₀ × S_i ⊂ Q₁ forces S₀ = {0} whenever S₀ and some S_i are both nonempty. So the forests split into three disjoint classes:
- **0 is a Podge leaf (S₀ = ∅).** Attach 0, as a leaf, to R₂ or to a Hodge vertex of a forest on the other n − 1 points. By the Hodge ↔ Podge symmetry, the average number of Hodge vertices is (n−1)/2. This gives 2(n+1)^{n−2}·(1 + (n−1)/2) = (n+1)^{n−1}.
- **S₀ = {0}, i.e. 0 is a Hodge leaf, and no S_i ∋ 0.** This means 0 has no proper ancestor, so it hangs directly from R₁: 2(n+1)^{n−2} forests.
- **S₀ ∉ {∅, {0}}.** Then every leaf i has S_i = ∅, i.e. is a Podge leaf. So 0 is Hodge and each i hangs from R₂ or from 0, with at least one child at 0: 2^{n−1} − 1 forests. ∎

**The fraction tends to a constant for stars.** Dividing by all 2(n+2)^{n−1} splittings gives **1/(2e) ≈ 0.184** as n → ∞. Chains instead lose a factor of about 0.40 per vertex (v2.12). Branching keeps compatibility alive; long lines exhaust it.

**Gloss.** *A star of histories radiating from one instant keeps about one split in 2e compatible, however many rays it has. A single long line of history loses more than half its compatible splits with each new instant.*

---

## 3. Proposed amendments (v2.15 → v2.16)

| amendment | content |
|---|---|
| **Q1¹⁶ Uniqueness** | **Theorem**: on acyclic quivers without parallel arrows, folds = canonical extensions of arrow-rule vertex folds |
| **Q1¹⁶′ Rigidity** | Theorem in the absence of parallel arrows. Conjecture 3′ (finitely many Aut_Δ-orbits) remains for parallel arrows |
| **P2¹⁶ Counts** | Fold counts = arrow-rule forest counts. Stars: closed form, limiting fraction 1/(2e) |

---

## 4. Open questions

1. **Parallel arrows (Conjecture 3′).** There ker Δ ⊋ A₀ (differences of parallel arrows), so P(A₀) ⊂ ker Δ, and the vertex reduction needs a "kernel reduction". Is every fold an Aut_Δ-translate of a canonical extension from ker Δ?
2. **Cycles without loops.** The proof uses acyclicity only in the endpoint step. The 2- and 3-cycles are rigid on truncations (v2.14). Is the theorem true for all loopless quivers?
3. **Looped catalogue** (v2.15): block collapse, spectators.
4. **Asymptotics**: the zigzag and chain constants. Is there a closed form for zigzags like the one for stars?
5. Carried: dihedral = C_L ∪ C_R.

---

## 5. Sources

- M. Aguiar, *On the associative analog of Lie bialgebras*, J. Algebra 244 (2001).
- I. Assem, D. Simson and A. Skowroński, *Elements of the Representation Theory of Associative Algebras*, Vol. 1, CUP, 2006.
- R. P. Stanley, *Enumerative Combinatorics*, Vol. 2, CUP, 1999 (rooted forests; Cayley's formula).
- W. Decker, G.-M. Greuel, G. Pfister and H. Schönemann, *Singular* (via passagemath-singular).

*"On a branching history without loops or parallel steps, a split of order from chaos is decided entirely at the instants. A star of histories keeps one split in 2e compatible; a single line keeps ever fewer."*
