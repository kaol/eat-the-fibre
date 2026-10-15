# Discordianism: Open Questions Explored (v2.14)

*Answers questions from `discordianism-open-questions-explored-v2.13.md` §4: unique extension (Q2) and quivers with loops (Q3). It also places the original k⟨a,u⟩ framework inside the quiver picture.*

*Tags **[Un]** are machine checks in `verify_v2_14.py` (folder `v2_14/`, zipped as `verify_v2_14.zip`; `pip install -r requirements.txt`, and `SLOW=1` for the larger cases). Labels: **proved**, **verified (exact)** or **conjectured**. Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

*Status of v2.13: your `SLOW=1` run of `verify_v2_13.py` gave 10/10, so the slow checks are now independently confirmed on a second machine.*

---

## 0. Summary

| # | Question | Answer |
|---|---|---|
| v2.13 Q2 | Unique extension (Conjecture 4) | **Half proved.** The fold's correction vanishes on the Hodge part H_m for every quiver without loops (Theorem 1). The Podge half is verified, and the uniqueness system is shown to split into exactly these two stages **[U1]**. Random acyclic quivers: certified on all 37 tested **[U2]** |
| v2.13 Q3 | Quivers with loops | **The original framework is a quiver.** k⟨a,u⟩ with a-deletion ≅ the path algebra of the one-vertex 2-loop quiver, via u = α − β **[U3]**. **Cycles without loops stay rigid** **[U4]**. **Loops alone create families**: the *loop-collapse folds* (Theorem 3, proved; negative controls **[U6]**), whose dimension matches every looped fibre computed **[U5]**. Loop translations are Δ-automorphisms, and modulo them these families are mostly single orbits. One exception, the **2-cycle with a loop, carries a genuine 1-dimensional modulus** **[U7]** |
| — | Dihedral; ρ; framed folds | Carried |

---

## 1. Unique extension: the Hodge half

**Setting** (v2.13). Q has no parallel arrows. A fold P has vertex part P₀, and by induction on length P = T + a on A_m, where T = P₀^{⊗•} is the canonical extension and a: A_m → A₀. Write H_m = T(A_m) and Q_m = (1 − T)(A_m). Write F1 and F2 for the Hodge and Podge blocks of the Cayley forest of P₀, so that H₀ = span{1_B : B ∈ F1}, with the blocks in F1 pairwise disjoint.

**Theorem 1 (proved; structure checked [U1]).** If Q has no loops (cycles are allowed), then a|_{H_m} = 0 for every m.

*Proof.*
- **Length 1.** H₁ = (P₀ ⊗ P₀)(A₁) ⊂ H₀ ⊗ H₀. Since the F1-blocks are disjoint, H₁ is spanned by the elements z = 1_B ⊗ 1_{B′} with B × B′ ⊂ Q₁. With no loops, B ∩ B′ = ∅, so B ≠ B′. Let w = a(z). The H₀-multiple equations are a(h·Tz) = π_Q(h·a(z)) and a(Tz·h) = π_Q(a(z)·h). Apply them with h = 1_{B′}:
  - on the left, 1_{B′}·z = 0, so 0 = π_Q(1_{B′} w);
  - on the right, z·1_{B′} = z, so w = π_Q(w·1_{B′}).

  Since A₀ is commutative, w = 0.
- **Length m ≥ 2.** H_m is spanned by the elements 1_{B₀} ⊗ ⋯ ⊗ 1_{B_m} with each B_i × B_{i+1} ⊂ Q₁. Each of these is a glued product (1_{B₀} ⊗ ⋯ ⊗ 1_{B_i})·(1_{B_i} ⊗ ⋯ ⊗ 1_{B_m}) of two elements of H_i and H_{m−i}, because 1_{B_i}² = 1_{B_i}. The H·H equation a(Tx·Ty) = 0 kills it. ∎

**The Podge half** (verified, not yet proved). Once a|_H = 0, the remaining unknown is ψ = a|_{Q_m} : Q_m → A₀. The following equations force ψ = 0 on every quiver tested **[U1]**:
- the H₀-conditions π_Q(h·ψ(u)) = 0;
- the Q₀-multiples ψ(q·u) = π_H(q·ψ(u)), with right-hand versions;
- the Q·Q products ψ(u′u″) = 0.

Equation-family experiments locate the work:
- at length 1, the vertex multiples alone suffice;
- at length ≥ 2, the Q₀-multiples together with Q·Q suffice;
- dropping idempotency changes nothing;
- dropping the path products leaves nullity 2–8.

**Random quivers [U2].** 25 random acyclic 4-vertex quivers and 12 random 5-vertex quivers without parallel arrows are all certified, at every length.

*Correction during this round.* A first random run reported nullity 3 on the orientation 0→1→2←3. The cause was a default argument in my code that silently dropped the path-product equations. With the full system that orientation is certified, and its full polynomial fibre is a single reduced point over each of its 109 vertex folds.

**Gloss.** *Where history has no loops, the orderly half of a split is pinned at the instants. The chaotic half is pinned too, in every case the Church has examined, but nobody has yet written down why.*

---

## 2. Loops, cycles, and where the old framework lives

### 2.1 The original framework is a quiver

**Proposition 2 (proved; checked [U3]).** Let K₂° be the quiver with one vertex and two loops α, β, with arrow deletion as Δ. The map

  **φ : k⟨a, u⟩ → kK₂°, a ↦ α, u ↦ α − β**

is an isomorphism of ε-bialgebras from k⟨a,u⟩ with a-deletion (the setting of v2.1–v2.10). Indeed φ is multiplicative and unitriangular, and Δφ(u) = 1⊗1 − 1⊗1 = 0 = (φ ⊗ φ)Δ_a(u). Two ε-derivations that agree on generators agree everywhere.

So the infinite-dimensional moduli found earlier are the folds of a quiver that has **both** loops **and** a Δ-primitive difference of parallel arrows:
- the identification classification R ∪ L (v2.8);
- the one-sided dihedral families (v2.9–v2.10).

v2.12–v2.13 showed that each ingredient alone is tame for acyclic quivers. The sections below take loops without parallel arrows.

### 2.2 Cycles without loops are rigid

- Theorem 1 covers all loopless quivers, cycles included.
- The full truncated fibres are single points: the 2-cycle at L = 2…5 and the 3-cycle at L = 3, 4 **[U4]**.

*Truncation note.* A fibre computed to length L must include Δ-compatibility for paths of length L + 1, since that is where a_L first appears. Without it, the linear system shows nullity 1 at every even length on the 2-cycle. That spurious nullity is exactly the kind of window artifact v2.9 warned about.

### 2.3 Loops: the loop-collapse folds

**Theorem 3 (proved; checked with symbolic parameters, with negative controls [U6]).** Let s be a vertex with a loop u, and let s′ ≠ s be another vertex. Define:
- σ(s) = 1, σ(s′) = −1, and σ = 0 elsewhere;
- a weight w(b) = u − c_b ∈ k[u] for each arrow b between vertices of {s, s′}, with w(u) = u;
- for each path p, χ(p) = Π_{v ∈ p} σ(v) · Π_{b ∈ p} w(b) ∈ k[u], where vertices are counted with multiplicity.

Then **P(x) = χ(x)·e_s**, where k[u] ≅ k[u]e_s via the loop, is a Δ-compatible fold. Its Hodge part is k[u]e_s.

*Proof.*
1. **Idempotent.** χ(u^k e_s) = u^k.
2. **Δ-compatible.** Let p = prefix·b·suffix. The vertices of p are those of the prefix and those of the suffix, so χ(p) = χ(prefix)·w(b)·χ(suffix). The divided-difference identity for monic linear factors,

   (F(x) − F(y))/(x − y) = Σ_i Π_{j<i}(x − c_j) Π_{j>i}(y − c_j),

   turns Δ(χ(p)e_s) into Σ χ(prefix) ⊗ χ(suffix) = (P ⊗ P)Δp. Paths through a vertex with σ = 0 give 0 on both sides.
3. **Kernel is a subalgebra.** χ(xy) = Σ_v σ(v)·χ(x e_v)·χ(e_v y). With σ supported on {s, s′} with opposite signs, the conditions χ(x) = χ(y) = 0 give χ(xy) = (σ_s + σ_{s′})·X·Y = 0. ∎

**Each hypothesis is needed [U6].** The construction fails with:
- three supported vertices;
- equal signs;
- a non-monic weight 2u − c;
- a translated loop at s.

**The fibres are exactly these families [U5].** In every looped example, the positive-dimensional fibre has dimension equal to the number of parameters c_b, the arrows among {s, s′} other than u. This holds at L = 3, 4, and 5 or 6 where computed. Where it was examined, the fibre has a single component, and the closed form satisfies all fold equations with symbolic parameters.

| quiver | family | fibre dim (stable over every L computed, 3–6) |
|---|---|---|
| loop at 0, arrow 0→1 | α ↦ c − u | 1 |
| arrow 0→1, loop at 1 | α ↦ c − u | 1 |
| loops at 0 and 1, arrow 0→1 | α ↦ c − u, v ↦ u − c′ | 2 |
| 2-cycle 0⇄1 with loop at 0 | α ↦ c₁ − u, β ↦ c₂ − u | 2 |

**Orbits [U7].** Translating a loop, u ↦ u − c·e_s, is a Δ-compatible automorphism: Δ(u − c·e_s) = e_s ⊗ e_s = Δu. So Aut_Δ is bigger than v2.13's formula Σ m(m−1), which counted only acyclic quivers: each loop adds a translation. The computed filtered Δ-derivations are exactly these.

Conjugating a collapse fold by translations:
- translating u shifts every c_b together;
- translating a loop at s′ shifts its own parameter alone.

Hence the number of **moduli modulo Aut_Δ** is n_{s→s′} + n_{s′→s} − 1:
- loop + arrow: 0, a single orbit (orbit dimension 1 = fibre dimension 1);
- two loops + arrow: 0 (2 = 2);
- **2-cycle + loop: 1.** The invariant is c₁ − c₂, with fibre dimension 2 against orbit dimension 1.

The grading rescaling (x ↦ λ^{|x|}x) preserves Δ only up to a scalar, so it maps folds to folds without lying in Aut_Δ. It acts by c ↦ λ⁻¹c, which leaves two orbits: c₁ = c₂ and c₁ ≠ c₂.

**Gloss.** *A loop is a place where history can repeat. Loops let chaos fold a neighbouring vertex onto the repeating one, and identify each connecting step with "a constant minus one more repetition". Shifting the loop's origin absorbs most of those constants. A 2-cycle through a loop keeps one genuine difference: two routes out and back whose offsets no symmetry can reconcile.*

---

## 3. Proposed amendments (v2.13 → v2.14)

| amendment | content |
|---|---|
| **Q1¹⁴ Uniqueness** | Hodge half proved for all loopless quivers. Podge half verified; Conjecture 4 certified on 37 random quivers |
| **Q2¹⁴ Loops** | Loop-collapse folds (Theorem 3). Aut_Δ contains loop translations. Without parallel arrows the only genuine modulus found modulo Aut_Δ is the 2-cycle-with-loop offset, which becomes finite modulo grading |
| **M1¹⁴ Framework** | k⟨a,u⟩ = path algebra of the 2-loop quiver. Its rich moduli need both loops and parallel arrows |
| **Revised Conjecture 3″** | For quivers without parallel arrows, folds form finitely many orbits of the group of algebra automorphisms preserving Δ up to scalar (Aut_Δ ⋊ grading) |

---

## 4. Open questions

1. **Podge half of Conjecture 4.** ψ: Q_m → A₀ with the H₀-conditions, the Q₀-multiples and Q·Q. The obstacle is that Q_m is not spanned by glued products inside A, since the path constraint breaks the vertex-sequence factorization used for H_m.
2. **Classification with loops.** Are all folds of a loopy quiver without parallel arrows either canonical extensions or loop-collapse folds? Is Conjecture 3″ true?
3. **Parallel loops.** In the 2-loop quiver the moduli are infinite-dimensional (v2.8). Which minimal configuration first produces infinitely many moduli modulo Aut_[Δ]: two parallel arrows plus one loop, or a loop with a parallel pair at the same vertex?
4. Carried: dihedral = C_L ∪ C_R; ρ; framed folds. The loop-collapse folds have 1 ∈ Q, the same "unit in the Podge" feature as the framed folds of v2.4–v2.5, so they may be their quiver shadow.

---

## 5. Sources

- M. Aguiar, *On the associative analog of Lie bialgebras*, J. Algebra 244 (2001) (ε-bialgebras).
- I. Assem, D. Simson and A. Skowroński, *Elements of the Representation Theory of Associative Algebras*, Vol. 1, CUP, 2006 (path algebras).
- W. Decker, G.-M. Greuel, G. Pfister and H. Schönemann, *Singular* (via passagemath-singular).
- The divided-difference (Newton) identity: any text on interpolation, e.g. C. de Boor, *Divided differences*, Surveys in Approximation Theory 1 (2005).

*"Where history cannot loop, the split is pinned at the instants. A loop lets chaos fold a neighbour onto the repeating vertex, and a 2-cycle through the loop leaves one difference that no symmetry absorbs."*
