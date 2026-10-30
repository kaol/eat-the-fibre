# Discordianism: Open Questions Explored (v2.76): Theorem 5 for paths; Conjecture 3ᶜ′ refuted, 3ᶜ″

*v2.75 ended with a quiver its Theorem 5 did not cover: 1⇛2, 2→3→5, 2→4 grows (3417, 3423), but only one of the arrows
leaving 2 ends at a sink. This round finds the reason.*
- *Theorem 5's proof never needs the tails to be arrows. Any two distinct paths from j to sinks work, provided each
  path runs through arrows that are alone in their class. That is Theorem 5′.*
- *One consequence: 1⇛2→3, 3→4, 3→5 has infinitely many orbits, although only one arrow leaves the class's target.
  That refutes v2.75's Conjecture 3ᶜ′.*
- *The vertex part of the family is not unique. When i is a source, the vertex part can be the idempotent of
  {i} ∪ F, or of its complement, where F is any set of intermediate vertices on the two paths. That gives
  2^{|F|+1} disjoint families.*
- *The revision is Conjecture 3ᶜ″: the descendants of the target, and the ancestors of the source, of every class of
  size ≥ 3 must each form a chain of vertices.*

---

## 0. Results

The predictions are in logs/v2_76/PREDICTIONS.md. 127–130 were committed before any run. Addendum 1 (131–132) was
committed before any run of this round. Addendum 2 (133–134) was committed after the path-family checks and before
the censuses it is about.

| # | Statement | Status |
|---|---|---|
| 127 | 1⇉2 with three tails is flat: 3714 at q = 3 | **held** (3714, 3714) |
| 128 | the 6 extra classes of 1⇛2, 2→3→5, 2→4 sit in the identity block | **not tested** (needs a 66-minute census rerun with saved representatives). §2 accounts for at least 4 of the 6 there |
| 129 | a component with generic orbit rank < dimension, invariant pairing a3 and a4 | **superseded** by 131: the invariant pairs the paths a3·a5 and a4 |
| 130 | 1⇛2, 2→3, 2→4, 3→4 is not flat | **held**: 681 → 689 |
| 131 | the path version is a fold, and φ is invariant, on three quivers | **held** on all three, symbolically |
| 132 | 1⇛2→3, 3→4, 3→5 jumps from q = 2 to q = 3 | **pending**: 4899 at q = 2. The q = 3 run was killed for lack of memory (another census ran beside it); it is rerun after 133's run. Theorem 5′ already gives at least 4(q − 1) classes |
| 133 | 1⇛2→3⇉4 is flat | **pending**: 2236 at q = 2. At q = 3 it had not finished after 76 minutes (stuck in fibre group 3 of 6) |
| 134 | the jumps are at least the number of vertex-part families (4 and 4) | **held** for 1⇛2, 2→3, 2→4, 3→4 (+8 ≥ 4); pending for 1⇛2→3, 3→4, 3→5 |

## 1. Theorem 5′ (paths)

Q is acyclic, and C = {a₁, …, a_n} is a class i→j with n ≥ 3. Let π ≠ π′ be paths from j to sinks, and let N(u·a₁·ρ) =
u·(λ_ρ e_i + m_ρ ℓ₀) for ρ ∈ {π, π′}, with N = 0 on every other path. P = Π₀ + N, exactly as in Theorem 5.

**Theorem 5′.**
1. P is a fold over every commutative ring.
2. If every arrow of π and π′ is alone in its class, φ (the cross-ratio of δ_π, δ_π′ against the vertex and V lines)
   is invariant under the census symmetries, and equals λm′/(λ′m) on the family.
3. Hence Q has infinitely many orbits over every infinite field, and at least q − 1 census classes in one block.

*Proof.* v2.75 §2, word for word, with x replaced by π. Two checks need a sentence more:
- **Kernel products.** A path through a₁ meets C once and cannot return from j to i, so the same three cases occur.
  A product N(p)·κ(u·a₁·r′) with r′ ∉ {π, π′} is λ·κ(u·a₁·r′), because N(u·a₁·r′) = 0.
- **Δ, deleting an arrow inside π.** Write π = π₁βπ₂. The term is u·a₁·π₁ ⊗ π₂, and N(u·a₁·π₁) = 0 because π₁ is a
  proper prefix of π. A proper prefix cannot be π′, since π′ ends at a sink.

**Invariance.** S(π) = t_π·π, where t_π is the product of the weights of π's arrows. That needs each arrow to be alone
in its class. The rest of the argument is unchanged. ∎

**Checks** (`analysis/v2_76_pathfamily.py`, INV=1, generic symbolic S). In every case the fold identities hold
identically, and φ(SPS⁻¹) = φ(P) = λm′/(λ′m):

| quiver | π, π′ |
|---|---|
| 1⇛2, 2→3, 2→4 | b, b′ (Theorem 5 itself) |
| 1⇛2, 2→3→5, 2→4 | b·a5, b′ |
| 1⇛2, 2→3, 2→4, 3→4 | b′, b·a5 (two parallel paths, both ending at 4) |
| 1⇛2→3, 3→4, 3→5 | b·c, b·c′ (sharing their first arrow) |

**The single-class hypothesis matters.** On 1⇛2→3⇉4, take π = b·c₁ and π′ = b·c₂ through the 2-class c₁, c₂. P is
still a fold, but φ is not invariant: the class substitution mixes the two paths, and φ(SPS⁻¹) depends on g. A
dimension count suggests the family then sits in finitely many orbits. 133 tests this.

## 2. How many families: the vertex part

`analysis/v2_76_vertexpart.py` replaces λe_i by a general vertex vector c_ρ = Σ_v x_{ρ,v} e_v for each path. It
solves the fold conditions, then saturates by m·m′ and decomposes (`v2_76_vparts_primes.py`, Singular minAssGTZ).
Every component with m·m′ ≠ 0 has the same idempotent for both paths, c = λ·e_X:

| quiver | i, j, sinks | X (and its complement) | families |
|---|---|---|---|
| 1⇛2, 2→3, 2→4 | 1, 2; 3, 4 | {1}; {2, 3, 4} | 2 (v2.74's components 12, 13) |
| 1⇛2, 2→3→5, 2→4 | 1, 2; 4, 5 | {1}, {1, 3}; {2, 3, 4, 5}, {2, 4, 5} | 4 |
| 1⇛2, 2→3, 2→4, 3→4 | 1, 2; 4 | {1}, {1, 3}; complements | 4 |
| 1⇛2→3, 3→4, 3→5 | 1, 2; 4, 5 | {1}, {1, 3}; complements | 4 |

In each case X = {i} ∪ F, or the complement of such a set, where F is any subset of the intermediate vertices of the
paths (strictly between j and the sinks).

Why, in outline (i is a source in all four quivers):
- the kernel is spanned by κ(a₁·r);
- a sink s ∈ X would put (a₁ − a₂)·π into the kernel, so s ∉ X;
- j ∈ X would put ℓ₀ there, so j ∉ X;
- an intermediate vertex v only meets κ(a₁·r′) with r′ ending at v, and that product is again in the kernel, so v is
  free.

**When i has an incoming arrow** u : h→i, Δ forces the vertex part to be e_i. The term u ⊗ c from deleting u must
match λ·u ⊗ e_i. This is why v2.75 found only the src family on 5→1⇛2, 2→3, 2→4.

The families are pairwise disjoint over any field: the vertex part of δ_π, as a line, is H-invariant, because S fixes
the vertices. So the identity block has at least 2^{|F|+1}·(q − 1) classes, where F is now the set of all
intermediate vertices. For 1⇛2, 2→3, 2→4 that is 2(q − 1), matching v2.74's jump of 2. For 1⇛2, 2→3→5, 2→4 it is
4(q − 1), against an observed jump of 6. Two classes are still unaccounted for.

## 3. Conjecture 3ᶜ″

v2.75's Conjecture 3ᶜ′ allowed at most one arrow leaving the target of a size-≥3 class. 1⇛2→3, 3→4, 3→5 has one such
arrow and is infinite by Theorem 5′, so 3ᶜ′ is false.

**Conjecture 3ᶜ″.** An acyclic quiver has finitely many census orbits if and only if:
- it has no consecutive classes;
- it has no converging or diverging pair of classes with unjoined ends;
- for every class i⇛j of size ≥ 3, j and the vertices reachable from it form a chain j → v₁ → v₂ → ⋯: each of them has
  at most one successor vertex, though parallel arrows are allowed. Dually, i and the vertices from which i is
  reachable form a chain.

**What it covers.** In an unobstructed quiver, two distinct paths from j to sinks exist exactly when the descendants
of j branch, unless the only branching is a parallel class. Theorem 5′ proves the "only if" half in the branching
case. 3ᶜ″ guesses that a branching through a 2-class (1⇛2→3⇉4) does not count; 133 tests that.

**Evidence.** All 19 multiplicity-3 quivers examined fit (one still pending):
- **6 proved finite:** 1⇛2, 1⇛2→3, 1⇛2 + 1→3, their opposites, and 1⇛2→3→4;
- **flat at q = 2, 3:** 1→2⇛3→4, 1⇛2 + 1→3 + 4→2, 1⇛2→3 + 4→2, 1⇛2→4 + 1→3, and 1⇛2→3⇉4 (pending: q = 2 only so far);
- **infinite by Theorem 5 or 5′:**
  - 1⇛2, 2→3, 2→4 and four of its v2.75 neighbours (1⇛2 with three tails, 5→1⇛2, 2→3, 2→4, the extra arrow 5→2, and
    1→3);
  - 1⇛2, 2→3→5, 2→4; 1⇛2, 2→3, 2→4, 3→4; 1⇛2→3, 3→4, 3→5.

Every multiplicity-2 census-list quiver is consistent, since the condition only concerns classes of size ≥ 3.

## 4. Census

| quiver | q = 2 | q = 3 | |
|---|---|---|---|
| 1⇉2 with three tails | 3714 | 3714 | flat (127) |
| 1⇛2, 2→3, 2→4, 3→4 | 681 | 689 | +8 ≥ 4 families (130, 134) |
| 1⇛2→3, 3→4, 3→5 | 4899 | killed for lack of memory; rerun pending | infinite by Theorem 5′ |
| 1⇛2→3⇉4 | 2236 | not finished (76 min, fibre group 3 of 6) | the test of 3ᶜ″'s parallel-class clause |

## 5. Next

1. **Lean.** Formalise Theorem 5′ on an acyclic quiver, generalising v2.30's `Path.lean` set-up. On 1⇛2→3, 3→4, 3→5
   the generated check would be d = 22; the d = 15 file took about 10 minutes and 2.3 GB.
2. **The remaining classes.** Find the 2 unexplained classes of 1⇛2, 2→3→5, 2→4. Rerun its census with saved
   representatives at q = 2 and 3, and run the pipeline on the blocks that grow.
3. **The other half of 3ᶜ″.** Prove that "chain" quivers are finite at multiplicity 3: 1→2⇛3→4, 1⇛2 + 1→3 + 4→2 and
   the two v2.75 neighbours are only flat so far. Finish the v2.66 chain on them.
4. **Exact counts** for 1⇛2, 2→3, 2→4 (706 + 2(q − 1)?) at q = 5. That needs a bigger machine than this container.

## Files

- **analysis/:**
  - v2_76_pathfamily.py (the path version, with optional invariance check);
  - v2_76_vertexpart.py and v2_76_vparts_primes.py (the vertex-part families);
  - v2_75_phi.py now accepts paths.
- **pcensus.py:** K3_tgtfork_sq, K3_after_fork, K3_after_pair.
- **logs/v2_76/:**
  - PREDICTIONS.md;
  - pathfam_*.json, vertexpart_*.{json,sing};
  - census_*;
  - run_census*.sh.
