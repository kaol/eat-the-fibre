# Discordianism: Open Questions Explored (v2.27)

*Takes up v2.26's mechanism claim ("the moduli come from the product of parallel differences") and tests it. It is half
right. Consecutive parallel classes always carry a family (Theorem 2, the general form of Theorem 1). But a census
designed to confirm the converse failed: 1⇉2⇇3 has no such product, and it still has infinitely many orbits. The cause
is a second, simpler mechanism (Theorem 3). There the modulus is a vertex coefficient of P, which no symmetry can move.*

*Scripts (in `analysis/`): `theorem2.py`, `theorem3.py`, `fixedlocus.py`, `family.py`, `family_invariant.py` (patched:
rational substitutions keep numerators). Logs: `logs/theorem2.log`, `logs/theorem3.log`,
`logs/fixedlocus_pairs_in_3.log`, `logs/family_pairs_in_3_2.log`, `logs/results_container.jsonl`. Labels: **proved**,
**verified (exact)**, **numerical**, **conjectured**. Predictions were committed to git before each census ran
(commits 2d85cef, 2e088c2, e398159). Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Summary

| # | Question | Answer |
|---|---|---|
| v2.26 Q3 | Do consecutive parallel classes always give moduli? | **Theorem 2:** yes. For classes A: u→v, B: v→w (|A|, |B| ≥ 2, u ≠ w, no arrow between u and w), P_t is a fold for every t, and t is invariant **[proved; symbolic on 4 quivers, fails on the control]** |
| — | Converse: no consecutive classes ⇒ finitely many orbits? | **False.** The prediction for 1⇉2⇇3 (committed first) failed. Local-torus orbits: 453, 461, 469, 485 over 𝔽₃, 𝔽₅, 𝔽₇, 𝔽₁₁ **[verified (exact)]** |
| — | The second mechanism | **Theorem 3:** for converging classes A: u→v, B: w→v (or diverging ones), there are three explicit one-parameter families of folds. Their invariant is the **vertex block** of P (P's coefficients of vertices on vertices). G, the class weightings and (on the quivers checked) Aut_Δ° fix it **[proved; symbolic on 7 quivers, fails on the control]** |
| — | Is that all of 1⇉2⇇3? | Yes over 𝔽₃ and 𝔽₅. Orbits = 449 + 4(q − 2): every vertex block carrying exactly one orbit is a generic member of a family, and the rest total 449 at both primes. Vertex blocks: 4q + 19 (q = 3, 5, 7) **[verified (exact)]** |
| — | Contrast | 1⇉2⇉3 has 23 vertex blocks over 𝔽₃ and over 𝔽₅. Its moduli sit in the radical (Theorem 1). The two mechanisms differ |
| — | Conjectures | **3ᵛ** (finite iff no two classes of size ≥ 2 share a vertex) was refuted by its own recorded test: 1⇉2⇇3 + 1→3 has 366 local-torus orbits over 𝔽₃ and 𝔽₅. **3ʰ:** infinitely many iff Theorem 2 or 3 applies. Tests: §5 |

---

## 1. Theorem 2: consecutive classes

Let A: u → v and B: v → w be parallel classes with |A|, |B| ≥ 2, u ≠ w, no arrow u → w or w → u, and no loop at u or w.
Take a₁, a₂ ∈ A, b₁, b₂ ∈ B, δ_a = a₁ − a₂, δ_b = b₁ − b₂. For t ∈ k put

  P_t(e_u) = e_u − δ_a − (1 + t)·δ_aδ_b,  P_t(e_w) = e_w − δ_b + t·δ_aδ_b,  P_t = 0 on every other basis path.

**Theorem 2.** P_t is a Δ-compatible fold for every t. Write P(e_u) = e_u − x − y and P(e_w) = e_w − z + w′ (x, z of
length 1; y, w′ of length 2). Then t = w′/(x·z) is preserved by G and by every class weighting. Hence the P_t lie in
pairwise distinct orbits.

*Proof.*
- **Image.** The image is span{P(e_u), P(e_w)}. It is closed under multiplication because P(e_u)P(e_w) = 0 = P(e_w)P(e_u),
  and both are idempotents.
- **Kernel.** The kernel is spanned by the other basis paths together with e_u − P(e_u) and e_w − P(e_w). Products of
  kernel elements never produce e_u or e_w with a nonzero coefficient.
- **Δ-compatibility.** Δ kills δ_a, δ_b and δ_aδ_b, so ΔP = 0. For a basis path b, (P⊗P)Δ(b) is nonzero only if b has a
  prefix and a suffix in {e_u, e_w}, with one arrow deleted between them. That requires an arrow between u and w, which
  the hypothesis excludes.
- **Invariant.** G and the class weightings fix vertices. They scale δ_a, δ_b and δ_aδ_b compatibly (as in Theorem 1),
  so w′/(xz) does not change. ∎

`theorem2.py` checks the fold equations symbolically over ℚ(t):
- **true** on 1⇉2⇉3, 1⇉2⇛3, 1→2⇉3⇉4 and 1⇉2⇉3→4;
- **false** on 1⇉2⇉3 + 1→3, where the hypothesis fails.

Theorem 1 (v2.26) is the case 1⇉2⇉3.

---

## 2. The failed prediction

v2.26 attributed the moduli to δ_aδ_b ∈ ker Δ. 1⇉2⇇3 has two parallel classes but no composable pair, so no such product.
Prediction (commit 2d85cef): the local-torus count would be the same over 𝔽₃ and 𝔽₅. Census (container, torch engine):

| q | N | G-orbits | G × grading = G × T = local torus |
|---|---|---|---|
| 3 | 6 376 | 485 | 453 |
| 5 | 59 972 | 557 | 461 |
| 7 | 281 568 | 629 | 469 |
| 11 | 2 385 368 | 773 | 485 |

**The prediction failed:** 445 + 4(q − 1). The growth sits in one orbit type, size q²(q − 1)²: 63, 71, 79, 95 orbits.
Only the grading preserves the whole fold set (torus_preserving = q − 1 of (q − 1)²). Some folds mix the two classes, so
a class weighting is not a symmetry there.

**Locating the family.** `fixedlocus.py` (v2.26 method) over 𝔽₃: 453 orbits, 4 with no member lifting to ℚ, 39
candidates (tangent of Folds^S > dim N·P). 14 candidates carry a family numerically, all of size 36 = q²(q − 1)².

Caveat: `fixedlocus.py` uses the full class torus in Lie H. Where a weighting is not a symmetry, dim N·P overcounts the
orbit inside Folds; rows with tangent < dim N·P show this. The family verdicts stay valid, because the overcount only
makes the test stricter. The "no family" verdicts may not be.

Singular (`family.py`) on two candidates found, exactly, one component of Folds^S through P of dimension 3 > n = 2.
`family_invariant.py` parametrized it, checked the fold equations for all parameter values, and found the N-invariant.
It is **y₂₆, the e₂-coefficient of P(e₁)**: a vertex coefficient.

---

## 3. Theorem 3: converging classes

**Lemma (vertex block).** Let Q be acyclic, and let π₀: kQ → kQ₀ be the projection along the radical (the paths of
positive length). Let h be an algebra automorphism with h(e_i) = e_i for all i. Then π₀(hPh⁻¹)|_{kQ₀} = π₀P|_{kQ₀}
for every linear P. This applies to every element of G, to every class weighting, and to the identity component of
Aut_Δ: every Δ-derivation of kQ kills the vertices (`theorem3.py`, check 2: dim Lie(Aut_Δ) = dim Lie G = 4 on 1⇉2⇇3,
1⇇2⇉3, 1⇉2⇉3 and 1⇉2⇇3→4).

*Proof.* For an arrow c: s → t, h(c) = h(e_s c e_t) = e_s h(c) e_t ∈ e_s·kQ·e_t, which lies in rad because s ≠ t
(Q is acyclic, so it has no loops). So h(rad) ⊆ rad and h fixes kQ₀; hence π₀h = π₀ and h⁻¹|_{kQ₀} = id. Then
π₀ hPh⁻¹ |_{kQ₀} = π₀ P |_{kQ₀}. For Aut_Δ° it suffices that its Lie algebra kills the vertices; that is checked, not
proved in general. ∎

So the vertex block B(P) = π₀P|_{kQ₀} is a matrix-valued invariant. If it takes infinitely many values on folds, there
are infinitely many orbits.

**Setting.** Let A: u → v and B: w → v be parallel classes with |A|, |B| ≥ 2, u ≠ w, and no arrow between u and w.
Diverging classes (A: v → u, B: v → w) work by op-duality: the fold conditions are invariant under kQ ↦ (kQ)^op with Δ
flipped. P vanishes on every basis path except e_u, e_w, a₂ and b₂. The families are:

| family | P(e_u) | P(e_w) | P(a₂) | P(b₂) | vertex block |
|---|---|---|---|---|---|
| F1(s) | e_u − s e_v + (s + 1)δ_a + δ_b | −P(e_u) | z = −s e_v + s δ_a + δ_b | −(1 + s) z | P(e_u) ≡ e_u − s e_v |
| F2(c) | c e_u + (c − 1)e_w + c(c − 1)(δ_a − δ_b) | −P(e_u) | z = e_u + e_w + (c − 1)δ_a − c δ_b | z | P(e_u) ≡ c e_u + (c − 1) e_w |
| F4(x) | e_u + x e_v + x(1 − x)(δ_a − δ_b) | e_w + (1 − x)e_v − x(1 − x)(δ_a − δ_b) | z = e_v − x δ_a + (x − 1)δ_b | z | P(e_u) ≡ e_u + x e_v, P(e_w) ≡ e_w + (1 − x)e_v |

F3 is F1 with u and w exchanged. F4 is unital on 1⇉2⇇3 (P(1) = 1). F1 and F2 have P(1) = 0.

**Theorem 3.** Under the hypothesis, F1(s), F2(c) and F4(x) are Δ-compatible folds for every value of the parameter.
Distinct parameter values give distinct vertex blocks. So, over an infinite field, the folds of kQ fall into
infinitely many orbits under G and the class weightings, and also under Aut_Δ° wherever its Lie algebra kills the
vertices (all quivers checked).

*Proof.* Put W = span{e_u, e_v, e_w} ⊕ kA ⊕ kB. This is the path algebra of the full subquiver on the two classes. P
maps into V = span{e_u, e_v, e_w, δ_a, δ_b} ⊆ W.
1. Idempotency and closure of the image are computations in V. They do not depend on Q.
2. Kernel closure. Every basis path outside W is killed, so ker P = (span of the paths outside W) ⊕ (ker P ∩ W).
   - Two paths outside W multiply to 0, to one of them (when the other is a vertex e_x with x ∉ {u, v, w}), or to a
     path of length ≥ 2. All of these lie outside W.
   - For a path π outside W and k ∈ ker P ∩ W: kπ and πk are a vertex coefficient of k times π, plus arrow·π or π·arrow
     (length ≥ 2). All of these lie outside W, so in the kernel.
   - (ker P ∩ W)² ⊆ ker P is a computation in W.
3. Δ-compatibility. P(·) ∈ V ⊆ ker Δ, so ΔP = 0. (P⊗P)Δ(b) ≠ 0 needs b = prefix·c·suffix with prefix and suffix in
   {e_u, e_w, a₂, b₂}. A prefix in {e_u, e_w} with a suffix in {e_u, e_w, a₂, b₂} needs an arrow c between u and w
   (excluded). A prefix in {a₂, b₂} needs c to leave v and the suffix to start at u or w: a cycle (excluded, Q acyclic).
   This is the same mechanism as Theorem 2.

So the claims reduce to the fold equations inside W = k(u ⇉ v ⇇ w), which are verified symbolically over ℚ(s), ℚ(c)
and ℚ(x) (`theorem3.py`, check 1: classes of sizes (2, 2) and (3, 2)). Larger classes add arrows that P kills, and step
2 applies to them verbatim. The invariance is the Lemma. ∎

Check 1 also runs the full equations on more quivers:
- **fold for all parameter values:** 1⇉2⇇3, 1⇇2⇉3, 1⇉2⇇3 with 2→4, with 4→2, with 4→1, with 1→4, and 1⇛2⇇3, for
  both choices of (u, w);
- **not a fold:** 1⇉2⇇3 + 1→3, where the hypothesis fails.

**Census check over 𝔽₃ and 𝔽₅** (`theorem3.py`, check 3; H = G × the weightings that preserve the fold set):

| | 𝔽₃ | 𝔽₅ | 𝔽₇ |
|---|---|---|---|
| members of F1–F4 in the census | all | all | all |
| distinct vertex blocks | 31 | 39 | 47 |
| H-orbits | 453 | 461 | (469, census) |
| blocks with exactly one orbit | 4 | 12 | |
| … equal to the generic family members (q − 2 per family) | yes | yes | |
| orbits on the other blocks | 449 | 449 | |

The histogram of orbits per block agrees at 𝔽₃ and 𝔽₅ except for the single-orbit blocks. So the vertex block accounts
for all the growth: **orbits = 449 + 4(q − 2)** (= 469 at q = 7 and 485 at q = 11, matching the census). The degenerate
parameter values are s ∈ {0, −1}, c ∈ {0, 1} and x ∈ {0, 1}. There the block carries more orbits. One member each of F1 and F3 has the
smaller orbit size q(q − 1)² (12 at q = 3, 80 at q = 5); every other member has size q²(q − 1)².

For contrast, 1⇉2⇉3 has 23 vertex blocks over 𝔽₃ and over 𝔽₅. Its moduli live in the radical.

**Gloss.** *When twins arrive at the same place from different homes, a split can decide how much of each home it hands
to the meeting place. No relabelling of twins, and no clock, touches a home: whatever share it chose, it keeps.*

---

## 4. Proposed amendments (v2.26 → v2.27)

| amendment | content |
|---|---|
| **Q1²⁷ Theorem 2** | Consecutive parallel classes (outer ends distinct and non-adjacent) ⇒ a one-parameter family with invariant w′/(xz) |
| **Q2²⁷ Vertex-block lemma** | Every vertex-fixing automorphism preserves the vertex block of P (π₀P restricted to kQ₀); on the quivers checked this includes Aut_Δ° |
| **Q3²⁷ Theorem 3** | Converging or diverging parallel classes (outer ends distinct and non-adjacent) ⇒ three explicit families with distinct vertex blocks |
| **Q4²⁷ Retraction** | v2.26 Q3²⁶ ("the moduli come from the product of parallel differences") is withdrawn as a characterization. It describes Theorem 2's mechanism only |
| **Q5²⁷ Conjecture 3ʰ** | Infinitely many local-torus orbits ⟺ two parallel classes of size ≥ 2 share a vertex with outer ends distinct and non-adjacent (the hypotheses of Theorems 2 and 3). Conjecture 3ᵖ (one class) is a special case. (3ᵛ, without the adjacency clause, is refuted: §5) |

---

## 5. Predictions and outcomes

The predictions as committed before the censuses (verbatim from the draft):

*Predictions recorded before the censuses run (this file is completed when they finish).*

- **1⇉2⇇3** (two parallel classes, not consecutive): Conjecture 3ᶜ predicts finitely many classes. Prediction: the
  orbit count under the local class torus (`orbits_GxTloc`) is the same over 𝔽₃ and 𝔽₅. (G-orbit counts may grow: the
  quiver has a same-orientation fork at vertex 2, which v2.22 showed gives grading/torus-absorbed families.)
- **1⇉2⇉3 + 1→3** (consecutive classes, but Theorem 2's hypothesis fails): no prediction.

**Outcome (recorded before any analysis of it): the 1⇉2⇇3 prediction failed.** Local-torus orbit counts over
𝔽₃, 𝔽₅, 𝔽₇, 𝔽₁₁: 453, 461, 469, 485 = 445 + 4(q − 1). G-orbits 485, 557, 629, 773. G × grading = G × T = local torus
at every prime. The only orbit type whose count grows is the one of size q²(q − 1)² (63, 71, 79, 95). Raw lines:
`logs/results_container.jsonl`.

**Second set of predictions (recorded after Theorem 3 below was found and checked, before these censuses run):**

- **1⇇2⇉3** (opposite quiver): identical numbers to 1⇉2⇇3 over 𝔽₃ and 𝔽₅ (N, and the orbit count under every group).
  This follows from op-duality, so it tests the engine, not the mathematics.
- **1⇉2⇇3, 2→4**: the four families of Theorem 3 persist (symbolic check), and no other mechanism is present. Prediction:
  the local-torus orbit count grows by exactly 8 from 𝔽₃ to 𝔽₅, and every vertex block that carries exactly one orbit is a
  generic member of one of the four families.

**1⇇2⇉3 outcome:** N = 6376 / 59 972 and 485 / 557 G-orbits, 453 / 461 under every larger group, over 𝔽₃ / 𝔽₅: identical
to 1⇉2⇇3, as op-duality requires.

**Third prediction (recorded before the census):** Conjecture 3ᵛ (below) says the local-torus orbit count is finite iff no
two parallel classes of size ≥ 2 share a vertex. **1⇉2⇇3 + 1→3** violates the hypothesis of Theorem 3 (the families
P_s are not folds there), but its classes share vertex 2. Prediction under 3ᵛ: its local-torus count grows from 𝔽₃ to
𝔽₅. A constant count refutes 3ᵛ as stated.

**Outcomes so far (appended in order; nothing above is edited):**

1. 1⇇2⇉3: as predicted (identical to 1⇉2⇇3).
2. 1⇉2⇇3, 2→4: **both parts of the prediction failed.** The local-torus count is 1919 over 𝔽₃ and 1955 over 𝔽₅, which
   grows by 36, not 8. Single-orbit vertex blocks: 21 and 57, not only generic members of F1–F4. What survived: the
   vertex block still accounts for all of the growth. The orbits-per-block histogram is identical at both primes except
   for the single-orbit blocks. The other blocks carry 1898 orbits at both primes. The single-orbit blocks are:
   - 3 constant ones;
   - the generic members of 18 lines of blocks. On {1, 2, 3}, each line restricts to one of F1–F4: 5 lines for F1, 5
     for F3, 3 for F2 and 5 for F4. The lines differ in how e₄ enters: P(e₄) ∈ {0, e₄}, and the e₄-coefficients of
     P(e₁) and P(e₃).

   So orbits = 1901 + 18(q − 2) at q = 3, 5.
3. 1⇉2⇇3 + 1→3: **the prediction failed.** The local-torus count is 366 over 𝔽₃ and 366 over 𝔽₅. (G × T grows,
   376 → 384, but local class-torus moves absorb the growth.) Conjecture 3ᵛ as stated is refuted by the criterion
   recorded in advance.

**Revised Conjecture 3ʰ** (recorded before the next census): Folds(kQ) has infinitely many local-torus orbits iff two
parallel classes of size ≥ 2 share a vertex and their outer ends are distinct and non-adjacent, i.e. iff Theorem 2 or
Theorem 3 applies. **Prediction for 1⇉2⇉3 + 1→3** (consecutive classes with adjacent outer ends): the local-torus count
is the same over 𝔽₃ and 𝔽₅.

---

## 6. Open questions

1. **Conjecture 3ʰ**: why does an arrow between the outer ends remove the moduli (1⇉2⇇3 + 1→3: 366 orbits at 𝔽₃ and
   𝔽₅)? A proof needs the "only if" direction, i.e. a finiteness argument.
2. Finiteness when classes are pairwise vertex-disjoint (1⇉2→3⇉4 is finite; v2.26). Is there a proof through the
   kernel fibration, by showing that the vertex block and the radical invariants are both discrete?
3. Are F1–F4 the complete list of vertex-block families for converging classes, on every quiver satisfying the
   hypothesis?
4. Carried: classify the 19 families of 1⇉2⇉3; Conjecture 3ᵖ; unobstructed violations; Kronecker plus a loop; zigzag;
   dihedral; ρ.

---

## 7. Sources

- M. Aguiar, *On the associative analog of Lie bialgebras*, J. Algebra 244 (2001).
- W. Decker, G.-M. Greuel, G. Pfister and H. Schönemann, *Singular* (via passagemath-singular): minimal associated primes.

*"When twins arrive at the same place from different homes, a split can decide how much of each home it hands to the
meeting place, and no relabelling and no clock can take that share back."*
