# Discordianism: Open Questions Explored (v2.75): the fork family (Theorem 5)

*v2.74 found that 1⇛2, 2→3, 2→4 has infinitely many orbits. It certified two invariants on two 7-dimensional
components, but only as ratios of stratum coordinates. This round finds what the invariant is and turns the example
into a theorem.*
- *The invariant is a cross-ratio of four lines in a plane.*
- *There is an explicit normal form P_t. It is a fold over every commutative ring, and distinct t give distinct
  census orbits. Both facts are kernel-checked in Lean.*
- *The family exists on every quiver with a class of at least three arrows i→j and two single arrows from j to
  sinks (Theorem 5), and dually at the source.*
- *The sink condition is not the boundary: 1⇛2, 2→3→5, 2→4 also grows (3417, 3423). Conjecture 3ᶜ′ replaces 3ᶜ:
  finite iff 3ᶜ's conditions hold and no class of size ≥ 3 has two arrows leaving its target, or two entering its
  source.*
- *Why multiplicity 2 is finite: the family needs a nonzero vector in ⟨a₂, …, a_n⟩ ∩ ker σ, which exists only for
  n ≥ 3.*

---

## 0. Results

The predictions are in logs/v2_75/PREDICTIONS.md, committed before the runs. Two facts were seen before that file was
written, and it records them as observations rather than predictions:
- the structure of the generic points (§1);
- that the normal form P_t is a fold on 1⇛2, 2→3, 2→4 (one run).

| # | Statement | Status |
|---|---|---|
| 121 | the variant P⁻ (vertex part e₂ + e₃ + e₄) is also a fold identically | **held** |
| 122 | φ (the cross-ratio) is H-invariant, symbolically; it equals v2.74's φ₁₂ and φ₁₃ on the rational points of components 12 and 13 | **held**: symbolic for generic S ∈ H on 7 quivers (random S on K₄'s variant); equal at 4 + 4 points |
| 123 | over every field the P_{1,1,1,t} are pairwise non-conjugate, as are the P⁻; the identity block has at least 2(q − 1) classes | **held**, proved twice: on paper (§2) and in Lean over every commutative ring (§3) |
| 124 | the obstruction at multiplicity 2 is ℓ₀ ∈ ker σ: with ℓ₀ = a2 the equations force m = m′ = 0 | **held**, on both 1⇛2, 2→3, 2→4 and 1⇉2, 2→3, 2→4 (the failing condition is Δ) |
| 125 | the projection family extends to six neighbouring quivers | **failed as stated** on 3 of 6: the "naive" extension moves only a0, a0·b, a0·b′. The corrected, multiplicative extension (§2) is a fold on 5 of 6. On 1⇛2, 2→3→5, 2→4 the tail that continues loses its parameters, which is why Theorem 5 asks for sinks |
| 126 | census flat at q = 2, 3 for the four non-family quivers; not flat for three tails | **mixed** (§5): flat for 1⇛2→3 + 4→2 (783) and 1⇛2→4 + 1→3 (640). For 1⇛2⇉3 it was a slip: that quiver has consecutive classes and is infinite by Theorem 4. For 1⇉2 with three tails the q = 3 run was killed by memory. For 1⇛2 with three tails only q = 2 was run (6692); it is infinite by Theorem 5 anyway |

## 1. What the generic points look like

`analysis/v2_75_show.py` prints rational points of a component as folds, labelled by paths. Arrows: a0, a1, a2 : 1→2
form the class; b = a3 : 2→3 and b′ = a4 : 2→4 are the tails. V = ⟨a0, a1, a2⟩, and σ : V → k is the sum of
coordinates. The class group is G = {g : σg = σ}, times weightings.

At every sampled point of component 12:
- **The image.** Im P = ⟨vertices, U, b, b′, U·b, U·b′⟩, where U ⊂ V is a plane. It has dimension 12.
- **The kernel of P on V.** P maps V to V, and ker(P|V) is a line w inside ker σ.
- **The tails.** P(v·y) = P(v)·y + f(v)·(λ_y e₁ + m_y ℓ₀) for y ∈ {b, b′}. Here f is the w-coordinate (it vanishes
  on U), and ℓ₀ spans the line U ∩ ker σ.

Component 13 is the same, with e₂ + e₃ + e₄ in place of e₁.

So each tail y contributes one vector δ_y = P(w·y) = λ_y e₁ + m_y ℓ₀, in the plane E = ⟨e₁, ℓ₀⟩. The torus scales
the two coordinates of δ_y by different characters, and each tail weight scales δ_y as a whole. That leaves one
modulus: the cross-ratio

  φ(P) = (v_b / v_b′) / (u_b / u_b′),  where v_y is the vertex part and u_y the V part of δ_y.

Equivalently, φ is the cross-ratio of the four lines E ∩ ⟨vertices⟩, E ∩ V, kδ_b and kδ_b′ in E.

On the 4 + 4 sampled points, the intrinsic φ (`analysis/v2_75_phi.py`) equals v2.74's φ₁₂ = y0_7·y6_8/(y0_8·y6_7)
on component 12, and φ₁₃ on component 13. Those coordinate formulas are this cross-ratio, read in one chart.

## 2. Theorem 5

**Setting.** Q is a finite acyclic quiver, and C = {a₁, …, a_n} is a class of n ≥ 3 parallel arrows i→j. Further:
- x : j→k and x′ : j→k′ are arrows that are each alone in their class;
- k and k′ are sinks;
- V = ⟨C⟩, σ(a_r) = 1, and ℓ₀ = a₂ − a₃, so ℓ₀ ∈ ker σ and ℓ₀ avoids a₁.

For λ, m, λ′, m′ in a commutative ring, define:
- **Π₀**, the algebra endomorphism a₁ ↦ a₂ that fixes the other arrows and the vertices;
- **N**, given by N(u·a₁·y) = u·(λ_y e_i + m_y ℓ₀) for y ∈ {x, x′} and u any path ending at i (u = e_i allowed),
  with (λ_x, m_x) = (λ, m) and (λ_x′, m_x′) = (λ′, m′); N = 0 on every other path;
- **P = Π₀ + N.**

So P(a₁) = a₂, P(a₁x) = a₂x + λe_i + mℓ₀, P(u·a₁·x) = u·a₂·x + λu + m·uℓ₀, and P fixes every path that avoids a₁.

**Theorem 5.**
1. P is a fold of kQ for every commutative ring k and all parameters.
2. Over a field F, let H be the census symmetries (class substitutions and weightings). Then φ is H(F)-invariant
   wherever it is defined, and on the family φ(P) = λm′/(λ′m).
3. Hence the folds P_t (λ = m = λ′ = 1, m′ = t ∈ F^×) lie in pairwise distinct H(F)-orbits. Q has infinitely many
   orbits over every infinite field, and its census has at least q − 1 classes in the block of P_t, for every q.

The dual statement, with two single arrows from sources into i, follows by the duality of v2.66.

*Proof of 1.* The four fold conditions:

(i) **Idempotent, and the image is closed.** Every path avoiding a₁ is fixed, and P(A) lies in their span A′. So
P² = P, and Im P = A′ is a subalgebra.

(ii) **The kernel is closed.** Ker P is spanned by κ(p) = p − P(p) for the paths p = u·a₁·r through a₁. In an acyclic
quiver a path meets C at most once, and no path returns from j to i. So in a product κ(p)κ(p′):
- the parts u(a₁−a₂)r · (anything from p′) vanish;
- m·uℓ₀ · (anything) vanishes;
- what survives is λu · u′(…), which needs u′ = e_i.

In that case κ(p)κ(p′) = −λ·κ(u·a₁·r′), and in every other case it is 0.

(iii) **Δ-compatibility.** Δ deletes one arrow, giving prefix ⊗ suffix. Substitution commutes with deletion, so Π₀ is
Δ-compatible, and it remains to check ΔN = (P⊗P)Δ − (Π₀⊗Π₀)Δ on paths.

For p = u·a₁·x, the only deletions that see N are those of an arrow inside u = u₁αu₂. They give
Σ u₁ ⊗ u₂(λe_i + mℓ₀). Deleting a₁ or x gives nothing new: x ends at a sink, so no path continues past it.

On the other side, Δ(λu + m·uℓ₀) is the same sum, plus the term m·u ⊗ Δ-part of ℓ₀, which is m·σ(ℓ₀)·u ⊗ e_j = 0.
For u = e_i, both sides are 0 for the same reason.

(iv) **Where the hypotheses enter.** σ(ℓ₀) = 0 is used exactly once, in (iii). The sink condition is used once too: a
path a₁·x·q would put (λe_i + mℓ₀) ⊗ q into (P⊗P)Δ, and no N can match it. ∎

*Proof of 2.* Let S ∈ H(F). Then S fixes the vertices, maps V to V, and S(y) = t_y·y for y = x, x′ (each is alone in
its class). Let P′ = SPS⁻¹. Then:
- ker(P′|V) = S·ker(P|V), so take w′ = Sw;
- δ′_y = P′(Sw·y) = S·P(w·S⁻¹y) = t_y⁻¹·S·δ_y;
- the vertex parts therefore scale by t_y⁻¹, and the V parts are mapped by t_y⁻¹·S, which is linear.

Both ratios in φ change by the same factor t_x′/t_x, so φ(P′) = φ(P).

On the family, ker(P|V) = ⟨a₁ − a₂⟩, and δ_y = P(a₁y) − P(a₂y) = λ_y e_i + m_y ℓ₀. ∎

**Why multiplicity 2 is different.** The correction ℓ₀ must lie in the image (so it avoids a₁) and in ker σ (for Δ).
⟨a₂, …, a_n⟩ ∩ ker σ is nonzero exactly when n ≥ 3. With n = 2 the V-parts u_y vanish and φ degenerates.

The computation agrees (`analysis/v2_75_family.py`). With ℓ₀ = a₂ (not in ker σ), the Δ equations force m = m′ = 0
on both 1⇛2, 2→3, 2→4 and 1⇉2, 2→3, 2→4. In Lean, `forkBad_check_fails` shows the same failure at Δ.

## 3. Checks

**Fold identities** (`analysis/v2_75_family.py`, FAMILY_MODE=mult). These are identities in ℚ[λ, m, …] with integer
coefficients, so they hold over ℤ. Vertex part c: e_i ("src"), all other vertices ("others"), or e_j + the tail
targets ("local").

| quiver | src | others | local |
|---|---|---|---|
| 1⇛2, 2→3, 2→4 | fold | fold | fold |
| K₄ variant 1⇶2, 2→3, 2→4 | fold | fold | fold |
| 1⇛2 with three tails (all three get parameters) | fold | fold | fold |
| 5→1⇛2, 2→3, 2→4 | fold | λ forced 0 | λ forced 0 |
| 1⇛2, 2→3, 2→4, 5→2 | fold | fold | fold |
| 1⇛2, 2→3, 2→4, 1→3 | fold | fold | fold |
| 1⇛2, 2→3→5, 2→4 (the continuing tail is a3) | λ₃ = m₃ = 0 | λ₃ = m₃ = 0 | λ₃ = m₃ = 0 |

The vertex part e_{t(y)} ("tail") fails on every quiver (it forces λ_y·m_y = 0, or more).

**Invariance** (`analysis/v2_75_invariance.py`). Take S generic in H, with all g, s symbolic. Then φ(SPS⁻¹) = φ(P) =
λm′/(λ′m) after cancellation:
- on 1⇛2, 2→3, 2→4 (src, others, local);
- on 1⇛2 with three tails, for all three tail pairs;
- on 5→1⇛2, 2→3, 2→4;
- on 1⇛2, 2→3, 2→4, 5→2 (src, others);
- on 1⇛2, 2→3, 2→4, 1→3 (src, others).

The symbolic 4×4 inverse for the K₄ variant did not finish; three random rational S gave equality.

**Lean** (`lean/EatTheFibre/Gen/K3_tgtfork.lean`, `K3_tgtfork_rigid.lean`, generated by `lean/gen/gen_fork.py`).
Every result below is kernel-checked and uses only the axioms propext, Classical.choice and Quot.sound
(`logs/v2_75/lean_audit_fork.txt`).
- **`fork_isFold`, `forkM_isFold`:** P_t and P⁻_t are folds over every commutative ring, for every t (per-row
  `decide +kernel`).
- **`forkBad_check_fails`:** with ℓ₀ = a2 the Δ condition fails, as a polynomial identity.
- **`fork_rigid_GT`, `forkM_rigid_GT`:** substK·P_t = P_t′·substK implies t = t′. Here substK is the census symmetry:
  the class substituted by a 3×3 matrix with column sums 1, the class weight wa, the tail weights wb, wc. The proof
  needs no hypotheses: four entries E₀…E₃ of the difference satisfy t·E₀ − t′·E₁ − t·E₂ + E₃ = t − t′ in
  ℤ[g, w, t, t′] (a Singular lift over ℤ, `analysis/v2_75_lift.py`).
- **`fork_forkM_disjoint`:** for the two families, one entry of substK·P_t − P⁻_t′·substK is the constant 1. So no
  census symmetry joins them over a nontrivial ring.

Hence, over 𝔽_q, the identity block of 1⇛2, 2→3, 2→4 has at least 2(q − 1) census classes. This is consistent with
v2.74's 111 at q = 2 and 113 at q = 3.

## 4. Conjectures 3ᶜ and 3ᵖ

**Conjecture 3ᵖ is false too.** It is v2.25's Conjecture 3\*ˡ for one parallel class (finitely many orbits under Aut_Δ
and local class-torus moves). 1⇛2, 2→3, 2→4 has one parallel class, and Theorem 5 applies to it.

**Unobstructed quivers.** In a quiver without consecutive classes, every arrow leaving the target of C is
automatically alone in its class. So the only extra hypothesis of Theorem 5 there is two arrows from j to sinks (or,
dually, two arrows from sources into i).

**The sink condition is not the boundary.** 1⇛2, 2→3→5, 2→4 has only one arrow from j to a sink, so Theorem 5 does
not apply to it, and the multiplicative family loses the continuing tail's parameters (§3). Yet its census is not
flat: 3417 at q = 2 and 3423 at q = 3 (§5). Whatever grows there is a second mechanism, not yet found.

So the natural revision drops the sinks.

**Conjecture 3ᶜ′.** An acyclic quiver has finitely many census orbits if and only if:
- it has no consecutive classes;
- it has no converging or diverging pair of classes with unjoined ends;
- for each class of size ≥ 3, at most one arrow leaves its target and at most one arrow enters its source.

**Evidence.** Every multiplicity-3 quiver computed so far fits:
- **proved finite (6):** 1⇛2, 1⇛2→3, 1⇛2 + 1→3, their two opposites, and 1⇛2→3→4 (v2.66, v2.74);
- **flat at q = 2, 3 (4):** 1→2⇛3→4 and 1⇛2 + 1→3 + 4→2 (v2.74); 1⇛2→3 + 4→2 and 1⇛2→4 + 1→3 (§5);
- **infinite by Theorem 5 (5):** 1⇛2, 2→3, 2→4 and its four neighbours in §3 other than the long tail;
- **not flat (1):** 1⇛2, 2→3→5, 2→4.

Arrows leaving the source or entering the target of a size-3 class do not matter in this data (1⇛2 + 1→3,
1⇛2→3 + 4→2). Classes of size 2 never trigger the condition, which is why 3ᶜ held on every multiplicity-2 quiver.

## 5. Census of the neighbours (torch engine)

| quiver | q = 2 | q = 3 | |
|---|---|---|---|
| 1⇛2→3 + 4→2 (K3_tgtinout) | 783 | 783 | flat, as predicted |
| 1⇛2→4 + 1→3 (K3_srcout_after) | 640 | 640 | flat, as predicted |
| 1⇛2⇉3 (K3_tgtpair) | 7020 | — | consecutive classes: infinite by Theorem 4; q = 3 not feasible (1.3·10⁸ kernel subspaces) |
| 1⇉2 with three tails (pair_tgtfork3) | 3714 | killed at 20480/23294 fibres by the container's memory limit (a Lean build ran beside it) | rerun not finished this round |
| 1⇛2 with three tails (K3_tgtfork3) | 6692 | not run (q = 2 took 21 minutes) | infinite by Theorem 5 |
| 1⇛2, 2→3→5, 2→4 (K3_tgtfork_long) | 3417 | 3423 | **not flat** (66 minutes at q = 3); not covered by Theorem 5 |

## 6. Next

1. **The second mechanism.** Find the family on 1⇛2, 2→3→5, 2→4:
   - rerun the census with saved representatives at q = 2 and 3, to find the block that grows (6 classes);
   - run the v2.55 pipeline and the exact generic rank on that block;
   - generalise Theorem 5 to any two arrows out of j, as Conjecture 3ᶜ′ predicts.
2. **Theorem 5 in Lean on every acyclic quiver,** as v2.30 did for Theorems 3 and 4. The paper proof (§2) is short
   and uses only path combinatorics.
3. **Exact counts.** Find the census of 1⇛2, 2→3, 2→4 (predicted 706 + 2(q − 1)), and the q = 5 run that tests it.
   That run did not fit this container (5160 H-orbit fibres; killed at the turn end in v2.74).
4. **Finish the chain** on 1⇛2→3→4, 1→2⇛3→4 and 1⇛2 + 1→3 + 4→2, and do K₄ at p = 2.

## Files

- **analysis/:**
  - v2_75_show.py (points of a component as folds; intrinsic φ against a coordinate formula);
  - v2_75_phi.py (the intrinsic cross-ratio);
  - v2_75_family.py (the families, naive and multiplicative, checked against the fold equations);
  - v2_75_invariance.py (φ(SPS⁻¹) = φ(P), symbolic or at random S);
  - v2_75_rigid.py, v2_75_lift.py (rigidity ideals and their integer certificates).
- **lean/:**
  - gen/gen_fork.py;
  - EatTheFibre/Gen/K3_tgtfork.lean, K3_tgtfork_rigid.lean;
  - Audit.lean and build.sh updated.
- **pcensus.py:** ten new quivers (v2.75 block in QUIVERS).
- **logs/v2_75/:**
  - PREDICTIONS.md;
  - family_*.json, familym_*.json, inv_*.json;
  - rigid/ (ideals, certificates.json);
  - show12/13 output;
  - census_*;
  - lean_build.log, lean_audit_fork.txt.
