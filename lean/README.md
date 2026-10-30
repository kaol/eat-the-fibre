# Lean 4 formalization (v2.29–v2.44)

Lean 4 (v4.35.0) with Mathlib (v4.35.0, pinned in `lake-manifest.json`). Mathlib is used from its prebuilt cache; `.lake/`
(packages, oleans, ~7 GB) is git-ignored and never committed.

```bash
curl -sSfL https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y --default-toolchain none
cd lean && ./build.sh        # lake exe cache get, then each module in turn (the kernel checks need ~5 GB each)
python3 gen/gen.py           # regenerate EatTheFibre/Gen/*.lean from the census engine (pcensus.build)
python3 gen/gen_orbit.py     # regenerate EatTheFibre/Gen/orbits_*.lean (v2.33)
python3 gen/gen_families.py NF.pkl REPS.pkl   # regenerate Gen/two_pairs_fam*.lean (v2.36; inputs in logs/v2_36, logs/v2_32)
python3 gen/gen_rigid.py CERT.pkl              # regenerate Gen/two_pairs_rigid.lean (v2.37; input logs/v2_37/certificates_two_pairs_7.pkl)
python3 gen/gen_disjoint.py CERT.pkl           # regenerate Gen/two_pairs_disjoint.lean (v2.38; same input)
python3 gen/gen_behaviour.py CERT.pkl          # regenerate Gen/{two_pairs,pairs_in}_behaviour.lean (v2.39; same input)
python3 gen/gen_swap.py                       # regenerate Gen/pairs_in_aut.lean (v2.41)
python3 gen/gen_line.py REPS.pkl 5107 LOG 0    # regenerate Gen/two_pairs_short_line.lean (v2.44; logs/v2_44/lines_5107.txt, merged F_5 reps)
python3 gen/gen_short_families.py ../logs/v2_42/curves_two_pairs_short_5_merged.txt REPS_5_merged.pkl   # Gen/two_pairs_short_fam{1..18}.lean (v2.49)
python3 gen/gen_curve.py ../logs/v2_50/curve_A.json A   # Gen/two_pairs_short_curve{A,B}.lean (v2.50; B with curve_B.json)
python3 gen/gen_vblocks.py ../logs/v2_52/realise_7.json   # Gen/pairs_in_short_vblocks.lean (v2.52)
python3 gen/gen_fork.py                                    # Gen/K3_tgtfork{,_rigid}.lean (v2.75; needs Singular for the lift)
```

## What is proved

| file | content |
|---|---|
| `Basic.lean` | Coordinate algebras (structure constants `c`, arrow deletion `D`), `IsFold` (idempotent, image and kernel closed under products, `Δ(Pv) = P Δ(v) Pᵀ`). **Basis criterion** `isFold_iff_basis`: the conditions hold for all vectors iff they hold on basis vectors. That is the method of every Python check. **`IsFold.conj`**: conjugating by a Δ-compatible automorphism gives a fold. **`vertex_block`**: an automorphism fixing the vertices and preserving the radical does not change the vertex-to-vertex entries |
| `Poly.lean` | Computable integer polynomials; `eval` is a ring map into any commutative ring; `isZero p → eval t p = 0` |
| `Reflect.lean` | `check c D P : Bool` evaluates the four basis conditions on polynomial matrices. **`check_sound`**: `check = true` gives `IsFold` over every commutative ring at every parameter value. The check is split into rows, one `decide +kernel` each |
| `Gen/*.lean` | Generated from `pcensus.build`, i.e. the census engine's own basis, products and Δ. Families proved to be folds (`*_isFold`, with readable images `*_a2`, `*_kills` …), plus controls proved to fail (`*_check_fails`) |
| `Invariants.lean` | Distinct parameters give distinct orbits: `thm3_rigid_GT`, `thm4_rigid_GT` (and `thm4_core`) |
| `Reduction.lean` | **Reduction to a sub-basis** (`isFold_of_support`); `Embeds` (a small algebra inside a big one), `isFold_push` (folds push along embeddings, given the Δ condition off the image), `Embeds.trans`, `conj_delta_eq_zero` |
| `Path.lean` | Path algebras of finite quivers (`FinQuiver`, `Path`, `pathAlg`, census conventions); finiteness from an acyclicity potential; injective quiver morphisms give embeddings (`Hom.embeds`) |
| `General.lean` | **Theorems 3 and 4 on every acyclic quiver** (`thm3_general`, `thm4_general`) and their rigidity (`thm3_general_rigid`, `thm4_general_rigid`); the census tables of 1⇉2⇉3, 1⇉2⇇3 are their path algebras (`decide +kernel`) |
| `Examples.lean` | The general theorems applied to 1⇉2⇉3 + 1→3, 1⇉2⇇3 + 2→4, 4→1⇉2⇇3 |
| `Discordia/Fold.lean` | The v2 algebra layer derived: M4, RB(−1), δ-L A1–A5, δ-L1, E1 from M3 (idempotent ring endomorphism); δ-D, δ-D1 from RB(−1); E1⋆ from idempotence |
| `Discordia/Completion.lean` | δ-D2 derived and refined: dendriform completions of a fold by class (any fold / `ker·im ⊆ ker` / bimodule / `ker·im·ker = 0`); every idempotent endomorphism is in the bimodule class (`endo_bimod`) |
| `Discordia/Exact.lean` | δ-D2 exactly (v2.32): exact defects of every completion family for an arbitrary fold; in a unital algebra each stratum is present iff an explicit condition holds (`point0_iff`, `line1_iff`, `pointB_iff`, `line2_iff`, `line0_iff`, `lineB_iff`); `mult_iff` (fold multiplicative ⟺ kernel an image-bimodule) |
| `Discordia/Bridge.lean` | **`isFold_iff_rb`**: folds = idempotent RB(−1) operators commuting with Δ; `fold_star` (δ-D1 for folds); `thm3_isRB`, `thm4_isRB`; δ-L fails for Theorem 4's fold (`thm4_not_dialgebra`, `thm4_not_multiplicative`, `thm4_not_module`) |
| `Discordia/Logic.lean` | 𝔼: seven values, collapse homomorphism to S_fde (L2), L3 (`Γ ⊨ φ ⟺ Γ ⊨_LP φ ∧ Var φ ⊆ Var Γ`), L4, the L5 reduction, Meta-1 (4 of 49) |
| `Discordia/Ops.lean` | M2, N5 (zero counit ⇒ trivial coalgebra), N6 on ℤ/5, β4, B4 (Eckmann–Hilton, and B4 ⇒ no neutral chaos), B1 (countable vs continuum; grid-named set null under i.i.d. Bernoulli products), B3 periods, §9 arithmetic |
| `OrbitReflect.lean`, `Gen/orbits_*.lean` | v2.33: `orbit_sound` (kernel checks ⇒ `Φ_s` is a multiplicative bijection with `Φ_s P_0 = P_s Φ_s`); Theorem 2, F2, F4 are single inner orbits (`thm2_one_orbit`, `thm3F2_one_orbit`, `thm3F4_one_orbit`) |
| `Gen/two_pairs_fam{1..19}.lean` | v2.36: the 19 families of 1⇉2⇉3 in closed form; `fam{k}_isFold` (fold over every commutative ring, every t) and `fam{k}_one_orbit` (one inner orbit) |
| `Gen/two_pairs_rigid.lean` | v2.37: `fam{k}_rigid` (an intertwiner acting on V by `ScalesV` and fixing e₃ forces t = t′) and `fam{k}_rigid_GT` (the same for every class substitution with column sums 1 and class weights), by integer `linear_combination` certificates |
| `Gen/two_pairs_disjoint.lean` | v2.38: `fam{k}_fam{k'}_disjoint` for all 342 ordered pairs (one constant entry of `h·P − P'·h`; no hypotheses on `h = subst (wa • ga) (wb • gb)`), and `fam_orbits_injective` ((family, t) ↦ (G × T²)-orbit is injective) |
| `Gen/two_pairs_behaviour.lean`, `Gen/pairs_in_behaviour.lean` | v2.39: `Filt` (ε∘Φ = ε, Φ(radⁿ) ⊆ radⁿ), closed under products; `subst_filt`; explicit `innM x` with `innM_conj` (it is y ↦ u y u⁻¹) and `innM_filt`; `row_zero`; `fam_behaviours_distinct` (19 behaviours) and `thm3_behaviours_distinct` (F1–F4, t ≠ 0) |
| `AutFilt.lean` | v2.40 (handwritten): `Mult`, `Mult.inv`, `Mult.filt` (every multiplicative invertible Φ on k(1⇉2⇉3), k reduced, is `Filt`), `fam_behaviours_distinct_aut` |
| `AutFiltIn.lean` | v2.41 (handwritten): automorphisms of k(1⇉2⇇3) over a domain: `Mult.filt_or_swap`, `Mult.cases`; the source swap `S` (`S_mult`, `S_sq`, `S_not_filt`) |
| `Gen/pairs_in_aut.lean` | v2.41: mirror tables `S·F_j·S`, swap-case cross lemmas, `thm3_F1_F3`, `thm3_behaviours_aut` |
| `Thm3Behaviours.lean` | v2.41 (handwritten): `Conj` (conjugacy by an automorphism), the F1 witness `conj_F1`, `thm3_exactly_three` |
| `Gen/two_pairs_short_line.lean` | v2.44: a kernel-tilting line P₀ + t·N on 1⇉2⇉3 + 1→3, `line_isFold` (fold over every commutative ring, every t); N is not a derivation direction (computed, `analysis/v2_44_exact.py`) |
| `Gen/two_pairs_short_fam{1..18}.lean` | v2.49: the 18 inner families of 1⇉2⇉3 + 1→3, P_t = Inn(1 + t·x)·P₀·Inn⁻¹; `famk_isFold`, `famk_one_orbit` (all built; audit 184 results, standard axioms only) |
| `Gen/two_pairs_short_curve{A,B}.lean` | v2.50: the slice line C = p + s·v inside the divisor component L_A / L_B; `curveA_isFold`, `curveB_isFold` |
| `Gen/pairs_in_short_vblocks.lean` | v2.52: an integer fold for each of the 21 vertex blocks of 1⇉2⇇3 + 1→3; `vbk_isFold`, k = 1..21 (that there are no others is a Singular elimination, `analysis/v2_52_vertex_ideal.py`) |
| `Gen/K3_tgtfork.lean`, `Gen/K3_tgtfork_rigid.lean` | v2.75 (Theorem 5 on its base quiver 1⇛2, 2→3, 2→4): `fork_isFold`, `forkM_isFold` (the fork families P_t, P⁻_t are folds over every commutative ring, every t); `forkBad_check_fails` (with ℓ₀ = a2 ∉ ker σ, Δ fails); `fork_rigid_GT`, `forkM_rigid_GT` (substK·P_t = P_t'·substK ⇒ t = t', for every census symmetry substK, no hypotheses: a four-term integer certificate); `fork_forkM_disjoint` (no census symmetry joins the two families over a nontrivial ring) |
| `ForkGeneral.lean` | v2.77: **Theorem 5 on every acyclic quiver with a full fork** (vertices i, j, k ≠ k′, three arrows i→j, x : j→k, x′ : j→k′, and no other arrow among these four vertices): `fork_general` (both fork families pushed into k Q are folds, every commutative ring, every t), `fork_general_rigid`, `forkM_general_rigid` (an h acting on the subquiver's span as substK and intertwining P_t, P_t′ forces t = t′), `fork_forkM_general_disjoint` |
| `Transport.lean` | v2.34: `twist` (Δ transported along an algebra automorphism), `twist_delta`, **`transport`**: φPφ⁻¹ is a fold iff P is a fold for Δ^φ |
| `Audit.lean` | `#print axioms`: only `propext`, `Classical.choice`, `Quot.sound` (no `sorry`, no `native_decide`) |

| result (docs) | Lean |
|---|---|
| Theorem 1 / Theorem 2 normal form on 1⇉2⇉3 (v2.26, v2.27) | `Q_two_pairs.thm2_isFold` |
| Theorem 2 fails with the arrow 1→3 | `Q_two_pairs_short.thm2_check_fails` (fails at Δ of the arrow 1→3) |
| Theorem 3, families F1–F4 on 1⇉2⇇3 (v2.27) | `Q_pairs_in.thm3F{1,2,3,4}_isFold`, rigidity `thm3_rigid_GT` |
| Theorem 3 fails with the arrow 1→3 | `Q_pairs_in_short.thm3F{1,2,4}_check_fails` |
| Theorem 4 on 1⇉2⇉3 and on 1⇉2⇉3 + 1→3 (v2.28) | `Q_two_pairs.thm4_isFold`, `Q_two_pairs_short.thm4_isFold`, rigidity `Q_two_pairs.thm4_rigid_GT` |
| Theorem 4 on every acyclic quiver with consecutive classes (v2.28) | `thm4_general`, rigidity `thm4_general_rigid` |
| Theorem 3 on every acyclic quiver with converging classes, `u ≠ w`, no arrow between `u` and `w` (v2.27) | `thm3_general`, rigidity `thm3_general_rigid` |
| Folds are idempotent RB(−1) operators commuting with Δ (v2.31) | `CoordAlg.isFold_iff_rb` |
| All 19 families of 1⇉2⇉3 (v2.36) | `Q_two_pairs.fam{1..19}_isFold`, `fam{1..19}_one_orbit` |
| Rigidity of the 19 families (v2.37) | `Q_two_pairs.fam{1..19}_rigid_GT` |
| Families inject into census orbits (v2.38) | `Q_two_pairs.fam_orbits_injective` |
| Families are distinct behaviours (v2.39) | `Q_two_pairs.fam_behaviours_distinct`, `Q_pairs_in.thm3_behaviours_distinct`, `innM_conj` |
| Every automorphism is filtered; 19 Aut(kQ)-orbits (v2.40) | `Q_two_pairs.Mult.filt`, `Q_two_pairs.fam_behaviours_distinct_aut` |
| Theorem 3's families: exactly three behaviours (v2.41) | `Q_pairs_in.Mult.filt_or_swap`, `Q_pairs_in.thm3_behaviours_aut`, `Q_pairs_in.thm3_exactly_three` |
| A behaviour-varying line of folds (v2.44) | `Q_two_pairs_short.line_isFold` |
| The 18 inner families of 1⇉2⇉3 + 1→3 (v2.49) | `Q_two_pairs_short.famk_isFold`, `Q_two_pairs_short.famk_one_orbit`, k = 1..18 |
| The slice lines in the divisor components (v2.50) | `Q_two_pairs_short.curveA_isFold`, `Q_two_pairs_short.curveB_isFold` |
| The 21 vertex blocks of 1⇉2⇇3 + 1→3 are realised (v2.52) | `Q_pairs_in_short.vb1_isFold` … `vb21_isFold` |
| Transport: fibres over behaviours are coproduct twists (v2.34) | `CoordAlg.transport` |
| Theorem 2, F2, F4 families are single inner orbits (v2.33) | `Q_two_pairs.thm2_one_orbit`, `Q_pairs_in.thm3F{2,4}_one_orbit` |
| δ-D2 exact (v2.32) | `Discordia.Completion.*_iff`, `mult_iff` |
| v2 axioms M2–M4, N5, N6, β4, δ-L, δ-D, δ-D1, δ-D2, E1, E1⋆, B1, B3, B4, L2–L5, Meta-1 (v2.31) | `Discordia/*.lean` |

Not formalized (v2 axioms): see v2.31 §6. Not formalized: Theorem 2 in general (only on 1⇉2⇉3); that the class substitutions of an arbitrary quiver satisfy the
hypotheses of `thm4_general_rigid` (they are stated as hypotheses on `h`); census counts; Conjecture 3ᶜ.
