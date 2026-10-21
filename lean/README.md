# Lean 4 formalization (v2.29–v2.36)

Lean 4 (v4.35.0) with Mathlib (v4.35.0, pinned in `lake-manifest.json`). Mathlib is used from its prebuilt cache; `.lake/`
(packages, oleans, ~7 GB) is git-ignored and never committed.

```bash
curl -sSfL https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y --default-toolchain none
cd lean && ./build.sh        # lake exe cache get, then each module in turn (the kernel checks need ~5 GB each)
python3 gen/gen.py           # regenerate EatTheFibre/Gen/*.lean from the census engine (pcensus.build)
python3 gen/gen_orbit.py     # regenerate EatTheFibre/Gen/orbits_*.lean (v2.33)
python3 gen/gen_families.py NF.pkl REPS.pkl   # regenerate Gen/two_pairs_fam*.lean (v2.36; inputs in logs/v2_36, logs/v2_32)
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
| Transport: fibres over behaviours are coproduct twists (v2.34) | `CoordAlg.transport` |
| Theorem 2, F2, F4 families are single inner orbits (v2.33) | `Q_two_pairs.thm2_one_orbit`, `Q_pairs_in.thm3F{2,4}_one_orbit` |
| δ-D2 exact (v2.32) | `Discordia.Completion.*_iff`, `mult_iff` |
| v2 axioms M2–M4, N5, N6, β4, δ-L, δ-D, δ-D1, δ-D2, E1, E1⋆, B1, B3, B4, L2–L5, Meta-1 (v2.31) | `Discordia/*.lean` |

Not formalized (v2 axioms): see v2.31 §6. Not formalized: Theorem 2 in general (only on 1⇉2⇉3); that the class substitutions of an arbitrary quiver satisfy the
hypotheses of `thm4_general_rigid` (they are stated as hypotheses on `h`); census counts; Conjecture 3ᶜ.
