# Documents

The project began with a Mu–Nu axiom framework for "Discordianism" drafted by Le Chat (v0). Its open questions became the study of
**folds**: idempotent operators P on a path algebra kQ whose image and kernel are subalgebras and which commute with
arrow deletion Δ (ΔP = (P⊗P)Δ). The census engine at the repository root counts them over finite fields.

## Layout

| folder | content |
|---|---|
| `v0-le-chat/` | Le Chat's v0 axioms, its explorations v1.1–v1.3, and its citation correction |
| `axioms/` | v1 axioms (refinement of v0, machine-checked), the review of v1.1–v1.3, and the v2 axioms |
| `explorations/` | `discordianism-open-questions-explored-v2.N.md`, one per round (v2.1–v2.54) |
| `../verify/` | the machine checks cited by the documents: `axioms/` (v1, v2), `v2_1-5/`, then one folder `v2_N/` per round with `verify_v2_N.py` |
| `../data/` | raw data of the v2.23 batch: the census checkpoints and their results (see `data/README.md`) |
| `../lean/` | Lean 4 + Mathlib formalization (v2.29–v2.44): basis criterion, conjugation and vertex-block lemmas, Theorems 1–4 on their base quivers, Theorems 3 and 4 on every acyclic quiver, rigidity, controls; the v2 axioms derived (`Discordia/`); see `lean/README.md` |
| `../analysis/` | analysis scripts used from v2.23 on (tangent spaces, local dimension, class torus, orbit checks); they run against the root engine. `completion_*.py` (v2.31) read the generated Lean tables; `module_classes.py`, `completion_words.py` (v2.32) classify census orbits from `--save-reps` dumps (`logs/v2_32/`); `behaviour_census.py`, `family_behaviour.py` (v2.33); `fibre_rigidity.py`, `v2_34/` (v2.34); `family_signatures.py`, `v2_35/` (v2.35); `family_normal_forms.py`, `family_curves.py`, `v2_36_*.py` (v2.36); `v2_37_certificates.py` (v2.37); `v2_38_*.py` (v2.38); `v2_42_*.py` (v2.42); `v2_44_*.py` (v2.44); `v2_45_*.py` (v2.45); `v2_46_*.py` (v2.46); `v2_47_divisor.py` (v2.47); `v2_48_divisor_components.py` (v2.48; env `T` v2.49); `v2_50_curve_slice.py`, `v2_50_invariants.py` (v2.50); `v2_51_profile.py` (v2.51); `v2_52_vertex_blocks.py`, `v2_52_vertex_ideal.py`, `v2_52_realise.py` (v2.52); `v2_53_block_census.py`, `v2_53_stratum.py`, `v2_53_orbits.py` (v2.53); `v2_54_orbit_strata.py` (v2.54) |

Tags such as **[Zn]** in a document refer to checks in the `verify_v2_N.py` of the same round. Labels used throughout:
**proved**, **verified (exact)**, **numerical**, **conjectured**. Discordian glosses are descriptive labels: nothing in the
mathematics forces them.

## Where the fold story stands (v2.54)

| result | status | round |
|---|---|---|
| acyclic quivers without parallel arrows: every fold is canonical, so finitely many folds | proved | v2.16 |
| Kronecker K_m: 14m − 8 orbits of Aut_Δ, point count, dimension ⌊m²/2⌋ | proved | v2.18 |
| Conjecture 3′ (finitely many Aut_Δ-orbits, acyclic Q) | **false**: 1⇉2, 1→3 | v2.22 |
| Conjecture 3″ (Aut_Δ ⋊ grading) | **false**: 1⇉2⇉3 (exact, 6-dimensional smooth locus) | v2.23 |
| class-weight torus lemma (split-compatible folds are moved to folds by every class weighting) | proved | v2.23 |
| Conjecture 3\* (Aut_Δ ⋊ T_Q, global torus) | **false**: 1⇉2, 2→3, 2→4 (456 orbits over 𝔽₅, 428 over 𝔽₃) | v2.24 |
| Conjecture 3\*ˡ (Aut_Δ and local class-torus moves) | **false** (theorem): 1⇉2⇉3 has the explicit family P_t with H-invariant t; counts 678 + 19(q − 1) | v2.25–v2.26 |
| consecutive parallel classes u⇉v⇉w (u, w distinct, non-adjacent): explicit family with invariant w′/(xz) | proved (Theorem 2) | v2.27 |
| converging or diverging classes u⇉v⇇w: three families whose invariant is the vertex block (π₀P on kQ₀) | proved (Theorem 3); 1⇉2⇇3 has 449 + 4(q − 2) orbits at q = 3, 5, 7, 11 | v2.27 |
| Conjecture 3ᵖ (3\*ˡ for one parallel class) | supported by every one-class quiver tested | v2.25 |
| Conjecture 3ᵛ (finite iff no two classes of size ≥ 2 share a vertex) | **false**: 1⇉2⇇3 + 1→3 (366 orbits over 𝔽₃ and 𝔽₅) | v2.27 |
| Conjecture 3ʰ (infinite iff Theorem 2 or 3 applies) | **false**: 1⇉2⇉3 + 1→3 grows (by Theorem 4; with all of Aut_Δ°, 939 → 983 over 𝔽₃, 𝔽₅) | v2.27, v2.28 |
| consecutive classes on any acyclic quiver: family vanishing on all vertices, invariant read off the image | proved (Theorem 4) | v2.28 |
| Theorems 1–4 on their base quivers, with rigidity under class substitutions and weightings; basis criterion; vertex-block lemma | **Lean 4 (kernel-checked)** | v2.29 |
| Theorems 3 and 4 on every acyclic quiver (path algebras, reduction to the two-class subalgebra), with rigidity | **Lean 4 (kernel-checked)** | v2.30 |
| folds = idempotent Rota–Baxter(−1) operators commuting with Δ (v2 open question 2); δ-D, δ-D1, E1⋆ hold for folds, δ-L fails; δ-D2 completions classified by how ker sits over im | **Lean 4** (sufficiency); necessity generic (computed) | v2.31 |
| δ-D2 exact: each completion stratum present iff an explicit condition (point 0 ⟺ ker·im ⊆ ker; line 2/point B ⟺ P multiplicative ⟺ ker an im-bimodule); census orbits classified by module class, growth split by class | **Lean 4** (iff); census computed | v2.32 |
| behaviour census (folds up to all algebra automorphisms): every Theorem 2–4 family is one behaviour on its generic range (explicit witnesses); behaviour counts independent of q on 1⇉2, 1→3 (91) and 1⇉2⇇3 (319) | Theorem 2, F2, F4: **Lean 4**; others verified exact / computed | v2.33 |
| transport: the fibre over a behaviour is the set of coproduct twists Δ^φ compatible with P (Conjecture 3ᶜ becomes a statement about twist varieties); first-order rigidity refuted, second-order fibre excess discriminates the obstruction on 3 quivers | **Lean 4** (transport); computed | v2.34 |
| the 19 families of 1⇉2⇉3 located by Aut(kQ)-invariant signatures (9 groups: 9 + 3 unclassed, 3 L, 3 R, 1 LRKE), all carried by the canonical derivations ad δ_a, ad δ_b, ad δ_aδ_b and the class scalings; joined-ends quiver has no growing behaviour type | computed (partly pre-registered) | v2.35 |
| all 19 families of 1⇉2⇉3 in closed form P_t = Inn(1 + t·x) P₀ Inn(1 + t·x)⁻¹, x ∈ {δ_a, δ_b, δ_aδ_b}; local-torus count 659 + 19q at q = 5, 7; joined-ends quiver has no family-carrying inner curve | **Lean 4** (each family a fold for all t, one inner orbit); distinctness computed | v2.36 |
| each of the 19 families is rigid: distinct t give distinct (G × T²)-orbits, over every commutative ring (integer certificates of degree ≤ 2) | **Lean 4** | v2.37 |
| the 19 families are pairwise disjoint over every nontrivial commutative ring, and (family, t) ↦ census orbit is injective (≥ 19q orbits over 𝔽_q) | **Lean 4** | v2.38 |
| Aut(kQ)-invariant depth vector of the vertex rows εP: the 19 families are 19 distinct behaviours; on 1⇉2⇇3, F2 is its own behaviour and F1, F3 are mirror images under the source swap; τ-strata: census growth is family growth stratum by stratum | computed (partly pre-registered) | v2.38 |
| the 19 families of 1⇉2⇉3 are 19 distinct behaviours for automorphisms generated by Inn(1 + rad) and class substitutions (both proved filtered; Inn(1 + x) proved to be conjugation); F1–F4 of 1⇉2⇇3 likewise modulo Inn ⋊ Sub, for t ≠ 0 | **Lean 4** | v2.39 |
| every algebra automorphism of k(1⇉2⇉3) is filtered (k reduced), so the 19 families are 19 distinct Aut(kQ)-orbits with no generation hypothesis | **Lean 4** | v2.40 |
| every automorphism of k(1⇉2⇇3) over a domain is filtered or filtered after the source swap; Theorem 3's F1–F4 are exactly three behaviours (F1 ≅ F3 via the swap) for t ∉ {0, −1} | **Lean 4** | v2.41 |
| 1⇉2⇉3 + 1→3 (full Δ-compatible group, slope 22): 18 inner-curve families found and cross-checked at q = 3, 5; 4 missing, in two large-image (dim im 7) groups, not inner or shear curves; q = 7 handed to the GPU (job 16) | computed (5 predictions failed or partly) | v2.42 |
| GPU: 1⇉2⇇3 + 1→3 constant 366 at q = 3, 5, 7, 11; 1⇛2⇇3 + 1→3 constant 752 at q = 3, 5; 1⇉2⇉3 + 1→3: unmerged 1668, 1840, 2028 not linear, but merged with the shear 939, 983, 1027 = 873 + 22q: 18 inner families + 4 open large-image families | computed (pre-registered; two held, one failed) | v2.43 |
| 1⇉2⇉3 + 1→3: the 4 large-image families are kernel-tilting lines P₀ + t·u·wᵀ (image fixed), not automorphism orbits: N ∉ [Der, P₀] over ℚ, so in characteristic 0 they carry infinitely many behaviours; one line a fold over every ring | **Lean 4** (fold) + exact over ℚ (non-derivation) | v2.44 |
| the folds with that fixed image form an irreducible 7-dimensional variety V_A over ℚ (Singular), smooth at the base point; automorphisms preserving A cover 6 of 7 tangent directions generically, so the behaviours with image A form a ≥ 1-dimensional family (= 1 to first order) | computed exactly / mod p (7 predictions: 3 held, 1 in part, 3 failed) | v2.45 |
| dominant slice: V_A geometrically irreducible of dim 7; Δ-compatible automorphisms and admissible weightings preserving A, applied to the line, give a dominant map onto V_A; each Aut(kQ)-orbit meets the line finitely, so generically the behaviours with image A form exactly a one-parameter family (char 0); one image per signature group | **proved** (char 0, from exact ℚ computations + Lean lemmas) | v2.46 |
| mirror theorem for the R image B; the rank-drop locus of the stabiliser on V_A is a divisor (dim 6) with 5-dimensional orbits, giving the second one-parameter family per group (2q growth) | Singular over ℚ / exact ranks | v2.47 |
| the divisor D has a unique divisorial component; at a smooth rational point one tangent direction of D is not a derivation direction, so D also carries a ≥ 1-dimensional family of behaviours | Singular / exact ranks over ℚ | v2.48 |
| the R divisor D_B computed directly: unique 6-dimensional component, smooth at the line's point t = −1 with one non-derivation tangent direction; the 18 inner families of 1⇉2⇉3 + 1→3 as closed forms over ℤ[t], each a fold over every ring and one inner orbit | Singular / exact over ℚ; **Lean 4** (all 18 built and audited) | v2.49 |
| V_A, V_B are affine-linear 7-spaces and their divisor components affine 6-spaces L_A, L_B; on each L the behaviours form **exactly** a one-parameter family (curve slice: Stab(L) has Lie algebra 𝔥′, dominant onto L from a line C; the generic Der-orbit dimension 21 gives the lower bound; the invariant κ restricts to 1 − s on C) | **proved** (char 0: exact ℚ computations + Lean `curveA_isFold`, `curveB_isFold` + paper algebraic-group steps) | v2.50 |
| 1⇉2⇇3 + 1→3 + 2→4 (P4, job 14): local-torus count 1220 at q = 3 and q = 5. Job 14's OOM was the column order: a per-fibre order portfolio (fibre_torch.portfolio) finishes it in the container | computed (pre-registered; job 14 reproduces job 13's 1220) | v2.51 |
| joined ends freeze Theorem 3's modulus: on 1⇉2⇇3 + 1→3 (char 0) the vertex block of a fold takes exactly 21 values, all 0/1, each realised; on 1⇉2⇇3 the same conditions leave a 1-dimensional set (census: 4q + 19 blocks) | **proved** (Singular elimination over ℚ + Lean realisations) | v2.52 |
| Folds(1⇉2⇇3 + 1→3) in char 0: 21 vertex-block strata, 118 irreducible components (top ones absolutely irreducible), each with a dense orbit of the connected census group; the 366 local-torus orbits split by block identically at q = 3, 5, 7. The boundary recursion remains for lemma (L) | Singular over ℚ + exact ranks | v2.53 |
| **Folds(k̄(1⇉2⇇3 + 1→3)) is a finite union of census orbits (char 0)**: on each of the 118 components the orbit-dimension loci satisfy dim Z_r ≤ r. Lemma (L) in strong form for the joined-ends configuration; uniformity in q for large p by spreading out (sketched) | **proved** (char 0; Singular + exact minors over ℚ) | v2.54 |
| Conjecture 3ᶜ (finite iff no consecutive classes and no converging/diverging pair with unjoined ends; the "if" half of the converse is Theorems 3 and 4) | conjectured; 1⇉2⇇3 + 1→3 has 366 orbits at q = 3, 5, 7; two tests handed over | v2.28 |
