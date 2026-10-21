# Documents

The project began with a Mu–Nu axiom framework for "Discordianism" drafted by Le Chat (v0). Its open questions became the study of
**folds**: idempotent operators P on a path algebra kQ whose image and kernel are subalgebras and which commute with
arrow deletion Δ (ΔP = (P⊗P)Δ). The census engine at the repository root counts them over finite fields.

## Layout

| folder | content |
|---|---|
| `v0-le-chat/` | Le Chat's v0 axioms, its explorations v1.1–v1.3, and its citation correction |
| `axioms/` | v1 axioms (refinement of v0, machine-checked), the review of v1.1–v1.3, and the v2 axioms |
| `explorations/` | `discordianism-open-questions-explored-v2.N.md`, one per round (v2.1–v2.35) |
| `../verify/` | the machine checks cited by the documents: `axioms/` (v1, v2), `v2_1-5/`, then one folder `v2_N/` per round with `verify_v2_N.py` |
| `../data/` | raw data of the v2.23 batch: the census checkpoints and their results (see `data/README.md`) |
| `../lean/` | Lean 4 + Mathlib formalization (v2.29–v2.34): basis criterion, conjugation and vertex-block lemmas, Theorems 1–4 on their base quivers, Theorems 3 and 4 on every acyclic quiver, rigidity, controls; the v2 axioms derived (`Discordia/`); see `lean/README.md` |
| `../analysis/` | analysis scripts used from v2.23 on (tangent spaces, local dimension, class torus, orbit checks); they run against the root engine. `completion_*.py` (v2.31) read the generated Lean tables; `module_classes.py`, `completion_words.py` (v2.32) classify census orbits from `--save-reps` dumps (`logs/v2_32/`); `behaviour_census.py`, `family_behaviour.py` (v2.33); `fibre_rigidity.py`, `v2_34/` (v2.34); `family_signatures.py`, `v2_35/` (v2.35) |

Tags such as **[Zn]** in a document refer to checks in the `verify_v2_N.py` of the same round. Labels used throughout:
**proved**, **verified (exact)**, **numerical**, **conjectured**. Discordian glosses are descriptive labels: nothing in the
mathematics forces them.

## Where the fold story stands (v2.35)

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
| Conjecture 3ᶜ (finite iff no consecutive classes and no converging/diverging pair with unjoined ends; the "if" half of the converse is Theorems 3 and 4) | conjectured; 1⇉2⇇3 + 1→3 has 366 orbits at q = 3, 5, 7; two tests handed over | v2.28 |
