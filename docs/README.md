# Documents

The project began with a Mu–Nu axiom framework for "Discordianism" drafted by Le Chat (v0). Its open questions became the study of
**folds**: idempotent operators P on a path algebra kQ whose image and kernel are subalgebras and which commute with
arrow deletion Δ (ΔP = (P⊗P)Δ). The census engine at the repository root counts them over finite fields.

## Layout

| folder | content |
|---|---|
| `v0-le-chat/` | Le Chat's v0 axioms, its explorations v1.1–v1.3, and its citation correction |
| `axioms/` | v1 axioms (refinement of v0, machine-checked), the review of v1.1–v1.3, and the v2 axioms |
| `explorations/` | `discordianism-open-questions-explored-v2.N.md`, one per round (v2.1–v2.24) |
| `../verify/` | the machine checks cited by the documents: `axioms/` (v1, v2), `v2_1-5/`, then one folder `v2_N/` per round with `verify_v2_N.py` |
| `../data/` | raw data of the v2.23 batch: the census checkpoints and their results (see `data/README.md`) |
| `../analysis/` | analysis scripts used from v2.23 on (tangent spaces, local dimension, class torus, orbit checks); they run against the root engine |

Tags such as **[Zn]** in a document refer to checks in the `verify_v2_N.py` of the same round. Labels used throughout:
**proved**, **verified (exact)**, **numerical**, **conjectured**. Discordian glosses are descriptive labels: nothing in the
mathematics forces them.

## Where the fold story stands (v2.24)

| result | status | round |
|---|---|---|
| acyclic quivers without parallel arrows: every fold is canonical, so finitely many folds | proved | v2.16 |
| Kronecker K_m: 14m − 8 orbits of Aut_Δ, point count, dimension ⌊m²/2⌋ | proved | v2.18 |
| Conjecture 3′ (finitely many Aut_Δ-orbits, acyclic Q) | **false**: 1⇉2, 1→3 | v2.22 |
| Conjecture 3″ (Aut_Δ ⋊ grading) | **false**: 1⇉2⇉3 (exact, 6-dimensional smooth locus) | v2.23 |
| class-weight torus lemma (split-compatible folds are moved to folds by every class weighting) | proved | v2.23 |
| Conjecture 3\* (Aut_Δ ⋊ T_Q, global torus) | **false**: 1⇉2, 2→3, 2→4 (456 orbits over 𝔽₅, 428 over 𝔽₃) | v2.24 |
| Conjecture 3\*ˡ (Aut_Δ and local class-torus moves) | holds in every count except 1⇉2⇉3 (716 over 𝔽₃, 754 over 𝔽₅); 𝔽₇ decides between ℚ(i)-orbits (716) and families (792) | v2.24 |
