# Discordianism: Open Questions Explored (v2.23)

*Takes up question 2 of `discordianism-open-questions-explored-v2.22.md` §5: Conjecture 3″ beyond one parallel pair. The first quiver with two pairs, **1⇉2⇉3**, refutes Conjecture 3″ exactly. A larger symmetry group than v2.22 recognized repairs it in every case tested. Also new: a parallel, resumable census (`pcensus.py`) and a batch of long runs prepared for the user's machine (`run_batch.py`).*

*Tags **[Zn]** are machine checks in `verify_v2_23.py` (folder `v2_23/`, zipped as `verify_v2_23.zip`; `pip install -r requirements.txt`; `JOBS=n` sets worker processes). The default run is Z1–Z4; `SLOW=1` adds Z5–Z8. Labels: **proved**, **verified (exact)**, **numerical** or **conjectured**. Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Summary

| # | Question | Answer |
|---|---|---|
| v2.14 Conj. 3″ | Finitely many orbits of Aut_Δ ⋊ grading? | **False (exact).** For 1⇉2⇉3 the fold variety has smooth ℚ-points where it is **6-dimensional** (tangent dimension 6, equal to the dimension of an orbit inside the variety). Aut_Δ ⋊ 𝔾_m has dimension 4 + 1 = 5, from the exact symmetry algebra **[Z6–Z8]** |
| — | What was missing | **The class-weight torus.** Conjugation by φ_w, which scales each parallel class c by w_c, maps Δ-folds to folds compatible with Δ_w = Σ w_c Δ_c. If a fold is compatible with each partial deletion Δ_c separately, every φ_w maps it to a fold. All folds of 1⇉2→3, 1⇛2→3 and 1⇉2⇉3 over 𝔽₃ have this property **[Z5]**, so the whole torus T_Q ≅ 𝔾_m^{#classes} acts. These maps are not Δ-conformal, so the derivation computation of v2.22 (X3) could not see them. v2.22's sentence "nothing hidden enlarges it" is **withdrawn** |
| — | Repair | **Conjecture 3\*: finitely many orbits of Aut_Δ ⋊ T_Q**, where T_Q is the subtorus of class weightings preserving the fold set. Supported for 1⇉2⇉3, where at generic points near every transversal ℚ-fold the local dimension equals the orbit rank of Aut_Δ ⋊ T_Q (numerical) **[Z8]**, and by every orbit count so far |
| v2.22 fork rule | Families iff a same-orientation arrow at an endpoint of a pair? | **Holds for one pair** (also for K₃: 1⇛2→3 has none, with 120 orbits under G, G × grading and G × T_Q; 1⇛2, 1→3 has 8 families, with 191 G-orbits and 183 under both larger groups **[Z5]**). **Fails for two pairs**: 1⇉2⇉3 has Aut_Δ-families with no such arrow. No replacement rule yet |
| — | Tools | `pcensus.py` splits single fibres across processes (task graph), checkpoints, and counts orbits under G, G × grading and G × T_Q. It reproduces every earlier census **[Z1–Z3]** |

---

## 1. The class-weight torus

Let the arrows of Q fall into parallel classes c (single arrows are classes of size 1). Let Δ_c delete only arrows of class c, so Δ = Σ_c Δ_c. For weights w ∈ (k^×)^{classes}, let φ_w scale each path by the product of its arrows' weights. φ_w is an algebra automorphism, and deleting one arrow divides the weight by that arrow's w_c, so

  Δ φ_w = (φ_w ⊗ φ_w) Δ_w,  Δ_w = Σ_c w_c Δ_c.

**Lemma (proved).** If a fold P satisfies Δ_c P = (P ⊗ P) Δ_c for every class c, then φ_w P φ_w⁻¹ is a fold for every w.
*Proof.* Linearity in Δ gives Δ_w P = (P ⊗ P) Δ_w, and Δ_w φ_w⁻¹ = (φ_w⁻¹ ⊗ φ_w⁻¹) Δ. Then Δ(φPφ⁻¹) = (φ ⊗ φ) Δ_w P φ⁻¹ = (φPφ⁻¹ ⊗ φPφ⁻¹) Δ. φ_w maps subalgebras to subalgebras. ∎

The uniform weighting w = (λ, …, λ) is the grading, which preserves every fold (Δ_w = λΔ). T_Q is the subgroup of weightings that preserve the fold set; it always contains the grading.

**Split compatibility, by census over 𝔽₃ [Z5].**

| quiver | folds compatible with each Δ_c | T_Q (preserving weightings over 𝔽₅, where computed) |
|---|---|---|
| 1⇉2→3 | 219 / 219 | all 16 of 𝔾_m²(𝔽₅) |
| 1⇛2→3 | 5439 / 5439 | (𝔽₃ only) |
| 1⇉2⇉3 | 16376 / 16376 | (𝔽₃ only) |
| 1⇉2, 1→3 | **346 / 378**. A fold fails **iff it mixes the two classes beyond ker Δ**: the α-part of P(γ) has nonzero coefficient sum (e.g. P(γ) = α₁ + α₂, P(e₃) = −e₂), or γ occurs in some P(α_i). All 48 family folds are split compatible | **4 of 16**: only the grading, since separate scalings of γ and the pair break the mixing folds |

Single-arrow deletions within a pair are *not* respected (e.g. 23 / 219 for 1⇉2→3), as they must not be, since G mixes the arrows of a class.

For 1⇉2→3, then, T² acts although the derivation algebra of the fold equations is only Lie(G) ⊕ k·E. It changes nothing there (71 orbits under G, G × grading and G × T alike, at q = 5).

---

## 2. 1⇉2⇉3 refutes Conjecture 3″

Census over 𝔽₃: **16376 folds; 814 G-orbits; 728 (G × grading)-orbits; 716 (G × T_Q)-orbits** **[Z6]**. Here Aut_Δ° = G ≅ Aff(1)² has dimension 4 (exact symmetry algebra: Lie Aut_Δ is 4-dimensional, with the grading 5).

**Transversality.** For one representative per (G × grading)-orbit, compare r_G = rank{[D, P] : D ∈ Lie G} with r_H (E added). r_H > r_G means the grading moves P across G-orbits, so a one-parameter family. This holds for 82 representatives mod 3, and for 23 exactly over ℚ among the representatives that lift (which representatives lift depends on enumeration order; the mod-3 count does not). So families exist although neither pair carries a same-orientation arrow.

**The exact refutation [Z7, Z8].** At the transversal ℚ-folds (all of them split compatible):

| (r_G, r_H, r_T, tangent dim) | folds |
|---|---|
| (4, 5, 6, 6) | **576** |
| (4, 5, 5, 5) | 96 |
| (3, 4, 4, 4) | 16 |
| others (tangent > r_T) | 808 |

At the 576 folds of type (4, 5, 6, 6), the orbit of Aut_Δ ⋊ T² (inside the variety, by the Lemma) has dimension 6, and so does the Zariski tangent space. So the fold variety is smooth of dimension 6 there. A 5-dimensional group has orbits of dimension ≤ 5, so **Aut_Δ ⋊ grading has infinitely many orbits on Folds(k(1⇉2⇉3))**. All inputs are exact ranks over ℚ.

**Concretely**, these folds put a vertex constant on each pair: P(b₁) = −e₂ + (b₁ − b₂) and P(a₂) = −e₁ − (a₁ − a₂), and so on. The grading scales both constants together, so their ratio survives as a modulus. The class torus scales them independently and removes it.

**Everything else is explained by T_Q (numerical) [Z8].** At each class with tangent > r_T, a generic nearby fold has local dimension = r_T (5 or 6), and its split residual is below 10⁻⁹. So at every point tested, Aut_Δ ⋊ T_Q accounts for the whole local family.

**Gloss.** *When twin steps follow twin steps, a split can lean on both instants at once, by amounts that a single change of the clock's scale cannot reconcile. Giving each pair of twins its own clock reconciles them.*

---

## 3. The parallel census and the batch

`pcensus.py` reorganizes v2.22's census as a task graph. A task extends one partial fold by one column; workers return survivors, which become new tasks, depth-first. So a single huge fibre (the cause of v2.22's stalled 𝔽₅ run) spreads over all cores. Completed fibres are checkpointed atomically, and an interrupted run resumes (tested: killed at 1043/2100 fibres, resumed to the same 1280 / 428 result). Orbit counting is vectorized over the group (needed for K₃ over 𝔽₅, |G| = 12000). Orbits are counted under G, G × grading and G × T_Q, and the transversality test is built in.

`run_batch.py` lists the long jobs with predictions made before they run:

| job | prediction |
|---|---|
| 1⇉2, 2→3, 2→4 over 𝔽₅ | 428 orbits under G and under G × grading |
| 1⇉2, 1→3, 4→2 over 𝔽₇ | 14928 folds; 615 G-orbits; 485 with grading |
| 1⇉2→3→4 over 𝔽₇ | 5257 folds, 334 orbits |
| 1→2⇉3→4 over 𝔽₇ | 3276 folds, 257 orbits |
| 1⇉2, 1→3 over 𝔽₁₁ | 8474 folds; 137 G-orbits; 101 with grading |
| K₃ + arrow over 𝔽₅; 1⇉2→3⇉4, 1⇉2→3→4→5 over 𝔽₃ | new data |
| 1⇉2⇉3 over 𝔽₅ (heaviest) | whether the (G × T_Q)-orbit count stays at its 𝔽₃ value |

---

## 4. Proposed amendments (v2.22 → v2.23)

| amendment | content |
|---|---|
| **Q1²³ Conjecture 3″** | **Withdrawn**: 1⇉2⇉3 has a 6-dimensional smooth locus against a 5-dimensional Aut_Δ ⋊ 𝔾_m (exact) |
| **Q2²³ Class torus (proved lemma)** | Folds compatible with every class deletion Δ_c are moved to folds by every class weighting. T_Q = weightings preserving Folds; T_Q ⊇ grading |
| **Q3²³ Conjecture 3\*** | Folds(kQ) is a finite union of orbits of Aut_Δ ⋊ T_Q (Q acyclic) |
| **Q4²³ Fork rule** | Restricted to quivers with one parallel class; it fails for 1⇉2⇉3 |
| **Q5²³ Correction to v2.22** | "Nothing hidden enlarges [the symmetry group]" is withdrawn: T_Q can exceed the grading, e.g. T² on 1⇉2→3 |

---

## 5. Open questions

1. **Characterize T_Q.** When is every fold split compatible? It holds for 1⇉2→3, 1⇛2→3 and 1⇉2⇉3. It fails for 1⇉2, 1→3 exactly at the folds mixing the two classes beyond ker Δ (classes with a common source), which suggests T_Q = 𝔾_m^{#classes} exactly when no two classes share an endpoint in the same orientation.
2. **Prove Conjecture 3\*** in the one-pair case, using the fibration over kernel folds.
3. The batch in §3, especially 1⇉2⇉3 over 𝔽₅.
4. Carried: unobstructed violations; Kronecker plus a loop; zigzag; dihedral; ρ.

---

## 6. Sources

- M. Aguiar, *On the associative analog of Lie bialgebras*, J. Algebra 244 (2001).
- I. Assem, D. Simson and A. Skowroński, *Elements of the Representation Theory of Associative Algebras*, Vol. 1, CUP, 2006.
- S. Lang, *Algebraic groups over finite fields*, Amer. J. Math. 78 (1956).

*"When twin steps follow twin steps, a split can lean on both instants at once by amounts no single change of the clock's scale reconciles; give each pair of twins its own clock, and they are reconciled."*
