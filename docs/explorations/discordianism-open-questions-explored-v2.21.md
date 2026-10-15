# Discordianism: Open Questions Explored (v2.21)

*Takes up questions 2 and 3 of `discordianism-open-questions-explored-v2.20.md` §4 for the quiver 1⇉2→3: whether Aut_Δ alone is enough, and the point count. The user's exhaustive 𝔽₅ run supplied the first out-of-sample data point.*

*Tags **[Wn]** are machine checks in `verify_v2_21.py` (folder `v2_21/`, zipped as `verify_v2_21.zip`; `pip install -r requirements.txt`). The default run (W1–W7, about 15 minutes) includes the 𝔽₅ census; `SLOW=1` adds 𝔽₇ and `SLOW=2` adds 𝔽₁₁. Labels: **proved**, **verified (exact)** or **conjectured**. Discordian glosses are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Summary

| # | Question | Answer |
|---|---|---|
| v2.20 Q3 | Point count of the fold variety of 1⇉2→3 | **N(q) = 21q² + 6q + 12.** The formula was derived from 𝔽₃ data alone and then confirmed by exhaustive censuses: <br>• 𝔽₅: 567 folds, by the user's 9.7-hour run and by an independent column-staged census **[W5]**; <br>• 𝔽₇: 1083 folds **[W8]**; 𝔽₁₁: 2619 folds **[W9]** |
| v2.20 Q2 | Is Aut_Δ alone enough, or is grading needed? | **Aut_Δ alone (strong evidence, not proof).** The number of G-orbits is exactly 71 at q = 3, 5, 7 and 11. A family of orbits would make that number grow with q |
| — | Orbit structure | **Conjecture 3‴ (1⇉2→3).** The fold variety is a disjoint union of 71 orbits of G ≅ Aff(1): <br>• 23 fixed points; <br>• 11 orbits ≅ G/U ≅ 𝔾_m (unipotent stabilizer); <br>• 16 orbits ≅ G/T ≅ 𝔸¹ (torus stabilizer); <br>• 21 free orbits ≅ G. <br>All stabilizers are connected, so by Lang's theorem each orbit has \|G(𝔽_q)/H(𝔽_q)\| points, and 23 + 11(q − 1) + 16q + 21q(q − 1) = 21q² + 6q + 12. The orbit sizes match this exactly at q = 3, 5, 7 and 11 **[W1, W5, W8, W9]**, and over ℚ every orbit has a representative with the predicted stabilizer **[W3]** |
| — | Base of the fibration | The algebra folds of K = ker Δ number **24q² + 46q + 30** (checked at q = 3, 5, 7, 11, 13) **[W4]** |

---

## 1. Why the 𝔽₅ run was decisive

The user's 𝔽₅ run (ChatGPT's resumable extension of `q123v.py`, about 9.7 hours) found **567 folds in 71 orbits of G(𝔽₅), |G| = 20**. I checked its code: it re-chunks and checkpoints the enumeration without changing what is enumerated or checked, and it reproduces 219 / 71 at p = 3. Its `[F5] PASS` line means only that the run completed.

To read 567, look at G. The column-sum-1 substitutions act on D = kδ by λ = a − b and trivially on V/D, so **G ≅ Aff(1) = {x ↦ λx + c}**, of order q(q − 1). Its algebraic subgroups up to conjugacy are 1, U (translations), T (a torus) and G. The connected ones give orbit sizes q(q − 1), q − 1, q and 1. Over 𝔽₃ the 71 orbits have sizes 1, 2, 3 and 6 (23, 11, 16 and 21 orbits). Size 3 is ambiguous at q = 3: it could be a torus stabilizer (q) or μ₂ (q(q − 1)/2). Two checks settle it:
- **Over ℚ [W3].** Every one of the 71 orbits has a member that lifts to an exact ℚ-fold. Its stabilizer in Lie(G) is 2-dimensional, a unipotent line, a torus line or 0, exactly as its 𝔽₃ size requires. So the size-3 orbits have torus stabilizers, not μ₂.
- **Arithmetic.** Of all mixes of torus and μ₂ readings for the 16 size-3 orbits, only the all-torus one gives 567 at q = 5.

So the formula was a genuine prediction. The 𝔽₅ run tested it and passed. The orbit-size distribution {1: 23, 4: 11, 5: 16, 20: 21}, also predicted, was confirmed by my independent census **[W5]**. (`f5_orbit_sizes.py` checks the same distribution directly from the user's checkpoint.)

---

## 2. A faster census

`q123v.py` enumerated the kernel components of (α₁, β) jointly and then those of (α₁β, α₂β) jointly: p⁸ candidates for each survivor of each stage. Every condition is a masked necessary condition on the columns known so far, and only the last stage checks everything. So the stages can be split **one column at a time** (p⁴ candidates each) without losing completeness (`q123f.py`). A second reduction uses the symmetry: translations fix K pointwise and the torus acts on K by δ ↦ λδ. So only one kernel fold per torus orbit needs a fibre, and G-orbits of folds are counted as Stab_G(P_K)-orbits on that fibre (`q123s.py`). Finally, `q123s2.py` computes products of the fixed columns once per batch and drops failed candidates after each check. All three reproduce 𝔽₃ exactly **[W1, W2]**; the times below are for `q123s2.py`.

| field | folds | orbits | orbit sizes | non-filtered | time |
|---|---|---|---|---|---|
| 𝔽₃ | 219 | 71 | 23×1, 11×2, 16×3, 21×6 | 2 | 2 s |
| 𝔽₅ | 567 | 71 | 23×1, 11×4, 16×5, 21×20 | 4 | 14 s (original method: 9.7 h) |
| 𝔽₇ | 1083 | 71 | 23×1, 11×6, 16×7, 21×42 | 6 | 2 min |
| 𝔽₁₁ | 2619 | 71 | 23×1, 11×10, 16×11, 21×110 | 10 | 76 min |

The two folds of v2.20 that are not length-filtered (β ↦ β + δβ) form **one orbit with unipotent stabilizer**: there are q − 1 of them (2, 4, 6) **[W6]**.

**Fibres [W2, W5, W8, W9].** Over the 24q² + 46q + 30 kernel folds, the fibre sizes are, exactly at q = 3, 5, 7 and 11:

| fibre size | number of kernel folds | contributes |
|---|---|---|
| 0 | 24q² + 19q + 18 | 0 |
| 1 | 11q + 12 | 11q + 12 |
| q | 16q − 2 | 16q² − 2q |
| q² | 1 (P_K: e₁ ↦ 0, the rest fixed) | q² |
| 4q² − 3q | 1 (P_K: vertices fixed, δ ↦ 0) | 4q² − 3q |

Together: 21q² + 6q + 12. The special fibre is **four free G-orbits and one orbit of size q** (4q(q − 1) + q), checked at q = 5 and 7 **[W7, W8]**. It is not a union of affine planes: I tested "four planes through a common line", and no point of the fibre is joined to all the others by lines inside it.

---

## 3. What is and isn't established

- **Exact:** the counts, orbit-size distributions and fibre histograms at q = 3, 5, 7 and 11; the stabilizer types over ℚ of one representative per 𝔽₃ orbit; the base counts at five primes.
- **Conjectured (3‴):** the 71-orbit decomposition as a statement about the fold *variety* over any field of characteristic ≠ 2. Constant orbit counts and exact polynomial point counts at four primes (and fibre histograms matching five separate polynomials) are strong evidence. They are not a proof. A proof would classify the fibres over the 24q² + 46q + 30 kernel folds by hand, as v2.18 did for Kronecker.
- **Consequence for Conjecture 3′:** the data say that for 1⇉2→3 the symmetries respecting Δ act with finitely many orbits by themselves; rescaling by path length is not needed. After the Kronecker theorem (v2.18), this is the second parallel-arrow case where Conjecture 3′ is supported, here by evidence rather than proof.

**Gloss.** *When twin steps are followed by a further step, the splits fall into 71 kinds, and every kind is visited by the twins' relabellings alone: 23 kinds the relabelling cannot move, 11 it can only rescale, 16 it can only shift, and 21 it moves freely.*

---

## 4. Proposed amendments (v2.20 → v2.21)

| amendment | content |
|---|---|
| **Q1²¹ Point count** | 1⇉2→3 has N(q) = 21q² + 6q + 12 folds (verified at q = 3, 5, 7, 11) in 71 Aut_Δ-orbits |
| **Q2²¹ Conjecture 3‴** | The fold variety of 1⇉2→3 is a union of 71 Aff(1)-orbits with connected stabilizers (23 / 11 / 16 / 21 of types G / U / T / 1). This supersedes v2.20's "modulo grading" |

---

## 5. Open questions

1. **Prove Conjecture 3‴** by classifying the fibres over the kernel folds; the base count 24q² + 46q + 30 suggests that the base also has a finite orbit decomposition, under Aut(K).
2. **General parallel-arrow quivers.** Kronecker K₂ and 1⇉2→3 share G ≅ Aff(1), with 20 and 71 orbits. Is there a uniform orbit count for "a parallel pair followed by a path of length n"?
3. Carried: the combinatorial criterion for unobstructed violations; Kronecker plus a loop; zigzag; dihedral; ρ.

---

## 6. Sources

- M. Aguiar, *On the associative analog of Lie bialgebras*, J. Algebra 244 (2001).
- S. Lang, *Algebraic groups over finite fields*, Amer. J. Math. 78 (1956) (Lang's theorem: for connected H, G(𝔽_q) acts transitively on (G/H)(𝔽_q)).
- W. Decker, G.-M. Greuel, G. Pfister and H. Schönemann, *Singular* (via passagemath-singular).

*"When twin steps are followed by a further step, the splits fall into 71 kinds: 23 the relabelling of twins cannot move, 11 it only rescales, 16 it only shifts, and 21 it moves freely."*
