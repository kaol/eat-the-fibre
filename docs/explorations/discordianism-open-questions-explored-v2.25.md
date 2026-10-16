# Discordianism: Open Questions Explored (v2.25)

*Decides the question v2.24 left open: 1⇉2⇉3 over 𝔽₇ (batch job 11, computed on the user's machine). Conjecture 3\*ˡ
fails. The census engine was optimized in parallel by another instance (kernel stage on the GPU, symmetry transport;
see `DESIGN.md` §6–8); job 11 was computed with the earlier engine and re-finished here with the new one.*

*Raw data: `data/two_pairs_7.pkl.gz` and the last line of `data/results.jsonl`. Scripts: `analysis/stabtype.py`,
`analysis/spansym.py` (new). Labels: **proved**, **verified (exact)**, **numerical**, **conjectured**. Discordian glosses
are descriptive labels (L3): nothing in the mathematics forces them.*

---

## 0. Summary

| # | Question | Answer |
|---|---|---|
| v2.24 §2a | 1⇉2⇉3: ℚ(i)-orbits (𝔽₇ → 716) or families (𝔽₇ → 792)? | **Families. 792** orbits of H = G × T² over 𝔽₇, with exactly the per-type counts predicted from 𝔽₃ and 𝔽₅ before the run (24, 70, 235) **[verified (exact)]** |
| v2.24 Conj. 3\*ˡ | Finitely many classes under Aut_Δ and local class-torus moves? | **False for 1⇉2⇉3** (exact counts at three primes, predicted in advance; no explicit family over ℚ(t) yet). Every fold of this quiver admits the whole class torus, so the local and global tori coincide here |
| — | A hidden symmetry? | **None of the natural kind.** The derivations preserving span{Δ_c} (whose flows map split-compatible folds to folds) form exactly Lie(G) ⊕ ⟨E_c⟩: dimension 6 for 1⇉2⇉3 **[verified (exact)]**, matching 3, 4, 4 for K₂, 1⇉2→3, 1⇉2, 1→3 |
| — | Where the moduli live | In **lower strata**: 19 one-parameter families of H-orbits of dimension 2, 3 and 4 (stabilizers of dimension 4, 3, 2). The open stratum is not affected: its orbit counts are constant |

---

## 1. The counts

H = G × T² (dim 4 + 2). Orbit types by orbit size, the only ones that change:

| orbit size | orbit dim | 𝔽₃ | 𝔽₅ | 𝔽₇ | form |
|---|---|---|---|---|---|
| (q − 1)² | 2 | 20 | 22 | 24 | 18 + (q − 1) |
| q(q − 1)² | 3 | 54 | 62 | 70 | 46 + 4(q − 1) |
| q²(q − 1)² | 4 | 179 | 207 | 235 | 151 + 14(q − 1) |
| all others | | 463 | 463 | 463 | constant |
| **H-orbits** | | **716** | **754** | **792** | **678 + 19(q − 1)** |

For comparison, the other counts at q = 3, 5, 7:
- G-orbits: 814, 1112, 1506. The second difference is 96, so the count is quadratic in q: two-parameter families of G-orbits, matching the two torus directions.
- G × grading orbits: 728, 790, 852, i.e. 31(q − 1) + 666.
- Folds: N = 16376, 251550, 1682084.

Each of the 19 families contributes q − 1 orbits over 𝔽_q, consistent with a 𝔾_m-parameter (a one-parameter family
whose parameter is a nonzero scalar). The stabilizer point counts of these types are those of connected groups at all
three primes (v2.24 §2a for 3 and 5; the 𝔽₇ sizes 36, 252 and 1764 against |H| = 63504 again give |Stab| = q²(q − 1)²,
q(q − 1)² and (q − 1)²).

**Why the Galois reading fails.** Orbits defined over ℚ(i) but not ℚ have 𝔽_q-points only when 4 ∣ q − 1. 7 ≡ 3 (mod
4), like 3, so that reading predicted 716 at q = 7. The count is 792.

---

## 2. No further symmetry of the natural kind

The class torus (v2.23) exists because every fold of these quivers is compatible with each partial deletion Δ_c, so
anything that preserves span{Δ_c} (which contains Δ) maps folds to folds. Its infinitesimal version: derivations D with

  Δ_c D − (D ⊗ 1 + 1 ⊗ D) Δ_c ∈ span{Δ_d}  for every class c.

`analysis/spansym.py` solves this exactly (sympy, over ℚ). Results (`logs/spansym.log`):

| quiver | dim | = Lie Aut_Δ + #classes |
|---|---|---|
| K₂ | 3 | 2 + 1 |
| 1⇉2→3 | 4 | 2 + 2 |
| 1⇉2, 1→3 | 4 | 2 + 2 |
| 1⇉2⇉3 | 6 | 4 + 2 |

So for 1⇉2⇉3 the connected symmetry group of this kind is exactly H, and the 19 families are moduli for it. (A
symmetry of some other kind, not induced by a derivation preserving span{Δ_c}, is not excluded.)

**Gloss.** *When twin steps follow twin steps, giving each pair of twins its own clock is not enough: in the
degenerate splits, where some symmetry already fixes the split, a further continuous parameter remains that no
relabelling and no clock can reach.*

---

## 3. Status of the finiteness conjectures

| conjecture | symmetry | status |
|---|---|---|
| 3 (v1) | fold types finite | superseded (v2.13) |
| 3′ (v2.13) | Aut_Δ | false: 1⇉2, 1→3 (v2.22) |
| 3″ (v2.14) | Aut_Δ ⋊ grading | false: 1⇉2⇉3, exact 6-dimensional smooth locus (v2.23) |
| 3\* (v2.23) | Aut_Δ ⋊ T_Q (global class torus) | false: 1⇉2, 2→3, 2→4 (v2.24) |
| 3\*ˡ (v2.24) | Aut_Δ and local class-torus moves | **false: 1⇉2⇉3 (v2.25)** |

What survives:
- every proved theorem: acyclic quivers without parallel arrows (v2.16), Kronecker (v2.18), and the torus lemmas
  (v2.23, v2.24);
- every *single parallel pair* quiver tested, all of which have finitely many classes under the right torus.

**Conjecture 3ᵖ (one parallel class).** If Q is acyclic with one parallel class, Folds(kQ) has finitely many classes
under Aut_Δ and local class-torus moves. Supported by 1⇉2→3; 1→2⇉3; 1⇉2, 1→3; 1⇉2, 1→3, 4→2; 1⇉2, 2→3, 2→4; the 4-vertex chains
1⇉2→3→4 and 1→2⇉3→4; and K₃ plus an arrow, each at two to four primes.

**Two pairs.** For 1⇉2⇉3 the fold variety has infinitely many orbits under every continuous symmetry found. The natural
next object is the **moduli of the lower strata**: describe the 19 families (one 𝔾_m-parameter each) and find the
invariant that separates their members.

---

## 4. Open questions

1. **Exhibit one of the 19 families over ℚ(t)** and an H-invariant that varies with t (turning the count evidence into a
   proof). Start from the 14 families of 4-dimensional orbits (stabilizer a 2-torus): fix the stabilizer torus S and
   study the fixed-point locus Folds^S under N_H(S).
2. **Conjecture 3ᵖ** for one parallel class: prove it, or find the first counterexample.
3. 1⇉2→3⇉4 (two pairs separated by an arrow): 1158 orbits over 𝔽₃ under every group. 𝔽₃ alone cannot see families
   (they showed up in 1⇉2⇉3 only from 𝔽₅ on), so its 𝔽₅ count decides whether separation by an arrow removes the moduli.
4. Carried: unobstructed violations; Kronecker plus a loop; zigzag; dihedral; ρ.

---

## 5. Sources

- M. Aguiar, *On the associative analog of Lie bialgebras*, J. Algebra 244 (2001).
- S. Lang, *Algebraic groups over finite fields*, Amer. J. Math. 78 (1956) (point counts of homogeneous spaces with
  connected stabilizers).

*"When twin steps follow twin steps, separate clocks for each pair are not enough: in the degenerate splits a
parameter remains that no relabelling and no clock can reach."*
