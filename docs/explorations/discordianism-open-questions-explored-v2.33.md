# Discordianism: Open Questions Explored (v2.33): the behaviour census

*Le Chat proposed, by analogy with observational type theory (TTobs), to separate the **realiser** (a fold matrix P)
from its **behaviour** (what P exhibits). The census would then be taken on the space of behaviours, and the
parametric families would become wide fibres over single behaviours. This round makes the proposal precise, corrects
two points in it, and tests it. Predictions were committed before the census was run
(`logs/v2_33/PREDICTIONS.md`).*

---

## 0. Summary

| # | Statement | Status |
|---|---|---|
| 1 | P is determined by the *subspaces* im P, ker P, so behaviour must forget the embedding. The natural choice: behaviour = P up to conjugation by all algebra automorphisms of kQ, not only the Δ-compatible group G | definition |
| 2 | Le Chat's pair (im P ↪ kQ, kQ ↠ kQ/ker P) is mostly redundant. kQ/ker P is an algebra only when ker P is an ideal (class E). In general it is a quotient *coalgebra*, which P identifies with im P | proved (one line, §1) |
| 3 | **Every Theorem 2–4 family is a single behaviour on its generic range.** Explicit automorphisms φ = Inn_u ∘ T carry P_t to P_t′ (table in §2). Theorem 2, F2 and F4 need only inner automorphisms, along a one-parameter unipotent group whose time is t | verified (exact, ℚ(t, t′)), `analysis/family_behaviour.py`; reachability over 𝔽₅ by exhaustive search |
| 4 | The degenerate parameters are exactly the poles of φ: t = 1 for Theorem 4, t ∈ {0, −1} for F1 and F3. At t = 1 (Theorem 4) and t = 0 (F1, F3) the image is not semisimple; t = −1 (F1, F3) is a different embedding of the same algebra. The census ranges q − 1 (Theorem 4) and q − 2 (Theorem 3) are exactly the generic ranges | computed |
| 5 | Behaviour census: 1⇉2, 1→3 has **91 behaviours over both 𝔽₃ and 𝔽₁₁**; its G-orbit growth (105 → 137) lies entirely in the fibres. 1⇉2⇇3 has **319 behaviours (170 with the 1↔3 quiver automorphism) over both 𝔽₃ and 𝔽₅**, with 485 → 557 G-orbits | computed, `analysis/behaviour_census.py`; all four pre-registered predictions hold |
| 6 | Module class (v2.32) is constant on behaviours. Theorem 3's vertex block is **not**: it is the coordinate along the four wide fibres of 1⇉2⇇3 | computed |
| 7 | Theorem 2, F2 and F4 are single orbits of `s ↦ Inn(1 + s·w)` over every commutative ring and every s | **proved (Lean)**, `thm2_one_orbit`, `thm3F2_one_orbit`, `thm3F4_one_orbit` (`Gen/orbits_*.lean`, `OrbitReflect.lean`); audit 77 results, standard axioms only |

## 1. What a behaviour is

An idempotent P is the projection onto im P along ker P, so P ↦ (im P, ker P) is injective. Le Chat wrote "an
idempotent is exactly a choice of section for this data, so the fold space fibers over behaviour space by
construction". That holds only after the pair is taken *up to isomorphism*: the fibres are the different ways of
placing the same data inside kQ.

Which isomorphisms? The pair's own structure answers it:
- im P is a subalgebra, and also a subcoalgebra for Δ, because Δ(Px) = (P⊗P)Δx.
- ker P is a subalgebra and a coideal, but not in general an ideal (that is class E, the endomorphisms of v2.32).
  So kQ/ker P is a quotient coalgebra, not an algebra, and P induces kQ/ker P ≅ im P as coalgebras. The second
  component of Le Chat's pair adds nothing to the first.
- **Behaviour** therefore means the splitting (im P, ker P) up to ambient isomorphism, i.e. the Aut(kQ)-orbit of P.
  The census already quotients by G = Aut(kQ) ∩ {Δ-compatible} (the column-sum-1 arrow substitutions). The behaviour
  census quotients by all of Aut(kQ) = Inn ⋊ (Sub ⋊ Aut(Q)), where:
  - Inn = conjugation by 1 + r, r ∈ rad kQ;
  - Sub = arrow substitutions within each parallel class, with any invertible matrix;
  - Aut(Q) = quiver automorphisms (reported separately).

  This decomposition holds when no path of length ≥ 2 runs parallel to an arrow, which is true for every quiver used
  here.
- Coarser choices (abstract isomorphism classes of im P and ker P separately) were also computed for the families.
  Their invariants are constant along each family except at the degenerate points (§2), so they give the same
  picture with less resolution.

In TTobs terms: realiser = P, behaviour = Aut(kQ)·P, and the fibre is the set of Δ-compatible positions of one
splitting. **The parameter t of every family is a modulus of position relative to Δ, invisible to the algebra.**

## 2. The families collapse, with explicit witnesses

With δ_a = a₁ − a₂ and δ_b = b₁ − b₂ (the two parallel classes), φ = Inn_u ∘ T with φ P_t φ⁻¹ = P_t′:

| family | arrow torus T | u | poles |
|---|---|---|---|
| Theorem 2 (1⇉2⇉3) | 1 | 1 + (t′ − t) δ_a δ_b | none |
| Theorem 4 (1⇉2⇉3) | a_i ↦ (t−1)/(t′−1) · a_i | 1 − (t′−t)/(t′−1) · δ_a δ_b | t, t′ = 1 |
| Theorem 3 F2, F4 (1⇉2⇇3) | 1 | 1 + (t′ − t)(δ_b − δ_a) | none |
| Theorem 3 F1 | a₂ ↦ (t/t′) a₂, b₁ ↦ (t/t′) b₁, b₂ ↦ t(t+1)/(t′(t′+1)) · b₂ | 1 + (t′−t)/t′ · a₂ + (t′−t)/(t′(t′+1)) · b₂ | t, t′ ∈ {0, −1} |
| Theorem 3 F3 | F1 with classes a and b exchanged | | t, t′ ∈ {0, −1} |

- **Verification.** The identity φ P_t = P_t′ φ is checked in ℚ(t, t′), with P_t taken from the kernel-checked Lean
  tables (`analysis/family_behaviour.py`).
- **How the witnesses were found.** First, exhaustive reachability over 𝔽₅: the inner automorphisms alone (5⁸ of them
  on 1⇉2⇉3), then all of Inn ⋊ Sub (230 400 × 390 625). Then the general solution was computed symbolically, and the
  free parameters were set to zero.
- **Why the census never saw the collapse.** For Theorem 2, F2 and F4 the family is one orbit of a one-parameter
  unipotent group of inner automorphisms, so t is literally its time. For Theorem 4, F1 and F3, neither the inner
  automorphisms nor the substitutions suffice alone: the composite is needed. The torus factor has non-unit column
  sums, so it is not Δ-compatible. That is why t survived every group the earlier rounds tried (G, the class torus,
  local torus moves).
- **The degenerate parameters are the poles.** At t = 1 (Theorem 4) the image is k[x]/(x²) instead of k × k. At
  t = 0 (F1, F3) it is k ⊕ (square-zero). At t = −1 (F1, F3) the image is still k × k, but its embedding is a different
  Aut(kQ)-orbit.

## 3. The behaviour census

`analysis/behaviour_census.py` computes the behaviours of a census from its fold set (`pcensus.py --save-folds`).
- It forms the G-orbits (unions along generators of G).
- It joins them whenever Inn_u ∘ s carries a member to a fold, with s running over right-coset representatives of G in
  Sub (Sub = ⋃ G·s, the cosets labelled by the column-sum vector 1ᵀs).
- This suffices because Inn is normal in Aut(kQ): g·Inn_u·s = Inn_{g(u)}·g·s.
- Quiver automorphisms are Δ-compatible, so they are merged in afterwards.

| census | G-orbits | behaviours (Aut⁰) | behaviours (Aut) | inner automorphisms | coset reps of G in Sub | G-orbits per behaviour (count) |
|---|---|---|---|---|---|---|
| 1⇉2, 1→3 / 𝔽₃ | 105 | 91 | 91 | 27 | 16 | 1 (77), 2 (14) |
| 1⇉2, 1→3 / 𝔽₁₁ | 137 | 91 | 91 | 1331 | 1200 | 1 (77), 2 (10), 10 (4) |
| 1⇉2⇇3 / 𝔽₃ | 485 | 319 | 170 | 81 | 64 | 1 (213), 2 (90), 4 (12), 11 (4) |
| 1⇉2⇇3 / 𝔽₅ | 557 | 319 | 170 | 625 | 576 | 1 (213), 2 (74), 4 (28), 21 (4) |
| 1⇇2⇉3 / 𝔽₃ | 485 | 319 | 170 | 81 | 64 | 1 (213), 2 (90), 4 (12), 11 (4) |
| 1⇇2⇉3 / 𝔽₅ | 557 | 319 | 170 | 625 | 576 | 1 (213), 2 (74), 4 (28), 21 (4) |

Aut⁰ = Inn ⋊ Sub; Aut adds the quiver automorphisms. Checks built into the runs: the G-orbit counts reproduce the
census; module class is constant on every behaviour.

Reading the table:
- **The q-dependence lives entirely in the fibres.**
  - On 1⇉2⇇3 the four widest behaviours (11 G-orbits at q = 3, 21 at q = 5) are the four Theorem 3 families: three
    have no module condition (F1, F2, F3) and one is in L (F4), the 3 + 1 split of v2.32.
  - Inside each, the vertex block takes q distinct values: it is the coordinate along the fibre.
  - On 1⇉2, 1→3 the four behaviours of width 10 at q = 11 carry the growth 105 → 137.
- **1⇉2⇇3 and its opposite agree** behaviour for behaviour in the counts, as duality requires.
- **Not run:** 𝔽₇ for 1⇉2⇇3. The method costs |cosets| × |Inn| = 2304 × 2401 conjugations per representative, which
  puts it out of reach in this container. 1⇉2⇉3 (Inn = q⁸) is out of reach for the same reason.

**Pre-registered predictions** (`logs/v2_33/PREDICTIONS.md`):

| # | prediction | outcome |
|---|---|---|
| 1 | behaviour count independent of q | holds: 91 (1⇉2, 1→3, q = 3, 11); 319 / 170 (1⇉2⇇3, q = 3, 5) |
| 2 | the 4(q−2) growth of 1⇉2⇇3 disappears; the Theorem 3 families are 4 wide behaviours | holds: 4 behaviours of width 11 → 21 |
| 3 | module class constant on behaviours | holds |
| 4 | pairs_in and pairs_out agree | holds, at both primes |

The run over 𝔽₇ named in the pre-registration was not done (cost, §3).

## 4. Le Chat's three sub-questions

1. **Finiteness** (behaviours independent of q). Holds on every census tested: 1⇉2, 1→3 (91 behaviours at q = 3 and 11) and 1⇉2⇇3 (319 / 170 at q = 3 and 5). The growth terms of the orbit counts are fibre widths.
2. **Width as a theorem about the type.** Inside each family, the members sharing the generic behaviour are exactly
   the parameters off the poles of the connecting automorphism:
   - all t for Theorem 2, F2 and F4;
   - t ≠ 1 for Theorem 4;
   - t ∉ {0, −1} for F1 and F3.

   So the width is forced by a closed-form witness, not measured. Whether every positive-degree term of every census
   polynomial arises this way is the open half.
3. **The lift filter.**
   - *Factor through behaviour:* module class (v2.32), the abstract isomorphism classes of im P and ker P, the
     completion variety (determined by the vanishing pattern, itself Aut(kQ)-invariant).
   - *Do not factor:* the family parameters, including Theorem 2's invariant w′/(xz) (it is the Δ-position), and,
     against Le Chat's expectation, Theorem 3's vertex block, which varies inside behaviours because inner
     automorphisms move the vertex idempotents.

## 5. Next

1. Lean for Theorem 4, F1 and F3. Their witnesses are rational in (t, t′). With a fixed generic base point, the
   denominators can be cleared, or the statement can be made over a field.
2. The behaviour census for 1⇉2⇉3 (inner group 5⁸ over 𝔽₅) needs the GPU engine or a smarter orbit method (a
   canonical form under Inn).
3. Characterise the fibres as Δ-compatible positions: the set {φ ∈ Aut(kQ) : φ P φ⁻¹ is a fold} modulo G.

Files: `analysis/behaviour_census.py`, `analysis/family_behaviour.py`, `lean/gen/gen_orbit.py`,
`lean/EatTheFibre/OrbitReflect.lean`, `lean/EatTheFibre/Gen/orbits_*.lean`, `logs/v2_33/`.
