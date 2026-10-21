# Discordianism: Open Questions Explored (v2.36): all nineteen families of 1⇉2⇉3 in closed form

*v2.35 located the 19 families behind 1⇉2⇉3's local-torus count 678 + 19(q−1) by behaviour type. This round gives every
one of them a closed form, P_t = Inn(1 + t·x)·P₀·Inn(1 + t·x)⁻¹, with P₀ an integer fold and x ∈ {δ_a, δ_b, δ_aδ_b}. It
checks the forms at two primes, and kernel-checks them in Lean. Theorems 2 and 4 are two of the nineteen; seventeen
are new. This answers v2.28 open question 2 for 1⇉2⇉3.*

---

## 0. Summary

| # | Statement | Status |
|---|---|---|
| 1 | Each of the 19 growing families of 1⇉2⇉3 is the orbit of a fold P₀ under a one-parameter group of inner automorphisms t ↦ Inn(1 + t·x), x ∈ {δ_a, δ_b, δ_aδ_b} | found by `analysis/family_normal_forms.py` (𝔽₇ census); fold for all t verified over ℚ[t] |
| 2 | Each family is a fold over every commutative ring for every t, and is one inner orbit | **proved (Lean)**, `fam{k}_isFold`, `fam{k}_one_orbit` (`Gen/two_pairs_fam*.lean`) |
| 3 | The 19 closed forms meet q pairwise distinct local-torus orbits each, disjoint across families, at q = 7 and q = 5 (where they were not fitted). The remaining 659 orbits are the same number at both primes: the local-torus count is **659 + 19q** | computed (`analysis/v2_36_crossprime.py`) |
| 5 | The same inner-curve search finds Theorem 3's four families on 1⇉2⇇3 (curves along δ_a + δ_b meeting q orbits at q = 5, 7), and **no** such curve on the joined-ends quiver 1⇉2⇇3 + 1→3 (q = 3, 5) | computed; pre-registered thresholds failed, corrected criterion post hoc (§4) |
| 4 | Theorem 4 is family 3 (δ_aδ_b, image k×k). Theorem 2 is family 13 (class LRKE). Each shares 6 of its 7 local orbits over 𝔽₇ with the closed form; Theorem 4's seventh is its degenerate t = 1 | computed (`analysis/v2_36_named.py`) |

## 1. The families

The table below lists, for each family:
- the signature group of v2.35;
- the direction x of its inner curve;
- the census representative P₀, a G-orbit representative over 𝔽₇, lifted to ℤ by symmetric residues.

| family | class | im P | ker P (dim, powers, ann_L, ann_R, centre) | x | named |
|---|---|---|---|---|---|
| 1, 4, 7 | none | k×k | (9, 9, 0, 0, 1) | δ_a | |
| 2, 5, 8 | none | k×k | (9, 9, 0, 0, 1) | δ_b | |
| 3, 6, 9 | none | k×k | (9, 9, 0, 0, 1) | δ_aδ_b | family 3 = **Theorem 4** |
| 10 | none | k³ | (8, 8, 0, 0, 0) | δ_a | |
| 11 | none | k³ | (8, 8, 0, 0, 0) | δ_b | |
| 12 | none | k³ | (8, 8, 0, 0, 0) | δ_aδ_b | |
| 13 | LRKE | k×k | (9, 9, 6, 6, 4) | δ_a | **Theorem 2** |
| 14 | L | k×k | (9, 9, 0, 0, 0) | δ_a | |
| 15 | L | k×k | (9, 9, 0, 5, 0) | δ_b | |
| 16 | L | k³ | (8, 8, 0, 4, 0) | δ_b | |
| 17 | R | k×k | (9, 9, 0, 0, 0) | δ_b | |
| 18 | R | k×k | (9, 9, 5, 0, 0) | δ_a | |
| 19 | R | k³ | (8, 8, 4, 0, 0) | δ_a | |

The groups have a visible structure:
- **Unclassed families come in triples.** Each triple is the same behaviour type along δ_a, δ_b and δ_aδ_b. The k×k
  group holds three such triples, the k³ group one.
- **L and R mirror each other** under opposite duality, which swaps δ_a with δ_b and the annihilators: 14 ↔ 17,
  15 ↔ 18, 16 ↔ 19.
- **The δ-direction of a curve is not a family invariant.** Theorem 2 appears here along δ_a; in v2.33 it was the curve
  along δ_aδ_b from the base point P_0 of the Lean table. A family can contain several inner curves through different
  base points.

Three of the closed forms (from `analysis/v2_36_show.py`; nonzero columns only; δ_a = a1 − a2, δ_b = b1 − b2,
δ_aδ_b = a1b1 − a1b2 − a2b1 + a2b2):

```
family 12 (no class, im k³, along δ_aδ_b):
   P(e2)   = e1 + e2 + e3
   P(a2)   = e1 − δ_a − t·δ_aδ_b
   P(b2)   = −e3 − δ_b − (t + 1)·δ_aδ_b
family 13 (LRKE = Theorem 2, along δ_a):
   P(e1)   = e1 − t·δ_a + δ_aδ_b
   P(e3)   = e3 − δ_b − (t + 1)·δ_aδ_b
family 16 (L, im k³, along δ_b):
   P(e1)   = e1 + e2 + e3
   P(a2)   = −e2 − (t + 1)·e3 − δ_a − t(t + 1)·δ_b
   P(a2b2) = e3 + (t + 1)·δ_b + δ_aδ_b
```

The full matrices are in `lean/EatTheFibre/Gen/two_pairs_fam*.lean`.

## 2. How they were found and checked

1. **Search** (`analysis/family_normal_forms.py`). For every local-torus orbit in a growing signature group, take a
   member whose symmetric residues form a fold over ℤ. 777 of the 792 local orbits over 𝔽₇ have such a member.
   - For x ∈ {δ_a, δ_b, δ_aδ_b}, form P_t = Inn(1 + t·x)·P₀·Inn(1 + t·x)⁻¹. Image and kernel closure are automatic by
     transport (v2.34). Δ-compatibility is checked as a polynomial identity in t.
   - Locate the 𝔽₇-points of the curve in the census. On 1⇉2⇉3 every fold admits every class weighting (v2.24), so a
     local-torus orbit is a (G × T²)-orbit; it is identified by a minimum hash over the 63 504 elements of G × T².
   - A curve that meets at least 3 local orbits is a family. Curves meeting 2 are degenerations: a generic orbit plus a
     special point.
2. **The search finds exactly the expected families.** It finds 9, 3, 1, 1, 1, 1, 1, 1, 1 families in the nine groups:
   the slopes of v2.35, with no extra and no missing family.
3. **Cross-prime check** (`analysis/v2_36_crossprime.py`). The same integer closed forms, reduced mod 5, were located
   in the 𝔽₅ census. Each meets exactly 5 local orbits, in its own module class, and the 19 are pairwise disjoint.
   With 754 local orbits over 𝔽₅ and 792 over 𝔽₇:
   - 754 − 19·5 = 659;
   - 792 − 19·7 = 659.

   The families account for the entire growth, and the 659 others are constant.
4. **Lean** (`lean/gen/gen_families.py`, `Gen/two_pairs_fam{1..19}.lean`). Each family is checked by the reflection
   of v2.29. `fam{k}_check` runs one `decide +kernel` per row, and `check_sound` turns it into `fam{k}_isFold` over
   every commutative ring, for every t. `fam{k}_one_orbit` (v2.33's `orbit_sound`) shows that
   Φ_t = Inn(1 + t·x) is a multiplicative bijection with Φ_t P₀ = P_t Φ_t.

## 3. What is and is not proved

- **Proved (Lean):** existence. Each closed form is a fold for every t over every commutative ring, and is a single
  inner orbit.
- **Computed, not proved:** distinctness. Different t give different local-torus orbits; this is checked over 𝔽₅ and
  𝔽₇, all t. The census-exhaustive statement that no other family exists, and that the remaining 659 orbits do not
  move with q, holds at these two primes.
- **The rigidity theorems of v2.29–v2.30** prove distinctness for Theorems 3 and 4 by exhibiting H-invariants. The
  same method should work for each family, using a scaling invariant of the δ-coefficients, as for Theorem 4.

## 4. What this says about Conjecture 3ᶜ

Every modulus of 1⇉2⇉3 is now explicit, and all are of one kind: an inner automorphism along δ_a, δ_b or δ_aδ_b,
applied to a fold. The single-class elements δ_a and δ_b already occur in quivers with one parallel class, which have no
families (Conjecture 3ᵖ). So the obstruction cannot be "the δ's exist". It must be that Δ-compatibility survives
conjugation along δ for some fold.

**The direct test** (`analysis/family_curves.py`; directions α δ_a + β δ_b (+ γ δ_aδ_b), |α|, |β|, |γ| ≤ 2):

| quiver | q | inner fold-curves | curves meeting q orbits | curves meeting exactly 3 orbits |
|---|---|---|---|---|
| 1⇉2⇇3 | 5 | 472 | 4 (F1–F4; classes −, −, −, L) | 4 (LRE, LRKE) |
| 1⇉2⇇3 | 7 | 472 | 4 (the same) | 4 (the same) |
| 1⇉2⇇3 + 1→3 | 3 | 260 | — (3 = q is ambiguous) | 3 (LRE, LRKE) |
| 1⇉2⇇3 + 1→3 | 5 | 260 | **0** | 3 (the same) |

- The joining arrow removes nearly half of the inner fold-curves, and every family among them.
- What survives meets a bounded number of orbits: 3 at every q, a generic orbit plus two special points.
- F1 and F3, which v2.33 could connect only with a torus composite, are pure inner curves from a suitable base point.
- The pre-registered criterion (≥ 3 orbits) was too crude, and both its predictions failed as written
  (`logs/v2_36/PREDICTIONS.md`). The criterion that separates, "the number of orbits met grows with q", is post hoc.

That points to a sharper form of lemma (L) (v2.34): **on an unobstructed quiver, for every fold P and every δ in the
span of parallel-class differences and their products, an inner curve Inn(1 + t·δ)·P of folds meets only boundedly
many census orbits.** On 1⇉2⇉3 the 19 families are exactly the failures of that statement. On 1⇉2⇇3 + 1→3 it holds
for every curve tested, at two primes.

## 5. Next

1. Rigidity (distinctness) proofs in Lean for the 17 new families, following `thm4_rigid_GT`.
2. A proof of the sharpened (L) for joined converging ends. The data say the arrow u → w destroys exactly the curves
   that would carry a modulus. The three surviving bounded curves (classes LRE, LRKE) are the cases a proof must
   handle.
3. The same normal-form search on 1⇉2⇉3 + 1→3, where Theorem 2 fails but Theorem 4 survives.

Files: `analysis/family_normal_forms.py`, `analysis/family_curves.py`, `analysis/v2_36_*.py`, `lean/gen/gen_families.py`,
`lean/EatTheFibre/Gen/two_pairs_fam*.lean`, `logs/v2_36/`.
