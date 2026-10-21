# Discordianism: Open Questions Explored (v2.37): the nineteen families of 1⇉2⇉3 are rigid

*v2.36 gave every family of 1⇉2⇉3 a closed form and proved in Lean that each is a fold and one inner orbit.
Distinctness, meaning that different t give different census orbits, was computed at q = 5 and 7 but not proved.
This round proves it for all 19 families, over every commutative ring, under the full census symmetry group G × T².
Each proof is a linear_combination certificate with integer coefficients. The audit now has 136 results, all on the
standard axioms (one on none).*

---

## 0. Summary

| # | Statement | Status |
|---|---|---|
| 1 | For every family k and every h acting on V = ⟨e₁, e₂, e₃, δ_a, δ_b, δ_aδ_b⟩ by fixing the e_i and scaling δ_a, δ_b, δ_aδ_b by l, m, lm: if h·P_t = P_t'·h then t = t' | **proved (Lean)**, `fam{k}_rigid` (`Gen/two_pairs_rigid.lean`) |
| 2 | Hence, for every class substitution with column sums 1 and any class weights, h·P_t = P_t'·h forces t = t'. So distinct parameters lie in distinct (G × T²)-orbits, which are the local-torus orbits of 1⇉2⇉3 | **proved (Lean)**, `fam{k}_rigid_GT` |
| 3 | Every certificate has integer coefficients and degree ≤ 2. It has degree 1 exactly when a2b2 is an image column of P_t (families 1, 2, 4, 5, 7, 8, 10, 11, 16, 19); otherwise it has degree 2 | computed (`analysis/v2_37_certificates.py`), then checked by Lean |

## 1. The argument

Every image column of each P_t lies in V, by `fam{k}_col{j}`, which Lean proves for each nonzero column j. Under h
the e-part of V is fixed and the δ-part is scaled (`ScalesV`, v2.30). The G × T² case comes from `subst_scalesV`
and `colsum_eigen`, together with the new `subst_ve3`: a class substitution fixes every vertex idempotent.

Apply h·P_t = P_t'·h to e_j for each nonzero column j, and write w_j = h·e_j. This gives

  S_{l,m}(P_t e_j) = P_t' w_j = Σ_{i ∈ nz} (w_j)_i · P_t' e_i,

because P_t' vanishes off the nonzero columns. The coordinates of this identity are polynomial equations in
(w, l, m, t, t′), bilinear in w and the scalings. A certificate is a set of cofactors c_q with Σ c_q·E_q = t′ − t:
- the search tries cofactor degree D = 1, then D = 2;
- it solves the coefficient system modulo a 31-bit prime and reconstructs rationals;
- sympy checks each certificate exactly;
- all denominators come out as 1, so `linear_combination` closes the goal over any commutative ring.

Family 16 is the smallest example. Its column a2 is P(a2) = −e₂ − (t + 1)·e₃ − δ_a − t(t + 1)·δ_b. With w = h·a2, four
coordinates of the identity give:

```
e1 :  −w_e1 = 0
e2 :  −w_e1 + w_a2 − 1 = 0
e3 :  −t + t′·w_a2 − w_e1 − w_a2b2 + w_a2 − 1 = 0
a1b1: −w_a2b2 = 0
```

Here w_a2 = 1, and the e₃-coefficient then reads t′ − t = 0. The parameter is the e₃-component of a column whose
δ-free part h cannot move.

The degree-2 certificates (the δ_aδ_b triples 3, 6, 9 and 12, the k×k families of classes L, R and LRKE) need one
product of a scaling with an unknown, such as m·E or t′·E. Without an a2b2 column, nothing pins the δ_aδ_b-coefficient
linearly. Theorem 4 is family 3, one of these; its hand proof (`thm4_core`, v2.30) is likewise a chain of
`linear_combination` steps with scaling-times-unknown terms.

## 2. Status of the family picture for 1⇉2⇉3

| claim | status |
|---|---|
| each family is a fold for every t, over every commutative ring | Lean (v2.36) |
| each family is one inner orbit, so one behaviour | Lean (v2.36) |
| **different t give different local-torus orbits** | **Lean (this round)** |
| the 19 families are pairwise disjoint | computed, q = 5, 7 |
| no further families; the other 659 local orbits do not move with q | computed, q = 5, 7 |

So every family contributes at least |k| local-torus orbits over any field k. Over 𝔽_q this gives the lower bound
19q on the growing part of the census. The equality 659 + 19q still rests on the census at two primes.

Disjointness across families could be proved the same way, with h·P^{(k)}_t = P^{(k′)}_t′·h and a certificate for
1 = 0. For families in different signature groups it follows from the Aut-invariants of v2.35, which this round has
not formalized.

## 3. Next

1. Disjointness certificates (1 ∈ ideal) for pairs of families in the same signature group: the three triples of k×k
   unclassed families, and the members of each triple.
2. Proof of the sharpened (L) for joined converging ends (v2.36 §4).
3. Normal-form search on 1⇉2⇉3 + 1→3.

Files: `analysis/v2_37_certificates.py`, `lean/gen/gen_rigid.py`, `lean/EatTheFibre/Gen/two_pairs_rigid.lean`,
`logs/v2_37/` (certificates and their log; `PREDICTIONS.md` notes this round was not pre-registered).
