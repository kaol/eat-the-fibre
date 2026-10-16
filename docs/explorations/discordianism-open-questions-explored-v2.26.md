# Discordianism: Open Questions Explored (v2.26)

*Takes up question 1 of v2.25: turn the count evidence for 1⇉2⇉3 into an explicit family over ℚ(t). Done. The
refutation of Conjecture 3\*ˡ is now a theorem with a two-line proof, and the families of all three growing orbit types
are located exactly, through fixed-point loci of stabilizers.*

*Scripts (all in `analysis/`): `fixedlocus.py` (exact stage and numerical stage), `family.py` (Singular decomposition of
a fixed-point locus), `family_invariant.py` (parametrization, symbolic fold check, invariant), `family_explicit.py` (the
theorem's checks). Logs: `logs/fixedlocus_two_pairs_3*.log`, `logs/family_*.log`. Labels: **proved**,
**verified (exact)**, **numerical**, **conjectured**. Discordian glosses are descriptive labels (L3): nothing in the
mathematics forces them.*

---

## 0. Summary

| # | Question | Answer |
|---|---|---|
| v2.25 Q1 | An explicit family for 1⇉2⇉3 | **Theorem 1:** P_t (below) is a Δ-compatible fold for every t, and t is invariant under H = G × T² (Aut_Δ° times the class torus). So Folds(k(1⇉2⇉3)) has infinitely many H-orbits, over any field with infinitely many elements **[proved; checked symbolically and over 𝔽₃, 𝔽₅]** |
| — | Where the families are | Fixed-point loci locate families of non-generic orbits exactly. Over 𝔽₃, **82 H-orbits** lie on families, all in the three orbit types whose counts grow with q (sizes (q − 1)², q(q − 1)², q²(q − 1)²). No family occurs in any type whose count is constant **[exact candidates; numerical verdicts; one exact Gröbner confirmation each way]** |
| — | Why two pairs | The family uses δ_αδ_β, the product of the two parallel differences. It lies in ker Δ, and a quiver with one parallel class has no such element. This supports Conjecture 3ᵖ (v2.25) and suggests where a proof should look |

---

## 1. Theorem 1

Arrows α₁, α₂: 1 → 2 and β₁, β₂: 2 → 3. Write δ_α = α₁ − α₂, δ_β = β₁ − β₂ and
δ_αδ_β = α₁β₁ − α₁β₂ − α₂β₁ + α₂β₂. For a, b, c ∈ k, define P = P_{a,b,c} by

  P(e₁) = e₁ − a·δ_α − (ab + c)·δ_αδ_β,  P(e₃) = e₃ − b·δ_β + c·δ_αδ_β,

and P = 0 on every other basis path (e₂, the arrows, the paths of length 2). Put P_t = P_{1,1,t}.

**Theorem 1.** (i) Every P_{a,b,c} is a Δ-compatible fold. (ii) For h ∈ H = G × T², hP_{a,b,c}h⁻¹ = P_{λa, μb, λμc}
for scalars λ, μ depending on h. Hence (iii) t = c/(ab) is H-invariant on {ab ≠ 0}, and the P_t (t ∈ k) lie in pairwise
distinct H-orbits.

*Proof.* (i) Symbolic verification of idempotency, of both subalgebra conditions and of ΔP = (P⊗P)Δ over ℚ(a, b, c)
(`family_explicit.py`, check 1). (ii) G substitutes arrows within each class with column sums 1, so g(δ_α) = λ_g δ_α
and g(δ_β) = μ_g δ_β with g fixing vertices. The class torus scales δ_α, δ_β by its weights. So every h ∈ H fixes e₁, e₃
and scales δ_α, δ_β and δ_αδ_β by λ, μ and λμ. Also h preserves path length, so h⁻¹ maps paths of positive length to
combinations of such paths, which P kills. Then hPh⁻¹ has the stated form. (iii) follows. ∎

Checks over 𝔽₃ and 𝔽₅ (`logs/family_explicit.log`): every element of H preserves t on the family, every P_t is a census
fold, the p members lie in p distinct H-orbits, and each orbit has size (q − 1)². That is exactly the orbit type whose
count grows (18 + (q − 1)).

Generic members have r = 2 (Lie H acts with rank 2, so the stabilizer has dimension 4). The family is a 3-dimensional
rational component of the fixed-point locus of that stabilizer (§2), swept by 2-dimensional orbits.

**Gloss.** *When twin steps follow twin steps, the difference of the first twins times the difference of the second is
invisible to deletion. A split can lean on it by an amount that no relabelling and no clock can reach: the ratio of
how far it leans on the product to how far it leans on each difference.*

---

## 2. Locating families: fixed-point loci

Families in lower strata cannot be found by walking to generic points (v2.23–v2.25). They live inside the fixed-point
locus of their stabilizer. For an H-orbit representative P with stabilizer S (positive dimension), the orbits with
stabilizer ⊇ S form Folds^S = Folds ∩ C_S, where C_S is the commutant of Lie S. The orbit of P in it is N·P,
N = N_H(S). So there is a family through P iff dim_P Folds^S > dim N·P.

`fixedlocus.py`, exact stage (over ℚ, DomainMatrix): for each of the 689 lifted H-orbits of 1⇉2⇉3 over 𝔽₃ (27 do not
lift), it computes Lie S, C_S, Lie N, n = dim N·P and t_S, the tangent dimension of Folds^S at P. **t_S = n proves that no
family passes through P.** That holds for 553 orbits. 12 more are free (no stabilizer, open orbits), and 124 are candidates:

| orbit size (r) | candidates | on a family (numerical local dim > n) | count over 𝔽_q |
|---|---|---|---|
| 4 (r = 2) | 4 | **4** | 18 + (q − 1) |
| 6 (r = 2) | 10 | 0 | constant |
| 12 (r = 3) | 26 | **18** | 46 + 4(q − 1) |
| 18 (r = 3) | 20 | 0 | constant |
| 36 (r = 4) | 64 | **60** | 151 + 14(q − 1) |

Families occur exactly in the three types whose counts grow. Two exact Gröbner computations (Singular over ℚ, `family.py`)
anchor the numerical verdicts:
- a size-6 candidate: the component of Folds^S through P has dimension 1 = n, so no family;
- a size-4 candidate: dimension 3 > n = 2, the family of Theorem 1.

Caveats:
- The numerical local dimension measures the span of the tangent cone, so where two family components cross it can
  read high (some size-36 orbits read 4). The family verdict is unaffected: any component through P other than N·P has
  dimension > n.
- One reading (local dimension 0 < n = 1) is a numerical failure.
- The 27 orbits with no member lifting to ℚ are untested.
- The member counts over 𝔽₃ (4, 18, 60) exceed what 𝔾_m-families would give (2, 8, 28). So some families are complete
  or have extra special members (Theorem 1's family has t ∈ 𝔸¹ plus degenerations with ab = 0).

---

## 3. Proposed amendments (v2.25 → v2.26)

| amendment | content |
|---|---|
| **Q1²⁶ Theorem 1** | Explicit family P_t on 1⇉2⇉3 with H-invariant t: Conjecture 3\*ˡ is false (proved) |
| **Q2²⁶ Families via fixed loci** | Families of non-generic orbits are detected by dim Folds^S > dim N_H(S)·P; t_S = n excludes them exactly |
| **Q3²⁶ Mechanism** | The moduli come from the product of parallel differences, which lies in ker Δ. Conjecture 3ᵖ (one parallel class) stands |

---

## 4. Open questions

1. **Conjecture 3ᵖ**: with one parallel class, ker Δ = A₀ ⊕ D has no products of differences. Is that enough for a proof
   through the kernel fibration?
2. Classify the 19 families of 1⇉2⇉3 (one is Theorem 1) and their parameter curves.
3. 1⇉2→3⇉4: are the two pairs, separated by an arrow, still coupled through ker Δ (δ_α β δ_γ)? Its 𝔽₅ census is
   running.
4. Carried: unobstructed violations; Kronecker plus a loop; zigzag; dihedral; ρ.

---

## 5. Sources

- M. Aguiar, *On the associative analog of Lie bialgebras*, J. Algebra 244 (2001).
- W. Decker, G.-M. Greuel, G. Pfister and H. Schönemann, *Singular* (via passagemath-singular): minimal associated primes.

*"When twin steps follow twin steps, the difference of the first twins times the difference of the second is invisible
to deletion. A split can lean on it, and the ratio of that lean to the leans on each difference is a modulus no
relabelling and no clock can reach."*
