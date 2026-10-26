# Discordianism: Open Questions Explored (v2.58): a uniform bound in q for 1⇉2⇇3 + 1→3, in every characteristic

*Program 2 of ChatGPT's Phase 2 brief: turn the characteristic-0 finiteness theorem of Q₀ = 1⇉2⇇3 + 1→3 (v2.54)
into a bound on the 𝔽_q census count that holds for every prime power q. Exceptional characteristics are not to be
waved away. The result: Q₀'s certificates specialise to every prime. The integral structure is unimodular, every
presolve pivot is ±1, and the only candidate bad prime, 2, passes a direct check. Lang's theorem and Bézout counts
then give an explicit bound, the same for all q. So Conjecture 3ᶜ holds for Q₀, as a theorem modulo the trusted
Singular computations listed in §4. For Q₁ the integral structure is checked too; its certificates need one more
rerun.*

---

## 0. Results

The predictions are in logs/v2_58/PREDICTIONS.md, committed before the runs.

| # | Statement | Status |
|---|---|---|
| 61 | Every specialisation claim of Q₀'s 21 strata holds at every prime; no exceptional prime | **held**: 1441 claims. Candidate primes (from leading coefficients and Rabinowitsch integers): only 2, in two strata, and the direct check over 𝔽₂ passes. All 21 strata were also checked directly over 𝔽₂ and 𝔽₃ |
| 62 | The 21 vertex blocks are the only ones in every characteristic | **held**: elimination over ℤ (every leading coefficient 1) and a radical certificate. Cross-checked by direct elimination over 𝔽₂, 𝔽₃, 𝔽₅: 21 points each |

**Theorem (Conjecture 3ᶜ for Q₀).** There is a constant B such that for every prime power q, the census count of
Q₀ over 𝔽_q is at most B. The census count is the number of classes of folds of 𝔽_qQ₀ under the census group and the
local class-weight moves. The proof gives B ≈ 10^51.5, which is absurd next to the observed 366 at q = 3, 5, 7, 11,
but it is uniform.

## 1. The integral structure (`analysis/v2_58_integral.py`, `logs/v2_58/integral_q0.txt`)

All of the following are exact computations over ℤ.

| check | Q₀ | Q₁ |
|---|---|---|
| Smith normal form of Δ (as a d² × d integer matrix) | rank d − kdim, elementary divisors all 1 | same |
| Kmat (the K-basis used by the strata) | in ker Δ over ℤ; gcd of maximal minors 1 (saturated) | same |
| det [Kmat \| e_C] | −1 | −1 |
| census group G (column-sum-1 substitutions per parallel class), generic parameters | Φ multiplicative; Δ ∘ Φ = (Φ⊗Φ) ∘ Δ and Δ_c ∘ Φ = (Φ⊗Φ) ∘ Δ_c for every class, as identities over ℤ[g] | same |
| Der_Δ basis of the pipeline | integral combinations of the tangents of G at 1 | same |
| class weightings φ_w | algebra automorphisms over ℤ[w^±1]; Δ ∘ φ_w = (φ_w⊗φ_w) ∘ Δ_w; Φ φ_w = φ_w Φ | same |
| presolve pivots of every stratum (now recorded in the certificates) | all ±1 (21 strata) | (to rerun) |

Consequences, for every prime p:
- **K_p is the reduction of K.** The kernel of Δ over 𝔽_p is spanned by Kmat mod p.
- **The stratum coordinates are valid over 𝔽_p.** The coordinates (P_K, P(C)) are an isomorphism over 𝔽_p.
- **The presolve is valid over 𝔽_p.** It eliminates variables by x = −rest/c with c = ±1, so the presolved equations
  `rem` are the stratum's equations mod p. X_{V,p} ≅ V(rem mod p) by the same substitutions.
- **The class-weight lemma holds in every characteristic.** If P is compatible with Δ_w, then P′ = φ_w P φ_w⁻¹ is a
  fold, since ΔP′ = (φ⊗φ)Δ_w P φ⁻¹ = (φ⊗φ)(P⊗P)Δ_w φ⁻¹ = (P′⊗P′)Δ.
- **G is smooth and connected.** G is the open subscheme det ≠ 0 of an affine space over ℤ. It preserves each Δ_c, so
  the compatible weightings W(P) are constant on G-orbits, and it commutes with the weightings. So H̃_W = G × T_W is a
  connected smooth group over ℤ, where T_W is the subtorus of weightings constant on the blocks of the partition W.
  H̃_W acts on X^W, the folds compatible with W, and preserves vertex blocks. Its Lie algebra mod p contains the
  pipeline's integral basis of 𝔥 = Der_Δ ⊕ E_W.

## 2. The vertex blocks in every characteristic (`analysis/v2_58_vertex_integral.py`)

v2.52's system has 389 equations in 55 unknowns: P_K idempotent, products closed on K, and the Δ-conditions of the
arrows a, b, c. It is rewritten over ℤ, with the integral rows of B⁻¹ in place of v2.52's rational left inverse, and
eliminated to the 9 vertex-block entries in Singular over `integer`.
- **The elimination ideal.** It has a strong Gröbner basis E_ℤ of 44 elements, **all with leading coefficient 1**.
  Over ℚ its vdim is 29: 21 points, with multiplicity.
- **The radical certificate.** The ideal R of the 21 points over ℚ has 35 integral generators, all quadratic with
  coefficients ±1, among them v² − v for each entry. Each generator ρ satisfies ρ^k ∈ E_ℤ with k ≤ 2. This is checked
  by exact reduction modulo E_ℤ, which over ℤ divides by leading coefficient 1 only.
- **R cuts out exactly the 21 points over every field.** The generators v² − v force 0/1 entries. Of the 512
  0/1-matrices, exactly the 21 are common zeros. Every other one has a generator taking the value ±1 there, so it
  stays excluded mod every p.

So over every field, a fold's vertex block is one of the 21. As a cross-check, direct elimination over 𝔽₂, 𝔽₃ and 𝔽₅
gives dimension 0, vdim 29 and 21 points each time (`vertex_char{2,3,5}_q0.txt`).

A fully checkable cofactor certificate, ρ^k = Σ a_i f_i over ℤ via Singular `lift` in 55 variables, was attempted
and stopped after 50 min without finishing. This step therefore rests on Singular's elimination over ℤ.

## 3. The char 0 certificates specialise (`analysis/v2_58_goodprimes.py`, `logs/v2_58/goodprimes_q0.txt`)

The Q₀ pipeline was rerun with certificates and presolve pivots (`certs_q0.jsonl`, same 118 components as v2.54).
For each stratum, with gens_Z the recorded generators of a component's prime and Z_p := V(gens_Z mod p), three
claims are needed at p:

| claim | statement mod p | how the candidate bad primes are found |
|---|---|---|
| DIM | dim V(gens_Z + F_r) ≤ r for each cut; dim V(gens_Z) ≤ dim Z | strong Gröbner basis over ℤ: for p dividing no leading coefficient, the leading ideal mod p contains the one over ℚ |
| DEF | the weightings W_Z are compatible on Z_p | the defects lie in (gens_Z) over ℚ; dividing by the ℤ-basis needs only its leading coefficients |
| COV | V(rem mod p) ⊆ ∪_Z Z_p | chain K_i = K_{i−1} ∩ (gens_i) with checked products K_{i−1}·(gens_i) ⊆ K_i; then each generator of K_n lies in √rem, by Rabinowitsch over ℤ: the ideal rem + (1 − t·k) contains an integer n_k |

The minors are used exactly as computed, with denominators cleared and no content divided out, so that mod p they are
still minors. Every candidate prime is then checked directly in characteristic p, with dimensions, membership or
Rabinowitsch, and coverage by intersection and Rabinowitsch. Independently of the candidates, every stratum is checked
directly over 𝔽₂ and 𝔽₃.

**Result.**
- 1441 claims in 21 strata.
- Every leading coefficient and every n_k is ±1, except in strata 100,000,000 and 100,010,001, where 2 occurs.
- Every direct check passes, at p = 2 everywhere and at p = 3 everywhere. **No exceptional prime.**

**The checker has been tested.** A certificate with every minor multiplied by 5 produces the candidate 5, and the
direct check fails there ("BAD DIM" at 9 cuts), as it should.

## 4. The proof of the theorem

Fix q = pⁿ. All geometry is over 𝔽̄_p.

1. **Vertex blocks (§2).** X_p = ⊔_V X_{V,p} over the 21 blocks, and X_{V,p} ≅ V(rem mod p) (§1).
2. **Components (COV).** Each irreducible component Y of X_{V,p} lies in Z_p for some component Z of the char 0
   stratum. On Z_p, W_Z is compatible (DEF), so Y ⊆ X^{W_Z}. As a component of X_{V,p}, Y is also a component of
   the closed subset X^{W_Z}_{V,p}. The connected group H̃_Z = G × T_{W_Z} acts there, so it preserves Y.
3. **Rank bounds orbit dimension in any characteristic.** For P ∈ Y, dim H̃_Z·P ≥ rank [𝔥_Z, P], because the
   differential of h ↦ hPh⁻¹ maps Lie H̃_Z ⊇ 𝔥_Z into T_P(H̃_Z·P) by ξ ↦ [ξ, P]. So an orbit of dimension r < dim Z
   lies in V(J_{Z,r}) = Z_p ∩ {all (r+1)-minors = 0}, which has dimension ≤ r (DIM).
   - Its closure is therefore an irreducible component of V(J_{Z,r}).
   - Two different orbits cannot share a closure, because orbits are open in their closures.
   - Orbits of dimension ≥ dim Z can only be dense in Y, since dim Y ≤ dim Z_p ≤ dim Z (DIM, top).

   So H̃_Z has at most o_Z = 1 + Σ_{r<dim Z} #comp V(J_{Z,r},p) orbits on Y.
4. **From geometric orbits to 𝔽_q-orbits (Lang).** Let P ∈ X(𝔽_q), in a component Y. Its orbit O = H̃_Z·P is
   Frobenius-stable, because H̃_Z is defined over 𝔽_p. H̃_Z is connected, so O(𝔽_q) is a union of at most
   |π₀(Stab P)| orbits of H̃_Z(𝔽_q) = G(𝔽_q) × T_{W_Z}(𝔽_q).
   - Each element of that group is a census move: G(𝔽_q) is the census group, and t ∈ T_{W_Z}(𝔽_q) is compatible at
     every point of Y (DEF; Y is H̃_Z-stable).
   - So every census class contains a whole H̃_Z(𝔽_q)-orbit, for some Y and some orbit O ⊆ Y.

   Hence c(q) ≤ Σ_V #comp(X_{V,p}) · max_Z o_Z · max_P |π₀(Stab P)|.
5. **Bézout counts, uniform in p and q (Heintz).** A variety in Aⁿ cut out by equations of degree ≤ D has at most
   Dⁿ irreducible components.
   - #comp(X_{V,p}) ≤ D_V^{n_V}, with n_V ≤ 30 presolved coordinates.
   - #comp V(J_{Z,r},p) ≤ D_{Z,r}^{n_V}.
   - Stab P is cut out of H̃ ⊂ A¹² by equations of degree ≤ 4. These are Φ φ_w P = P Φ φ_w, where paths have length
     ≤ 2, together with the group equations. So |π₀| ≤ 4¹² ≈ 1.7·10⁷.

   From the certificates (`analysis/v2_58_bound.py`, `logs/v2_58/bound_q0.txt`), Σ_V D_V^{n_V} · max_Z o_Z ≤ 10^44.3.
   So c(q) ≤ B ≈ 10^51.5 for every q. ∎

**What the proof trusts.** Exact computations, but run inside Singular and not independently re-certified:
- `minAssGTZ` over ℚ (the components; v2.57's coverage check confirmed them on Q₁, though not on one stratum);
- `std` over ℤ (the leading coefficients);
- `eliminate` over ℤ (§2);
- `intersect` over ℚ.

Each of these is cross-checked by direct computation at p = 2, 3 (and 5 for §2). A step that would remove one of
them: cofactor certificates (`lift`) for each Gröbner basis element, which are checkable by expansion. The one
attempted, in 55 variables, did not finish.

**Exceptional primes.** None. The brief asked for every prime outside the good-reduction set to be treated
separately. The candidate set is {2}, and the claims hold there by direct computation. Characteristic 2 is not
special for Q₀.

## 5. Q₁ (next)

The integral structure of Q₁ is the same: Smith normal form 1, a saturated K-basis, det B = −1, and the group
identities over ℤ (`integral_q1.txt`). Two things remain:
1. **Rerun the certified pipeline with pivots recorded.** The v2.57 certificates predate the pivot field.
2. **Treat block 0000,0000,0000,0001 in characteristic p.** It was decomposed by case splits. Its rewrites are
   char-free when the factorisation constants are p-units: splitting f₁f₂ = 0, dividing by a known-nonzero factor,
   and the radical of c·∏f_i^{e_i}. Recording those constants gives that block's candidate primes.

Then `v2_58_goodprimes.py` and `v2_58_bound.py` apply unchanged. Program 3 (the leaf extension) comes after.

## Files

- analysis/v2_58_integral.py, v2_58_vertex_integral.py, v2_58_goodprimes.py, v2_58_bound.py (new). In
  fold_ideal.py the presolve records pivots (`PIVOTS`), and in v2_55_pipeline.py the stratum certificates now carry
  them.
- logs/v2_58/:
  - PREDICTIONS.md;
  - integral_q0.txt, integral_q1.txt;
  - pipeline_q0_certs.txt and certs_q0.jsonl;
  - goodprimes_q0.txt;
  - vertex_integral_q0.txt, vertex_char{2,3,5}_q0.txt, vertex_lift_q0.txt (stopped);
  - bound_q0.txt.
