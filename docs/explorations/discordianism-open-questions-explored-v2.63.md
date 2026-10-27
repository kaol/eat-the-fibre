# Discordianism: Open Questions Explored (v2.63): the λ = 0 lifts by hand, and the mixed kind is rank one

*Three lemmas for a sink leaf ℓ: u → v (kQ⁺ = kQ ⊕ J), and a classification in characteristic 0 of every fold of
Q₁ = Q₀ + 2→4.
- **Lemma A (proved, char-free).** In the J-preserving kind with λ = 0, the lift is determined by a character of the
  algebra im P̄. So there are finitely many over each fold P̄ of Q.
- **Lemma B (proved, char-free).** In the kQ-preserving kind with λ = 0, the lift is determined by data in im P̄,
  including the idempotent −P(e_v) of im P̄ ∩ K.
- **Lemma C (proved, char-free).** A rank-one fold is x ⊗ φ with x ∈ ker Δ and (φ⊗φ)∘Δ = 0.
- **Classification (exact over ℚ).** Every one of Q₁'s 426 components is J-preserving or kQ-preserving on all of it,
  except 6 components in the 6 vertex blocks whose corner is not a Q₀ block. Those 6 are of the mixed kind and have
  rank ≤ 1 identically.

Together with v2.60–v2.62 this accounts for every fold of Q₁ in terms of Q₀.*

---

## 0. Results

| # | Statement | Status |
|---|---|---|
| 76 | In char 0, every component of Q₁ outside the 6 cornerless blocks is J- or kQ-preserving on all of it; the components in those blocks are mixed, with rank ≤ 1 | **held**. Of the 426 components, 234 are split (including 1 point stratum), 98 J-preserving, 88 kQ-preserving (including 3 point strata), and 6 mixed with rank ≤ 1, exactly in the 6 blocks (`kinds_q1.txt`) |

The census observations behind the lemmas (`lemmas_AB_census.txt`, 7 sink-leaf extensions) were made first, then
proved.

## 1. Notation

The leaf is a sink: ℓ: u → v. B is the span of the paths of Q ending at u, so J = k e_v ⊕ Bℓ. For p ∈ kQ, write
Δp = Σ pre ⊗ suf over the deletions of one arrow.

Two facts about products:
- kQ·e_v = 0 = e_v·kQ;
- Bℓ·kQ = 0, and Bℓ·e_v = Bℓ.

As in v2.60, Δ(pℓ) = (1⊗R_ℓ)Δp + p ⊗ e_v, and K ∩ J = k e_v.

In both one-sided kinds P(e_v) has e_v-coefficient λ ∈ {0, 1}.

**Step 0 (both kinds, λ = 0): P_J = 0.**
- Write the J-component of P(pℓ) as φ(p)ℓ + μ(p)e_v.
- The kQ ⊗ e_v component of ΔP(pℓ) = (P⊗P)Δ(pℓ) gives φ(p) = Σ μ(suf)·P̄(pre). Any J-components of P on kQ, or kQ-components of P on J, land in other tensor components.
- The J-part of P² = P gives μ∘φ = μ.
- Corollary 4 of v2.60 then applies verbatim: induct on the length of p. The suffixes are shorter paths ending at u, so μ vanishes on them, then φ(p) = 0, then μ(p) = μ(φ(p)) = 0.

So P maps J into kQ, which is 0 in the J-preserving kind. ∎

## 2. Lemma A: J-preserving, λ = 0

**Lemma A.** Let P be a fold of kQ⁺ with P(J) ⊆ J and λ = 0. Then:
- P = 0 on J;
- for x ∈ kQ, P(x) = P̄(x) + ξ(x)e_v + ψ(x)ℓ, where ξ = ξ∘P̄ and ψ = (P̄ ⊗ ξ)∘Δ (that is, ψ(x) = Σ ξ(suf)·P̄(pre));
- the restriction of ξ to im P̄ is an algebra homomorphism im P̄ → k, possibly zero.

Hence over each P̄ there are at most 1 + #Hom_alg(im P̄, k) such folds, a finite number.

*Proof.*
- **P = 0 on J.** This is Step 0.
- **The form on kQ.** Write the J-part of P(x) as ξ(x)e_v + ψ(x)ℓ. The kQ ⊗ e_v component of ΔP(x) = (P⊗P)Δx gives
  ψ = (P̄⊗ξ)Δ.
- **ξ = ξ∘P̄.** Idempotency gives P(P(x)) = P̄x + ξ(P̄x)e_v + ψ(P̄x)ℓ, using P = 0 on J. This must equal P(x), so
  ξ∘P̄ = ξ and ψ∘P̄ = ψ.
- **Multiplicativity.** In particular the element of im P with kQ-part w ∈ im P̄ is unique: P(w) = w + ξ(w)e_v + ψ(w)ℓ.
  The product of P(x) and P(y) is P̄xP̄y + ξ(x)ξ(y)e_v + (…)ℓ, using the products listed in §1. It lies in im P and
  has kQ-part P̄xP̄y, so ξ(P̄xP̄y) = ξ(x)ξ(y) = ξ(P̄x)ξ(P̄y).
- **Finiteness.** A homomorphism from a finite-dimensional algebra A to k kills the nilpotent radical. It then factors
  through the semisimple quotient, a product of matrix algebras, and only the 1 × 1 factors map onto k. ∎

**Census check.** This holds on every J-preserving λ = 0 representative of the 7 sink-leaf extensions at q = 3. For
example, Q₀ → Q₁ has 276 such representatives: P_J = 0, ξ∘P̄ = ξ, and ξ is multiplicative on im P̄.

## 3. Lemma B: kQ-preserving, λ = 0

**Lemma B.** Let P be a fold of kQ⁺ with P(kQ) ⊆ kQ and λ = 0. Then:
- P(e_v) = y ∈ ker Δ ∩ im P̄, and ε = −y is an idempotent;
- P(pℓ) = Y(pℓ) ∈ im P̄, with ΔY(pℓ) = Σ P̄(pre) ⊗ Y(suf ℓ) + P̄(p) ⊗ y.

In particular ΔY(ℓ) = P̄(e_u) ⊗ y.

*Proof.*
- **The J-part vanishes.** By Step 0 the J-components vanish.
- **y ∈ ker Δ.** ΔP(e_v) = (P⊗P)Δe_v = 0.
- **The images lie in im P̄.** Idempotency gives P̄y = y and P̄Y = Y.
- **ε is idempotent.** e_v − y ∈ ker P, a subalgebra. By the products of §1, (e_v − y)² = e_v + y², whose image
  y + y² must vanish.
- **The Δ-identity.** It is the kQ ⊗ kQ component of ΔP(pℓ) = (P⊗P)Δ(pℓ). ∎

The idempotents of ker Δ come in families: e₁ + t(a − a′) is idempotent for every t. So Lemma B alone does not give
finiteness over a fixed P̄. Over the 366 classes of Q₀ (v2.62), every kQ-preserving λ = 0 fibre is finite: 345 single
points and 21 fibres with 48 points in all. The families are cut down by the Δ-identity and by membership in im P̄.
This is consistent with the census on the 7 extensions, where every representative satisfies Lemma B. A hand proof of
finiteness is still missing.

## 4. Lemma C: rank one

**Lemma C.** A rank-one fold P = x ⊗ φ (with φ(x) = 1) satisfies:
- x ∈ ker Δ;
- (φ⊗φ)∘Δ = 0;
- k x and ker φ are subalgebras.

Conversely, those conditions make x ⊗ φ a fold.

*Proof.* Δx = ΔP(x) = (P⊗P)Δx = c·x⊗x with c = (φ⊗φ)(Δx). Δ lowers the total path length by one, and x ⊗ x has a
component of total length twice the top degree of x, so c ≠ 0 is impossible. So c = 0 and Δx = 0. Then
ΔP(y) = φ(y)Δx = 0, so (P⊗P)Δy = (φ⊗φ)(Δy)·x⊗x = 0 for every y. The subalgebra conditions are the image and kernel
conditions of a fold. ∎

**On Q₀ → Q₁, all 16 mixed census classes have rank one**, at q = 3, 5, 7. Examples:
- P(e₃) = −e₁ − e₄, P(e₄) = e₁ + e₄: x = e₁ + e₄, φ = e₄* − e₃*;
- P(e₃) = −e₁ − e₄ + a − a′, P(e₄) = e₁ + e₄ − a + a′.

In each, x ∈ ker Δ = K₀ ⊕ k e₄ is an idempotent mixing e₄ with Q₀'s idempotents. In the other eight extensions the
mixed kind also has ranks 2, 3, 4 and 9 (`mixed_ranks.txt`), so "mixed means rank one" is special to Q₀ → Q₁.

## 5. The classification of Folds(Q₁) in characteristic 0 (`analysis/v2_63_kinds.py`)

**Method.** For each of the 426 components of Q₁'s strata (v2.59 certificates; the case-split block from its integral
leaves), the off-diagonal blocks of P = Pr(vs) are reduced modulo a Gröbner basis of the component's prime. These are
X, the J-rows of the kQ₀-columns, and Y, the kQ₀-rows of the J-columns. The kinds are closed conditions, so a
component is of one kind identically or is generically mixed.

| kind on the whole component | components |
|---|---|
| split (X ≡ 0 and Y ≡ 0) | 234 |
| J-preserving only (Y ≡ 0) | 98 |
| kQ₀-preserving only (X ≡ 0) | 88 |
| mixed, with all 2 × 2 minors of P ≡ 0 (rank ≤ 1) | 6: one in each of 00m1,00m1,0000,00m1, 00m1,0000,0000,00m1, 0000,m001,m001,m001, 0000,m001,0000,m001, 0000,00m1,0000,00m1, 0000,0000,m001,m001 |

(The split and kQ₀-preserving counts include 1 and 3 single-point strata.)

**So in characteristic 0 every fold of Q₁ is one of:**
- one-sided or split, with a compression that is a fold of Q₀;
- a rank-one fold in one of 6 components.

## 6. The leaf extension Q₀ → Q₁, all parts (characteristic 0)

| part | finiteness over each Q₀ class | how |
|---|---|---|
| split, λ = 0 | the trivial lift | proved (v2.60) |
| J-preserving, λ = 0 | ≤ 1 + #characters of im P̄ | **proved** (Lemma A) |
| kQ₀-preserving, λ = 0 | finite | structure proved (Lemma B); finiteness computed over the 366 classes (v2.62) |
| λ = 1, both kinds | finitely many Stab(P̄)-orbits | relative criterion, computed over the 366 classes (v2.62); the split part is an affine-linear bundle (v2.61) |
| mixed | 6 components, rank one | Lemma C; finitely many orbits by v2.57's certificates on those 6 blocks |

The census classes of Q₁ are thus finitely many, explained by Q₀'s 366 classes, the fibres over them, and the
rank-one mixed folds. The explanation has two computational parts. The first is the representatives of v2.62, whose
completeness rests on the per-block agreement with the census. The second is the fibre computations over them.

**What generalises to every sink leaf, proved:**
- the four kinds and the compression lemma (v2.60);
- the trivial lift, Corollary 4, and the identity lift (v2.60);
- Lemma A (characters), Lemma B (structure) and Lemma C (rank one).

**What does not:** the mixed kind has higher rank elsewhere.

**Next for Program 3:**
1. Prove finiteness in Lemma B: the idempotent and the Δ-identity over a fixed P̄.
2. Understand the mixed kind in general, by its kernel-level shape: K⁺ = K × k e_v, so mixed folds mix the idempotent
   e_v with K's.
3. The λ = 1 relative bound in general.

## Files

- analysis/v2_63_kinds.py, v2_63_lemmas_census.py, v2_63_mixed_ranks.py
- logs/v2_63/: PREDICTIONS.md, kinds_q1.txt, lemmas_AB_census.txt, mixed_ranks.txt
