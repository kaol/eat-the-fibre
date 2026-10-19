/-
# Folds of coordinate algebras

A finite-dimensional algebra is given by structure constants on the coordinate space `ι → k` (`ι` a finite basis):
`(v * w) m = ∑ i, ∑ j, v i * w j * c i j m`. Arrow deletion `Δ : A → A ⊗ A` is given by
`Δ(v) a b = ∑ l, v l * D l a b` (a `d × d` matrix per vector).

A *fold* is a linear `P` (a matrix acting by `mulVec`) that is idempotent, whose image and kernel are closed under
multiplication, and that commutes with `Δ`: `Δ (P v) = P · Δ v · Pᵀ`.

The image of `P` is `{P v}` and the kernel of an idempotent `P` is `{(1 - P) v}`, so closure is stated for all
`v, w` directly. `isFold_iff_basis` reduces every condition to the standard basis. That is exactly what the Python
checks (`theorem*.py`, `family_*.py`, `q_fold` in `pcensus.py`) verify, so it certifies their method.
-/
import Mathlib

open Matrix

namespace EatTheFibre

variable {k : Type*} [CommRing k] {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- structure constants of an algebra on `ι → k` (`ι` a finite basis), and of arrow deletion `Δ` -/
structure CoordAlg (k : Type*) (ι : Type*) where
  c : ι → ι → ι → k
  D : ι → ι → ι → k

namespace CoordAlg

variable (A : CoordAlg k ι)

/-- multiplication from structure constants -/
def mul (v w : ι → k) : ι → k := fun m => ∑ i, ∑ j, v i * w j * A.c i j m

/-- arrow deletion, as a matrix: `Δ v a b = ∑ l, v l * D l a b` -/
def delta (v : ι → k) : Matrix (ι) (ι) k := fun a b => ∑ l, v l * A.D l a b

/-- the defining conditions of a fold, for all vectors -/
structure IsFold (P : Matrix (ι) (ι) k) : Prop where
  idem : P * P = P
  image : ∀ v w, P *ᵥ A.mul (P *ᵥ v) (P *ᵥ w) = A.mul (P *ᵥ v) (P *ᵥ w)
  kernel : ∀ v w, P *ᵥ A.mul ((1 - P) *ᵥ v) ((1 - P) *ᵥ w) = 0
  delta : ∀ v, A.delta (P *ᵥ v) = P * A.delta v * Pᵀ

/-- the same conditions on the standard basis only (what the computer checks) -/
structure IsFoldBasis (P : Matrix (ι) (ι) k) : Prop where
  idem : P * P = P
  image : ∀ i j, P *ᵥ A.mul (P *ᵥ Pi.single i 1) (P *ᵥ Pi.single j 1) =
    A.mul (P *ᵥ Pi.single i 1) (P *ᵥ Pi.single j 1)
  kernel : ∀ i j, P *ᵥ A.mul ((1 - P) *ᵥ Pi.single i 1) ((1 - P) *ᵥ Pi.single j 1) = 0
  delta : ∀ l, A.delta (P *ᵥ Pi.single l 1) = P * A.delta (Pi.single l 1) * Pᵀ

/-- multiplication as a bilinear map -/
def mulₗ : (ι → k) →ₗ[k] (ι → k) →ₗ[k] (ι → k) :=
  LinearMap.mk₂ k A.mul
    (fun u v w => by ext m; simp [mul, add_mul, Finset.sum_add_distrib])
    (fun a v w => by ext m; simp [mul, Finset.mul_sum, mul_assoc])
    (fun u v w => by ext m; simp [mul, mul_add, add_mul, Finset.sum_add_distrib])
    (fun a v w => by
      ext m; simp only [mul, Pi.smul_apply, smul_eq_mul, Finset.mul_sum]
      exact Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => by ring)

@[simp] lemma mulₗ_apply (v w : ι → k) : A.mulₗ v w = A.mul v w := rfl

/-- arrow deletion as a linear map -/
def deltaₗ : (ι → k) →ₗ[k] Matrix (ι) (ι) k where
  toFun := A.delta
  map_add' u v := by ext a b; simp [delta, add_mul, Finset.sum_add_distrib]
  map_smul' c v := by ext a b; simp [delta, Finset.mul_sum, mul_assoc]

@[simp] lemma deltaₗ_apply (v : ι → k) : A.deltaₗ v = A.delta v := rfl

/-- **Basis criterion.** The fold conditions hold for all vectors iff they hold on the standard basis. -/
theorem isFold_iff_basis (P : Matrix (ι) (ι) k) : A.IsFold P ↔ A.IsFoldBasis P := by
  constructor
  · intro h; exact ⟨h.idem, fun i j => h.image _ _, fun i j => h.kernel _ _, fun l => h.delta _⟩
  · intro h
    have hb : ∀ i, (Pi.basisFun k (ι)) i = Pi.single i 1 := fun i => by simp
    refine ⟨h.idem, ?_, ?_, ?_⟩
    · -- both sides are bilinear in (v, w)
      have := LinearMap.ext_basis (Pi.basisFun k (ι)) (Pi.basisFun k (ι))
        (B := (A.mulₗ.compl₁₂ P.mulVecLin P.mulVecLin).compr₂ P.mulVecLin)
        (B' := A.mulₗ.compl₁₂ P.mulVecLin P.mulVecLin) (fun i j => by simpa [hb] using h.image i j)
      intro v w; simpa using LinearMap.congr_fun₂ this v w
    · have := LinearMap.ext_basis (Pi.basisFun k (ι)) (Pi.basisFun k (ι))
        (B := (A.mulₗ.compl₁₂ (1 - P).mulVecLin (1 - P).mulVecLin).compr₂ P.mulVecLin) (B' := 0)
        (fun i j => by
          simpa only [LinearMap.compr₂_apply, LinearMap.compl₁₂_apply, Matrix.mulVecLin_apply, mulₗ_apply, hb,
            LinearMap.zero_apply] using h.kernel i j)
      intro v w
      simpa only [LinearMap.compr₂_apply, LinearMap.compl₁₂_apply, Matrix.mulVecLin_apply, mulₗ_apply,
        LinearMap.zero_apply] using LinearMap.congr_fun₂ this v w
    · -- both sides are linear in v
      let R : Matrix (ι) (ι) k →ₗ[k] Matrix (ι) (ι) k :=
        { toFun := fun M => P * M * Pᵀ
          map_add' := fun M N => by simp [Matrix.mul_add, Matrix.add_mul]
          map_smul' := fun c M => by simp }
      have := (Pi.basisFun k (ι)).ext (f₁ := A.deltaₗ ∘ₗ P.mulVecLin) (f₂ := R ∘ₗ A.deltaₗ)
        (fun l => by simpa [hb, R] using h.delta l)
      intro v; simpa [R] using LinearMap.congr_fun this v

end CoordAlg

end EatTheFibre

namespace EatTheFibre.CoordAlg

variable {k : Type*} [CommRing k] {ι : Type*} [Fintype ι] [DecidableEq ι] (A : CoordAlg k ι)

/-- `h` is a Δ-compatible algebra automorphism with inverse `hi` -/
structure IsDeltaAut (h hi : Matrix (ι) (ι) k) : Prop where
  left : h * hi = 1
  right : hi * h = 1
  mul : ∀ v w, h *ᵥ A.mul v w = A.mul (h *ᵥ v) (h *ᵥ w)
  delta : ∀ v, A.delta (h *ᵥ v) = h * A.delta v * hᵀ

/-- **Orbits.** Conjugating a fold by a Δ-compatible algebra automorphism gives a fold. -/
theorem IsFold.conj {P h hi : Matrix (ι) (ι) k} (hP : A.IsFold P) (hh : A.IsDeltaAut h hi) :
    A.IsFold (h * P * hi) := by
  have hhi : ∀ v, h *ᵥ (hi *ᵥ v) = v := fun v => by rw [Matrix.mulVec_mulVec, hh.left, Matrix.one_mulVec]
  have hih : ∀ v, hi *ᵥ (h *ᵥ v) = v := fun v => by rw [Matrix.mulVec_mulVec, hh.right, Matrix.one_mulVec]
  have conjv : ∀ v, (h * P * hi) *ᵥ v = h *ᵥ (P *ᵥ (hi *ᵥ v)) := fun v => by
    simp [Matrix.mulVec_mulVec, Matrix.mul_assoc]
  have conjq : ∀ v, (1 - h * P * hi) *ᵥ v = h *ᵥ ((1 - P) *ᵥ (hi *ᵥ v)) := fun v => by
    rw [Matrix.sub_mulVec, Matrix.one_mulVec, conjv, Matrix.sub_mulVec, Matrix.one_mulVec, Matrix.mulVec_sub, hhi]
  refine ⟨?_, ?_, ?_, ?_⟩
  · calc h * P * hi * (h * P * hi) = h * P * (hi * h) * P * hi := by simp only [Matrix.mul_assoc]
      _ = h * P * hi := by rw [hh.right, Matrix.mul_one, Matrix.mul_assoc h P P, hP.idem]
  · intro v w
    simp only [conjv]; rw [← hh.mul, hih, hP.image]
  · intro v w
    simp only [conjq, conjv]; rw [← hh.mul, hih, hP.kernel, Matrix.mulVec_zero]
  · intro v
    -- Δ(hi v) = hi Δ(v) hiᵀ, from compatibility of h applied at hi v
    have dhi : A.delta (hi *ᵥ v) = hi * A.delta v * hiᵀ := by
      have e := hh.delta (hi *ᵥ v); rw [hhi] at e
      rw [e]
      calc A.delta (hi *ᵥ v) = (hi * h) * A.delta (hi *ᵥ v) * (hi * h)ᵀ := by rw [hh.right]; simp
        _ = hi * (h * A.delta (hi *ᵥ v) * hᵀ) * hiᵀ := by
          rw [Matrix.transpose_mul]; simp only [Matrix.mul_assoc]
    rw [conjv, hh.delta, hP.delta, dhi, Matrix.transpose_mul, Matrix.transpose_mul]
    simp only [Matrix.mul_assoc]

/-- **Vertex-block lemma.** Let `isV` mark the vertex coordinates. If `h` fixes every vertex basis vector and maps
non-vertex basis vectors into the span of non-vertex ones, and `hi` is its inverse, then conjugation by `h` does not
change the vertex-to-vertex entries of any matrix `P`. -/
theorem vertex_block (isV : ι → Prop) [DecidablePred isV] {h hi : Matrix (ι) (ι) k}
    (hinv : hi * h = 1) (hfix : ∀ a, isV a → ∀ i, h i a = if i = a then 1 else 0)
    (hrad : ∀ a, ¬ isV a → ∀ i, isV i → h i a = 0) (P : Matrix (ι) (ι) k)
    {i v : ι} (hi_ : isV i) (hv : isV v) : (h * P * hi) i v = P i v := by
  -- the column of `hi` at a vertex is the basis vector: hi = hi * (h * e_v) = e_v
  have hcol : ∀ b, hi b v = if b = v then 1 else 0 := by
    intro b
    have : (hi * h) b v = (1 : Matrix (ι) (ι) k) b v := by rw [hinv]
    rw [Matrix.mul_apply] at this
    rw [Finset.sum_eq_single v (fun a _ hav => by simp [hfix v hv a, hav])
      (by simp)] at this
    simp only [hfix v hv v, ite_true, mul_one] at this
    rw [this, Matrix.one_apply]
  -- the row of `h` at a vertex is the basis covector
  have hrow : ∀ a, h i a = if i = a then 1 else 0 := by
    intro a
    by_cases ha : isV a
    · exact hfix a ha i
    · rw [hrad a ha i hi_]; split_ifs with hia
      · exact absurd (hia ▸ hi_) ha
      · rfl
  simp only [Matrix.mul_apply, hcol, hrow]
  simp

end EatTheFibre.CoordAlg
