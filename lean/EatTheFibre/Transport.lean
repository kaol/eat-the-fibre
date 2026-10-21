/-
# Transport: the fibre over a behaviour is a set of coproduct twists (v2.34)

Let `φ` be an algebra automorphism of a coordinate algebra `A` (multiplicative, with two-sided inverse `ψ`), not
necessarily compatible with Δ. Then `φ P φ⁻¹` is a fold of `A` iff `P` is a fold of `A.twist φ ψ`: the same algebra,
with Δ replaced by the transported coproduct `Δ^φ = (ψ ⊗ ψ) ∘ Δ ∘ φ` (`transport`).

Consequences, for the realisation / behaviour / fibre split of v2.33:
- image and kernel closure are invariant under all of Aut(A); only Δ-compatibility varies along a behaviour;
- the folds in the behaviour of `P` are `{φ P φ⁻¹ : P is Δ^φ-compatible}`, so the fibre over the behaviour of `P` is
  `X_P / Stab(P)` with `X_P = {φ ∈ Aut(A) : P is a fold of A.twist φ φ⁻¹}`, and its census orbits are the double cosets
  `G \ X_P / Stab(P)` (orbit–stabiliser, since G ⊆ Aut(A));
- a census-orbit count therefore decomposes as a sum over behaviours of double-coset counts. Growth in q can come from
  more behaviours or from larger `G \ X_P / Stab(P)`, and v2.33 found it is always the latter on the quivers tested.
-/
import EatTheFibre.Basic

open Matrix

namespace EatTheFibre.CoordAlg

variable {k : Type*} [CommRing k] {ι : Type*} [Fintype ι] [DecidableEq ι] (A : CoordAlg k ι)

/-- the same product, with Δ transported along `φ` (inverse `ψ`): `Δ' v = ψ · Δ(φ v) · ψᵀ` -/
def twist (φ ψ : Matrix ι ι k) : CoordAlg k ι :=
  ⟨A.c, fun l a b => ∑ m, φ m l * ∑ x, ∑ y, ψ a x * A.D m x y * ψ b y⟩

lemma twist_mul (φ ψ : Matrix ι ι k) (v w : ι → k) : (A.twist φ ψ).mul v w = A.mul v w := rfl

lemma twist_delta (φ ψ : Matrix ι ι k) (v : ι → k) : (A.twist φ ψ).delta v = ψ * A.delta (φ *ᵥ v) * ψᵀ := by
  ext a b
  simp only [delta, twist, Matrix.mul_apply, Matrix.transpose_apply, Matrix.mulVec, dotProduct, Finset.mul_sum,
    Finset.sum_mul]
  set f : ι → ι → ι → ι → k := fun l m x y => v l * (φ m l * (ψ a x * A.D m x y * ψ b y)) with hf
  calc ∑ l, ∑ m, ∑ x, ∑ y, f l m x y
      = ∑ m, ∑ l, ∑ x, ∑ y, f l m x y := Finset.sum_comm
    _ = ∑ m, ∑ x, ∑ l, ∑ y, f l m x y := Finset.sum_congr rfl fun m _ => Finset.sum_comm
    _ = ∑ m, ∑ x, ∑ y, ∑ l, f l m x y :=
        Finset.sum_congr rfl fun m _ => Finset.sum_congr rfl fun x _ => Finset.sum_comm
    _ = ∑ m, ∑ y, ∑ x, ∑ l, f l m x y := Finset.sum_congr rfl fun m _ => Finset.sum_comm
    _ = ∑ y, ∑ m, ∑ x, ∑ l, f l m x y := Finset.sum_comm
    _ = ∑ y, ∑ x, ∑ m, ∑ l, f l m x y := Finset.sum_congr rfl fun y _ => Finset.sum_comm
    _ = _ := by
        refine Finset.sum_congr rfl fun y _ => Finset.sum_congr rfl fun x _ => Finset.sum_congr rfl fun m _ =>
          Finset.sum_congr rfl fun l _ => ?_
        simp only [hf]; ring

variable {A}

/-- **Transport.** Conjugating by an algebra automorphism moves folds of the twisted algebra to folds of `A`. -/
theorem transport {φ ψ P : Matrix ι ι k} (hmul : ∀ v w, φ *ᵥ A.mul v w = A.mul (φ *ᵥ v) (φ *ᵥ w))
    (hψφ : ψ * φ = 1) (hφψ : φ * ψ = 1) : A.IsFold (φ * P * ψ) ↔ (A.twist φ ψ).IsFold P := by
  -- normal forms (everything right-associated)
  have L1 : ∀ M : Matrix ι ι k, ψ * (φ * M) = M := fun M => by rw [← Matrix.mul_assoc, hψφ, Matrix.one_mul]
  have L2 : ∀ M : Matrix ι ι k, φ * (ψ * M) = M := fun M => by rw [← Matrix.mul_assoc, hφψ, Matrix.one_mul]
  have L3 : ∀ M : Matrix ι ι k, φᵀ * (ψᵀ * M) = M := fun M => by
    rw [← Matrix.mul_assoc, ← Matrix.transpose_mul, hψφ, Matrix.transpose_one, Matrix.one_mul]
  have L4 : ∀ M : Matrix ι ι k, ψᵀ * (φᵀ * M) = M := fun M => by
    rw [← Matrix.mul_assoc, ← Matrix.transpose_mul, hφψ, Matrix.transpose_one, Matrix.one_mul]
  have L5 : φᵀ * ψᵀ = 1 := by rw [← Matrix.transpose_mul, hψφ, Matrix.transpose_one]
  have L6 : ψᵀ * φᵀ = 1 := by rw [← Matrix.transpose_mul, hφψ, Matrix.transpose_one]
  have cancelL : ∀ u u' : ι → k, φ *ᵥ u = φ *ᵥ u' → u = u' := fun u u' h => by
    have := congrArg (ψ *ᵥ ·) h; simpa [Matrix.mulVec_mulVec, hψφ] using this
  have conjv : ∀ u, (φ * P * ψ) *ᵥ (φ *ᵥ u) = φ *ᵥ (P *ᵥ u) := fun u => by
    simp only [Matrix.mulVec_mulVec, Matrix.mul_assoc, hψφ, Matrix.mul_one]
  have conjQ : 1 - φ * P * ψ = φ * (1 - P) * ψ := by rw [Matrix.mul_sub, Matrix.sub_mul, Matrix.mul_one, hφψ]
  have conjvQ : ∀ u, (1 - φ * P * ψ) *ᵥ (φ *ᵥ u) = φ *ᵥ ((1 - P) *ᵥ u) := fun u => by
    rw [conjQ]; simp only [Matrix.mulVec_mulVec, Matrix.mul_assoc, hψφ, Matrix.mul_one]
  have surj : ∀ v : ι → k, φ *ᵥ (ψ *ᵥ v) = v := fun v => by rw [Matrix.mulVec_mulVec, hφψ, Matrix.one_mulVec]
  constructor
  · intro h
    refine ⟨?_, fun v w => ?_, fun v w => ?_, fun v => ?_⟩
    · have e := congrArg (fun M => ψ * (M * φ)) h.idem
      simpa only [Matrix.mul_assoc, L1, hψφ, Matrix.mul_one] using e
    · rw [twist_mul]
      apply cancelL
      have := h.image (φ *ᵥ v) (φ *ᵥ w)
      rw [conjv, conjv, ← hmul, conjv] at this
      rw [this]
    · rw [twist_mul]
      apply cancelL
      have := h.kernel (φ *ᵥ v) (φ *ᵥ w)
      rw [conjvQ, conjvQ, ← hmul, conjv] at this
      rw [Matrix.mulVec_zero]; exact this
    · rw [twist_delta, twist_delta]
      have := h.delta (φ *ᵥ v)
      rw [conjv] at this
      rw [this]
      simp only [Matrix.transpose_mul, Matrix.mul_assoc, L1, L5, Matrix.mul_one]
  · intro h
    refine ⟨?_, fun v w => ?_, fun v w => ?_, fun v => ?_⟩
    · simp only [Matrix.mul_assoc, L1]
      rw [← Matrix.mul_assoc P P ψ, h.idem]
    · have := h.image (ψ *ᵥ v) (ψ *ᵥ w)
      rw [twist_mul] at this
      rw [← surj v, ← surj w, conjv, conjv, ← hmul, conjv, this]
    · have := h.kernel (ψ *ᵥ v) (ψ *ᵥ w)
      rw [twist_mul] at this
      rw [← surj v, ← surj w, conjvQ, conjvQ, ← hmul, conjv, this, Matrix.mulVec_zero]
    · have := h.delta (ψ *ᵥ v)
      rw [twist_delta, twist_delta, surj] at this
      have e := congrArg (fun M => φ * (M * φᵀ)) this
      simp only [Matrix.mul_assoc, L2, L6, Matrix.mul_one] at e
      rw [← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec, e]
      simp only [Matrix.transpose_mul, Matrix.mul_assoc]

end EatTheFibre.CoordAlg
