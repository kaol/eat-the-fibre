/-
# From the Discordianism axioms to the research project

The project began from open question 2 of the v2 axioms: "is there an idempotent Rota–Baxter(−1) fold on the grid
algebra that is also a morphism for the ε-coproduct Δ?" (v2.1 restated it as a question about subalgebra
splittings). This file makes that identification a theorem. Then it sorts the axioms into those that transfer to the
project's folds and those that do not.

* `isFold_iff_rb`: for any coordinate algebra (associativity is not needed), `P` is a fold iff it is idempotent, a
  Rota–Baxter operator of weight −1, and commutes with Δ. So the project's folds are exactly the "folds that
  deconstruct" of the axioms' open question 2.
* `fold_star`: δ-D1 for folds. `P` maps the Baxter star to the product, and `1 − P` maps it to minus the product.
  So (A, ⋆) ≅ im P × ker P for every fold.
* `thm4_isRB`, `thm3_isRB`: the general families of Theorems 3 and 4 (every acyclic quiver) are idempotent
  Rota–Baxter(−1) operators commuting with Δ.
* `thm4_not_dialgebra`: δ-L does **not** transfer. For Theorem 4's fold on 1⇉2⇉3, Loday's axiom A1 fails at every
  parameter value over every nontrivial ring. δ-L needs `f` to be an algebra endomorphism, and a fold need not be one.
-/
import EatTheFibre.General
import EatTheFibre.Discordia.Fold

open Matrix

namespace EatTheFibre

namespace CoordAlg

variable {k : Type*} [CommRing k] {ι : Type*} [Fintype ι] [DecidableEq ι] (A : CoordAlg k ι)

/-- `P` is a Rota–Baxter operator of weight −1 for the product of `A` -/
def IsRB (P : Matrix ι ι k) : Prop :=
  ∀ v w, A.mul (P *ᵥ v) (P *ᵥ w) = P *ᵥ (A.mul (P *ᵥ v) w + A.mul v (P *ᵥ w) - A.mul v w)

/-- the bilinear identity behind everything: `P x · y + x · P y − x y = P x · P y − (x − P x)(y − P y)` -/
lemma baxter_split (P : Matrix ι ι k) (v w : ι → k) :
    A.mul (P *ᵥ v) w + A.mul v (P *ᵥ w) - A.mul v w =
      A.mul (P *ᵥ v) (P *ᵥ w) - A.mul ((1 - P) *ᵥ v) ((1 - P) *ᵥ w) := by
  simp only [← mulₗ_apply, Matrix.sub_mulVec, Matrix.one_mulVec, map_sub, LinearMap.sub_apply]
  abel

lemma mul_zero' (v : ι → k) : A.mul v 0 = 0 := by rw [← mulₗ_apply, map_zero]
lemma zero_mul' (w : ι → k) : A.mul 0 w = 0 := by rw [← mulₗ_apply, map_zero, LinearMap.zero_apply]

/-- **Folds are idempotent Rota–Baxter(−1) operators commuting with Δ, and conversely.** -/
theorem isFold_iff_rb (P : Matrix ι ι k) :
    A.IsFold P ↔ P * P = P ∧ A.IsRB P ∧ ∀ v, A.delta (P *ᵥ v) = P * A.delta v * Pᵀ := by
  constructor
  · intro h
    refine ⟨h.idem, fun v w => ?_, h.delta⟩
    rw [baxter_split, Matrix.mulVec_sub, h.image, h.kernel, sub_zero]
  · rintro ⟨hid, hrb, hdel⟩
    have hPP : ∀ v, P *ᵥ (P *ᵥ v) = P *ᵥ v := fun v => by rw [Matrix.mulVec_mulVec, hid]
    have hz : ∀ u, (1 - P) *ᵥ (P *ᵥ u) = 0 := fun u => by
      rw [Matrix.sub_mulVec, Matrix.one_mulVec, hPP, sub_self]
    have hz' : ∀ u, P *ᵥ ((1 - P) *ᵥ u) = 0 := fun u => by
      rw [Matrix.sub_mulVec, Matrix.one_mulVec, Matrix.mulVec_sub, hPP, sub_self]
    refine ⟨hid, fun v w => ?_, fun v w => ?_, hdel⟩
    · have := hrb (P *ᵥ v) (P *ᵥ w)
      rw [baxter_split, hPP, hPP, hz, zero_mul', sub_zero] at this
      exact this.symm
    · have := hrb ((1 - P) *ᵥ v) ((1 - P) *ᵥ w)
      rw [hz' v, hz' w] at this
      simp only [zero_mul', mul_zero', zero_add, add_zero, zero_sub, Matrix.mulVec_neg] at this
      exact neg_eq_zero.mp this.symm

/-- the Baxter star of a matrix `P` -/
def starP (P : Matrix ι ι k) (v w : ι → k) : ι → k := A.mul (P *ᵥ v) w + A.mul v (P *ᵥ w) - A.mul v w

/-- **δ-D1 for folds.** `P` sends `⋆` to the product, and `1 − P` to minus the product: (A, ⋆) ≅ im P × ker P. -/
theorem fold_star {P : Matrix ι ι k} (h : A.IsFold P) (v w : ι → k) :
    P *ᵥ A.starP P v w = A.mul (P *ᵥ v) (P *ᵥ w) ∧
      (1 - P) *ᵥ A.starP P v w = - A.mul ((1 - P) *ᵥ v) ((1 - P) *ᵥ w) := by
  have hrb := ((A.isFold_iff_rb P).mp h).2.1
  have h1 : P *ᵥ A.starP P v w = A.mul (P *ᵥ v) (P *ᵥ w) := (hrb v w).symm
  refine ⟨h1, ?_⟩
  rw [Matrix.sub_mulVec, Matrix.one_mulVec, h1, starP, baxter_split]; abel

end CoordAlg

/-! ## The project's general families are idempotent Rota–Baxter(−1) operators -/

theorem thm4_isRB {Q : FinQuiver} [Fintype Q.Path] (C : Consecutive Q) {k : Type*} [CommRing k] (t : k) :
    (Q.pathAlg k).IsRB (CoordAlg.pushM C.ret (Reflect.evalM Q_two_pairs.thm4P t)) :=
  (((Q.pathAlg k).isFold_iff_rb _).mp (thm4_general C t)).2.1

theorem thm3_isRB {Q : FinQuiver} [Fintype Q.Path] (C : Converging Q) {k : Type*} [CommRing k] (t : k) :
    (Q.pathAlg k).IsRB (CoordAlg.pushM C.ret (Reflect.evalM Q_pairs_in.thm3F1P t)) :=
  (((Q.pathAlg k).isFold_iff_rb _).mp (thm3_general C t).1).2.1

/-! ## δ-L does not transfer to folds -/

namespace Q_two_pairs

open EatTheFibre.Reflect EatTheFibre.Poly

/-- basis vectors as constant polynomial vectors -/
def eP (i : Fin 11) : Fin 11 → Poly := fun j => if j = i then [1] else []

/-- `(x ⊣ y) ⊣ z` with `x ⊣ y = x · P y`, symbolically -/
def lhsV (x y z : Fin 11) : Fin 11 → Poly := algMulP c (algMulP c (eP x) (mulVecP thm4P (eP y))) (mulVecP thm4P (eP z))
/-- `x ⊣ (y ⊣ z)` -/
def rhsV (x y z : Fin 11) : Fin 11 → Poly := algMulP c (eP x) (mulVecP thm4P (algMulP c (eP y) (mulVecP thm4P (eP z))))

/-- at `(e₂, e₂, b₂)`, coordinate `e₂`: the left side is `0`, the right side is `1` (independent of `t`) -/
lemma lhs_val : isZero (lhsV 1 1 6 1) = true := by decide +kernel
lemma rhs_val : isZero (sub (rhsV 1 1 6 1) [1]) = true := by decide +kernel

section
variable {k : Type*} [CommRing k] (t : k)

/-- evaluate a polynomial vector -/
def ev (v : Fin 11 → Poly) : Fin 11 → k := fun i => eval t (v i)

lemma ev_eP (i : Fin 11) : ev t (eP i) = Pi.single i 1 := by
  funext j; by_cases h : j = i <;> simp [ev, eP, h]

lemma ev_mulVecP (Q : Fin 11 → Fin 11 → Poly) (v : Fin 11 → Poly) :
    ev t (mulVecP Q v) = evalM Q t *ᵥ ev t v := by
  funext i; exact mulVecP_eval Q t v i

lemma ev_algMulP (u w : Fin 11 → Poly) : ev t (algMulP c u w) = (alg k).mul (ev t u) (ev t w) := by
  funext m; exact algMulP_eval c D t u w m

end

/-- **δ-L fails for project folds.** For Theorem 4's fold on 1⇉2⇉3, Loday's A1 `(x ⊣ y) ⊣ z = x ⊣ (y ⊣ z)`
(with `x ⊣ y = x · P y`) fails at `(e₂, e₂, b₂)`, for every `t` and over every nontrivial commutative ring.
(The full list of failures, computed in `analysis/`, is 11 basis triples. All lie along the arrows `a₂`, `b₂` that
Theorem 4's fold twists.) -/
theorem thm4_not_dialgebra {k : Type*} [CommRing k] [Nontrivial k] (t : k) :
    (alg k).mul ((alg k).mul (Pi.single 1 1) (evalM thm4P t *ᵥ Pi.single 1 1)) (evalM thm4P t *ᵥ Pi.single 6 1) ≠
      (alg k).mul (Pi.single 1 1)
        (evalM thm4P t *ᵥ (alg k).mul (Pi.single 1 1) (evalM thm4P t *ᵥ Pi.single 6 1)) := by
  intro heq
  have eL : ev t (lhsV 1 1 6) = (alg k).mul ((alg k).mul (Pi.single 1 1) (evalM thm4P t *ᵥ Pi.single 1 1))
      (evalM thm4P t *ᵥ Pi.single 6 1) := by
    simp only [lhsV, ev_algMulP, ev_mulVecP, ev_eP]
  have eR : ev t (rhsV 1 1 6) = (alg k).mul (Pi.single 1 1)
      (evalM thm4P t *ᵥ (alg k).mul (Pi.single 1 1) (evalM thm4P t *ᵥ Pi.single 6 1)) := by
    simp only [rhsV, ev_algMulP, ev_mulVecP, ev_eP]
  have h1 := congrFun (eL.trans (heq.trans eR.symm)) 1
  simp only [ev] at h1
  rw [eval_eq_zero_of_isZero t _ lhs_val, eval_eq_of_isZero_sub t _ _ rhs_val] at h1
  simp at h1

/-- …and the fold is not multiplicative: `P (e₁ a₂) ≠ P e₁ · P a₂` (coordinate `a₂`: `1` vs `0`). So it is not an
M3 endomorphism, and the hypothesis of δ-L is genuinely absent. -/
lemma multL_val : isZero (sub (mulVecP thm4P (algMulP c (eP 0) (eP 4)) 4) [1]) = true := by decide +kernel
lemma multR_val : isZero (algMulP c (mulVecP thm4P (eP 0)) (mulVecP thm4P (eP 4)) 4) = true := by decide +kernel

theorem thm4_not_multiplicative {k : Type*} [CommRing k] [Nontrivial k] (t : k) :
    evalM thm4P t *ᵥ (alg k).mul (Pi.single 0 1) (Pi.single 4 1) ≠
      (alg k).mul (evalM thm4P t *ᵥ Pi.single 0 1) (evalM thm4P t *ᵥ Pi.single 4 1) := by
  intro heq
  have eL : ev t (mulVecP thm4P (algMulP c (eP 0) (eP 4))) =
      evalM thm4P t *ᵥ (alg k).mul (Pi.single 0 1) (Pi.single 4 1) := by
    simp only [ev_algMulP, ev_mulVecP, ev_eP]
  have eR : ev t (algMulP c (mulVecP thm4P (eP 0)) (mulVecP thm4P (eP 4))) =
      (alg k).mul (evalM thm4P t *ᵥ Pi.single 0 1) (evalM thm4P t *ᵥ Pi.single 4 1) := by
    simp only [ev_algMulP, ev_mulVecP, ev_eP]
  have h1 := congrFun (eL.trans (heq.trans eR.symm)) 4
  simp only [ev] at h1
  rw [eval_eq_of_isZero_sub t _ _ multL_val, eval_eq_zero_of_isZero t _ multR_val] at h1
  simp at h1

/-- **Theorem 4's fold is in the bottom class of `Completion`.** Its kernel is neither a right nor a left module over
its image: `P((e₁ − P e₁) · P a₂)` and `P(P a₂ · (e₂ − P e₂))` both have `a₂`-coordinate `1`. Generically such a fold
has the Baxter point as its only dendriform completion (`analysis/completion_generic.py`). For this fold that is
confirmed at symbolic `t`. -/
lemma leftMod_val : isZero (sub (mulVecP thm4P (algMulP c (mulVecP (oneSubP thm4P) (eP 0)) (mulVecP thm4P (eP 4))) 4)
    [1]) = true := by decide +kernel
lemma rightMod_val : isZero (sub (mulVecP thm4P (algMulP c (mulVecP thm4P (eP 4)) (mulVecP (oneSubP thm4P) (eP 1))) 4)
    [1]) = true := by decide +kernel

theorem thm4_not_module {k : Type*} [CommRing k] [Nontrivial k] (t : k) :
    let P := evalM thm4P t
    (P *ᵥ (alg k).mul ((1 - P) *ᵥ Pi.single 0 1) (P *ᵥ Pi.single 4 1)) ≠ 0 ∧
      (P *ᵥ (alg k).mul (P *ᵥ Pi.single 4 1) ((1 - P) *ᵥ Pi.single 1 1)) ≠ 0 := by
  intro P
  constructor
  · intro heq
    have e : ev t (mulVecP thm4P (algMulP c (mulVecP (oneSubP thm4P) (eP 0)) (mulVecP thm4P (eP 4)))) =
        P *ᵥ (alg k).mul ((1 - P) *ᵥ Pi.single 0 1) (P *ᵥ Pi.single 4 1) := by
      simp only [ev_algMulP, ev_mulVecP, ev_eP, oneSub_eval, P]
    have h1 := congrFun (e.trans heq) 4
    simp only [ev] at h1
    rw [eval_eq_of_isZero_sub t _ _ leftMod_val] at h1
    simp at h1
  · intro heq
    have e : ev t (mulVecP thm4P (algMulP c (mulVecP thm4P (eP 4)) (mulVecP (oneSubP thm4P) (eP 1)))) =
        P *ᵥ (alg k).mul (P *ᵥ Pi.single 4 1) ((1 - P) *ᵥ Pi.single 1 1) := by
      simp only [ev_algMulP, ev_mulVecP, ev_eP, oneSub_eval, P]
    have h1 := congrFun (e.trans heq) 4
    simp only [ev] at h1
    rw [eval_eq_of_isZero_sub t _ _ rightMod_val] at h1
    simp at h1

end Q_two_pairs

end EatTheFibre
