/-
# Fold checks by reflection

An algebra with integer structure constants `c`, `D` gives `CoordAlg.ofInt c D` over every commutative ring. A matrix
family with polynomial entries `P : Fin d → Fin d → Poly` evaluates at `t : k` to `evalM t P`.

`check c D P : Bool` computes the four basis conditions of a fold (idempotency, image and kernel closure, Δ) on `Poly`
and tests each difference with `isZero`. The kernel evaluates it (`decide +kernel`). `check_sound` turns
`check c D P = true` into `IsFold` over every commutative ring and every `t`.
-/
import EatTheFibre.Basic
import EatTheFibre.Poly

open Matrix

namespace EatTheFibre

variable {d : ℕ}

/-- the coordinate algebra with integer structure constants, over any commutative ring -/
def CoordAlg.ofInt (k : Type*) [CommRing k] (c D : Fin d → Fin d → Fin d → ℤ) : CoordAlg k (Fin d) :=
  ⟨fun i j m => (c i j m : k), fun l a b => (D l a b : k)⟩

namespace Reflect

open EatTheFibre.Poly

variable (c D : Fin d → Fin d → Fin d → ℤ) (P : Fin d → Fin d → Poly)

/-- `(P *ᵥ v) i` -/
def mulVecP (v : Fin d → Poly) (i : Fin d) : Poly := sumF fun j => mul (P i j) (v j)
/-- `(P * Q) i j` -/
def matMulP (Q : Fin d → Fin d → Poly) (i j : Fin d) : Poly := sumF fun l => mul (P i l) (Q l j)
/-- algebra product -/
def algMulP (u w : Fin d → Poly) (m : Fin d) : Poly :=
  sumF fun i => sumF fun j => if c i j m = 0 then [] else smul (c i j m) (mul (u i) (w j))
/-- `1 - P` -/
def oneSubP (i j : Fin d) : Poly := sub (if i = j then [1] else []) (P i j)
/-- `Δ v` -/
def deltaP (v : Fin d → Poly) (a b : Fin d) : Poly := sumF fun l => if D l a b = 0 then [] else smul (D l a b) (v l)
/-- `(P · Δ(e_l) · Pᵀ) a b` -/
def conjDeltaP (l a b : Fin d) : Poly :=
  sumF fun x => sumF fun y => if D l x y = 0 then [] else smul (D l x y) (mul (P a x) (P b y))

/-- column `j` of a polynomial matrix -/
def colP (Q : Fin d → Fin d → Poly) (j : Fin d) : Fin d → Poly := fun i => Q i j

/-- idempotency -/
def idemOK : Bool :=
  (List.finRange d).all fun i => (List.finRange d).all fun j => isZero (sub (matMulP P P i j) (P i j))

/-- image closure, for the products `P e_i · P e_j`, all `j` -/
def imageOK (i : Fin d) : Bool :=
  (List.finRange d).all fun j => (List.finRange d).all fun m =>
    isZero (sub (mulVecP P (algMulP c (colP P i) (colP P j)) m) (algMulP c (colP P i) (colP P j) m))

/-- kernel closure, for the products `(1 - P) e_i · (1 - P) e_j`, all `j` -/
def kernelOK (i : Fin d) : Bool :=
  (List.finRange d).all fun j => (List.finRange d).all fun m =>
    isZero (mulVecP P (algMulP c (colP (oneSubP P) i) (colP (oneSubP P) j)) m)

/-- Δ-compatibility at the basis vector `e_l` -/
def deltaOK (l : Fin d) : Bool :=
  (List.finRange d).all fun a => (List.finRange d).all fun b =>
    isZero (sub (deltaP D (colP P l) a b) (conjDeltaP D P l a b))

/-- all four basis conditions (split so that each piece can be checked by the kernel on its own) -/
def check : Bool :=
  idemOK P && ((List.finRange d).all fun i => imageOK c P i) && ((List.finRange d).all fun i => kernelOK c P i) &&
    ((List.finRange d).all fun l => deltaOK D P l)

variable {k : Type*} [CommRing k] (t : k)

/-- the evaluated matrix -/
def evalM : Matrix (Fin d) (Fin d) k := fun i j => eval t (P i j)

lemma all_finRange {p : Fin d → Bool} : (List.finRange d).all p = true ↔ ∀ i, p i = true := by
  simp [List.all_eq_true, List.mem_finRange]

lemma col_eval (j : Fin d) : (evalM P t) *ᵥ Pi.single j 1 = fun i => eval t (P i j) := by
  funext i; simp [Matrix.mulVec, dotProduct, Pi.single_apply, evalM]

lemma oneSub_eval : 1 - evalM P t = evalM (oneSubP P) t := by
  ext i j; by_cases h : i = j <;> simp [oneSubP, evalM, h]

lemma mulVecP_eval (v : Fin d → Poly) (i : Fin d) :
    eval t (mulVecP P v i) = ((evalM P t) *ᵥ fun j => eval t (v j)) i := by
  simp [mulVecP, eval_sumF, Matrix.mulVec, dotProduct, evalM]

lemma algMulP_eval (u w : Fin d → Poly) (m : Fin d) :
    eval t (algMulP c u w m) = (CoordAlg.ofInt k c D).mul (fun i => eval t (u i)) (fun j => eval t (w j)) m := by
  simp only [algMulP, eval_sumF, CoordAlg.mul, CoordAlg.ofInt]
  exact Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => by
    split_ifs with h0 <;> simp [h0]; ring

/-- **Soundness.** A successful check proves the fold conditions over every commutative ring, at every `t`. -/
lemma check_of_parts (h0 : idemOK P = true) (h1 : ∀ i, imageOK c P i = true) (h2 : ∀ i, kernelOK c P i = true)
    (h3 : ∀ l, deltaOK D P l = true) : check c D P = true := by
  simp only [check, Bool.and_eq_true, all_finRange]; exact ⟨⟨⟨h0, h1⟩, h2⟩, h3⟩

theorem check_sound (h : check c D P = true) : (CoordAlg.ofInt k c D).IsFold (evalM P t) := by
  simp only [check, idemOK, imageOK, kernelOK, deltaOK, Bool.and_eq_true, all_finRange] at h
  obtain ⟨⟨⟨hid, him⟩, hker⟩, hdel⟩ := h
  rw [CoordAlg.isFold_iff_basis]
  refine ⟨?_, ?_, ?_, ?_⟩
  · ext i j
    have := eval_eq_of_isZero_sub t _ _ (hid i j)
    simpa [matMulP, eval_sumF, Matrix.mul_apply, evalM] using this
  · intro i j; funext m
    have := eval_eq_of_isZero_sub t _ _ (him i j m)
    rw [mulVecP_eval] at this
    have e : (fun x => eval t (algMulP c (colP P i) (colP P j) x)) =
        (CoordAlg.ofInt k c D).mul (evalM P t *ᵥ Pi.single i 1) (evalM P t *ᵥ Pi.single j 1) := by
      funext x; rw [algMulP_eval (D := D), col_eval, col_eval]; rfl
    rw [e, show eval t (algMulP c (colP P i) (colP P j) m) = _ from congrFun e m] at this
    exact this
  · intro i j; funext m
    have := eval_eq_zero_of_isZero t _ (hker i j m)
    rw [mulVecP_eval] at this
    have e : (fun x => eval t (algMulP c (colP (oneSubP P) i) (colP (oneSubP P) j) x)) =
        (CoordAlg.ofInt k c D).mul ((1 - evalM P t) *ᵥ Pi.single i 1) ((1 - evalM P t) *ᵥ Pi.single j 1) := by
      funext x; rw [algMulP_eval (D := D), oneSub_eval, col_eval, col_eval]; rfl
    rw [e] at this; simpa using this
  · intro l; ext a b
    have h1 := eval_eq_of_isZero_sub t _ _ (hdel l a b)
    simp only [deltaP, conjDeltaP, colP, eval_sumF, apply_ite (eval t), eval_nil, eval_smul, eval_mul] at h1
    have z1 : ∀ x, (if D x a b = 0 then (0 : k) else (D x a b : k) * eval t (P x l)) = (D x a b : k) * eval t (P x l) :=
      fun x => by split_ifs with h0 <;> simp [h0]
    have z2 : ∀ x y, (if D l x y = 0 then (0 : k) else (D l x y : k) * (eval t (P a x) * eval t (P b y))) =
        (D l x y : k) * (eval t (P a x) * eval t (P b y)) := fun x y => by split_ifs with h0 <;> simp [h0]
    simp only [z1, z2] at h1
    have hs : ∀ x y, (CoordAlg.ofInt k c D).delta (Pi.single l 1) x y = (D l x y : k) := by
      intro x y; simp [CoordAlg.delta, CoordAlg.ofInt, Pi.single_apply]
    rw [col_eval]
    simp only [Matrix.mul_apply, Matrix.transpose_apply, hs]
    simp only [CoordAlg.delta, CoordAlg.ofInt, evalM]
    calc ∑ x, eval t (P x l) * (D x a b : k) = ∑ x, (D x a b : k) * eval t (P x l) :=
          Finset.sum_congr rfl fun x _ => mul_comm _ _
      _ = ∑ x, ∑ y, (D l x y : k) * (eval t (P a x) * eval t (P b y)) := h1
      _ = ∑ y, (∑ x, eval t (P a x) * (D l x y : k)) * eval t (P b y) := by
          rw [Finset.sum_comm]
          exact Finset.sum_congr rfl fun y _ => by rw [Finset.sum_mul]; exact Finset.sum_congr rfl fun x _ => by ring

end Reflect

end EatTheFibre

namespace EatTheFibre.Reflect

variable {d : ℕ} (c D : Fin d → Fin d → Fin d → ℤ) (P : Fin d → Fin d → Poly)

/-! A failing part makes the whole check fail (used by the generated controls). -/

lemma check_false_of_idem (h : idemOK P = false) : check c D P = false := by simp [check, h]
lemma check_false_of_image (i : Fin d) (h : imageOK c P i = false) : check c D P = false := by
  simp only [check, Bool.and_eq_false_iff]; left; left; right
  rw [List.all_eq_false]; exact ⟨i, List.mem_finRange i, by simp [h]⟩
lemma check_false_of_kernel (i : Fin d) (h : kernelOK c P i = false) : check c D P = false := by
  simp only [check, Bool.and_eq_false_iff]; left; right
  rw [List.all_eq_false]; exact ⟨i, List.mem_finRange i, by simp [h]⟩
lemma check_false_of_delta (l : Fin d) (h : deltaOK D P l = false) : check c D P = false := by
  simp only [check, Bool.and_eq_false_iff]; right
  rw [List.all_eq_false]; exact ⟨l, List.mem_finRange l, by simp [h]⟩

end EatTheFibre.Reflect
