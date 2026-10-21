/-
# Reflection for conjugating automorphisms (v2.33)

A one-parameter family of folds `P_s` is a single orbit if there are matrices `Φ_s`, `Ψ_s` (polynomial in `s`) with
* `Φ_s` multiplicative on the algebra (`homOK`),
* `Ψ_s Φ_s = 1 = Φ_s Ψ_s` (`invOK`),
* `Φ_s P₀ = P_s Φ_s`, where `P₀` is the constant-term matrix of `P` (`interOK`).

These are univariate polynomial identities, checked by `decide +kernel` as in `Reflect.lean`. `orbit_sound` turns
them into the statement over every commutative ring at every `s`: `Φ_s` is an algebra automorphism (a bijective
multiplicative linear map) carrying `P_0` to `P_s`. It is not Δ-compatible in general, which is the point.
-/
import EatTheFibre.Reflect

open Matrix

namespace EatTheFibre.Reflect

open EatTheFibre.Poly

variable {d : ℕ} (c : Fin d → Fin d → Fin d → ℤ)

/-- the basis vector `e_i` as a constant polynomial vector -/
def eP (i : Fin d) : Fin d → Poly := fun j => if j = i then [1] else []

/-- the constant-term matrix -/
def constP (P : Fin d → Fin d → Poly) (i j : Fin d) : Poly := match P i j with | [] => [] | a :: _ => [a]

/-- multiplicativity on the products `e_i · e_j`, all `j` (one row, for a separate kernel check) -/
def homRow (Φ : Fin d → Fin d → Poly) (i : Fin d) : Bool :=
  (List.finRange d).all fun j => (List.finRange d).all fun m =>
    isZero (sub (mulVecP Φ (algMulP c (eP i) (eP j)) m) (algMulP c (colP Φ i) (colP Φ j) m))

def homOK (Φ : Fin d → Fin d → Poly) : Bool := (List.finRange d).all fun i => homRow c Φ i

lemma homOK_of_rows (Φ : Fin d → Fin d → Poly) (h : ∀ i, homRow c Φ i = true) : homOK c Φ = true := by
  simp only [homOK, all_finRange]; exact h

def idP (i j : Fin d) : Poly := if i = j then [1] else []

def invOK (Φ Ψ : Fin d → Fin d → Poly) : Bool :=
  (List.finRange d).all fun i => (List.finRange d).all fun j =>
    isZero (sub (matMulP Ψ Φ i j) (idP i j)) && isZero (sub (matMulP Φ Ψ i j) (idP i j))

def interOK (Φ P : Fin d → Fin d → Poly) : Bool :=
  (List.finRange d).all fun i => (List.finRange d).all fun j =>
    isZero (sub (matMulP Φ (constP P) i j) (matMulP P Φ i j))

variable {k : Type*} [CommRing k] (s : k)

lemma eval_eP (i : Fin d) : (fun j => eval s (eP i j)) = Pi.single i 1 := by
  funext j; by_cases h : j = i <;> simp [eP, h]

lemma evalM_constP (P : Fin d → Fin d → Poly) : evalM (constP P) s = evalM P 0 := by
  ext i j; simp only [evalM, constP]; cases P i j <;> simp

lemma evalM_matMulP (P Q : Fin d → Fin d → Poly) : evalM (matMulP P Q) s = evalM P s * evalM Q s := by
  ext i j; simp [evalM, matMulP, eval_sumF, Matrix.mul_apply]

lemma evalM_idP : evalM (idP (d := d)) s = 1 := by
  ext i j; by_cases h : i = j <;> simp [evalM, idP, h, Matrix.one_apply]

lemma evalM_of_check {P Q : Fin d → Fin d → Poly} (h : ∀ i j, isZero (sub (P i j) (Q i j)) = true) :
    evalM P s = evalM Q s := by
  ext i j; exact eval_eq_of_isZero_sub s _ _ (h i j)

/-- **Soundness.** -/
theorem orbit_sound (D : Fin d → Fin d → Fin d → ℤ) (Φ Ψ P : Fin d → Fin d → Poly)
    (hh : homOK c Φ = true) (hi : invOK Φ Ψ = true) (hc : interOK Φ P = true) :
    (∀ v w, evalM Φ s *ᵥ (CoordAlg.ofInt k c D).mul v w =
      (CoordAlg.ofInt k c D).mul (evalM Φ s *ᵥ v) (evalM Φ s *ᵥ w)) ∧
    evalM Ψ s * evalM Φ s = 1 ∧ evalM Φ s * evalM Ψ s = 1 ∧
    evalM Φ s * evalM P 0 = evalM P s * evalM Φ s := by
  simp only [homOK, homRow, invOK, interOK, Bool.and_eq_true, all_finRange] at hh hi hc
  refine ⟨?_, ?_, ?_, ?_⟩
  · -- on the basis, then by bilinearity
    set A := CoordAlg.ofInt k c D
    have hb : ∀ i, (Pi.basisFun k (Fin d)) i = Pi.single i 1 := fun i => by simp
    have base : ∀ i j, evalM Φ s *ᵥ A.mul (Pi.single i 1) (Pi.single j 1) =
        A.mul (evalM Φ s *ᵥ Pi.single i 1) (evalM Φ s *ᵥ Pi.single j 1) := by
      intro i j; funext m
      have := eval_eq_of_isZero_sub s _ _ (hh i j m)
      rw [mulVecP_eval, algMulP_eval (D := D)] at this
      rw [col_eval, col_eval]
      have e : (fun x => eval s (algMulP c (eP i) (eP j) x)) = A.mul (Pi.single i 1) (Pi.single j 1) := by
        funext x; rw [algMulP_eval (D := D), eval_eP, eval_eP]
      rw [← e]; exact this
    have := LinearMap.ext_basis (Pi.basisFun k (Fin d)) (Pi.basisFun k (Fin d))
      (B := A.mulₗ.compr₂ (evalM Φ s).mulVecLin)
      (B' := A.mulₗ.compl₁₂ (evalM Φ s).mulVecLin (evalM Φ s).mulVecLin) (fun i j => by simpa [hb] using base i j)
    intro v w; simpa using LinearMap.congr_fun₂ this v w
  · rw [← evalM_matMulP, evalM_of_check s (fun i j => (hi i j).1), evalM_idP]
  · rw [← evalM_matMulP, evalM_of_check s (fun i j => (hi i j).2), evalM_idP]
  · rw [← evalM_constP s P, ← evalM_matMulP, ← evalM_matMulP, evalM_of_check s (fun i j => hc i j)]

end EatTheFibre.Reflect
