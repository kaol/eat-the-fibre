/-
# Dense univariate integer polynomials, computable, with a certified evaluation map

`Poly = List ℤ`, lowest degree first. `add`, `smul`, `mul` are structural recursions that the kernel can evaluate
(`decide +kernel`). `eval t : Poly → k` is a ring map into any commutative ring (`eval_add`, `eval_mul`, …), and
`isZero p = true` implies `eval t p = 0` for every `t`. That is all reflection needs: an identity checked once on
`Poly` holds in every commutative ring at every parameter value.
-/
import Mathlib

namespace EatTheFibre

abbrev Poly := List ℤ

namespace Poly

def add : Poly → Poly → Poly
  | [], q => q
  | p, [] => p
  | a :: p, b :: q => (a + b) :: add p q

def smul (a : ℤ) (p : Poly) : Poly := p.map (a * ·)

def neg (p : Poly) : Poly := smul (-1) p

def sub (p q : Poly) : Poly := add p (neg q)

def mul : Poly → Poly → Poly
  | [], _ => []
  | a :: p, q => add (smul a q) (0 :: mul p q)

def isZero (p : Poly) : Bool := p.all (· == 0)

section eval

variable {k : Type*} [CommRing k] (t : k)

def eval : Poly → k
  | [] => 0
  | a :: p => (a : k) + t * eval p

@[simp] lemma eval_nil : eval t [] = 0 := rfl
@[simp] lemma eval_cons (a : ℤ) (p : Poly) : eval t (a :: p) = (a : k) + t * eval t p := rfl

@[simp] lemma eval_add : ∀ p q : Poly, eval t (add p q) = eval t p + eval t q
  | [], q => by simp [add]
  | a :: p, [] => by simp [add]
  | a :: p, b :: q => by simp only [add, eval_cons, eval_add p q]; push_cast; ring

@[simp] lemma eval_smul (a : ℤ) : ∀ p : Poly, eval t (smul a p) = (a : k) * eval t p
  | [] => by simp [smul]
  | b :: p => by
    have := eval_smul a p
    simp only [smul, List.map_cons, eval_cons] at this ⊢; rw [this]; push_cast; ring

@[simp] lemma eval_neg (p : Poly) : eval t (neg p) = - eval t p := by simp [neg]

@[simp] lemma eval_sub (p q : Poly) : eval t (sub p q) = eval t p - eval t q := by
  simp [sub, sub_eq_add_neg]

@[simp] lemma eval_mul : ∀ p q : Poly, eval t (mul p q) = eval t p * eval t q
  | [], q => by simp [mul]
  | a :: p, q => by simp only [mul, eval_add, eval_smul, eval_cons, eval_mul p q]; push_cast; ring

lemma eval_eq_zero_of_isZero : ∀ p : Poly, isZero p = true → eval t p = 0
  | [], _ => rfl
  | a :: p, h => by
    simp only [isZero, List.all_cons, Bool.and_eq_true, beq_iff_eq] at h
    simp [h.1, eval_eq_zero_of_isZero p (by simpa [isZero] using h.2)]

lemma eval_eq_of_isZero_sub (p q : Poly) (h : isZero (sub p q) = true) : eval t p = eval t q := by
  have := eval_eq_zero_of_isZero t _ h; rw [eval_sub] at this; exact sub_eq_zero.mp this

end eval

/-- finite sums over `Fin d`, computably -/
def sumF {d : ℕ} (f : Fin d → Poly) : Poly := (List.finRange d).foldr (fun i acc => add (f i) acc) []

lemma eval_sumF {k : Type*} [CommRing k] (t : k) {d : ℕ} (f : Fin d → Poly) :
    eval t (sumF f) = ∑ i, eval t (f i) := by
  rw [Fin.sum_univ_def, sumF]
  induction (List.finRange d) with
  | nil => simp
  | cons a l ih => simp [ih]

end Poly

end EatTheFibre
