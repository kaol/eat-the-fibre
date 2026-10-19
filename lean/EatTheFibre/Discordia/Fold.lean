/-
# Discordianism v2: the algebraic layer, derived

The v2 axioms (`docs/axioms/discordianism-axioms-v2.md`) state M3, M4, δ-L, δ-D, E1 and E1⋆ as axioms checked on
models. Here they are **derived**.
- M4, the Rota–Baxter identity, δ-L and E1 follow from M3 alone: `f` is an idempotent ring endomorphism.
- The dendriform layer δ-D and the star splitting δ-D1 need only that `P` is a Rota–Baxter operator of weight −1.
- The annihilator law E1⋆ needs only that `P` is idempotent.

So δ-D, δ-D1 and E1⋆ also hold for every idempotent Rota–Baxter(−1) operator. By `Bridge.lean`, these are exactly the
folds of the research project.
-/
import Mathlib

namespace EatTheFibre.Discordia

variable {A : Type*} [Ring A]

/-! ## M3 ⇒ M4: the splitting -/

section Endo

variable (f : A →+* A) (hf : ∀ x, f (f x) = f x)
include hf

/-- **M4.** Hodge = im f is a subring, Podge = ker f is a two-sided ideal, and they are complementary. -/
theorem M4_split : (f.range.toAddSubgroup ⊓ (RingHom.ker f).toAddSubgroup = ⊥) ∧
    (f.range.toAddSubgroup ⊔ (RingHom.ker f).toAddSubgroup = ⊤) := by
  constructor
  · rw [eq_bot_iff]; rintro x ⟨⟨y, rfl⟩, hk⟩
    change f (f y) = 0 at hk
    rw [hf] at hk; simp [hk]
  · rw [eq_top_iff]; intro x _
    rw [show x = f x + (x - f x) by abel]
    exact AddSubgroup.add_mem_sup ⟨x, rfl⟩ (by change f (x - f x) = 0; rw [map_sub, hf, sub_self])

/-- **M4 (Rota–Baxter).** An idempotent endomorphism is a Rota–Baxter operator of weight −1. -/
theorem M4_rotaBaxter (x y : A) : f x * f y = f (f x * y + x * f y - x * y) := by
  simp only [map_sub, map_add, map_mul, hf]; abel

/-! ## δ-L: the Loday dialgebra x ⊣ y = x f(y), x ⊢ y = f(x) y -/

theorem dL_A1 (x y z : A) : (x * f y) * f z = x * f (y * f z) := by simp [map_mul, hf, mul_assoc]
theorem dL_A2 (x y z : A) : x * f (y * f z) = x * f (f y * z) := by simp [map_mul, hf]
theorem dL_A3 (x y z : A) : (f x * y) * f z = f x * (y * f z) := mul_assoc _ _ _
theorem dL_A4 (x y z : A) : f (x * f y) * z = f (f x * y) * z := by simp [map_mul, hf]
theorem dL_A5 (x y z : A) : f (f x * y) * z = f x * (f y * z) := by simp [map_mul, hf, mul_assoc]

/-- the Leibniz bracket `[x, y] = x ⊣ y − y ⊢ x = x f(y) − f(y) x` -/
def bracket (x y : A) : A := x * f y - f y * x

/-- **δ-L1.** The bracket satisfies the right Leibniz identity. -/
theorem dL1_leibniz (x y z : A) :
    bracket f (bracket f x y) z = bracket f (bracket f x z) y + bracket f x (bracket f y z) := by
  simp only [bracket, map_sub, map_mul, hf]; noncomm_ring

/-- **E1, right half.** Podge annihilates on the right: `[x, y] = 0` whenever `f y = 0`. -/
theorem E1_right (x y : A) (hy : f y = 0) : bracket f x y = 0 := by simp [bracket, hy]

/-- **E1, left half.** The left annihilator is the centralizer of Hodge. It is exactly Hodge iff Hodge is
self-centralizing, as in the free model k⟨a,b⟩, where im f = k[a]. -/
theorem E1_left (x : A) : (∀ y, bracket f x y = 0) ↔ ∀ z ∈ f.range, x * z = z * x := by
  constructor
  · rintro h z ⟨y, rfl⟩; have := h y; rwa [bracket, sub_eq_zero] at this
  · intro h y; rw [bracket, sub_eq_zero]; exact h _ ⟨y, rfl⟩

end Endo

/-! ## δ-D: dendriform structure from any Rota–Baxter operator of weight −1 -/

section RB

/-- `P` is a Rota–Baxter operator of weight −1 -/
def IsRB (P : A →+ A) : Prop := ∀ x y, P x * P y = P (P x * y + x * P y - x * y)

variable (P : A →+ A)

/-- `x ≺ y = x P(y)`, `x ≻ y = P(x) y − x y`, `x ⋆ y = x ≺ y + x ≻ y` (Aguiar) -/
def prec (x y : A) : A := x * P y
def succ (x y : A) : A := P x * y - x * y
def star (x y : A) : A := prec P x y + succ P x y

lemma star_eq (x y : A) : star P x y = P x * y + x * P y - x * y := by simp [star, prec, succ]; abel

/-- the star is mapped to the product: `P (x ⋆ y) = P x · P y` (this is the Rota–Baxter identity) -/
lemma P_star (hP : IsRB P) (x y : A) : P (star P x y) = P x * P y := by rw [star_eq, hP]

/-- **δ-D, three dendriform axioms.** -/
theorem dD_1 (hP : IsRB P) (x y z : A) : prec P (prec P x y) z = prec P x (star P y z) := by
  simp only [prec]; rw [P_star P hP, mul_assoc]
theorem dD_2 (x y z : A) : prec P (succ P x y) z = succ P x (prec P y z) := by
  simp only [prec, succ]; noncomm_ring
theorem dD_3 (hP : IsRB P) (x y z : A) : succ P x (succ P y z) = succ P (star P x y) z := by
  simp only [succ]; rw [P_star P hP, star_eq]; noncomm_ring

/-- **δ-D1, star splitting** (a ring identity, no hypothesis): with `p = P x`, `n = x − P x`, … ,
`x ⋆ y = p q − n m`. -/
theorem dD1_split (x y : A) : star P x y = P x * P y - (x - P x) * (y - P y) := by rw [star_eq]; noncomm_ring

/-- **δ-D1, the product decomposition.** `x ↦ (P x, P x − x)` turns `⋆` into the componentwise product.
For idempotent `P` it is a bijection onto Hodge × Podge, so (A, ⋆) ≅ im P × ker P. -/
theorem dD1_hom (hP : IsRB P) (x y : A) :
    P (star P x y) = P x * P y ∧ P (star P x y) - star P x y = (P x - x) * (P y - y) := by
  refine ⟨P_star P hP x y, ?_⟩
  rw [P_star P hP, dD1_split]; noncomm_ring

theorem dD1_inj (x y : A) (h1 : P x = P y) (h2 : P x - x = P y - y) : x = y := by
  rw [h1] at h2; exact (sub_right_inj.mp h2)

/-- **E1⋆, the annihilator law.** For idempotent `P` (no other hypothesis), `x` ⋆-commutes with everything iff its
Hodge part is central in Hodge and its Podge part is central in Podge. -/
theorem E1star (hidem : ∀ x, P (P x) = P x) (x : A) :
    (∀ y, star P x y = star P y x) ↔
      (∀ q, P q = q → P x * q = q * P x) ∧ (∀ m, P m = 0 → (x - P x) * m = m * (x - P x)) := by
  constructor
  · intro h
    refine ⟨fun q hq => ?_, fun m hm => ?_⟩
    · have := h q; rw [dD1_split, dD1_split, hq, sub_self, mul_zero, zero_mul, sub_zero, sub_zero] at this
      exact this
    · have := h m; rw [dD1_split, dD1_split, hm, mul_zero, zero_mul, sub_zero] at this
      simpa using this
  · rintro ⟨h1, h2⟩ y
    rw [dD1_split, dD1_split, h1 (P y) (hidem y), h2 (y - P y) (by rw [map_sub, hidem, sub_self])]

end RB

/-- An idempotent ring endomorphism is a Rota–Baxter(−1) operator, so δ-D, δ-D1 and E1⋆ hold for M3's `f`. -/
theorem M3_isRB (f : A →+* A) (hf : ∀ x, f (f x) = f x) : IsRB (f : A →+ A) := fun x y => M4_rotaBaxter f hf x y

end EatTheFibre.Discordia
