/-
# δ-D2 derived: the dendriform completions of a fold

The v2 axioms fix `x ≺ y = x P(y)` and let `x ≻ y = b₁ xy + b₂ P(x)y + b₃ xP(y) + b₄ P(x)P(y)` range over k⁴.
Axiom δ-D2 claims the dendriform completions are two affine lines and two isolated points [R4]. That was computed in
the free model k⟨a,b⟩, for the endomorphism `f`.

Here `P` is any fold: an idempotent `R`-linear Rota–Baxter operator of weight −1. Then:

| hypothesis on the fold                      | completions proved dendriform               |
|---------------------------------------------|---------------------------------------------|
| none                                        | the Baxter point `b = (−1, 1, 0, 0)`        |
| `ker·im ⊆ ker` (`LeftMod`)                  | line 1 `(−1−s, 1+s, s, −s)`, point `0`      |
| `ker` an `im`-bimodule (e.g. endomorphisms) | line 1, line 2 `(−t, 1+t, t, −1−t)`, points `0`, `(0, 1, −1, 0)` |
| … and `ker·im·ker = 0` (`KIK`)              | the two points thicken to lines `(c, −c, −c, c)`, `(c, 1−c, −1−c, c)` |

The last row is v2's δ-D2, derived for every idempotent endomorphism instead of checked on one model. The converse
(nothing else is a completion, generically) is the computation `analysis/completion_generic.py`. The census folds sit
in different rows: Theorem 4's fold and three of Theorem 3's families are in the first row; Theorem 3's F4 is in the
second; Theorem 2 (an endomorphism, on a radical-cube-zero algebra) is in the last.

Proof method: write every variable as `p + n` with `P p = p`, `P n = 0`. Then `P` of any product of two pieces is
known from the class hypotheses, and `simp` followed by `module` closes each axiom.
-/
import Mathlib

namespace EatTheFibre.Discordia.Completion

variable {R A : Type*} [CommRing R] [Ring A] [Algebra R A] (P : A →ₗ[R] A)

/-- a fold: idempotent and Rota–Baxter of weight −1 (equivalently, by `Bridge.isFold_iff_rb`, image and kernel are
subalgebras) -/
structure Fold : Prop where
  idem : ∀ x, P (P x) = P x
  rb : ∀ x y, P x * P y = P (P x * y + x * P y - x * y)

/-- the kernel is a right module over the image -/
def LeftMod : Prop := ∀ x y, P ((x - P x) * P y) = 0
/-- the kernel is a left module over the image -/
def RightMod : Prop := ∀ x y, P (P x * (y - P y)) = 0
/-- `ker · im · ker = 0` -/
def KIK : Prop := ∀ x y z, (x - P x) * (P y * (z - P z)) = 0

def prec (x y : A) : A := x * P y
def succ (b₁ b₂ b₃ b₄ : R) (x y : A) : A := b₁ • (x * y) + b₂ • (P x * y) + b₃ • (x * P y) + b₄ • (P x * P y)

/-- the three dendriform axioms, with `x ⋆ y = x ≺ y + x ≻ y` -/
def Dend (pr sc : A → A → A) : Prop :=
  ∀ x y z, pr (pr x y) z = pr x (pr y z + sc y z) ∧ pr (sc x y) z = sc x (pr y z) ∧
    sc x (sc y z) = sc (pr x y + sc x y) z

variable {P}

section decomposition
variable (h : Fold P)
include h

lemma decomp (x : A) : ∃ p n, P p = p ∧ P n = 0 ∧ x = p + n :=
  ⟨P x, x - P x, h.idem x, by rw [map_sub, h.idem, sub_self], by abel⟩

lemma P_ii {a b : A} (ha : P a = a) (hb : P b = b) : P (a * b) = a * b := by
  have := h.rb a b; rw [ha, hb, add_sub_cancel_right] at this; exact this.symm

lemma P_kk {a b : A} (ha : P a = 0) (hb : P b = 0) : P (a * b) = 0 := by
  have := h.rb a b; rw [ha, hb, zero_mul, zero_mul, mul_zero, zero_add, zero_sub, map_neg] at this
  exact neg_eq_zero.mp this.symm

end decomposition

lemma P_ki (hL : LeftMod P) {a b : A} (ha : P a = 0) (hb : P b = b) : P (a * b) = 0 := by
  have := hL a b; rwa [ha, sub_zero, hb] at this

lemma P_ik (hR : RightMod P) {a b : A} (ha : P a = a) (hb : P b = 0) : P (a * b) = 0 := by
  have := hR a b; rwa [ha, hb, sub_zero] at this

lemma kik (hK : KIK P) {a b c : A} (ha : P a = 0) (hb : P b = b) (hc : P c = 0) : a * (b * c) = 0 := by
  have := hK a b c; rwa [ha, hb, hc, sub_zero, sub_zero] at this


/-- the simp set: distribute, pull scalars out, right-associate -/
macro "fold_simp" : tactic => `(tactic|
  simp only [prec, succ, map_add, map_sub, map_smul, map_neg, mul_add, add_mul, mul_sub, sub_mul, smul_add,
    smul_sub, smul_mul_assoc, mul_smul_comm, smul_smul, mul_assoc, neg_mul, mul_neg, smul_neg, map_zero,
    mul_zero, zero_mul, smul_zero, add_zero, zero_add, sub_zero, zero_sub, neg_zero, *] at *)

/-- **The Baxter point** `(−1, 1, 0, 0)`: every fold. -/
theorem baxter (h : Fold P) : Dend (prec P) (succ P (-1) 1 0 0) := by
  intro x y z
  obtain ⟨p₁, n₁, hp₁, hn₁, rfl⟩ := decomp h x
  obtain ⟨p₂, n₂, hp₂, hn₂, rfl⟩ := decomp h y
  obtain ⟨p₃, n₃, hp₃, hn₃, rfl⟩ := decomp h z
  have ii := fun {a b : A} => P_ii h (a := a) (b := b)
  have kk := fun {a b : A} => P_kk h (a := a) (b := b)
  refine ⟨?_, ?_, ?_⟩ <;> fold_simp <;> module

/-- the decomposition preamble: split `x y z`, collect the available product rules -/
macro "decomp3" h:term : tactic => `(tactic| (
  intro x y z
  obtain ⟨p₁, n₁, hp₁, hn₁, rfl⟩ := decomp $h x
  obtain ⟨p₂, n₂, hp₂, hn₂, rfl⟩ := decomp $h y
  obtain ⟨p₃, n₃, hp₃, hn₃, rfl⟩ := decomp $h z))

/-- **Line 1** `(−1−s, 1+s, s, −s)` and **the zero point**: folds whose kernel is a right module over the image. -/
theorem line1 (h : Fold P) (hL : LeftMod P) (s : R) : Dend (prec P) (succ P (-1 - s) (1 + s) s (-s)) := by
  decomp3 h
  have ii := fun {a b : A} => P_ii h (a := a) (b := b)
  have kk := fun {a b : A} => P_kk h (a := a) (b := b)
  have ki := fun {a b : A} => P_ki hL (a := a) (b := b)
  refine ⟨?_, ?_, ?_⟩ <;> fold_simp <;> module

theorem point0 (h : Fold P) (hL : LeftMod P) : Dend (prec P) (succ P 0 0 0 0) := by
  decomp3 h
  have ii := fun {a b : A} => P_ii h (a := a) (b := b)
  have kk := fun {a b : A} => P_kk h (a := a) (b := b)
  have ki := fun {a b : A} => P_ki hL (a := a) (b := b)
  refine ⟨?_, ?_, ?_⟩ <;> fold_simp <;> module

/-- **Line 2** `(−t, 1+t, t, −1−t)` and **the point** `(0, 1, −1, 0)`: folds whose kernel is an image-bimodule. -/
theorem line2 (h : Fold P) (hL : LeftMod P) (hR : RightMod P) (t : R) :
    Dend (prec P) (succ P (-t) (1 + t) t (-1 - t)) := by
  decomp3 h
  have ii := fun {a b : A} => P_ii h (a := a) (b := b)
  have kk := fun {a b : A} => P_kk h (a := a) (b := b)
  have ki := fun {a b : A} => P_ki hL (a := a) (b := b)
  have ik := fun {a b : A} => P_ik hR (a := a) (b := b)
  refine ⟨?_, ?_, ?_⟩ <;> fold_simp <;> module

theorem pointB (h : Fold P) (hL : LeftMod P) (hR : RightMod P) : Dend (prec P) (succ P 0 1 (-1) 0) := by
  decomp3 h
  have ii := fun {a b : A} => P_ii h (a := a) (b := b)
  have kk := fun {a b : A} => P_kk h (a := a) (b := b)
  have ki := fun {a b : A} => P_ki hL (a := a) (b := b)
  have ik := fun {a b : A} => P_ik hR (a := a) (b := b)
  refine ⟨?_, ?_, ?_⟩ <;> fold_simp <;> module

/-- **The points thicken** when `ker · im · ker = 0`: `(c, −c, −c, c)`, i.e. `x ≻ y = c (x − Px)(y − Py)`. -/
theorem line0 (h : Fold P) (hL : LeftMod P) (hK : KIK P) (c : R) : Dend (prec P) (succ P c (-c) (-c) c) := by
  decomp3 h
  have ii := fun {a b : A} => P_ii h (a := a) (b := b)
  have kk := fun {a b : A} => P_kk h (a := a) (b := b)
  have ki := fun {a b : A} => P_ki hL (a := a) (b := b)
  have k3 := fun {a b c : A} => kik hK (a := a) (b := b) (c := c)
  refine ⟨?_, ?_, ?_⟩ <;> fold_simp <;> module

theorem lineB (h : Fold P) (hL : LeftMod P) (hR : RightMod P) (hK : KIK P) (c : R) :
    Dend (prec P) (succ P c (1 - c) (-1 - c) c) := by
  decomp3 h
  have ii := fun {a b : A} => P_ii h (a := a) (b := b)
  have kk := fun {a b : A} => P_kk h (a := a) (b := b)
  have ki := fun {a b : A} => P_ki hL (a := a) (b := b)
  have ik := fun {a b : A} => P_ik hR (a := a) (b := b)
  have k3 := fun {a b c : A} => kik hK (a := a) (b := b) (c := c)
  refine ⟨?_, ?_, ?_⟩ <;> fold_simp <;> module

/-- **An idempotent algebra endomorphism is a fold whose kernel is an image-bimodule.** So v2's δ-D2 family
(lines 1, 2 and both points) is dendriform for every M3 projection, not just in the free model. -/
theorem endo_bimod (f : A →ₐ[R] A) (hf : ∀ x, f (f x) = f x) :
    Fold f.toLinearMap ∧ LeftMod f.toLinearMap ∧ RightMod f.toLinearMap := by
  refine ⟨⟨hf, fun x y => ?_⟩, fun x y => ?_, fun x y => ?_⟩ <;>
    simp only [AlgHom.toLinearMap_apply, map_sub, map_add, map_mul, hf, sub_self, mul_zero, zero_mul]
  abel

end EatTheFibre.Discordia.Completion
