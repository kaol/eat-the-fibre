/-
# δ-D2 exactly: which completions a given fold admits

`Completion.lean` proves that each stratum of dendriform completions is present under a class hypothesis. Here, in
any unital algebra, the hypotheses are shown to be **necessary** as well, fold by fold:

| completion                                | dendriform for the fold `P` iff                         |
|-------------------------------------------|----------------------------------------------------------|
| Baxter point                              | always                                                   |
| line 1, parameter `s`                     | `s · (x − Px) · P((y − Py) · Pz) = 0` (`line1_iff`)      |
| point 0                                   | `ker · im ⊆ ker` (`point0_iff`)                          |
| point B, or any single point of line 2    | `P` is multiplicative (`pointB_iff`, `line2_iff`)        |
| line 0, parameter `c`                     | `ker · im ⊆ ker` and `c · ker · im · ker = 0` (`line0_iff`) |
| line B, parameter `c`                     | `P` multiplicative and `c · ker · im · ker = 0` (`lineB_iff`) |

Also `mult_iff`: a fold is multiplicative iff its kernel is a two-sided module over its image (`LeftMod ∧ RightMod`).
So v2's δ-D2 family, the bimodule row of `Completion.lean`, is exactly the family of M3 projections.

Line 1 is the only stratum governed by a condition weaker than `LeftMod`. Its defect, in both D2 and D3, is the
single trilinear expression `s · (x − Px) · P((y − Py) · Pz)`.
-/
import EatTheFibre.Discordia.Completion

namespace EatTheFibre.Discordia.Completion

variable {R A : Type*} [CommRing R] [Ring A] [Algebra R A] {P : A →ₗ[R] A}

/-- `P` is multiplicative -/
def Mult (P : A →ₗ[R] A) : Prop := ∀ x y, P (x * y) = P x * P y

section
variable (h : Fold P)
include h

lemma P_PP (x y : A) : P (P x * P y) = P x * P y := P_ii h (h.idem x) (h.idem y)

/-- **A fold is multiplicative iff its kernel is an image-bimodule.** -/
theorem mult_iff : Mult P ↔ LeftMod P ∧ RightMod P := by
  constructor
  · intro hm
    refine ⟨fun x y => ?_, fun x y => ?_⟩
    · rw [hm, map_sub, h.idem, sub_self, zero_mul]
    · rw [hm, map_sub, h.idem y, sub_self, mul_zero]
  · rintro ⟨hL, hR⟩ x y
    obtain ⟨p₁, n₁, hp₁, hn₁, rfl⟩ := decomp h x
    obtain ⟨p₂, n₂, hp₂, hn₂, rfl⟩ := decomp h y
    simp only [mul_add, add_mul, map_add, hp₁, hn₁, hp₂, hn₂, P_ii h hp₁ hp₂, P_kk h hn₁ hn₂, P_ki hL hn₁ hp₂,
      P_ik hR hp₁ hn₂, add_zero]

/-- the decomposition preamble with the product rules of an arbitrary fold -/
macro "defect" h:term:max x:ident y:ident z:ident : tactic => `(tactic| (
  obtain ⟨p₁, n₁, hp₁, hn₁, rfl⟩ := decomp $h $x
  obtain ⟨p₂, n₂, hp₂, hn₂, rfl⟩ := decomp $h $y
  obtain ⟨p₃, n₃, hp₃, hn₃, rfl⟩ := decomp $h $z
  have ii := fun {a b : A} => P_ii $h (a := a) (b := b)
  have kk := fun {a b : A} => P_kk $h (a := a) (b := b)))

/-! ### Exact defects (any fold). `x − Px` is the Podge part, `Px` the Hodge part. -/

lemma point0_D1 (x y z : A) :
    prec P (prec P x y) z - prec P x (prec P y z + succ P 0 0 0 0 y z) = - (x * P ((y - P y) * P z)) := by
  defect h x y z; fold_simp; module

lemma pointB_D1 (x y z : A) :
    prec P (prec P x y) z - prec P x (prec P y z + succ P 0 1 (-1) 0 y z) = - (x * P (P y * (z - P z))) := by
  defect h x y z; fold_simp; module

lemma pointB_D2 (x y z : A) :
    prec P (succ P 0 1 (-1) 0 x y) z - succ P 0 1 (-1) 0 x (prec P y z) = x * P ((y - P y) * P z) := by
  defect h x y z; fold_simp; module

lemma line2_D1 (t : R) (x y z : A) :
    prec P (prec P x y) z - prec P x (prec P y z + succ P (-t) (1 + t) t (-1 - t) y z) =
      - (x * (P ((y - P y) * P z) + P (P y * (z - P z)))) := by
  defect h x y z; fold_simp; module

lemma line0_D1 (c : R) (x y z : A) :
    prec P (prec P x y) z - prec P x (prec P y z + succ P c (-c) (-c) c y z) = - (x * P ((y - P y) * P z)) := by
  defect h x y z; fold_simp; module

lemma line0_D2 (c : R) (x y z : A) :
    prec P (succ P c (-c) (-c) c x y) z - succ P c (-c) (-c) c x (prec P y z) =
      c • ((x - P x) * P ((y - P y) * P z)) := by
  defect h x y z; fold_simp; module

lemma line0_D3 (c : R) (x y z : A) :
    succ P c (-c) (-c) c x (succ P c (-c) (-c) c y z) - succ P c (-c) (-c) c (prec P x y + succ P c (-c) (-c) c x y) z =
      c • (P ((x - P x) * P y) * (z - P z)) - c • ((x - P x) * (P y * (z - P z))) := by
  defect h x y z; fold_simp; module

lemma lineB_D3 (hL : LeftMod P) (hR : RightMod P) (c : R) (x y z : A) :
    succ P c (1 - c) (-1 - c) c x (succ P c (1 - c) (-1 - c) c y z) -
      succ P c (1 - c) (-1 - c) c (prec P x y + succ P c (1 - c) (-1 - c) c x y) z =
      c • ((x - P x) * (P y * (z - P z))) := by
  defect h x y z
  have ki := fun {a b : A} => P_ki hL (a := a) (b := b)
  have ik := fun {a b : A} => P_ik hR (a := a) (b := b)
  fold_simp; module

/-- the line-1 defect: D1 holds, and D2, D3 fail by `−s · (x − Px) · P((y − Py) · Pz)` -/
lemma line1_defect (s : R) (x y z : A) :
    prec P (prec P x y) z = prec P x (prec P y z + succ P (-1 - s) (1 + s) s (-s) y z) ∧
    prec P (succ P (-1 - s) (1 + s) s (-s) x y) z - succ P (-1 - s) (1 + s) s (-s) x (prec P y z) =
      - s • ((x - P x) * P ((y - P y) * P z)) ∧
    succ P (-1 - s) (1 + s) s (-s) x (succ P (-1 - s) (1 + s) s (-s) y z) -
      succ P (-1 - s) (1 + s) s (-s) (prec P x y + succ P (-1 - s) (1 + s) s (-s) x y) z =
      - s • ((x - P x) * P ((y - P y) * P z)) := by
  defect h x y z
  refine ⟨?_, ?_, ?_⟩ <;> fold_simp <;> module

/-! ### The equivalences -/

/-- **Line 1** at parameter `s` is dendriform iff `s · (x − Px) · P((y − Py) · Pz) = 0`. -/
theorem line1_iff (s : R) : Dend (prec P) (succ P (-1 - s) (1 + s) s (-s)) ↔
    ∀ x y z : A, s • ((x - P x) * P ((y - P y) * P z)) = 0 := by
  constructor
  · intro hD x y z
    have := (line1_defect h s x y z).2.1
    rw [(hD x y z).2.1, sub_self, eq_comm, neg_smul, neg_eq_zero] at this; exact this
  · intro hN x y z
    obtain ⟨d1, d2, d3⟩ := line1_defect h s x y z
    rw [neg_smul, hN, neg_zero] at d2 d3
    exact ⟨d1, sub_eq_zero.mp d2, sub_eq_zero.mp d3⟩

/-- **Point 0** (`x ≻ y = 0`) is dendriform iff `ker · im ⊆ ker`. -/
theorem point0_iff : Dend (prec P) (succ P 0 0 0 0) ↔ LeftMod P := by
  refine ⟨fun hD x y => ?_, point0 h⟩
  have := point0_D1 h 1 x y
  rw [(hD 1 x y).1, sub_self, one_mul, eq_comm, neg_eq_zero] at this; exact this

/-- **Point B** (`x ≻ y = P(x)y − xP(y)`) is dendriform iff `P` is multiplicative. -/
theorem pointB_iff : Dend (prec P) (succ P 0 1 (-1) 0) ↔ Mult P := by
  refine ⟨fun hD => (mult_iff h).mpr ⟨fun x y => ?_, fun x y => ?_⟩, fun hm => ?_⟩
  · have := pointB_D2 h 1 x y
    rw [(hD 1 x y).2.1, sub_self, one_mul, eq_comm] at this; exact this
  · have := pointB_D1 h 1 x y
    rw [(hD 1 x y).1, sub_self, one_mul, eq_comm, neg_eq_zero] at this; exact this
  · obtain ⟨hL, hR⟩ := (mult_iff h).mp hm; exact pointB h hL hR

/-- **Any single point of line 2** is dendriform iff `P` is multiplicative. -/
theorem line2_iff (t : R) : Dend (prec P) (succ P (-t) (1 + t) t (-1 - t)) ↔ Mult P := by
  refine ⟨fun hD => (mult_iff h).mpr ⟨fun x y => ?_, fun x y => ?_⟩, fun hm => ?_⟩
  · have := line2_D1 h t 1 x (P y)
    rw [(hD 1 x (P y)).1, sub_self, one_mul, h.idem, sub_self, mul_zero, map_zero, add_zero, eq_comm,
      neg_eq_zero] at this
    exact this
  · have := line2_D1 h t 1 (P x) y
    rw [(hD 1 (P x) y).1, sub_self, one_mul, h.idem, sub_self, zero_mul, map_zero, zero_add, eq_comm,
      neg_eq_zero] at this
    exact this
  · obtain ⟨hL, hR⟩ := (mult_iff h).mp hm; exact line2 h hL hR t

/-- **Line 0** at parameter `c` (`x ≻ y = c (x − Px)(y − Py)`) is dendriform iff `ker · im ⊆ ker` and
`c · ker · im · ker = 0`. -/
theorem line0_iff (c : R) : Dend (prec P) (succ P c (-c) (-c) c) ↔
    LeftMod P ∧ ∀ x y z : A, c • ((x - P x) * (P y * (z - P z))) = 0 := by
  have L_of : Dend (prec P) (succ P c (-c) (-c) c) → LeftMod P := fun hD x y => by
    have := line0_D1 h c 1 x y
    rw [(hD 1 x y).1, sub_self, one_mul, eq_comm, neg_eq_zero] at this; exact this
  constructor
  · intro hD
    have hL := L_of hD
    refine ⟨hL, fun x y z => ?_⟩
    have := line0_D3 h c x y z
    rw [(hD x y z).2.2, sub_self, hL, zero_mul, smul_zero, zero_sub, eq_comm, neg_eq_zero] at this
    exact this
  · rintro ⟨hL, hK⟩ x y z
    have d1 := line0_D1 h c x y z
    have d2 := line0_D2 h c x y z
    have d3 := line0_D3 h c x y z
    rw [hL, mul_zero, neg_zero] at d1
    rw [hL, mul_zero, smul_zero] at d2
    rw [hL, zero_mul, smul_zero, hK, sub_zero] at d3
    exact ⟨sub_eq_zero.mp d1, sub_eq_zero.mp d2, sub_eq_zero.mp d3⟩

/-- **Line B as a whole** is dendriform iff `P` is multiplicative and `ker · im · ker = 0`. -/
theorem lineB_iff : (∀ c : R, Dend (prec P) (succ P c (1 - c) (-1 - c) c)) ↔ Mult P ∧ KIK P := by
  constructor
  · intro hD
    have hm : Mult P := (pointB_iff h).mp (by simpa using hD 0)
    obtain ⟨hL, hR⟩ := (mult_iff h).mp hm
    refine ⟨hm, fun x y z => ?_⟩
    have := lineB_D3 h hL hR 1 x y z
    rw [(hD 1 x y z).2.2, sub_self, one_smul, eq_comm] at this; exact this
  · rintro ⟨hm, hK⟩ c
    obtain ⟨hL, hR⟩ := (mult_iff h).mp hm
    exact lineB h hL hR hK c

end

end EatTheFibre.Discordia.Completion
