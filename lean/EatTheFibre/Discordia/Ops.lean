/-
# Discordianism v2: the operational, Nu and bound axioms, derived or model-checked

* **M2** (no opposites) holds in the free monoid.
* **N5** (a zero counit kills the coalgebra) holds for every Mathlib `Coalgebra`.
* **N6** (the chaos operation) is witnessed by `(ℤ/5, −)`.
* **β4** (equivariance) is an identity of abelian groups.
* **B4** is Eckmann–Hilton. It *forces* N6's "no unit": a non-associative chaos that interchanges with a unital
  order cannot be unital (`B4_no_neutral_chaos`).
* **B1** (duality): grids are countable, chaos has the cardinality of the continuum, and every countable set of
  chaoses (in particular the grid-named ones) is null under every non-degenerate i.i.d. Bernoulli product measure.
* **B3** (periods): the minimal-period distribution of the 256 binary words of length 8.
* **§9** worked examples: 1680 = 14 · 5!, 5 · 73 = 365, 73 ≡ 3 (mod 5).
-/
import Mathlib

namespace EatTheFibre.Discordia.Ops

open MeasureTheory

/-! ## M2: no opposites -/

theorem M2 {α : Type*} (x y : List α) (h : x ++ y = []) : x = [] ∧ y = [] := List.append_eq_nil_iff.mp h

/-! ## N5: chaos has a step map, not a counit -/

/-- **N5.** A coalgebra whose counit is zero is trivial. -/
theorem N5 {R C : Type*} [CommSemiring R] [AddCommMonoid C] [Module R C] [Coalgebra R C]
    (h : (Coalgebra.counit : C →ₗ[R] R) = 0) (x : C) : x = 0 := by
  have := Coalgebra.lTensor_counit_comul (R := R) x
  rw [h, LinearMap.lTensor_zero, LinearMap.zero_apply] at this
  simpa using (congrArg (TensorProduct.rid R C) this).symm

/-! ## N6: the chaos operation on ℤ/5 -/

def ν (a b : ZMod 5) : ZMod 5 := a - b

theorem N6_self (a : ZMod 5) : ν a a = 0 := sub_self a
theorem N6_nonassoc : ∃ a b c : ZMod 5, ν (ν a b) c ≠ ν a (ν b c) := by decide
theorem N6_no_unit : ¬ ∃ e : ZMod 5, ∀ a, ν e a = a ∧ ν a e = a := by decide

/-! ## β4: equivariance -/

theorem β4 {G : Type*} [AddCommGroup G] (a b c : G) : (a + b) - (a + c) = b - c := add_sub_add_left_eq_sub b c a

/-! ## B4: Eckmann–Hilton exclusion -/

/-- **B4.** Interchange plus units for both operations forces them to be equal. -/
theorem B4 {X : Type*} {μ ν : X → X → X} {e₁ e₂ : X} (h₁ : EckmannHilton.IsUnital μ e₁)
    (h₂ : EckmannHilton.IsUnital ν e₂) (interchange : ∀ a b c d, μ (ν a b) (ν c d) = ν (μ a c) (μ b d)) : μ = ν :=
  EckmannHilton.mul h₁ h₂ interchange

/-- **B4 ⇒ N6.** If a unital order interchanges with a non-associative chaos, the chaos has no unit. -/
theorem B4_no_neutral_chaos {X : Type*} {μ ν : X → X → X} {e₁ : X} (h₁ : EckmannHilton.IsUnital μ e₁)
    (interchange : ∀ a b c d, μ (ν a b) (ν c d) = ν (μ a c) (μ b d))
    (hna : ∃ a b c, ν (ν a b) c ≠ ν a (ν b c)) : ¬ ∃ e, EckmannHilton.IsUnital ν e := by
  rintro ⟨e₂, h₂⟩
  obtain ⟨a, b, c, h⟩ := hna
  exact h ((EckmannHilton.mul_assoc h₁ h₂ interchange).assoc a b c)

/-! ## B1: duality -/

/-- a hand-off grid `(u, v)` names the chaos `u · v^ω` -/
def handoff (g : List Bool × List Bool) : ℕ → Bool :=
  fun n => if h : n < g.1.length then g.1[n] else
    if hv : g.2 = [] then false else g.2[(n - g.1.length) % g.2.length]'(Nat.mod_lt _ (List.length_pos_iff.mpr hv))

theorem B1_grids_countable : Countable (List Bool × List Bool) := inferInstance

theorem B1_chaos_continuum : Cardinal.mk (ℕ → Bool) = Cardinal.continuum := by
  simp

theorem B1_chaos_uncountable : ¬ Countable (ℕ → Bool) := by
  rw [← Cardinal.mk_le_aleph0_iff, B1_chaos_continuum, not_le]; exact Cardinal.aleph0_lt_continuum

theorem B1_named_countable : (Set.range handoff).Countable := Set.countable_range _

/-- singletons are null under a non-degenerate i.i.d. product measure -/
theorem B1_singleton_null (β : Measure Bool) [IsProbabilityMeasure β] (hβ : ∀ b, β {b} < 1) (x : ℕ → Bool) :
    Measure.infinitePi (fun _ : ℕ => β) {x} = 0 := by
  set q := max (β {true}) (β {false}) with hq
  have hq1 : q < 1 := max_lt (hβ true) (hβ false)
  have hb : ∀ b, β {b} ≤ q := fun b => by cases b <;> simp [q]
  have bound : ∀ n : ℕ, Measure.infinitePi (fun _ : ℕ => β) {x} ≤ q ^ n := fun n => by
    calc Measure.infinitePi (fun _ : ℕ => β) {x}
        ≤ Measure.infinitePi (fun _ : ℕ => β) (Set.pi (Finset.range n : Set ℕ) fun i => {x i}) :=
          measure_mono fun y hy => by simp_all
      _ = ∏ i ∈ Finset.range n, β {x i} := Measure.infinitePi_pi _ fun _ _ => measurableSet_singleton _
      _ ≤ ∏ _i ∈ Finset.range n, q := Finset.prod_le_prod' fun i _ => hb _
      _ = q ^ n := by simp
  exact le_antisymm (ge_of_tendsto' (ENNReal.tendsto_pow_atTop_nhds_zero_of_lt_one hq1) bound) zero_le

/-- **B1 (measure).** Order covers a null set: the grid-named chaoses are null under every non-degenerate i.i.d.
Bernoulli product measure. -/
theorem B1_named_null (β : Measure Bool) [IsProbabilityMeasure β] (hβ : ∀ b, β {b} < 1) :
    Measure.infinitePi (fun _ : ℕ => β) (Set.range handoff) = 0 := by
  have : NullSingletonClass (Measure.infinitePi fun _ : ℕ => β) := ⟨B1_singleton_null β hβ⟩
  exact B1_named_countable.measure_zero _

/-! ## B3: minimal periods at depth 8 -/

/-- the minimal period of a word (its length if it has no smaller one) -/
def minPeriod (w : List Bool) : ℕ :=
  ((List.range w.length).map (· + 1)).find? (fun p => (List.range (w.length - p)).all fun i => w[i]? == w[i + p]?)
    |>.getD w.length

/-- all binary words of length `n` -/
def words : ℕ → List (List Bool)
  | 0 => [[]]
  | n + 1 => (words n).flatMap fun w => [false :: w, true :: w]

/-- **B3.** Minimal periods of the 256 words of length 8: {1:2, 2:2, 3:6, 4:12, 5:30, 6:48, 7:82, 8:74}. -/
theorem B3_periods :
    (List.range 8).map (fun p => ((words 8).filter fun w => minPeriod w = p + 1).length) =
      [2, 2, 6, 12, 30, 48, 82, 74] := by
  decide +kernel

/-! ## §9: worked examples -/

/-- the five-fingered hand: 14 bracketings × 5! orders -/
theorem hand : catalan 4 * Nat.factorial 5 = 1680 := by
  rw [catalan_eq_centralBinom_div]; decide

theorem calendar : 5 * 73 = 365 ∧ 73 % 5 = 3 := by decide

end EatTheFibre.Discordia.Ops
