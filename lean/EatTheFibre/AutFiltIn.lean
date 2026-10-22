/-
# Automorphisms of k(1⇉2⇇3): filtered, or filtered after the source swap (v2.41)

Over an integral domain `k`, let `Φ` be multiplicative for the path algebra of 1⇉2⇇3 (sources `e₀`, `e₂`, sink `e₁`)
with a two-sided matrix inverse. Then either `Filt Φ` or `Filt (S * Φ)`, where `S` is the source swap
(`e₀ ↔ e₂`, `aᵢ ↔ bᵢ`).

Proof:
- `rad` is preserved: `b² = 0` and `k` is reduced.
- The sink annihilates `rad` from the left and the sources from the right. This gives `Φ₀₁ = Φ₂₁ = Φ₁₀ = Φ₁₂ = 0`,
  and `Φ 1 = 1` gives `Φ₁₁ = 1`.
- The remaining block `[[Φ₀₀, Φ₀₂], [Φ₂₀, Φ₂₂]]` has idempotent entries, row sums `1` and orthogonal columns.
- No column vanishes: an idempotent in `rad` is `0`, since `rad² = 0`.
- In a domain this leaves the identity or the swap.
-/
import EatTheFibre.Gen.pairs_in_behaviour

open Matrix EatTheFibre.Reflect

namespace EatTheFibre.Q_pairs_in

variable {k : Type*} [CommRing k]

/-- `Φ` is multiplicative for the path algebra of 1⇉2⇇3 -/
def Mult (Φ : Matrix (Fin 7) (Fin 7) k) : Prop :=
  ∀ v w, Φ *ᵥ (alg k).mul v w = (alg k).mul (Φ *ᵥ v) (Φ *ᵥ w)

lemma Mult.inv {Φ Ψ : Matrix (Fin 7) (Fin 7) k} (hΦ : Mult Φ) (h1 : Φ * Ψ = 1) (h2 : Ψ * Φ = 1) : Mult Ψ := by
  intro v w
  have h := hΦ (Ψ *ᵥ v) (Ψ *ᵥ w)
  rw [Matrix.mulVec_mulVec, Matrix.mulVec_mulVec, h1, Matrix.one_mulVec, Matrix.one_mulVec] at h
  rw [← h, Matrix.mulVec_mulVec, h2, Matrix.one_mulVec]

lemma Mult.mul {A B : Matrix (Fin 7) (Fin 7) k} (hA : Mult A) (hB : Mult B) : Mult (A * B) := by
  intro v w; rw [← Matrix.mulVec_mulVec, hB, hA, Matrix.mulVec_mulVec, Matrix.mulVec_mulVec]

lemma col (Φ : Matrix (Fin 7) (Fin 7) k) (j : Fin 7) : Φ *ᵥ Pi.single j 1 = fun i => Φ i j := by
  funext i; simp [Matrix.mulVec, dotProduct, Pi.single_apply]

lemma mul_single (i j : Fin 7) : (alg k).mul (Pi.single i 1) (Pi.single j 1) = fun m => (c i j m : k) := by
  funext m; simp [CoordAlg.mul, CoordAlg.ofInt, Pi.single_apply]

lemma Mult.basis {Φ : Matrix (Fin 7) (Fin 7) k} (hΦ : Mult Φ) (i j : Fin 7) :
    Φ *ᵥ (fun m => (c i j m : k)) = (alg k).mul (fun r => Φ r i) (fun r => Φ r j) := by
  rw [← mul_single, hΦ, col, col]

/-- the unit `e₀ + e₁ + e₂` -/
def one7 : Fin 7 → k := ![1, 1, 1, 0, 0, 0, 0]

lemma one_mul7 (y : Fin 7 → k) : (alg k).mul one7 y = y := by
  rw [mul_explicit]; funext m; fin_cases m <;> simp [one7]

lemma mul_one7 (y : Fin 7 → k) : (alg k).mul y one7 = y := by
  rw [mul_explicit]; funext m; fin_cases m <;> simp [one7]

/-- the source swap `e₀ ↔ e₂`, `a₁ ↔ b₁`, `a₂ ↔ b₂` -/
def S : Matrix (Fin 7) (Fin 7) k :=
  !![0, 0, 1, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 1, 0;
     0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 1, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0]

lemma S_mult : Mult (S : Matrix (Fin 7) (Fin 7) k) := by
  intro v w; rw [mul_explicit, mul_explicit]; funext m; fin_cases m <;>
    simp [S, Matrix.mulVec, dotProduct, Fin.sum_univ_succ] <;> ring

lemma S_sq : (S : Matrix (Fin 7) (Fin 7) k) * S = 1 := by
  ext i j; fin_cases i <;> fin_cases j <;> simp [S, Matrix.mul_apply, Fin.sum_univ_succ, Matrix.one_apply]

lemma S_not_filt [Nontrivial k] : ¬ Filt (S : Matrix (Fin 7) (Fin 7) k) := by
  intro h; have := h.vert 0 0 rfl; simp [S] at this

lemma S_mul_apply (Φ : Matrix (Fin 7) (Fin 7) k) (i j : Fin 7) :
    ((S : Matrix (Fin 7) (Fin 7) k) * Φ) i j = Φ ((![2, 1, 0, 5, 6, 3, 4] : Fin 7 → Fin 7) i) j := by
  fin_cases i <;> simp [S, Matrix.mul_apply, Fin.sum_univ_succ]

section domain

variable [IsDomain k]

lemma Mult.vert_rad {Φ : Matrix (Fin 7) (Fin 7) k} (hΦ : Mult Φ) (i j : Fin 7) (hi : i.val < 3)
    (hj : 3 ≤ j.val) : Φ i j = 0 := by
  have h := congrFun (hΦ.basis j j) i
  have h0 : (fun m => (c j j m : k)) = 0 := by
    funext m; fin_cases j <;> simp at hj <;> fin_cases m <;> simp [c]
  rw [h0, Matrix.mulVec_zero, mul_explicit] at h
  have hsq : Φ i j ^ 2 = 0 := by
    fin_cases i <;> simp at hi <;> simp [pow_two] at h ⊢ <;> first | exact h | exact h.symm
  exact IsReduced.eq_zero _ ⟨2, hsq⟩

lemma idem01 {x : k} (h : x * x = x) : x = 0 ∨ x = 1 := by
  have : x * (x - 1) = 0 := by linear_combination h
  rcases mul_eq_zero.mp this with h | h
  · exact Or.inl h
  · exact Or.inr (by linear_combination h)

set_option maxHeartbeats 1000000 in
/-- **Every algebra automorphism of k(1⇉2⇇3) is filtered, or filtered after the source swap** (k a domain). -/
theorem Mult.filt_or_swap {Φ Ψ : Matrix (Fin 7) (Fin 7) k} (hΦ : Mult Φ) (h1 : Φ * Ψ = 1) (h2 : Ψ * Φ = 1) :
    Filt Φ ∨ Filt (S * Φ) := by
  have hΨ := hΦ.inv h1 h2
  have z : ∀ r s : Fin 7, r.val < 3 → 3 ≤ s.val → Φ r s = 0 := fun r s hr hs => hΦ.vert_rad r s hr hs
  have zΨ : ∀ s : Fin 7, 3 ≤ s.val → ∀ r : Fin 7, r.val < 3 → (Ψ *ᵥ Pi.single s 1) r = 0 := by
    intro s hs r hr; rw [col]; exact hΨ.vert_rad r s hr hs
  have hone : Φ *ᵥ one7 = one7 := by
    have h := hΦ one7 (Ψ *ᵥ one7)
    rw [one_mul7, Matrix.mulVec_mulVec, h1, Matrix.one_mulVec, mul_one7] at h
    exact h.symm
  have r0 : Φ 0 0 + Φ 0 1 + Φ 0 2 = 1 := by
    have h := congrFun hone 0
    simp [Matrix.mulVec, dotProduct, Fin.sum_univ_succ, one7] at h; linear_combination h
  have r1 : Φ 1 0 + Φ 1 1 + Φ 1 2 = 1 := by
    have h := congrFun hone 1
    simp [Matrix.mulVec, dotProduct, Fin.sum_univ_succ, one7] at h; linear_combination h
  have r2 : Φ 2 0 + Φ 2 1 + Φ 2 2 = 1 := by
    have h := congrFun hone 2
    simp [Matrix.mulVec, dotProduct, Fin.sum_univ_succ, one7] at h; linear_combination h
  -- sink e₁ kills rad from the left; sources e₀, e₂ kill it from the right
  have snk : ∀ s : Fin 7, 3 ≤ s.val → (alg k).mul (fun r => Φ r 1) (Pi.single s 1) = 0 := by
    intro s hs
    have h := hΦ (Pi.single 1 1) (Ψ *ᵥ Pi.single s 1)
    rw [Matrix.mulVec_mulVec, h1, Matrix.one_mulVec, col Φ 1] at h
    rw [← h]
    have h' := zΨ s hs 1 (by decide)
    have hm : (alg k).mul (Pi.single 1 1) (Ψ *ᵥ Pi.single s 1) = 0 := by
      rw [mul_explicit]; funext m; fin_cases m <;> simp [h']
    rw [hm, Matrix.mulVec_zero]
  have src : ∀ j : Fin 7, j = 0 ∨ j = 2 → ∀ s : Fin 7, 3 ≤ s.val →
      (alg k).mul (Pi.single s 1) (fun r => Φ r j) = 0 := by
    intro j hj s hs
    have h := hΦ (Ψ *ᵥ Pi.single s 1) (Pi.single j 1)
    rw [Matrix.mulVec_mulVec, h1, Matrix.one_mulVec, col Φ j] at h
    rw [← h]
    have h0 := zΨ s hs 0 (by decide)
    have h2' := zΨ s hs 2 (by decide)
    have hm : (alg k).mul (Ψ *ᵥ Pi.single s 1) (Pi.single j 1) = 0 := by
      rcases hj with rfl | rfl <;> rw [mul_explicit] <;> funext m <;> fin_cases m <;> simp [h0, h2']
    rw [hm, Matrix.mulVec_zero]
  have e01 : Φ 0 1 = 0 := by
    have h := congrFun (snk 3 (by decide)) 3; rw [mul_explicit] at h; simpa using h
  have e21 : Φ 2 1 = 0 := by
    have h := congrFun (snk 5 (by decide)) 5; rw [mul_explicit] at h; simpa using h
  have e10 : Φ 1 0 = 0 := by
    have h := congrFun (src 0 (Or.inl rfl) 3 (by decide)) 3; rw [mul_explicit] at h; simpa using h
  have e12 : Φ 1 2 = 0 := by
    have h := congrFun (src 2 (Or.inr rfl) 3 (by decide)) 3; rw [mul_explicit] at h; simpa using h
  have e11 : Φ 1 1 = 1 := by linear_combination r1 - e10 - e12
  -- idempotent and orthogonal columns
  have idem0 : (alg k).mul (fun r => Φ r 0) (fun r => Φ r 0) = fun r => Φ r 0 := by
    rw [← hΦ.basis 0 0]; funext r; simp [Matrix.mulVec, dotProduct, Fin.sum_univ_succ, c]
  have idem2 : (alg k).mul (fun r => Φ r 2) (fun r => Φ r 2) = fun r => Φ r 2 := by
    rw [← hΦ.basis 2 2]; funext r; simp [Matrix.mulVec, dotProduct, Fin.sum_univ_succ, c]
  have orth : (alg k).mul (fun r => Φ r 0) (fun r => Φ r 2) = 0 := by
    rw [← hΦ.basis 0 2]; funext r; simp [Matrix.mulVec, dotProduct, Fin.sum_univ_succ, c]
  have i00 : Φ 0 0 * Φ 0 0 = Φ 0 0 := by have := congrFun idem0 0; rw [mul_explicit] at this; exact this
  have i22 : Φ 2 2 * Φ 2 2 = Φ 2 2 := by have := congrFun idem2 2; rw [mul_explicit] at this; exact this
  have o0 : Φ 0 0 * Φ 0 2 = 0 := by have := congrFun orth 0; rw [mul_explicit] at this; exact this
  have o2 : Φ 2 0 * Φ 2 2 = 0 := by have := congrFun orth 2; rw [mul_explicit] at this; exact this
  -- no column is an idempotent of rad
  have nd : ∀ j : Fin 7, (alg k).mul (fun r => Φ r j) (fun r => Φ r j) = (fun r => Φ r j) →
      Φ 0 j = 0 → Φ 1 j = 0 → Φ 2 j = 0 → False := by
    intro j hid a0 a1 a2
    have hx : (fun r => Φ r j) = 0 := by
      rw [← hid, mul_explicit]; funext m; fin_cases m <;> simp [a0, a1, a2]
    have hv : Φ *ᵥ Pi.single j 1 = 0 := by rw [col]; exact hx
    have := congrArg (Ψ *ᵥ ·) hv
    simp only [Matrix.mulVec_mulVec, h2, Matrix.one_mulVec, Matrix.mulVec_zero] at this
    have := congrFun this j
    simp at this
  have filtOf : ∀ A : Matrix (Fin 7) (Fin 7) k, (∀ i j : Fin 7, i.val < 3 → 3 ≤ j.val → A i j = 0) →
      A 0 0 = 1 → A 0 1 = 0 → A 0 2 = 0 → A 1 0 = 0 → A 1 1 = 1 → A 1 2 = 0 →
      A 2 0 = 0 → A 2 1 = 0 → A 2 2 = 1 → Filt A := by
    intro A zA a00 a01 a02 a10 a11 a12 a20 a21 a22
    constructor
    · intro i j hi
      have hi' : i.val < 3 := by fin_cases i <;> simp_all [len]
      by_cases hj : 3 ≤ j.val
      · rw [zA i j hi' hj, if_neg (fun e => by rw [e] at hi'; omega)]
      · have hj' : j.val < 3 := by omega
        fin_cases i <;> fin_cases j <;> simp at hi' hj' ⊢ <;>
          simp [a00, a01, a02, a10, a11, a12, a20, a21, a22]
    · intro i j hij
      have hi' : i.val < 3 := by fin_cases i <;> fin_cases j <;> simp_all [len]
      have hj' : 3 ≤ j.val := by fin_cases i <;> fin_cases j <;> simp_all [len]
      exact zA i j hi' hj'
  rcases idem01 i00 with h00 | h00
  · -- swap
    have e02 : Φ 0 2 = 1 := by linear_combination r0 - h00 - e01
    have e20 : Φ 2 0 = 1 := by
      rcases idem01 (show Φ 2 0 * Φ 2 0 = Φ 2 0 by have := congrFun idem0 2; rw [mul_explicit] at this; exact this)
        with h | h
      · exact (nd 0 idem0 h00 e10 h).elim
      · exact h
    have e22 : Φ 2 2 = 0 := by linear_combination r2 - e20 - e21
    right
    apply filtOf
    · intro i j hi hj; rw [S_mul_apply]; fin_cases i <;> simp at hi <;> exact z _ j (by decide) hj
    all_goals simp [S_mul_apply, e01, e02, e10, e11, e12, e20, e21, e22, h00]
  · -- identity
    have e02 : Φ 0 2 = 0 := by linear_combination r0 - h00 - e01
    have e22 : Φ 2 2 = 1 := by
      rcases idem01 i22 with h | h
      · exact (nd 2 idem2 e02 e12 h).elim
      · exact h
    have e20 : Φ 2 0 = 0 := by linear_combination r2 - e21 - e22
    left
    exact filtOf Φ z h00 e01 e02 e10 e11 e12 e20 e21 e22

/-- An automorphism and its inverse are both filtered, or both become filtered after the swap. -/
theorem Mult.cases {Φ Ψ : Matrix (Fin 7) (Fin 7) k} (hΦ : Mult Φ) (h1 : Φ * Ψ = 1) (h2 : Ψ * Φ = 1) :
    (Filt Φ ∧ Filt Ψ) ∨ (Filt (S * Φ) ∧ Filt (Ψ * S)) := by
  have hΨ := hΦ.inv h1 h2
  rcases hΦ.filt_or_swap h1 h2 with hF | hF
  · rcases hΨ.filt_or_swap h2 h1 with hG | hG
    · exact Or.inl ⟨hF, hG⟩
    · exfalso; apply S_not_filt (k := k)
      have := hG.mul hF
      rwa [Matrix.mul_assoc, h2, Matrix.mul_one] at this
  · have h1' : (S * Φ) * (Ψ * S) = 1 := by
      rw [Matrix.mul_assoc, ← Matrix.mul_assoc Φ, h1, Matrix.one_mul, S_sq]
    have h2' : (Ψ * S) * (S * Φ) = 1 := by
      rw [Matrix.mul_assoc, ← Matrix.mul_assoc S, S_sq, Matrix.one_mul, h2]
    rcases (hΨ.mul S_mult).filt_or_swap h2' h1' with hG | hG
    · exact Or.inr ⟨hF, hG⟩
    · exfalso; apply S_not_filt (k := k)
      have e : S * (Ψ * S) * (S * Φ) = S := by
        calc S * (Ψ * S) * (S * Φ) = S * (Ψ * (S * S) * Φ) := by simp only [Matrix.mul_assoc]
          _ = S := by rw [S_sq, Matrix.mul_one, h2, Matrix.mul_one]
      have := hG.mul hF
      rwa [e] at this

end domain

end EatTheFibre.Q_pairs_in
