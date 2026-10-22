/-
# Every automorphism of k(1⇉2⇉3) is filtered (v2.40)

Over a reduced commutative ring `k`, let `Φ` be multiplicative for the path algebra of 1⇉2⇉3 and have a two-sided
matrix inverse `Ψ`. Then `Φ` fixes the vertex rows (`ε ∘ Φ = ε`) and preserves `rad` and `rad²`, i.e. `Filt Φ`
(v2.39). Hence the behaviour-distinctness of the 19 families (`fam_behaviours_distinct`) holds for every algebra
automorphism, with no appeal to a generation statement for `Aut(kQ)`.

Proof:
- **`rad` is preserved.** A basis element `b` of `rad` has `b² = 0`, so the vertex coordinates of `Φ b` square to
  zero, and vanish because `k` is reduced.
- **`rad²` is preserved.** A path of length 2 is a product of two arrows, whose images lie in `rad`.
- **The off-diagonal end entries vanish.** The sink `e₂` annihilates `rad` from the left and the source `e₀` from the
  right. Since `Ψ` also preserves `rad`, every arrow is `Φ` of an element of `rad`, which gives
  `Φ 0 2 = Φ 1 2 = Φ 1 0 = Φ 2 0 = 0`.
- **The end vertices are fixed.** For `c = Φ 0 0`, the vector `(1 - c) • Φ e₀` is an idempotent lying in `rad`. Since
  `rad³ = 0` it is zero, so `c = 1`; the same holds for `Φ 2 2`.
- **The vertex block is the identity.** `Φ 1 = 1` gives the row sums, which finish it.
-/
import EatTheFibre.Gen.two_pairs_behaviour

open Matrix EatTheFibre.Reflect

namespace EatTheFibre.Q_two_pairs

variable {k : Type*} [CommRing k]

/-- `Φ` is multiplicative for the path algebra of 1⇉2⇉3 -/
def Mult (Φ : Matrix (Fin 11) (Fin 11) k) : Prop :=
  ∀ v w, Φ *ᵥ (alg k).mul v w = (alg k).mul (Φ *ᵥ v) (Φ *ᵥ w)

lemma Mult.inv {Φ Ψ : Matrix (Fin 11) (Fin 11) k} (hΦ : Mult Φ) (h1 : Φ * Ψ = 1) (h2 : Ψ * Φ = 1) : Mult Ψ := by
  intro v w
  have h := hΦ (Ψ *ᵥ v) (Ψ *ᵥ w)
  rw [Matrix.mulVec_mulVec, Matrix.mulVec_mulVec, h1, Matrix.one_mulVec, Matrix.one_mulVec] at h
  rw [← h, Matrix.mulVec_mulVec, h2, Matrix.one_mulVec]

lemma col (Φ : Matrix (Fin 11) (Fin 11) k) (j : Fin 11) : Φ *ᵥ Pi.single j 1 = fun i => Φ i j := by
  funext i; simp [Matrix.mulVec, dotProduct, Pi.single_apply]

lemma mul_single (i j : Fin 11) : (alg k).mul (Pi.single i 1) (Pi.single j 1) = fun m => (c i j m : k) := by
  funext m; simp [CoordAlg.mul, CoordAlg.ofInt, Pi.single_apply]

lemma Mult.basis {Φ : Matrix (Fin 11) (Fin 11) k} (hΦ : Mult Φ) (i j : Fin 11) :
    Φ *ᵥ (fun m => (c i j m : k)) = (alg k).mul (fun r => Φ r i) (fun r => Φ r j) := by
  rw [← mul_single, hΦ, col, col]

lemma mul_smul' (a : k) (v w : Fin 11 → k) : (alg k).mul (a • v) (a • w) = (a * a) • (alg k).mul v w := by
  rw [mul_explicit, mul_explicit]; funext m; fin_cases m <;> simp <;> ring

lemma rad_cube (x : Fin 11 → k) (h0 : x 0 = 0) (h1 : x 1 = 0) (h2 : x 2 = 0) :
    (alg k).mul x ((alg k).mul x x) = 0 := by
  rw [mul_explicit, mul_explicit]; funext m; fin_cases m <;> simp [h0, h1, h2]

/-- the unit `e₀ + e₁ + e₂` -/
def one11 : Fin 11 → k := ![1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0]

lemma one_mul11 (y : Fin 11 → k) : (alg k).mul one11 y = y := by
  rw [mul_explicit]; funext m; fin_cases m <;> simp [one11]

lemma mul_one11 (y : Fin 11 → k) : (alg k).mul y one11 = y := by
  rw [mul_explicit]; funext m; fin_cases m <;> simp [one11]

section reduced

variable [IsReduced k]

/-- `Φ` maps `rad` into `rad`: vertex coordinates of images of arrows and paths vanish -/
lemma Mult.vert_rad {Φ : Matrix (Fin 11) (Fin 11) k} (hΦ : Mult Φ) (i j : Fin 11) (hi : i.val < 3)
    (hj : 3 ≤ j.val) : Φ i j = 0 := by
  have h := congrFun (hΦ.basis j j) i
  have h0 : (fun m => (c j j m : k)) = 0 := by
    funext m; fin_cases j <;> simp at hj <;> fin_cases m <;> simp [c]
  rw [h0, Matrix.mulVec_zero, mul_explicit] at h
  have hsq : Φ i j ^ 2 = 0 := by
    fin_cases i <;> simp at hi <;> simp [pow_two] at h ⊢ <;> exact h.symm
  exact IsReduced.eq_zero _ ⟨2, hsq⟩

/-- `Φ` maps `rad²` into `rad²` -/
lemma Mult.rad2 {Φ : Matrix (Fin 11) (Fin 11) k} (hΦ : Mult Φ) (i j : Fin 11) (hi : len i = 1)
    (hj : len j = 2) : Φ i j = 0 := by
  have z : ∀ r s : Fin 11, r.val < 3 → 3 ≤ s.val → Φ r s = 0 := fun r s hr hs => hΦ.vert_rad r s hr hs
  fin_cases j <;> simp [len] at hj
  · have h := congrFun (hΦ.basis 3 5) i
    fin_cases i <;> simp [len] at hi <;> simp [Matrix.mulVec, dotProduct, Fin.sum_univ_succ, c, mul_explicit] at h <;>
      exact h.trans (by simp [z])
  · have h := congrFun (hΦ.basis 3 6) i
    fin_cases i <;> simp [len] at hi <;> simp [Matrix.mulVec, dotProduct, Fin.sum_univ_succ, c, mul_explicit] at h <;>
      exact h.trans (by simp [z])
  · have h := congrFun (hΦ.basis 4 5) i
    fin_cases i <;> simp [len] at hi <;> simp [Matrix.mulVec, dotProduct, Fin.sum_univ_succ, c, mul_explicit] at h <;>
      exact h.trans (by simp [z])
  · have h := congrFun (hΦ.basis 4 6) i
    fin_cases i <;> simp [len] at hi <;> simp [Matrix.mulVec, dotProduct, Fin.sum_univ_succ, c, mul_explicit] at h <;>
      exact h.trans (by simp [z])

set_option maxHeartbeats 1000000 in
/-- **Every algebra automorphism of k(1⇉2⇉3) is filtered**, over a reduced commutative ring. -/
theorem Mult.filt {Φ Ψ : Matrix (Fin 11) (Fin 11) k} (hΦ : Mult Φ) (h1 : Φ * Ψ = 1) (h2 : Ψ * Φ = 1) :
    Filt Φ := by
  have hΨ := hΦ.inv h1 h2
  have z : ∀ r s : Fin 11, r.val < 3 → 3 ≤ s.val → Φ r s = 0 := fun r s hr hs => hΦ.vert_rad r s hr hs
  have zΨ : ∀ s : Fin 11, 3 ≤ s.val → ∀ r : Fin 11, r.val < 3 → (Ψ *ᵥ Pi.single s 1) r = 0 := by
    intro s hs r hr; rw [col]; exact hΨ.vert_rad r s hr hs
  -- Φ 1 = 1, row sums
  have hone : Φ *ᵥ one11 = one11 := by
    have h := hΦ one11 (Ψ *ᵥ one11)
    rw [one_mul11, Matrix.mulVec_mulVec, h1, Matrix.one_mulVec, mul_one11] at h
    exact h.symm
  have r0 : Φ 0 0 + Φ 0 1 + Φ 0 2 = 1 := by
    have h := congrFun hone 0
    simp [Matrix.mulVec, dotProduct, Fin.sum_univ_succ, one11] at h; linear_combination h
  have r1 : Φ 1 0 + Φ 1 1 + Φ 1 2 = 1 := by
    have h := congrFun hone 1
    simp [Matrix.mulVec, dotProduct, Fin.sum_univ_succ, one11] at h; linear_combination h
  have r2 : Φ 2 0 + Φ 2 1 + Φ 2 2 = 1 := by
    have h := congrFun hone 2
    simp [Matrix.mulVec, dotProduct, Fin.sum_univ_succ, one11] at h; linear_combination h
  -- the source e₀ annihilates rad from the right, the sink e₂ from the left
  have src : ∀ s : Fin 11, 3 ≤ s.val → (alg k).mul (Pi.single s 1) (fun r => Φ r 0) = 0 := by
    intro s hs
    have h := hΦ (Ψ *ᵥ Pi.single s 1) (Pi.single 0 1)
    rw [Matrix.mulVec_mulVec, h1, Matrix.one_mulVec, col Φ 0] at h
    rw [← h]
    have h0 := zΨ s hs 0 (by decide)
    have hm : (alg k).mul (Ψ *ᵥ Pi.single s 1) (Pi.single 0 1) = 0 := by
      rw [mul_explicit]; funext m; fin_cases m <;> simp [h0]
    rw [hm, Matrix.mulVec_zero]
  have snk : ∀ s : Fin 11, 3 ≤ s.val → (alg k).mul (fun r => Φ r 2) (Pi.single s 1) = 0 := by
    intro s hs
    have h := hΦ (Pi.single 2 1) (Ψ *ᵥ Pi.single s 1)
    rw [Matrix.mulVec_mulVec, h1, Matrix.one_mulVec, col Φ 2] at h
    rw [← h]
    have h2' := zΨ s hs 2 (by decide)
    have hm : (alg k).mul (Pi.single 2 1) (Ψ *ᵥ Pi.single s 1) = 0 := by
      rw [mul_explicit]; funext m; fin_cases m <;> simp [h2']
    rw [hm, Matrix.mulVec_zero]
  have e10 : Φ 1 0 = 0 := by
    have h := congrFun (src 3 (by decide)) 3; rw [mul_explicit] at h; simpa using h
  have e20 : Φ 2 0 = 0 := by
    have h := congrFun (src 5 (by decide)) 5; rw [mul_explicit] at h; simpa using h
  have e02 : Φ 0 2 = 0 := by
    have h := congrFun (snk 3 (by decide)) 3; rw [mul_explicit] at h; simpa using h
  have e12 : Φ 1 2 = 0 := by
    have h := congrFun (snk 5 (by decide)) 5; rw [mul_explicit] at h; simpa using h
  -- idempotents in rad vanish
  have kill : ∀ j : Fin 11, j.val < 3 → ∀ a : k,
      (a • fun r => Φ r j) 0 = 0 → (a • fun r => Φ r j) 1 = 0 → (a • fun r => Φ r j) 2 = 0 →
      (alg k).mul (a • fun r => Φ r j) (a • fun r => Φ r j) = (a • fun r => Φ r j) → a = 0 := by
    intro j hj a x0 x1 x2 hid
    have hx : (a • fun r => Φ r j) = 0 := by
      have h3 := rad_cube _ x0 x1 x2
      rw [hid] at h3; rw [hid] at h3; exact h3
    have hv : Φ *ᵥ (a • Pi.single j 1) = 0 := by rw [Matrix.mulVec_smul, col]; exact hx
    have := congrArg (Ψ *ᵥ ·) hv
    simp only [Matrix.mulVec_mulVec, h2, Matrix.one_mulVec, Matrix.mulVec_zero] at this
    have := congrFun this j
    simpa using this
  have idem0 : (alg k).mul (fun r => Φ r 0) (fun r => Φ r 0) = fun r => Φ r 0 := by
    rw [← hΦ.basis 0 0]; funext r; simp [Matrix.mulVec, dotProduct, Fin.sum_univ_succ, c]
  have idem2 : (alg k).mul (fun r => Φ r 2) (fun r => Φ r 2) = fun r => Φ r 2 := by
    rw [← hΦ.basis 2 2]; funext r; simp [Matrix.mulVec, dotProduct, Fin.sum_univ_succ, c]
  have hc0 : Φ 0 0 * Φ 0 0 = Φ 0 0 := by
    have := congrFun idem0 0; rw [mul_explicit] at this; exact this
  have hc2 : Φ 2 2 * Φ 2 2 = Φ 2 2 := by
    have := congrFun idem2 2; rw [mul_explicit] at this; exact this
  have e00 : Φ 0 0 = 1 := by
    have := kill 0 (by decide) (1 - Φ 0 0)
      (show (1 - Φ 0 0) * Φ 0 0 = 0 by linear_combination -hc0)
      (show (1 - Φ 0 0) * Φ 1 0 = 0 by rw [e10, mul_zero])
      (show (1 - Φ 0 0) * Φ 2 0 = 0 by rw [e20, mul_zero])
      (by rw [mul_smul', idem0]; congr 1; linear_combination hc0)
    linear_combination -this
  have e22 : Φ 2 2 = 1 := by
    have := kill 2 (by decide) (1 - Φ 2 2)
      (show (1 - Φ 2 2) * Φ 0 2 = 0 by rw [e02, mul_zero])
      (show (1 - Φ 2 2) * Φ 1 2 = 0 by rw [e12, mul_zero])
      (show (1 - Φ 2 2) * Φ 2 2 = 0 by linear_combination -hc2)
      (by rw [mul_smul', idem2]; congr 1; linear_combination hc2)
    linear_combination -this
  have e01 : Φ 0 1 = 0 := by linear_combination r0 - e00 - e02
  have e11 : Φ 1 1 = 1 := by linear_combination r1 - e10 - e12
  have e21 : Φ 2 1 = 0 := by linear_combination r2 - e20 - e22
  constructor
  · intro i j hi
    have hi' : i.val < 3 := by fin_cases i <;> simp_all [len]
    by_cases hj : 3 ≤ j.val
    · rw [z i j hi' hj, if_neg (fun e => by rw [e] at hi'; omega)]
    · have hj' : j.val < 3 := by omega
      fin_cases i <;> fin_cases j <;> simp at hi' hj' ⊢ <;>
        simp [e00, e01, e02, e10, e11, e12, e20, e21, e22]
  · intro i j hij
    have hl : ∀ r : Fin 11, len r = 0 ∨ len r = 1 ∨ len r = 2 := by intro r; fin_cases r <;> simp [len]
    rcases hl i with hi | hi | hi
    · have hi' : i.val < 3 := by fin_cases i <;> simp_all [len]
      have hj' : 3 ≤ j.val := by fin_cases j <;> simp_all [len]
      exact z i j hi' hj'
    · have hj2 : len j = 2 := by rcases hl j with h | h | h <;> omega
      exact hΦ.rad2 i j hi hj2
    · rcases hl j with h | h | h <;> omega

end reduced

/-- **The nineteen families of 1⇉2⇉3 are nineteen behaviours**, for every algebra automorphism, over a reduced
nontrivial commutative ring. -/
theorem fam_behaviours_distinct_aut [IsReduced k] [Nontrivial k] {Φ Ψ : Matrix (Fin 11) (Fin 11) k} (hΦ : Mult Φ)
    (h1 : Φ * Ψ = 1) (h2 : Ψ * Φ = 1) {t t' : k} (i j : Fin 19)
    (hc : Φ * evalM (famP i) t = evalM (famP j) t' * Φ) : i = j :=
  fam_behaviours_distinct (hΦ.filt h1 h2) ((hΦ.inv h1 h2).filt h2 h1) h1 h2 i j hc

end EatTheFibre.Q_two_pairs
