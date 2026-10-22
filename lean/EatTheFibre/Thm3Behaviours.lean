/-
# Theorem 3's four families are exactly three behaviours (v2.41)

Over a field, for parameters `t, t' ∉ {0, -1}`: `F_i(t)` and `F_j(t')` are conjugate by an algebra automorphism of
k(1⇉2⇇3) iff `i = j` or `{i, j} = {F1, F3}`.

- `→` is `thm3_behaviours_aut`.
- `←` uses the single-orbit theorems for F2 and F4 (v2.33).
- For F1 it uses the witness of v2.33, an arrow scaling followed by `Inn(1 + x)`, rational in `t, t'` with poles at
  `0, -1`.
- For F1 ↔ F3 it uses the source swap `S` (`thm3_F1_F3`).
-/
import EatTheFibre.Gen.pairs_in_aut
import EatTheFibre.Gen.orbits_thm3F2
import EatTheFibre.Gen.orbits_thm3F4

open Matrix EatTheFibre.Reflect

namespace EatTheFibre.Q_pairs_in

variable {k : Type*}

section ring

variable [CommRing k]

/-- `A` and `B` are conjugate by an algebra automorphism of k(1⇉2⇇3) -/
def Conj (A B : Matrix (Fin 7) (Fin 7) k) : Prop :=
  ∃ Φ Ψ : Matrix (Fin 7) (Fin 7) k, Mult Φ ∧ Φ * Ψ = 1 ∧ Ψ * Φ = 1 ∧ Φ * A = B * Φ

lemma Conj.refl (A : Matrix (Fin 7) (Fin 7) k) : Conj A A :=
  ⟨1, 1, fun v w => by simp, Matrix.mul_one 1, Matrix.mul_one 1, by rw [Matrix.one_mul, Matrix.mul_one]⟩

lemma Conj.symm {A B : Matrix (Fin 7) (Fin 7) k} (h : Conj A B) : Conj B A := by
  obtain ⟨Φ, Ψ, hΦ, h1, h2, hc⟩ := h
  refine ⟨Ψ, Φ, hΦ.inv h1 h2, h2, h1, ?_⟩
  calc Ψ * B = Ψ * B * (Φ * Ψ) := by rw [h1, Matrix.mul_one]
    _ = Ψ * (B * Φ) * Ψ := by simp only [Matrix.mul_assoc]
    _ = Ψ * (Φ * A) * Ψ := by rw [hc]
    _ = (Ψ * Φ) * A * Ψ := by simp only [Matrix.mul_assoc]
    _ = A * Ψ := by rw [h2, Matrix.one_mul]

lemma Conj.trans {A B C : Matrix (Fin 7) (Fin 7) k} (h : Conj A B) (h' : Conj B C) : Conj A C := by
  obtain ⟨Φ, Ψ, hΦ, h1, h2, hc⟩ := h
  obtain ⟨Φ', Ψ', hΦ', h1', h2', hc'⟩ := h'
  refine ⟨Φ' * Φ, Ψ * Ψ', hΦ'.mul hΦ, ?_, ?_, ?_⟩
  · rw [Matrix.mul_assoc, ← Matrix.mul_assoc Φ, h1, Matrix.one_mul, h1']
  · rw [Matrix.mul_assoc, ← Matrix.mul_assoc Ψ', h2', Matrix.one_mul, h2]
  · rw [Matrix.mul_assoc, hc, ← Matrix.mul_assoc, hc', Matrix.mul_assoc]

lemma conj_F2 (s s' : k) : Conj (evalM thm3F2P s) (evalM thm3F2P s') := by
  have orb : ∀ s : k, Conj (evalM thm3F2P 0) (evalM thm3F2P s) := fun s =>
    let ⟨hm, h2, h1, hc⟩ := thm3F2_one_orbit (k := k) s
    ⟨_, _, hm, h1, h2, hc⟩
  exact (orb s).symm.trans (orb s')

lemma conj_F4 (s s' : k) : Conj (evalM thm3F4P s) (evalM thm3F4P s') := by
  have orb : ∀ s : k, Conj (evalM thm3F4P 0) (evalM thm3F4P s) := fun s =>
    let ⟨hm, h2, h1, hc⟩ := thm3F4_one_orbit (k := k) s
    ⟨_, _, hm, h1, h2, hc⟩
  exact (orb s).symm.trans (orb s')

lemma conj_F1_F3 (s : k) : Conj (evalM thm3F1P s) (evalM thm3F3P s) := by
  refine ⟨S, S, S_mult, S_sq, S_sq, ?_⟩
  rw [← thm3F1P_mirror s, Matrix.mul_assoc, Matrix.mul_assoc, S_sq, Matrix.mul_one]

lemma innM_mult (x : Fin 7 → k) : Mult (innM x) := by
  intro v w; rw [mul_explicit, mul_explicit]; funext m; fin_cases m <;>
    simp [innM, Matrix.mulVec, dotProduct, Fin.sum_univ_succ] <;> ring

lemma innM_neg (x : Fin 7 → k) : innM x * innM (-x) = 1 := by
  ext i j; fin_cases i <;> fin_cases j <;> simp [innM, Matrix.mul_apply, Fin.sum_univ_succ, Matrix.one_apply] <;> ring

lemma diag_mult (p q r s : k) : Mult (Matrix.diagonal ![1, 1, 1, p, q, r, s]) := by
  intro v w; rw [mul_explicit, mul_explicit]; funext m; fin_cases m <;>
    simp [Matrix.mulVec, dotProduct, Fin.sum_univ_succ, Matrix.diagonal] <;> ring

end ring

section field

variable [Field k]

/-- the v2.33 witness for F1: scale `a₂, b₁` by `t/t'` and `b₂` by `t(t+1)/(t'(t'+1))`, then `Inn(1 + x)` -/
def f1T (t t' : k) : Matrix (Fin 7) (Fin 7) k :=
  Matrix.diagonal ![1, 1, 1, 1, t / t', t / t', t * (t + 1) / (t' * (t' + 1))]
def f1Ti (t t' : k) : Matrix (Fin 7) (Fin 7) k :=
  Matrix.diagonal ![1, 1, 1, 1, t' / t, t' / t, t' * (t' + 1) / (t * (t + 1))]
def f1x (t t' : k) : Fin 7 → k := ![0, 0, 0, 0, (t' - t) / t', 0, (t' - t) / (t' * (t' + 1))]

set_option maxHeartbeats 2000000 in
lemma conj_F1 {t t' : k} (ht : t ≠ 0) (ht1 : t + 1 ≠ 0) (ht' : t' ≠ 0) (ht1' : t' + 1 ≠ 0) :
    Conj (evalM thm3F1P t) (evalM thm3F1P t') := by
  refine ⟨innM (f1x t t') * f1T t t', f1Ti t t' * innM (-f1x t t'), (innM_mult _).mul (diag_mult _ _ _ _), ?_, ?_, ?_⟩
  · have hd : f1T t t' * f1Ti t t' = 1 := by
      rw [f1T, f1Ti, Matrix.diagonal_mul_diagonal, ← Matrix.diagonal_one]; congr 1; funext i
      fin_cases i <;> simp <;> field_simp
    rw [Matrix.mul_assoc, ← Matrix.mul_assoc (f1T t t'), hd, Matrix.one_mul, innM_neg]
  · have hd : f1Ti t t' * f1T t t' = 1 := by
      rw [f1T, f1Ti, Matrix.diagonal_mul_diagonal, ← Matrix.diagonal_one]; congr 1; funext i
      fin_cases i <;> simp <;> field_simp
    have hn : innM (-f1x t t') * innM (f1x t t') = 1 := by simpa using innM_neg (-f1x t t')
    rw [Matrix.mul_assoc, ← Matrix.mul_assoc (innM _), hn, Matrix.one_mul, hd]
  · ext i j
    simp only [Matrix.mul_apply, Fin.sum_univ_succ, Fin.sum_univ_zero, f1T, Matrix.mul_diagonal]
    fin_cases i <;> fin_cases j <;> simp [innM, f1x, evalM, thm3F1P] <;> field_simp <;> ring

/-- **Theorem 3's four families are exactly three behaviours.** Over a field, for `t, t' ∉ {0, -1}`, `F_i(t)` and
`F_j(t')` are conjugate by an algebra automorphism iff `i = j` or `{i, j} = {F1, F3}` (indices 0 and 2). -/
theorem thm3_exactly_three {t t' : k} (ht : t ≠ 0) (ht1 : t + 1 ≠ 0) (ht' : t' ≠ 0) (ht1' : t' + 1 ≠ 0)
    (i j : Fin 4) :
    Conj (evalM (thm3FP i) t) (evalM (thm3FP j) t') ↔ (i = j ∨ (i = 0 ∧ j = 2) ∨ (i = 2 ∧ j = 0)) := by
  constructor
  · rintro ⟨Φ, Ψ, hΦ, h1, h2, hc⟩
    exact thm3_behaviours_aut hΦ h1 h2 ht ht' i j hc
  · intro h
    have F3 : ∀ s s' : k, s ≠ 0 → s + 1 ≠ 0 → s' ≠ 0 → s' + 1 ≠ 0 →
        Conj (evalM thm3F3P s) (evalM thm3F3P s') := fun s s' a b c d =>
      (conj_F1_F3 s).symm.trans ((conj_F1 a b c d).trans (conj_F1_F3 s'))
    fin_cases i <;> fin_cases j <;> simp at h <;> simp only [thm3FP, Fin.isValue, Matrix.cons_val, Fin.reduceFinMk]
    · exact conj_F1 ht ht1 ht' ht1'
    · exact (conj_F1 ht ht1 ht' ht1').trans (conj_F1_F3 t')
    · exact conj_F2 t t'
    · exact (conj_F1_F3 t).symm.trans (conj_F1 ht ht1 ht' ht1')
    · exact F3 t t' ht ht1 ht' ht1'
    · exact conj_F4 t t'

end field

end EatTheFibre.Q_pairs_in
