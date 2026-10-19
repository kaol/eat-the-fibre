/-
# Invariants: distinct parameters give distinct orbits

**Theorem 3 (1⇉2⇇3).** Every automorphism we use (substitutions in a parallel class, class weightings, and on the
quivers checked the identity component of `Aut_Δ`) fixes each vertex and maps paths of positive length into their
span. By `CoordAlg.vertex_block`, conjugation by such an `h` does not change the vertex-to-vertex entries of a fold,
and each family of Theorem 3 carries its parameter in one such entry. So two members with different parameters are
never conjugate.

**Theorem 4 (1⇉2⇉3).** The family is zero on the vertices, so its vertex block carries nothing. Its parameter is read
off the image instead. If `h` fixes the vertices and scales `δ_a`, `δ_b`, `δ_aδ_b` by `λ`, `μ`, `λμ` (as every
substitution and weighting does), and `h · P_t = P_t' · h`, then `t = t'`.
-/
import EatTheFibre.Gen.pairs_in
import EatTheFibre.Gen.two_pairs

open Matrix EatTheFibre.Reflect

namespace EatTheFibre

/-- the basis is ordered vertices first, so the vertices of a quiver with `n` vertices are the indices `< n` -/
def isVertex {d : ℕ} (n : ℕ) (i : Fin d) : Prop := i.val < n

instance {d n : ℕ} : DecidablePred (isVertex (d := d) n) := fun i => inferInstanceAs (Decidable (i.val < n))

/-- automorphisms that fix the vertices and preserve the radical, with inverse -/
structure VertexFixing {k : Type*} [CommRing k] {d : ℕ} (n : ℕ) (h hi : Matrix (Fin d) (Fin d) k) : Prop where
  inv : hi * h = 1
  fix : ∀ a, isVertex n a → ∀ i, h i a = if i = a then 1 else 0
  rad : ∀ a, ¬ isVertex n a → ∀ i, isVertex n i → h i a = 0

namespace Q_pairs_in

variable {k : Type*} [CommRing k]

/-- the vertex entry that carries the parameter, for a vertex-fixing conjugation -/
lemma entry_conj {P P' : Matrix (Fin 7) (Fin 7) k} {h hi : Matrix (Fin 7) (Fin 7) k} (hh : VertexFixing 3 h hi)
    (hc : h * P * hi = P') (i v : Fin 7) (hi_ : i.val < 3) (hv : v.val < 3) : P' i v = P i v := by
  rw [← hc]; exact CoordAlg.vertex_block (isVertex 3) hh.inv hh.fix hh.rad P hi_ hv

/-- **Theorem 3, F1:** members with different parameters are not conjugate by any vertex-fixing automorphism
(the `e₂`-coefficient of `P(e₁)` is `-s`). -/
theorem thm3F1_rigid {s s' : k} {h hi : Matrix (Fin 7) (Fin 7) k} (hh : VertexFixing 3 h hi)
    (hc : h * evalM thm3F1P s * hi = evalM thm3F1P s') : s = s' := by
  have := entry_conj hh hc 1 0 (by decide) (by decide)
  simp [evalM, thm3F1P] at this; exact this.symm

/-- **Theorem 3, F2** (the `e₁`-coefficient of `P(e₁)` is the parameter). -/
theorem thm3F2_rigid {s s' : k} {h hi : Matrix (Fin 7) (Fin 7) k} (hh : VertexFixing 3 h hi)
    (hc : h * evalM thm3F2P s * hi = evalM thm3F2P s') : s = s' := by
  have := entry_conj hh hc 0 0 (by decide) (by decide)
  simp [evalM, thm3F2P] at this; exact this.symm

/-- **Theorem 3, F4** (the `e₂`-coefficient of `P(e₁)` is the parameter). -/
theorem thm3F4_rigid {s s' : k} {h hi : Matrix (Fin 7) (Fin 7) k} (hh : VertexFixing 3 h hi)
    (hc : h * evalM thm3F4P s * hi = evalM thm3F4P s') : s = s' := by
  have := entry_conj hh hc 1 0 (by decide) (by decide)
  simp [evalM, thm3F4P] at this; exact this.symm

/-- **Theorem 3, F3** (the `e₂`-coefficient of `P(e₃)` is `-s`). -/
theorem thm3F3_rigid {s s' : k} {h hi : Matrix (Fin 7) (Fin 7) k} (hh : VertexFixing 3 h hi)
    (hc : h * evalM thm3F3P s * hi = evalM thm3F3P s') : s = s' := by
  have := entry_conj hh hc 1 2 (by decide) (by decide)
  simp [evalM, thm3F3P] at this; exact this.symm

/-- the automorphism of `k(1⇉2⇇3)` that substitutes the class `{a₁, a₂}` by `ga` and `{b₁, b₂}` by `gb`
(column `i` of `ga` is the image of `aᵢ`). Class substitutions (column sums 1) and class weightings (`w • 1`) are
of this form, and so are their products. -/
def subst (ga gb : Matrix (Fin 2) (Fin 2) k) : Matrix (Fin 7) (Fin 7) k :=
  !![1, 0, 0, 0, 0, 0, 0;
     0, 1, 0, 0, 0, 0, 0;
     0, 0, 1, 0, 0, 0, 0;
     0, 0, 0, ga 0 0, ga 0 1, 0, 0;
     0, 0, 0, ga 1 0, ga 1 1, 0, 0;
     0, 0, 0, 0, 0, gb 0 0, gb 0 1;
     0, 0, 0, 0, 0, gb 1 0, gb 1 1]

lemma subst_mul (ga gb ga' gb' : Matrix (Fin 2) (Fin 2) k) :
    subst ga' gb' * subst ga gb = subst (ga' * ga) (gb' * gb) := by
  ext i j; fin_cases i <;> fin_cases j <;> simp [subst, Matrix.mul_apply, Fin.sum_univ_succ]

lemma subst_one : subst (1 : Matrix (Fin 2) (Fin 2) k) 1 = 1 := by
  ext i j; fin_cases i <;> fin_cases j <;> simp [subst]

/-- every invertible class substitution fixes the vertices and preserves the radical -/
lemma subst_vertexFixing {ga gb ga' gb' : Matrix (Fin 2) (Fin 2) k} (ha : ga' * ga = 1) (hb : gb' * gb = 1) :
    VertexFixing 3 (subst ga gb) (subst ga' gb') where
  inv := by rw [subst_mul, ha, hb, subst_one]
  fix a ha_ i := by
    simp only [isVertex] at ha_
    fin_cases a <;> simp at ha_ <;> fin_cases i <;> simp [subst]
  rad a ha_ i hi := by
    simp only [isVertex, not_lt] at ha_ hi
    fin_cases a <;> simp at ha_ <;> fin_cases i <;> simp at hi <;> simp [subst]

/-- **Theorem 3, all four families, under `G × T`:** no invertible class substitution or weighting conjugates a
member to a member with a different parameter. -/
theorem thm3_rigid_GT {ga gb ga' gb' : Matrix (Fin 2) (Fin 2) k} (ha : ga' * ga = 1) (hb : gb' * gb = 1)
    {s s' : k} :
    (subst ga gb * evalM thm3F1P s * subst ga' gb' = evalM thm3F1P s' → s = s') ∧
    (subst ga gb * evalM thm3F2P s * subst ga' gb' = evalM thm3F2P s' → s = s') ∧
    (subst ga gb * evalM thm3F3P s * subst ga' gb' = evalM thm3F3P s' → s = s') ∧
    (subst ga gb * evalM thm3F4P s * subst ga' gb' = evalM thm3F4P s' → s = s') :=
  have hh := subst_vertexFixing ha hb
  ⟨thm3F1_rigid hh, thm3F2_rigid hh, thm3F3_rigid hh, thm3F4_rigid hh⟩

end Q_pairs_in

namespace Q_two_pairs

variable {k : Type*} [CommRing k]

/-! Basis of 1⇉2⇉3: 0 e1, 1 e2, 2 e3, 3 a1, 4 a2, 5 b1, 6 b2, 7 a1b1, 8 a1b2, 9 a2b1, 10 a2b2. -/

/-- `δ_a`, `δ_b`, `δ_aδ_b` as coordinate vectors -/
def δa : Fin 11 → k := ![0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0]
def δb : Fin 11 → k := ![0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0]
def δab : Fin 11 → k := ![0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 1]
def ve1 : Fin 11 → k := ![1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
def ve2 : Fin 11 → k := ![0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]

/-- how every class substitution and every class weighting acts on `V = span{e₁, e₂, δ_a, δ_b, δ_aδ_b}` -/
structure ScalesV (h : Matrix (Fin 11) (Fin 11) k) (l m : k) : Prop where
  e1 : h *ᵥ ve1 = ve1
  e2 : h *ᵥ ve2 = ve2
  da : h *ᵥ δa = l • δa
  db : h *ᵥ δb = m • δb
  dab : h *ᵥ δab = (l * m) • δab

/-- `P_t` applied to any vector: only the `a₂`- and `b₂`-coordinates matter -/
lemma thm4_apply (t : k) (w : Fin 11 → k) (i : Fin 11) :
    (evalM thm4P t *ᵥ w) i = w 4 * (![t - 1, 0, 0, -1, 1, 0, 0, t, -t, -t, t] : Fin 11 → k) i +
      w 6 * (![1, 1, 0, 0, 0, -1, 1, 1, -1, -1, 1] : Fin 11 → k) i := by
  fin_cases i <;> simp [evalM, thm4P, Matrix.mulVec, dotProduct, Fin.sum_univ_succ] <;> ring

/-- **Theorem 4, invariance.** If `h` acts on `V` by the scalings of a substitution or weighting and intertwines
`P_t` with `P_t'`, then `t = t'`. (It also forces `λ = μ = 1` on the image.) -/
theorem thm4_rigid {t t' l m : k} {h : Matrix (Fin 11) (Fin 11) k} (hV : ScalesV h l m)
    (hc : h * evalM thm4P t = evalM thm4P t' * h) : t = t' := by
  -- the images of b₂ and a₂ under P_t, written in V
  have Yv : evalM thm4P t *ᵥ Pi.single 6 1 = ve1 + ve2 - δb + δab := by
    rw [thm4_b2]; funext i; fin_cases i <;> simp [ve1, ve2, δb, δab]
  have Xv : evalM thm4P t *ᵥ Pi.single 4 1 = (t - 1) • ve1 - δa + t • δab := by
    rw [thm4_a2]; funext i; fin_cases i <;> simp [ve1, δa, δab] <;> ring
  have hY := congrArg (· *ᵥ Pi.single (6 : Fin 11) (1 : k)) hc
  have hX := congrArg (· *ᵥ Pi.single (4 : Fin 11) (1 : k)) hc
  simp only [← Matrix.mulVec_mulVec] at hY hX
  rw [Yv, Matrix.mulVec_add, Matrix.mulVec_sub, Matrix.mulVec_add, hV.e1, hV.e2, hV.db, hV.dab] at hY
  rw [Xv, Matrix.mulVec_add, Matrix.mulVec_sub, Matrix.mulVec_smul, Matrix.mulVec_smul, hV.e1, hV.da, hV.dab] at hX
  set w := h *ᵥ Pi.single (6 : Fin 11) (1 : k)
  set u := h *ᵥ Pi.single (4 : Fin 11) (1 : k)
  have y1 := congrFun hY 1; have y3 := congrFun hY 3; have y5 := congrFun hY 5; have y7 := congrFun hY 7
  have x1 := congrFun hX 1; have x3 := congrFun hX 3; have x7 := congrFun hX 7
  simp only [thm4_apply] at y1 y3 y5 y7 x1 x3 x7
  simp [ve1, ve2, δa, δb, δab] at y1 y3 y5 y7 x1 x3 x7
  -- y1: 1 = w 6;  y3: 0 = -w 4;  y5: -m = -w 6;  y7: l*m = t'*w 4 + w 6
  -- x1: 0 = u 6;  x3: -l = -u 4;  x7: t*(l*m) = t'*u 4 + u 6
  have hw6 : w 6 = 1 := by linear_combination -y1
  have hw4 : w 4 = 0 := by linear_combination y3
  have hm : m = 1 := by linear_combination y5 + hw6
  have hl : l = 1 := by rw [hm, hw4, hw6] at y7; linear_combination y7
  have hu6 : u 6 = 0 := by linear_combination -x1
  have hu4 : u 4 = 1 := by rw [hl] at x3; linear_combination -x3
  rw [hl, hm, hu4, hu6] at x7; linear_combination x7

/-- the automorphism of `k(1⇉2⇉3)` substituting `{a₁, a₂}` by `ga` and `{b₁, b₂}` by `gb`; on paths
`a_i b_j ↦ Σ ga p i · gb q j · a_p b_q` -/
def subst (ga gb : Matrix (Fin 2) (Fin 2) k) : Matrix (Fin 11) (Fin 11) k :=
  !![1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, ga 0 0, ga 0 1, 0, 0, 0, 0, 0, 0;
     0, 0, 0, ga 1 0, ga 1 1, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, gb 0 0, gb 0 1, 0, 0, 0, 0;
     0, 0, 0, 0, 0, gb 1 0, gb 1 1, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0, 0, ga 0 0 * gb 0 0, ga 0 0 * gb 0 1, ga 0 1 * gb 0 0, ga 0 1 * gb 0 1;
     0, 0, 0, 0, 0, 0, 0, ga 0 0 * gb 1 0, ga 0 0 * gb 1 1, ga 0 1 * gb 1 0, ga 0 1 * gb 1 1;
     0, 0, 0, 0, 0, 0, 0, ga 1 0 * gb 0 0, ga 1 0 * gb 0 1, ga 1 1 * gb 0 0, ga 1 1 * gb 0 1;
     0, 0, 0, 0, 0, 0, 0, ga 1 0 * gb 1 0, ga 1 0 * gb 1 1, ga 1 1 * gb 1 0, ga 1 1 * gb 1 1]

/-- a class substitution (column sums 1), scaled by a class weight, has `(1, -1)` as an eigenvector with eigenvalue
`w · (g₀₀ - g₀₁)` -/
lemma colsum_eigen (g : Matrix (Fin 2) (Fin 2) k) (w : k) (h0 : g 0 0 + g 1 0 = 1) (h1 : g 0 1 + g 1 1 = 1) :
    (w • g) *ᵥ ![1, -1] = (w * (g 0 0 - g 0 1)) • ![1, -1] := by
  funext i; fin_cases i <;> simp [Matrix.mulVec, dotProduct, Fin.sum_univ_succ]
  · ring
  · linear_combination w * (h0 - h1)

/-- any `subst ga gb` with `(1, -1)` an eigenvector of `ga` and `gb` acts on `V` by `ScalesV` -/
lemma subst_scalesV {ga gb : Matrix (Fin 2) (Fin 2) k} {l m : k}
    (ha : ga *ᵥ ![1, -1] = l • ![1, -1]) (hb : gb *ᵥ ![1, -1] = m • ![1, -1]) : ScalesV (subst ga gb) l m := by
  have a0 := congrFun ha 0; have a1 := congrFun ha 1; have b0 := congrFun hb 0; have b1 := congrFun hb 1
  simp [Matrix.mulVec, dotProduct, Fin.sum_univ_succ] at a0 a1 b0 b1
  refine ⟨?_, ?_, ?_, ?_, ?_⟩ <;> funext i <;> fin_cases i <;>
    simp [subst, ve1, ve2, δa, δb, δab, Matrix.mulVec, dotProduct, Fin.sum_univ_succ]
  all_goals first
    | linear_combination a0 | linear_combination a1 | linear_combination b0 | linear_combination b1
    | linear_combination (gb 0 0 - gb 0 1) * a0 + l * b0
    | linear_combination (gb 1 0 - gb 1 1) * a0 + l * b1
    | linear_combination (gb 0 0 - gb 0 1) * a1 - l * b0
    | linear_combination (gb 1 0 - gb 1 1) * a1 - l * b1

/-- **Theorem 4 under `G × T`:** if a class substitution scaled by class weights intertwines `P_t` and `P_t'`
(in particular, if it conjugates one to the other), then `t = t'`. -/
theorem thm4_rigid_GT {ga gb : Matrix (Fin 2) (Fin 2) k} {wa wb t t' : k}
    (ha0 : ga 0 0 + ga 1 0 = 1) (ha1 : ga 0 1 + ga 1 1 = 1) (hb0 : gb 0 0 + gb 1 0 = 1) (hb1 : gb 0 1 + gb 1 1 = 1)
    (hc : subst (wa • ga) (wb • gb) * evalM thm4P t = evalM thm4P t' * subst (wa • ga) (wb • gb)) : t = t' :=
  thm4_rigid (subst_scalesV (colsum_eigen ga wa ha0 ha1) (colsum_eigen gb wb hb0 hb1)) hc

end Q_two_pairs

end EatTheFibre
