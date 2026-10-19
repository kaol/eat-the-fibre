/-
# Theorems 3 and 4 on every acyclic quiver

**Theorem 4 (v2.28).** Let `Q` be a finite quiver with an acyclicity potential and two consecutive parallel classes
`a₁, a₂ : u → v`, `b₁, b₂ : v → w` (`a₁ ≠ a₂`, `b₁ ≠ b₂`). For every commutative ring `k` and every `t : k`, the
normal form `P_t` of 1⇉2⇉3, pushed into `k Q` (zero on every other path), is a fold of `k Q`. No condition on the outer
ends is needed.

**Theorem 3 (v2.27).** Same with converging classes `a₁, a₂ : u → v`, `b₁, b₂ : w → v`, `u ≠ w`, and no arrow
with both ends in `{u, w}`. All four families of 1⇉2⇇3 push to folds of `k Q`.

Proof: the census tables of 1⇉2⇉3 (resp. 1⇉2⇇3) are the path algebras of those quivers (`decide +kernel`).
The classes give an injective quiver morphism into `Q`, hence an `Embeds` (`Hom.embeds`). `isFold_push` reduces the
claim to `P Δ(e_l) Pᵀ = 0` for the other paths `l`. That holds because the deleted arrow would have to run from the
end of one column path to the start of another, which the potential (and, for Theorem 3, the hypothesis on `u, w`)
excludes.
-/
import EatTheFibre.Path
import EatTheFibre.Invariants

open Matrix

namespace EatTheFibre

open FinQuiver CoordAlg Reflect

/-! ## The small quivers and their paths, in census order -/

/-- 1⇉2⇉3: arrows a₁ a₂ : 0 → 1, b₁ b₂ : 1 → 2 -/
abbrev W2 : FinQuiver := ⟨3, 4, ![0, 0, 1, 1], ![1, 1, 2, 2]⟩

/-- 1⇉2⇇3: arrows a₁ a₂ : 0 → 1, b₁ b₂ : 2 → 1 -/
abbrev W3 : FinQuiver := ⟨3, 4, ![0, 0, 2, 2], ![1, 1, 1, 1]⟩

/-- e1 e2 e3 a1 a2 b1 b2 a1b1 a1b2 a2b1 a2b2 (the basis of `Gen/two_pairs.lean`) -/
def enW2 : Fin 11 → W2.Path := ![⟨0, [], trivial⟩, ⟨1, [], trivial⟩, ⟨2, [], trivial⟩,
  ⟨0, [0], by decide⟩, ⟨0, [1], by decide⟩, ⟨1, [2], by decide⟩, ⟨1, [3], by decide⟩,
  ⟨0, [0, 2], by decide⟩, ⟨0, [0, 3], by decide⟩, ⟨0, [1, 2], by decide⟩, ⟨0, [1, 3], by decide⟩]

/-- e1 e2 e3 a1 a2 b1 b2 (the basis of `Gen/pairs_in.lean`) -/
def enW3 : Fin 7 → W3.Path := ![⟨0, [], trivial⟩, ⟨1, [], trivial⟩, ⟨2, [], trivial⟩,
  ⟨0, [0], by decide⟩, ⟨0, [1], by decide⟩, ⟨2, [2], by decide⟩, ⟨2, [3], by decide⟩]

lemma enW2_inj : Function.Injective enW2 := by
  intro i j h; revert i j; decide

lemma enW3_inj : Function.Injective enW3 := by
  intro i j h; revert i j; decide

/-- every path of 1⇉2⇉3 is enumerated -/
lemma enW2_surj (p : W2.Path) : ∃ i, enW2 i = p := by
  obtain ⟨x, l, hw⟩ := p
  rcases l with _ | ⟨a, _ | ⟨b, _ | ⟨c, l⟩⟩⟩
  · fin_cases x
    exacts [⟨0, rfl⟩, ⟨1, rfl⟩, ⟨2, rfl⟩]
  · fin_cases x <;> fin_cases a <;> simp [Walk, W2] at hw
    exacts [⟨3, rfl⟩, ⟨4, rfl⟩, ⟨5, rfl⟩, ⟨6, rfl⟩]
  · fin_cases x <;> fin_cases a <;> fin_cases b <;> simp [Walk, W2] at hw
    exacts [⟨7, rfl⟩, ⟨8, rfl⟩, ⟨9, rfl⟩, ⟨10, rfl⟩]
  · exfalso; fin_cases a <;> fin_cases b <;> fin_cases c <;> simp [Walk, W2] at hw

lemma enW3_surj (p : W3.Path) : ∃ i, enW3 i = p := by
  obtain ⟨x, l, hw⟩ := p
  rcases l with _ | ⟨a, _ | ⟨b, l⟩⟩
  · fin_cases x
    exacts [⟨0, rfl⟩, ⟨1, rfl⟩, ⟨2, rfl⟩]
  · fin_cases x <;> fin_cases a <;> simp [Walk, W3] at hw
    exacts [⟨3, rfl⟩, ⟨4, rfl⟩, ⟨5, rfl⟩, ⟨6, rfl⟩]
  · exfalso; fin_cases a <;> fin_cases b <;> simp [Walk, W3] at hw

instance : Fintype W2.Path := Fintype.ofSurjective enW2 enW2_surj
instance : Fintype W3.Path := Fintype.ofSurjective enW3 enW3_surj

/-- the census tables are the path algebras (checked by the kernel) -/
lemma W2_c : ∀ i j m, W2.cZ (enW2 i) (enW2 j) (enW2 m) = Q_two_pairs.c i j m := by decide +kernel
lemma W2_D : ∀ l a b, W2.DZ (enW2 l) (enW2 a) (enW2 b) = Q_two_pairs.D l a b := by decide +kernel
lemma W3_c : ∀ i j m, W3.cZ (enW3 i) (enW3 j) (enW3 m) = Q_pairs_in.c i j m := by decide +kernel
lemma W3_D : ∀ l a b, W3.DZ (enW3 l) (enW3 a) (enW3 b) = Q_pairs_in.D l a b := by decide +kernel

/-- the inverse enumeration, as a partial inverse that is never `none` -/
noncomputable def rW2 (p : W2.Path) : Option (Fin 11) := some ((Equiv.ofBijective enW2 ⟨enW2_inj, enW2_surj⟩).symm p)
noncomputable def rW3 (p : W3.Path) : Option (Fin 7) := some ((Equiv.ofBijective enW3 ⟨enW3_inj, enW3_surj⟩).symm p)

lemma embeds_W2 (k : Type*) [CommRing k] : Embeds (W2.pathAlg k) (Q_two_pairs.alg k) enW2 rW2 where
  re i := by simp only [rW2, Option.some.injEq]; exact (Equiv.ofBijective enW2 ⟨enW2_inj, enW2_surj⟩).symm_apply_apply i
  er x i h := by
    simp only [rW2, Option.some.injEq] at h
    rw [← h]; exact (Equiv.ofBijective enW2 ⟨enW2_inj, enW2_surj⟩).apply_symm_apply x
  c_emb i j m := by simp [pathAlg, CoordAlg.ofInt, W2_c]
  closed i j m hm := by simp [rW2] at hm
  ideal x y m hxy _ := by rcases hxy with h | h <;> simp [rW2] at h
  D_emb l a b := by simp [pathAlg, CoordAlg.ofInt, W2_D]
  D_local l x y hxy := by rcases hxy with h | h <;> simp [rW2] at h

lemma embeds_W3 (k : Type*) [CommRing k] : Embeds (W3.pathAlg k) (Q_pairs_in.alg k) enW3 rW3 where
  re i := by simp only [rW3, Option.some.injEq]; exact (Equiv.ofBijective enW3 ⟨enW3_inj, enW3_surj⟩).symm_apply_apply i
  er x i h := by
    simp only [rW3, Option.some.injEq] at h
    rw [← h]; exact (Equiv.ofBijective enW3 ⟨enW3_inj, enW3_surj⟩).apply_symm_apply x
  c_emb i j m := by simp [pathAlg, CoordAlg.ofInt, W3_c]
  closed i j m hm := by simp [rW3] at hm
  ideal x y m hxy _ := by rcases hxy with h | h <;> simp [rW3] at h
  D_emb l a b := by simp [pathAlg, CoordAlg.ofInt, W3_D]
  D_local l x y hxy := by rcases hxy with h | h <;> simp [rW3] at h

/-! ## Theorem 4 -/

/-- two consecutive parallel classes `a₁, a₂ : u → v`, `b₁, b₂ : v → w` in a quiver with a potential -/
structure Consecutive (Q : FinQuiver) where
  h : Fin Q.nv → ℕ
  pot : Q.Potential h
  u : Fin Q.nv
  v : Fin Q.nv
  w : Fin Q.nv
  a₁ : Fin Q.na
  a₂ : Fin Q.na
  b₁ : Fin Q.na
  b₂ : Fin Q.na
  ha₁ : Q.s a₁ = u ∧ Q.t a₁ = v
  ha₂ : Q.s a₂ = u ∧ Q.t a₂ = v
  hb₁ : Q.s b₁ = v ∧ Q.t b₁ = w
  hb₂ : Q.s b₂ = v ∧ Q.t b₂ = w
  na : a₁ ≠ a₂
  nb : b₁ ≠ b₂

namespace Consecutive

variable {Q : FinQuiver} (C : Consecutive Q)

lemma huv : C.h C.u < C.h C.v := by have := C.pot C.a₁; rwa [C.ha₁.1, C.ha₁.2] at this
lemma hvw : C.h C.v < C.h C.w := by have := C.pot C.b₁; rwa [C.hb₁.1, C.hb₁.2] at this

/-- the quiver morphism 1⇉2⇉3 → Q -/
def hom : Hom W2 Q where
  fv := ![C.u, C.v, C.w]
  fa := ![C.a₁, C.a₂, C.b₁, C.b₂]
  hs a := by fin_cases a <;> simp [W2, C.ha₁, C.ha₂, C.hb₁, C.hb₂]
  ht a := by fin_cases a <;> simp [W2, C.ha₁, C.ha₂, C.hb₁, C.hb₂]
  injv := by
    have h1 := C.huv; have h2 := C.hvw
    intro i j hij; fin_cases i <;> fin_cases j <;> simp at hij ⊢ <;> (have := congrArg C.h hij; omega)
  inja := by
    have h1 := C.huv; have h2 := C.hvw
    have sa : ∀ x y, Q.s x = C.u → Q.s y = C.v → x ≠ y := fun x y hx hy hxy => by
      subst hxy; rw [hx] at hy; rw [hy] at h1; omega
    intro i j hij; fin_cases i <;> fin_cases j <;> simp at hij ⊢
    all_goals first
      | exact C.na hij | exact C.na hij.symm | exact C.nb hij | exact C.nb hij.symm
      | exact sa _ _ C.ha₁.1 C.hb₁.1 hij | exact sa _ _ C.ha₁.1 C.hb₂.1 hij
      | exact sa _ _ C.ha₂.1 C.hb₁.1 hij | exact sa _ _ C.ha₂.1 C.hb₂.1 hij
      | exact sa _ _ C.ha₁.1 C.hb₁.1 hij.symm | exact sa _ _ C.ha₁.1 C.hb₂.1 hij.symm
      | exact sa _ _ C.ha₂.1 C.hb₁.1 hij.symm | exact sa _ _ C.ha₂.1 C.hb₂.1 hij.symm

/-- the embedding of the census algebra of 1⇉2⇉3 into `k Q` -/
noncomputable def emb : Fin 11 → Q.Path := C.hom.mapPath ∘ enW2
noncomputable def ret (x : Q.Path) : Option (Fin 11) := (C.hom.pre x).bind rW2

lemma embeds (k : Type*) [CommRing k] [Fintype Q.Path] : Embeds (Q.pathAlg k) (Q_two_pairs.alg k) C.emb C.ret :=
  (C.hom.embeds k).trans (embeds_W2 k)

end Consecutive

/-- `P_t` has nonzero columns only at `a₂` (4) and `b₂` (6) -/
lemma thm4P_cols : ∀ i j : Fin 11, j ≠ 4 → j ≠ 6 → Q_two_pairs.thm4P i j = [] := by decide

/-- **Theorem 4 on every acyclic quiver.** -/
theorem thm4_general {Q : FinQuiver} [Fintype Q.Path] (C : Consecutive Q) {k : Type*} [CommRing k] (t : k) :
    (Q.pathAlg k).IsFold (pushM C.ret (evalM Q_two_pairs.thm4P t)) := by
  have he := C.embeds k
  refine isFold_push he _ (Q_two_pairs.thm4_isFold t) fun l _ => ?_
  refine conj_delta_eq_zero he _ (fun j => j = 4 ∨ j = 6) (fun i j hj => ?_) l (fun i j hi hj => ?_)
  · push_neg at hj; simp [evalM, thm4P_cols i j hj.1 hj.2]
  · -- an arrow from the end of one column path to the start of another would contradict the potential
    have h1 := C.huv; have h2 := C.hvw
    show ((Q.DZ l (C.emb i) (C.emb j) : ℤ) : k) = 0
    rw [DZ, if_neg, Int.cast_zero]
    rintro ⟨c, hc1, hc2, -, -⟩
    have hp := C.pot c
    have tgt : ∀ i, (C.emb i).tgt = C.hom.fv (enW2 i).tgt := fun i => C.hom.mapPath_tgt _
    have src : ∀ i, (C.emb i).src = C.hom.fv (enW2 i).src := fun i => rfl
    rw [tgt] at hc1; rw [src] at hc2
    rcases hi with rfl | rfl <;> rcases hj with rfl | rfl <;>
      simp [enW2, Path.tgt, endv, W2, Consecutive.hom] at hc1 hc2 <;> rw [← hc1, hc2] at hp <;> omega

/-! ## Theorem 3 -/

/-- converging parallel classes `a₁, a₂ : u → v`, `b₁, b₂ : w → v`, `u ≠ w`, and no arrow with both ends in `{u, w}` -/
structure Converging (Q : FinQuiver) where
  h : Fin Q.nv → ℕ
  pot : Q.Potential h
  u : Fin Q.nv
  v : Fin Q.nv
  w : Fin Q.nv
  a₁ : Fin Q.na
  a₂ : Fin Q.na
  b₁ : Fin Q.na
  b₂ : Fin Q.na
  ha₁ : Q.s a₁ = u ∧ Q.t a₁ = v
  ha₂ : Q.s a₂ = u ∧ Q.t a₂ = v
  hb₁ : Q.s b₁ = w ∧ Q.t b₁ = v
  hb₂ : Q.s b₂ = w ∧ Q.t b₂ = v
  na : a₁ ≠ a₂
  nb : b₁ ≠ b₂
  huw : u ≠ w
  noArrow : ∀ c, (Q.s c = u ∨ Q.s c = w) → (Q.t c = u ∨ Q.t c = w) → False

namespace Converging

variable {Q : FinQuiver} (C : Converging Q)

lemma huv : C.h C.u < C.h C.v := by have := C.pot C.a₁; rwa [C.ha₁.1, C.ha₁.2] at this
lemma hwv : C.h C.w < C.h C.v := by have := C.pot C.b₁; rwa [C.hb₁.1, C.hb₁.2] at this

/-- the quiver morphism 1⇉2⇇3 → Q -/
def hom : Hom W3 Q where
  fv := ![C.u, C.v, C.w]
  fa := ![C.a₁, C.a₂, C.b₁, C.b₂]
  hs a := by fin_cases a <;> simp [W3, C.ha₁, C.ha₂, C.hb₁, C.hb₂]
  ht a := by fin_cases a <;> simp [W3, C.ha₁, C.ha₂, C.hb₁, C.hb₂]
  injv := by
    have h1 := C.huv; have h2 := C.hwv; have h3 := C.huw
    intro i j hij; fin_cases i <;> fin_cases j <;> simp at hij ⊢ <;> first | exact h3 hij | exact h3 hij.symm |
      (have := congrArg C.h hij; omega)
  inja := by
    have sa : ∀ x y, Q.s x = C.u → Q.s y = C.w → x ≠ y := fun x y hx hy hxy => by
      subst hxy; rw [hx] at hy; exact C.huw hy
    intro i j hij; fin_cases i <;> fin_cases j <;> simp at hij ⊢
    all_goals first
      | exact C.na hij | exact C.na hij.symm | exact C.nb hij | exact C.nb hij.symm
      | exact sa _ _ C.ha₁.1 C.hb₁.1 hij | exact sa _ _ C.ha₁.1 C.hb₂.1 hij
      | exact sa _ _ C.ha₂.1 C.hb₁.1 hij | exact sa _ _ C.ha₂.1 C.hb₂.1 hij
      | exact sa _ _ C.ha₁.1 C.hb₁.1 hij.symm | exact sa _ _ C.ha₁.1 C.hb₂.1 hij.symm
      | exact sa _ _ C.ha₂.1 C.hb₁.1 hij.symm | exact sa _ _ C.ha₂.1 C.hb₂.1 hij.symm

noncomputable def emb : Fin 7 → Q.Path := C.hom.mapPath ∘ enW3
noncomputable def ret (x : Q.Path) : Option (Fin 7) := (C.hom.pre x).bind rW3

lemma embeds (k : Type*) [CommRing k] [Fintype Q.Path] : Embeds (Q.pathAlg k) (Q_pairs_in.alg k) C.emb C.ret :=
  (C.hom.embeds k).trans (embeds_W3 k)

/-- the remaining Δ condition, for any matrix with columns at `e_u, e_w, a₂, b₂` -/
lemma delta_out [Fintype Q.Path] {k : Type*} [CommRing k] (P : Matrix (Fin 7) (Fin 7) k)
    (hcols : ∀ i j, ¬ (j = 0 ∨ j = 2 ∨ j = 4 ∨ j = 6) → P i j = 0) (l : Q.Path) :
    pushM C.ret P * (Q.pathAlg k).delta (Pi.single l 1) * (pushM C.ret P)ᵀ = 0 := by
  refine conj_delta_eq_zero (C.embeds k) _ (fun j => j = 0 ∨ j = 2 ∨ j = 4 ∨ j = 6) hcols l (fun i j hi hj => ?_)
  have h1 := C.huv; have h2 := C.hwv
  show ((Q.DZ l (C.emb i) (C.emb j) : ℤ) : k) = 0
  rw [DZ, if_neg, Int.cast_zero]
  rintro ⟨c, hc1, hc2, -, -⟩
  have hp := C.pot c
  have tgt : ∀ i, (C.emb i).tgt = C.hom.fv (enW3 i).tgt := fun i => C.hom.mapPath_tgt _
  have src : ∀ i, (C.emb i).src = C.hom.fv (enW3 i).src := fun i => rfl
  rw [tgt] at hc1; rw [src] at hc2
  rcases hi with rfl | rfl | rfl | rfl <;> rcases hj with rfl | rfl | rfl | rfl <;>
    simp [enW3, Path.tgt, endv, W3, Converging.hom] at hc1 hc2 <;>
    first
      | (rw [← hc1, hc2] at hp; omega)
      | exact C.noArrow c (by simp [← hc1]) (by simp [hc2])

end Converging

lemma thm3_cols (P : Fin 7 → Fin 7 → Poly) (hP : ∀ i j : Fin 7, ¬ (j = 0 ∨ j = 2 ∨ j = 4 ∨ j = 6) → P i j = [])
    {k : Type*} [CommRing k] (t : k) : ∀ i j, ¬ (j = 0 ∨ j = 2 ∨ j = 4 ∨ j = 6) → evalM P t i j = 0 :=
  fun i j hj => by simp [evalM, hP i j hj]

/-- **Theorem 3 on every acyclic quiver:** the four families of 1⇉2⇇3, pushed into `k Q`, are folds. -/
theorem thm3_general {Q : FinQuiver} [Fintype Q.Path] (C : Converging Q) {k : Type*} [CommRing k] (t : k) :
    (Q.pathAlg k).IsFold (pushM C.ret (evalM Q_pairs_in.thm3F1P t)) ∧
    (Q.pathAlg k).IsFold (pushM C.ret (evalM Q_pairs_in.thm3F2P t)) ∧
    (Q.pathAlg k).IsFold (pushM C.ret (evalM Q_pairs_in.thm3F3P t)) ∧
    (Q.pathAlg k).IsFold (pushM C.ret (evalM Q_pairs_in.thm3F4P t)) := by
  have he := C.embeds k
  refine ⟨?_, ?_, ?_, ?_⟩
  · exact isFold_push he _ (Q_pairs_in.thm3F1_isFold t)
      fun l _ => C.delta_out _ (thm3_cols _ (by decide) t) l
  · exact isFold_push he _ (Q_pairs_in.thm3F2_isFold t)
      fun l _ => C.delta_out _ (thm3_cols _ (by decide) t) l
  · exact isFold_push he _ (Q_pairs_in.thm3F3_isFold t)
      fun l _ => C.delta_out _ (thm3_cols _ (by decide) t) l
  · exact isFold_push he _ (Q_pairs_in.thm3F4_isFold t)
      fun l _ => C.delta_out _ (thm3_cols _ (by decide) t) l

/-! ## Rigidity on every acyclic quiver -/

/-- vertex paths -/
def isVertexPath {Q : FinQuiver} (p : Q.Path) : Prop := p.arrs = []
instance {Q : FinQuiver} : DecidablePred (isVertexPath (Q := Q)) := fun p => inferInstanceAs (Decidable (p.arrs = []))

/-- **Theorem 3, rigidity on every acyclic quiver.** No automorphism of `k Q` that fixes the vertices and preserves the
span of the paths of positive length conjugates a member of a family to a member with another parameter. -/
theorem thm3_general_rigid {Q : FinQuiver} [Fintype Q.Path] (C : Converging Q) {k : Type*} [CommRing k]
    {h hi : Matrix Q.Path Q.Path k} (hinv : hi * h = 1)
    (hfix : ∀ a, isVertexPath a → ∀ i, h i a = if i = a then 1 else 0)
    (hrad : ∀ a, ¬ isVertexPath a → ∀ i, isVertexPath i → h i a = 0) {s s' : k} :
    (h * pushM C.ret (evalM Q_pairs_in.thm3F1P s) * hi = pushM C.ret (evalM Q_pairs_in.thm3F1P s') → s = s') ∧
    (h * pushM C.ret (evalM Q_pairs_in.thm3F2P s) * hi = pushM C.ret (evalM Q_pairs_in.thm3F2P s') → s = s') ∧
    (h * pushM C.ret (evalM Q_pairs_in.thm3F3P s) * hi = pushM C.ret (evalM Q_pairs_in.thm3F3P s') → s = s') ∧
    (h * pushM C.ret (evalM Q_pairs_in.thm3F4P s) * hi = pushM C.ret (evalM Q_pairs_in.thm3F4P s') → s = s') := by
  have he := C.embeds k
  have vx : ∀ i : Fin 7, i.val < 3 → isVertexPath (C.emb i) := fun i hi => by
    fin_cases i <;> simp at hi <;> rfl
  -- conjugation keeps the entry at (emb i, emb j) for vertices i, j
  have key : ∀ (P : Matrix (Fin 7) (Fin 7) k) (P' : Matrix (Fin 7) (Fin 7) k) (i j : Fin 7), i.val < 3 → j.val < 3 →
      h * pushM C.ret P * hi = pushM C.ret P' → P' i j = P i j := fun P P' i j hi_ hj hc => by
    have := vertex_block isVertexPath hinv hfix hrad (pushM C.ret P) (vx i hi_) (vx j hj)
    rw [hc, pushM_e _ _ _ he, pushM_e _ _ _ he] at this; exact this
  refine ⟨fun hc => ?_, fun hc => ?_, fun hc => ?_, fun hc => ?_⟩
  · have := key _ _ 1 0 (by decide) (by decide) hc; simp [evalM, Q_pairs_in.thm3F1P] at this; exact this.symm
  · have := key _ _ 0 0 (by decide) (by decide) hc; simp [evalM, Q_pairs_in.thm3F2P] at this; exact this.symm
  · have := key _ _ 1 2 (by decide) (by decide) hc; simp [evalM, Q_pairs_in.thm3F3P] at this; exact this.symm
  · have := key _ _ 1 0 (by decide) (by decide) hc; simp [evalM, Q_pairs_in.thm3F4P] at this; exact this.symm

open Q_two_pairs in
/-- **Theorem 4, rigidity on every acyclic quiver.** If `h` fixes `e_u, e_v` and scales `δ_a, δ_b, δ_aδ_b` by
`λ, μ, λμ` (as every substitution within the two classes and every class weighting does), and intertwines `P_t`
with `P_t'`, then `t = t'`. -/
theorem thm4_general_rigid {Q : FinQuiver} [Fintype Q.Path] (C : Consecutive Q) {k : Type*} [CommRing k]
    {h : Matrix Q.Path Q.Path k} {l m t t' : k}
    (he1 : h *ᵥ pushV C.ret ve1 = pushV C.ret ve1) (he2 : h *ᵥ pushV C.ret ve2 = pushV C.ret ve2)
    (hda : h *ᵥ pushV C.ret δa = l • pushV C.ret δa) (hdb : h *ᵥ pushV C.ret δb = m • pushV C.ret δb)
    (hdab : h *ᵥ pushV C.ret δab = (l * m) • pushV C.ret δab)
    (hc : h * pushM C.ret (evalM thm4P t) = pushM C.ret (evalM thm4P t') * h) : t = t' := by
  have he := C.embeds k
  have col : ∀ j, pushM C.ret (evalM thm4P t) *ᵥ Pi.single (C.emb j) 1 = pushV C.ret (evalM thm4P t *ᵥ Pi.single j 1) :=
    fun j => by rw [push_single he, push_mulVec he]
  have hY := congrArg (· *ᵥ Pi.single (C.emb 6) (1 : k)) hc
  have hX := congrArg (· *ᵥ Pi.single (C.emb 4) (1 : k)) hc
  simp only [← Matrix.mulVec_mulVec] at hY hX
  rw [col, thm4_Y, pushV_add he, pushV_sub he, pushV_add he, Matrix.mulVec_add, Matrix.mulVec_sub, Matrix.mulVec_add,
    he1, he2, hdb, hdab, push_apply he] at hY
  rw [col, thm4_X, pushV_add he, pushV_sub he, pushV_smul he, pushV_smul he, Matrix.mulVec_add, Matrix.mulVec_sub,
    Matrix.mulVec_smul, Matrix.mulVec_smul, he1, hda, hdab, push_apply he] at hX
  refine thm4_core (l := l) (m := m) ((h *ᵥ Pi.single (C.emb 6) 1) ∘ C.emb) ((h *ᵥ Pi.single (C.emb 4) 1) ∘ C.emb)
    (pushV_inj he ?_) (pushV_inj he ?_)
  · rw [pushV_add he, pushV_sub he, pushV_add he, pushV_smul he, pushV_smul he]; exact hY
  · rw [pushV_add he, pushV_sub he, pushV_smul he, pushV_smul he, pushV_smul he, pushV_smul he]; exact hX

end EatTheFibre
