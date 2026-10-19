/-
# Reduction to a sub-basis

The hand proofs of Theorems 3 and 4 (v2.27, v2.28) reduce an arbitrary acyclic quiver to the path algebra `W` of the
two parallel classes. This file proves that reduction for coordinate algebras.

* `isFold_of_support`: let `S` be a set of basis indices whose complement spans a two-sided ideal. If `P` kills the
  basis vectors outside `S`, satisfies the fold conditions on basis vectors in `S`, and `P · Δ(e_l) · Pᵀ = 0` for every
  `l ∉ S`, then `P` is a fold.
* `Embeds`: a small algebra `B` sits inside `A` along an injection `e` (with partial inverse `r`): the structure constants
  and `Δ` agree, products of image vectors stay in the image, `Δ` of an image vector has both factors in the image,
  and the complement of the image spans an ideal.
* `isFold_push`: a fold `Q` of `B`, pushed into `A` (zero outside the image), is a fold of `A` as soon as
  `P · Δ(e_l) · Pᵀ = 0` for the basis vectors `l` outside the image. That last condition is the only place where the
  rest of the quiver matters.
-/
import EatTheFibre.Basic

open Matrix

namespace EatTheFibre.CoordAlg

variable {k : Type*} [CommRing k] {ι : Type*} [Fintype ι] [DecidableEq ι] (A : CoordAlg k ι)

lemma mulVec_single_one' (P : Matrix ι ι k) (j : ι) : P *ᵥ Pi.single j 1 = fun i => P i j := by
  funext i; simp [Matrix.mulVec, dotProduct, Pi.single_apply]

@[simp] lemma mul_zero_left (w : ι → k) : A.mul 0 w = 0 := by ext m; simp [mul]
@[simp] lemma mul_zero_right (w : ι → k) : A.mul w 0 = 0 := by ext m; simp [mul]
@[simp] lemma delta_zero : A.delta 0 = 0 := by ext a b; simp [delta]

/-- **Reduction to a sub-basis.** -/
theorem isFold_of_support (S : ι → Prop) [DecidablePred S]
    (ideal : ∀ i j m, (¬ S i ∨ ¬ S j) → A.c i j m ≠ 0 → ¬ S m)
    (P : Matrix ι ι k) (hcol : ∀ j, ¬ S j → ∀ i, P i j = 0) (hidem : P * P = P)
    (himg : ∀ i j, S i → S j → P *ᵥ A.mul (P *ᵥ Pi.single i 1) (P *ᵥ Pi.single j 1) =
      A.mul (P *ᵥ Pi.single i 1) (P *ᵥ Pi.single j 1))
    (hker : ∀ i j, S i → S j → P *ᵥ A.mul ((1 - P) *ᵥ Pi.single i 1) ((1 - P) *ᵥ Pi.single j 1) = 0)
    (hdel : ∀ l, S l → A.delta (P *ᵥ Pi.single l 1) = P * A.delta (Pi.single l 1) * Pᵀ)
    (hdelOut : ∀ l, ¬ S l → P * A.delta (Pi.single l 1) * Pᵀ = 0) : A.IsFold P := by
  have zcol : ∀ j, ¬ S j → P *ᵥ Pi.single j 1 = 0 := fun j hj => by
    rw [mulVec_single_one']; funext i; exact hcol j hj i
  have qcol : ∀ j, ¬ S j → (1 - P) *ᵥ Pi.single j 1 = Pi.single j 1 := fun j hj => by
    rw [Matrix.sub_mulVec, Matrix.one_mulVec, zcol j hj, sub_zero]
  -- P kills every vector supported outside S
  have kills : ∀ x : ι → k, (∀ m, S m → x m = 0) → P *ᵥ x = 0 := fun x hx => by
    funext r; simp only [Matrix.mulVec, dotProduct, Pi.zero_apply]
    exact Finset.sum_eq_zero fun m _ => by
      by_cases hm : S m
      · rw [hx m hm, mul_zero]
      · rw [hcol m hm r, zero_mul]
  -- a product with a basis vector outside S is supported outside S
  have outL : ∀ i w, ¬ S i → ∀ m, S m → A.mul (Pi.single i 1) w m = 0 := fun i w hi m hm => by
    simp only [mul]
    refine Finset.sum_eq_zero fun i' _ => Finset.sum_eq_zero fun j _ => ?_
    by_cases h : i' = i
    · subst h
      by_cases hc : A.c i' j m = 0
      · rw [hc, mul_zero]
      · exact absurd hm (ideal i' j m (Or.inl hi) hc)
    · simp [Pi.single_apply, h]
  have outR : ∀ j u, ¬ S j → ∀ m, S m → A.mul u (Pi.single j 1) m = 0 := fun j u hj m hm => by
    simp only [mul]
    refine Finset.sum_eq_zero fun i _ => Finset.sum_eq_zero fun j' _ => ?_
    by_cases h : j' = j
    · subst h
      by_cases hc : A.c i j' m = 0
      · rw [hc, mul_zero]
      · exact absurd hm (ideal i j' m (Or.inr hj) hc)
    · simp [Pi.single_apply, h]
  rw [isFold_iff_basis]
  refine ⟨hidem, ?_, ?_, ?_⟩
  · intro i j
    by_cases hi : S i
    · by_cases hj : S j
      · exact himg i j hi hj
      · simp [zcol j hj]
    · simp [zcol i hi]
  · intro i j
    by_cases hi : S i
    · by_cases hj : S j
      · exact hker i j hi hj
      · rw [qcol j hj]; exact kills _ (outR j _ hj)
    · rw [qcol i hi]; exact kills _ (outL i _ hi)
  · intro l
    by_cases hl : S l
    · exact hdel l hl
    · rw [zcol l hl, delta_zero, hdelOut l hl]

/-! ## Embedding a small algebra -/

variable {κ : Type*} [Fintype κ] [DecidableEq κ] (B : CoordAlg k κ)

/-- `B` sits inside `A` along `e`, with partial inverse `r` (`r (e i) = some i`, `r x = none` off the image) -/
structure Embeds (e : κ → ι) (r : ι → Option κ) : Prop where
  re : ∀ i, r (e i) = some i
  er : ∀ x i, r x = some i → e i = x
  c_emb : ∀ i j m, A.c (e i) (e j) (e m) = B.c i j m
  closed : ∀ i j m, r m = none → A.c (e i) (e j) m = 0
  ideal : ∀ x y m, (r x = none ∨ r y = none) → A.c x y m ≠ 0 → r m = none
  D_emb : ∀ l a b, A.D (e l) (e a) (e b) = B.D l a b
  D_local : ∀ l x y, (r x = none ∨ r y = none) → A.D (e l) x y = 0

variable {A B} {e : κ → ι} {r : ι → Option κ}

/-- extend a vector by zero -/
def pushV (r : ι → Option κ) (v : κ → k) : ι → k := fun x => match r x with | some i => v i | none => 0
/-- extend a matrix by zero -/
def pushM (r : ι → Option κ) (M : Matrix κ κ k) : Matrix ι ι k := fun x y =>
  match r x, r y with | some i, some j => M i j | _, _ => 0

lemma sum_push (h : Embeds A B e r) (f : ι → k) (hf : ∀ x, r x = none → f x = 0) : ∑ x, f x = ∑ i, f (e i) := by
  have inj : Function.Injective e := fun i j hij => by
    have := h.re i; rw [hij, h.re j] at this; exact (Option.some_injective _ this).symm
  have hm := Finset.sum_map Finset.univ ⟨e, inj⟩ f
  simp only [Function.Embedding.coeFn_mk] at hm
  rw [← hm]
  symm; apply Finset.sum_subset (Finset.subset_univ _)
  intro x _ hx
  apply hf
  cases hr : r x with
  | none => rfl
  | some i => exact absurd (Finset.mem_map.mpr ⟨i, Finset.mem_univ _, h.er x i hr⟩) hx

lemma pushV_e (v : κ → k) (i : κ) (h : Embeds A B e r) : pushV r v (e i) = v i := by simp [pushV, h.re]
lemma pushM_e (M : Matrix κ κ k) (i j : κ) (h : Embeds A B e r) : pushM r M (e i) (e j) = M i j := by
  simp [pushM, h.re]

lemma pushV_none (v : κ → k) {x : ι} (hx : r x = none) : pushV r v x = 0 := by simp [pushV, hx]
lemma pushM_none_left (M : Matrix κ κ k) {x : ι} (y : ι) (hx : r x = none) : pushM r M x y = 0 := by
  simp [pushM, hx]
lemma pushM_none_right (M : Matrix κ κ k) (x : ι) {y : ι} (hy : r y = none) : pushM r M x y = 0 := by
  unfold pushM; rw [hy]; cases r x <;> rfl

/-- evaluate a function on `ι` that is pushed: either on the image or zero -/
lemma cases_r (x : ι) (h : Embeds A B e r) : r x = none ∨ ∃ i, x = e i := by
  cases hr : r x with
  | none => exact Or.inl rfl
  | some i => exact Or.inr ⟨i, (h.er x i hr).symm⟩

lemma push_mulVec (h : Embeds A B e r) (Q : Matrix κ κ k) (v : κ → k) :
    pushM r Q *ᵥ pushV r v = pushV r (Q *ᵥ v) := by
  funext x
  rcases cases_r x h with hx | ⟨i, rfl⟩
  · rw [pushV_none _ hx]; simp [Matrix.mulVec, dotProduct, pushM_none_left _ _ hx]
  · simp only [Matrix.mulVec, dotProduct]
    rw [sum_push h _ (fun y hy => by rw [pushV_none _ hy, mul_zero])]
    simp [pushM_e _ _ _ h, pushV_e _ _ h, Matrix.mulVec, dotProduct]

lemma push_single (h : Embeds A B e r) (i : κ) : (Pi.single (e i) 1 : ι → k) = pushV r (Pi.single i 1) := by
  funext x
  rcases cases_r x h with hx | ⟨j, rfl⟩
  · rw [pushV_none _ hx]
    have : x ≠ e i := fun hxe => by rw [hxe, h.re] at hx; exact absurd hx (by simp)
    simp [Pi.single_apply, this]
  · rw [pushV_e _ _ h]
    by_cases hij : j = i
    · subst hij; simp
    · have : e j ≠ e i := fun hh => hij (by have := h.re j; rw [hh, h.re i] at this; exact (Option.some_injective _ this).symm)
      simp [Pi.single_apply, this, hij]

lemma push_one_sub (h : Embeds A B e r) (Q : Matrix κ κ k) (v : κ → k) :
    (1 - pushM r Q) *ᵥ pushV r v = pushV r ((1 - Q) *ᵥ v) := by
  rw [Matrix.sub_mulVec, Matrix.one_mulVec, push_mulVec h, Matrix.sub_mulVec, Matrix.one_mulVec]
  funext x; rcases cases_r x h with hx | ⟨i, rfl⟩
  · simp [pushV_none _ hx]
  · simp [pushV_e _ _ h]

lemma push_mul (h : Embeds A B e r) (u w : κ → k) : A.mul (pushV r u) (pushV r w) = pushV r (B.mul u w) := by
  funext m
  have inner : ∀ x, (∑ y, pushV r u x * pushV r w y * A.c x y m) = ∑ j, pushV r u x * w j * A.c x (e j) m :=
    fun x => by
      rw [sum_push h _ (fun y hy => by rw [pushV_none _ hy, mul_zero, zero_mul])]
      simp [pushV_e _ _ h]
  simp only [mul, inner]
  rw [sum_push h _ (fun x hx => by simp [pushV_none _ hx])]
  simp only [pushV_e _ _ h]
  rcases cases_r m h with hm | ⟨m', rfl⟩
  · rw [pushV_none _ hm]
    exact Finset.sum_eq_zero fun i _ => Finset.sum_eq_zero fun j _ => by rw [h.closed i j m hm, mul_zero]
  · rw [pushV_e _ _ h]; simp [h.c_emb]; rfl

lemma push_delta (h : Embeds A B e r) (v : κ → k) : A.delta (pushV r v) = pushM r (B.delta v) := by
  ext x y
  simp only [delta]
  rw [sum_push h _ (fun z hz => by rw [pushV_none _ hz, zero_mul])]
  simp only [pushV_e _ _ h]
  rcases cases_r x h with hx | ⟨a, rfl⟩
  · rw [pushM_none_left _ _ hx]; exact Finset.sum_eq_zero fun l _ => by rw [h.D_local l _ _ (Or.inl hx), mul_zero]
  rcases cases_r y h with hy | ⟨b, rfl⟩
  · rw [pushM_none_right _ _ hy]; exact Finset.sum_eq_zero fun l _ => by rw [h.D_local l _ _ (Or.inr hy), mul_zero]
  rw [pushM_e _ _ _ h]; simp [h.D_emb]; rfl

lemma push_matmul (h : Embeds A B e r) (Q M : Matrix κ κ k) : pushM r Q * pushM r M = pushM r (Q * M) := by
  ext x y
  simp only [Matrix.mul_apply]
  rw [sum_push h _ (fun z hz => by rw [pushM_none_left _ _ hz, mul_zero])]
  rcases cases_r x h with hx | ⟨a, rfl⟩
  · simp [pushM_none_left _ _ hx]
  rcases cases_r y h with hy | ⟨b, rfl⟩
  · simp [pushM_none_right _ _ hy]
  simp [pushM_e _ _ _ h, Matrix.mul_apply]

lemma push_transpose (h : Embeds A B e r) (Q : Matrix κ κ k) : (pushM r Q)ᵀ = pushM r Qᵀ := by
  ext x y
  simp only [Matrix.transpose_apply]
  rcases cases_r x h with hx | ⟨a, rfl⟩
  · rw [pushM_none_right _ _ hx, pushM_none_left _ _ hx]
  rcases cases_r y h with hy | ⟨b, rfl⟩
  · rw [pushM_none_left _ _ hy, pushM_none_right _ _ hy]
  simp [pushM_e _ _ _ h]

lemma push_conj (h : Embeds A B e r) (Q M : Matrix κ κ k) :
    pushM r Q * pushM r M * (pushM r Q)ᵀ = pushM r (Q * M * Qᵀ) := by
  rw [push_transpose h, push_matmul h, push_matmul h]

lemma ne_none_of_isSome {x : ι} (hx : (r x).isSome = true) : r x ≠ none := by
  intro h'; rw [h'] at hx; exact absurd hx (by simp)

/-- **Pushing a fold.** A fold of the small algebra, extended by zero, is a fold of the big one, provided
`P · Δ(e_l) · Pᵀ = 0` for the basis vectors `l` outside the image. -/
theorem isFold_push (h : Embeds A B e r) (Q : Matrix κ κ k) (hQ : B.IsFold Q)
    (hdelOut : ∀ l, r l = none → pushM r Q * A.delta (Pi.single l 1) * (pushM r Q)ᵀ = 0) :
    A.IsFold (pushM r Q) := by
  refine A.isFold_of_support (fun x => (r x).isSome)
    (fun i j m hij hc => ?_) _ (fun j hj i => ?_) ?_ ?_ ?_ ?_ (fun l hl => hdelOut l (by simpa using hl))
  · rcases hij with hi | hj
    · simpa using h.ideal i j m (Or.inl (by simpa using hi)) hc
    · simpa using h.ideal i j m (Or.inr (by simpa using hj)) hc
  · exact pushM_none_right _ _ (by simpa using hj)
  · ext x y
    simp only [Matrix.mul_apply]
    rw [sum_push h _ (fun z hz => by rw [pushM_none_right _ _ hz, zero_mul])]
    rcases cases_r x h with hx | ⟨a, rfl⟩
    · simp [pushM_none_left _ _ hx]
    rcases cases_r y h with hy | ⟨b, rfl⟩
    · simp [pushM_none_right _ _ hy]
    simp only [pushM_e _ _ _ h]
    rw [← Matrix.mul_apply, hQ.idem]
  · intro x y hx hy
    obtain ⟨i, rfl⟩ := (cases_r x h).resolve_left (ne_none_of_isSome hx)
    obtain ⟨j, rfl⟩ := (cases_r y h).resolve_left (ne_none_of_isSome hy)
    rw [push_single h, push_single h, push_mulVec h, push_mulVec h, push_mul h, push_mulVec h,
      hQ.image]
  · intro x y hx hy
    obtain ⟨i, rfl⟩ := (cases_r x h).resolve_left (ne_none_of_isSome hx)
    obtain ⟨j, rfl⟩ := (cases_r y h).resolve_left (ne_none_of_isSome hy)
    rw [push_single h, push_single h, push_one_sub h, push_one_sub h, push_mul h, push_mulVec h,
      hQ.kernel]
    funext x; rcases cases_r x h with hx | ⟨_, rfl⟩ <;> simp [pushV_none, pushV_e _ _ h, *]
  · intro x hx
    obtain ⟨l, rfl⟩ := (cases_r x h).resolve_left (ne_none_of_isSome hx)
    rw [push_single h, push_mulVec h, push_delta h, push_delta h, push_conj h, hQ.delta]

end EatTheFibre.CoordAlg

namespace EatTheFibre.CoordAlg

variable {k : Type*} [CommRing k] {ι κ μ : Type*} [Fintype ι] [DecidableEq ι] [Fintype κ] [DecidableEq κ]
  [Fintype μ] [DecidableEq μ] {A : CoordAlg k ι} {B : CoordAlg k κ} {C : CoordAlg k μ}

/-- embeddings compose -/
theorem Embeds.trans {e : κ → ι} {r : ι → Option κ} {e' : μ → κ} {r' : κ → Option μ}
    (h : Embeds A B e r) (h' : Embeds B C e' r') : Embeds A C (e ∘ e') (fun x => (r x).bind r') where
  re i := by simp [h.re, h'.re]
  er x i hx := by
    cases hr : r x with
    | none => simp [hr] at hx
    | some b => simp only [hr, Option.bind_some] at hx; simp [h'.er b i hx, h.er x b hr]
  c_emb i j m := by simp [h.c_emb, h'.c_emb]
  closed i j m hm := by
    cases hr : r m with
    | none => exact h.closed _ _ m hr
    | some b =>
      simp only [hr, Option.bind_some] at hm
      rw [← h.er m b hr, Function.comp_apply, Function.comp_apply, h.c_emb, h'.closed i j b hm]
  ideal x y m hxy hc := by
    cases hrm : r m with
    | none => rfl
    | some bm =>
      simp only [Option.bind_some]
      cases hrx : r x with
      | none => exact absurd hrm (by rw [h.ideal x y m (Or.inl hrx) hc]; simp)
      | some bx =>
        cases hry : r y with
        | none => exact absurd hrm (by rw [h.ideal x y m (Or.inr hry) hc]; simp)
        | some byy =>
          rw [← h.er x bx hrx, ← h.er y byy hry, ← h.er m bm hrm, h.c_emb] at hc
          refine h'.ideal bx byy bm ?_ hc
          simpa [hrx, hry] using hxy
  D_emb l a b := by simp [h.D_emb, h'.D_emb]
  D_local l x y hxy := by
    cases hrx : r x with
    | none => exact h.D_local _ x y (Or.inl hrx)
    | some bx =>
      cases hry : r y with
      | none => exact h.D_local _ x y (Or.inr hry)
      | some byy =>
        rw [← h.er x bx hrx, ← h.er y byy hry, Function.comp_apply, h.D_emb]
        exact h'.D_local l bx byy (by simpa [hrx, hry] using hxy)

/-- **The remaining Δ condition.** If the pushed matrix has its nonzero columns at indices in `T`, and `Δ(e_l)` has no
component on `e(T) ⊗ e(T)`, then `P · Δ(e_l) · Pᵀ = 0`. -/
theorem conj_delta_eq_zero {e : κ → ι} {r : ι → Option κ} (h : Embeds A B e r) (Q : Matrix κ κ k) (T : κ → Prop)
    (hT : ∀ i j, ¬ T j → Q i j = 0) (l : ι) (hD : ∀ i j, T i → T j → A.D l (e i) (e j) = 0) :
    pushM r Q * A.delta (Pi.single l 1) * (pushM r Q)ᵀ = 0 := by
  have hd : ∀ a b, A.delta (Pi.single l 1) a b = A.D l a b := fun a b => by simp [delta, Pi.single_apply]
  ext x y
  simp only [Matrix.mul_apply, Matrix.transpose_apply, hd, Matrix.zero_apply]
  rw [sum_push h _ (fun z hz => by rw [pushM_none_right _ _ hz, mul_zero])]
  refine Finset.sum_eq_zero fun j _ => ?_
  rw [sum_push h _ (fun z hz => by rw [pushM_none_right _ _ hz, zero_mul])]
  rcases cases_r x h with hx | ⟨a, rfl⟩
  · simp [pushM_none_left _ _ hx]
  rcases cases_r y h with hy | ⟨b, rfl⟩
  · simp [pushM_none_left _ _ hy]
  simp only [pushM_e _ _ _ h]
  by_cases hj : T j
  · rw [Finset.sum_eq_zero fun i _ => ?_, zero_mul]
    by_cases hi : T i
    · rw [hD i j hi hj, mul_zero]
    · rw [hT a i hi, zero_mul]
  · rw [hT b j hj, mul_zero]

end EatTheFibre.CoordAlg

namespace EatTheFibre.CoordAlg

variable {k : Type*} [CommRing k] {ι κ : Type*} [Fintype ι] [DecidableEq ι] [Fintype κ] [DecidableEq κ]
  {A : CoordAlg k ι} {B : CoordAlg k κ} {e : κ → ι} {r : ι → Option κ}

/-- a pushed matrix sees only the image coordinates of a vector -/
lemma push_apply (h : Embeds A B e r) (Q : Matrix κ κ k) (w : ι → k) : pushM r Q *ᵥ w = pushV r (Q *ᵥ (w ∘ e)) := by
  funext x
  rcases cases_r x h with hx | ⟨i, rfl⟩
  · rw [pushV_none _ hx]; simp [Matrix.mulVec, dotProduct, pushM_none_left _ _ hx]
  · simp only [Matrix.mulVec, dotProduct]
    rw [sum_push h _ (fun y hy => by rw [pushM_none_right _ _ hy, zero_mul])]
    simp [pushM_e _ _ _ h, pushV_e _ _ h, Matrix.mulVec, dotProduct]

lemma pushV_add (h : Embeds A B e r) (u v : κ → k) : pushV r (u + v) = pushV r u + pushV r v := by
  funext x; rcases cases_r x h with hx | ⟨i, rfl⟩ <;> simp [pushV_none, pushV_e _ _ h, *]
lemma pushV_sub (h : Embeds A B e r) (u v : κ → k) : pushV r (u - v) = pushV r u - pushV r v := by
  funext x; rcases cases_r x h with hx | ⟨i, rfl⟩ <;> simp [pushV_none, pushV_e _ _ h, *]
lemma pushV_smul (h : Embeds A B e r) (a : k) (u : κ → k) : pushV r (a • u) = a • pushV r u := by
  funext x; rcases cases_r x h with hx | ⟨i, rfl⟩ <;> simp [pushV_none, pushV_e _ _ h, *]
lemma pushV_inj (h : Embeds A B e r) {u v : κ → k} (huv : pushV r u = pushV r v) : u = v := by
  funext i; have := congrFun huv (e i); rwa [pushV_e _ _ h, pushV_e _ _ h] at this

end EatTheFibre.CoordAlg
