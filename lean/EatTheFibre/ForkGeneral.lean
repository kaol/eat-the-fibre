/-
# Theorem 5 on every acyclic quiver containing the fork as a full subquiver (v2.77)

**Setting.** `Q` is a finite quiver with an acyclicity potential. Inside it sit vertices `i, j, k, k'` (`k ≠ k'`), three
distinct arrows `a₁, a₂, a₃ : i → j` and arrows `x : j → k`, `x' : j → k'`, and these five are the only arrows of `Q`
with both ends in `{i, j, k, k'}` (a full subquiver isomorphic to 1⇛2, 2→3, 2→4). Nothing is assumed about the rest of
`Q`: more arrows into or out of these vertices, longer paths, other classes.

* `fork_general`: the fork families `P_t`, `P⁻_t` of `Gen/K3_tgtfork.lean`, pushed into `k Q` (zero on every other
  path), are folds of `k Q`, over every commutative ring and for every `t`.
* `fork_general_rigid`, `forkM_general_rigid`: if an endomorphism `h` of `k Q` acts on the span of the paths of the
  subquiver as the census symmetry `substK g wa wb wc` does (every class substitution and class weighting of `Q`
  does: by fullness the three classes of the subquiver are classes of `Q`) and intertwines `P_t` with `P_t'`, then
  `t = t'`. No invertibility is assumed.
* `fork_forkM_general_disjoint`: no such `h` intertwines a member of `fork` with a member of `forkM` over a nontrivial
  ring.

Proof: as `General.lean` (Theorems 3, 4). The census table of 1⇛2, 2→3, 2→4 is its path algebra (`decide +kernel`);
the subquiver gives an injective quiver morphism, hence an `Embeds`; `isFold_push` leaves the condition
`P Δ(e_l) Pᵀ = 0` for paths `l` off the image. A deletion `l = p·c·q` with `p, q` image paths has `c` joining two
vertices of the subquiver, so `c` is one of the five arrows by fullness, and then `l` is itself an image path.
-/
import EatTheFibre.Path
import EatTheFibre.Gen.K3_tgtfork_rigid

open Matrix

namespace EatTheFibre

open FinQuiver CoordAlg Reflect

/-- 1⇛2, 2→3, 2→4: arrows a₁ a₂ a₃ : 0 → 1, x : 1 → 2, x' : 1 → 3 -/
abbrev W5 : FinQuiver := ⟨4, 5, ![0, 0, 0, 1, 1], ![1, 1, 1, 2, 3]⟩

/-- the basis of `Gen/K3_tgtfork.lean`: e1 e2 e3 e4 a1 a2 a3 b c a1b a1c a2b a2c a3b a3c -/
def enW5 : Fin 15 → W5.Path := ![⟨0, [], trivial⟩, ⟨1, [], trivial⟩, ⟨2, [], trivial⟩, ⟨3, [], trivial⟩,
  ⟨0, [0], by decide⟩, ⟨0, [1], by decide⟩, ⟨0, [2], by decide⟩, ⟨1, [3], by decide⟩, ⟨1, [4], by decide⟩,
  ⟨0, [0, 3], by decide⟩, ⟨0, [0, 4], by decide⟩, ⟨0, [1, 3], by decide⟩, ⟨0, [1, 4], by decide⟩,
  ⟨0, [2, 3], by decide⟩, ⟨0, [2, 4], by decide⟩]

lemma enW5_inj : Function.Injective enW5 := by
  intro i j h; revert i j; decide

lemma enW5_surj (p : W5.Path) : ∃ i, enW5 i = p := by
  obtain ⟨x, l, hw⟩ := p
  rcases l with _ | ⟨a, _ | ⟨b, _ | ⟨c, l⟩⟩⟩
  · fin_cases x
    exacts [⟨0, rfl⟩, ⟨1, rfl⟩, ⟨2, rfl⟩, ⟨3, rfl⟩]
  · fin_cases x <;> fin_cases a <;> simp [Walk, W5] at hw
    exacts [⟨4, rfl⟩, ⟨5, rfl⟩, ⟨6, rfl⟩, ⟨7, rfl⟩, ⟨8, rfl⟩]
  · fin_cases x <;> fin_cases a <;> fin_cases b <;> simp [Walk, W5] at hw
    exacts [⟨9, rfl⟩, ⟨10, rfl⟩, ⟨11, rfl⟩, ⟨12, rfl⟩, ⟨13, rfl⟩, ⟨14, rfl⟩]
  · exfalso; fin_cases a <;> fin_cases b <;> fin_cases c <;> simp [Walk, W5] at hw

instance : Fintype W5.Path := Fintype.ofSurjective enW5 enW5_surj

/-- the census tables of 1⇛2, 2→3, 2→4 are its path algebra (checked by the kernel) -/
lemma W5_c : ∀ i j m, W5.cZ (enW5 i) (enW5 j) (enW5 m) = Q_K3_tgtfork.c i j m := by decide +kernel
lemma W5_D : ∀ l a b, W5.DZ (enW5 l) (enW5 a) (enW5 b) = Q_K3_tgtfork.D l a b := by decide +kernel

noncomputable def rW5 (p : W5.Path) : Option (Fin 15) := some ((Equiv.ofBijective enW5 ⟨enW5_inj, enW5_surj⟩).symm p)

lemma embeds_W5 (k : Type*) [CommRing k] : Embeds (W5.pathAlg k) (Q_K3_tgtfork.alg k) enW5 rW5 where
  re i := by simp only [rW5, Option.some.injEq]; exact (Equiv.ofBijective enW5 ⟨enW5_inj, enW5_surj⟩).symm_apply_apply i
  er x i h := by
    simp only [rW5, Option.some.injEq] at h
    rw [← h]; exact (Equiv.ofBijective enW5 ⟨enW5_inj, enW5_surj⟩).apply_symm_apply x
  c_emb i j m := by simp [pathAlg, CoordAlg.ofInt, W5_c]
  closed i j m hm := by simp [rW5] at hm
  ideal x y m hxy _ := by rcases hxy with h | h <;> simp [rW5] at h
  D_emb l a b := by simp [pathAlg, CoordAlg.ofInt, W5_D]
  D_local l x y hxy := by rcases hxy with h | h <;> simp [rW5] at h

/-- a full subquiver 1⇛2, 2→3, 2→4 of a quiver with a potential -/
structure Fork (Q : FinQuiver) where
  h : Fin Q.nv → ℕ
  pot : Q.Potential h
  i : Fin Q.nv
  j : Fin Q.nv
  k : Fin Q.nv
  k' : Fin Q.nv
  a₁ : Fin Q.na
  a₂ : Fin Q.na
  a₃ : Fin Q.na
  x : Fin Q.na
  x' : Fin Q.na
  ha₁ : Q.s a₁ = i ∧ Q.t a₁ = j
  ha₂ : Q.s a₂ = i ∧ Q.t a₂ = j
  ha₃ : Q.s a₃ = i ∧ Q.t a₃ = j
  hx : Q.s x = j ∧ Q.t x = k
  hx' : Q.s x' = j ∧ Q.t x' = k'
  n12 : a₁ ≠ a₂
  n13 : a₁ ≠ a₃
  n23 : a₂ ≠ a₃
  hkk : k ≠ k'
  full : ∀ c, (Q.s c = i ∨ Q.s c = j ∨ Q.s c = k ∨ Q.s c = k') → (Q.t c = i ∨ Q.t c = j ∨ Q.t c = k ∨ Q.t c = k') →
    c = a₁ ∨ c = a₂ ∨ c = a₃ ∨ c = x ∨ c = x'

namespace Fork

variable {Q : FinQuiver} (C : Fork Q)

lemma hij : C.h C.i < C.h C.j := by have := C.pot C.a₁; rwa [C.ha₁.1, C.ha₁.2] at this
lemma hjk : C.h C.j < C.h C.k := by have := C.pot C.x; rwa [C.hx.1, C.hx.2] at this
lemma hjk' : C.h C.j < C.h C.k' := by have := C.pot C.x'; rwa [C.hx'.1, C.hx'.2] at this

/-- the quiver morphism 1⇛2, 2→3, 2→4 → Q -/
def hom : Hom W5 Q where
  fv := ![C.i, C.j, C.k, C.k']
  fa := ![C.a₁, C.a₂, C.a₃, C.x, C.x']
  hs a := by fin_cases a <;> simp [W5, C.ha₁, C.ha₂, C.ha₃, C.hx, C.hx']
  ht a := by fin_cases a <;> simp [W5, C.ha₁, C.ha₂, C.ha₃, C.hx, C.hx']
  injv := by
    have h1 := C.hij; have h2 := C.hjk; have h3 := C.hjk'; have h4 := C.hkk
    intro u v huv; fin_cases u <;> fin_cases v <;> simp at huv ⊢ <;>
      first | exact h4 huv | exact h4 huv.symm | (have := congrArg C.h huv; omega)
  inja := by
    have hij : C.i ≠ C.j := fun e => by have := C.hij; rw [e] at this; omega
    have ns : ∀ a b, Q.s a ≠ Q.s b → a ≠ b := fun a b h e => h (by rw [e])
    have nt : ∀ a b, Q.t a ≠ Q.t b → a ≠ b := fun a b h e => h (by rw [e])
    have ax : ∀ a, Q.s a = C.i → a ≠ C.x := fun a ha => ns _ _ (by rw [ha, C.hx.1]; exact hij)
    have ax' : ∀ a, Q.s a = C.i → a ≠ C.x' := fun a ha => ns _ _ (by rw [ha, C.hx'.1]; exact hij)
    have xx : C.x ≠ C.x' := nt _ _ (by rw [C.hx.2, C.hx'.2]; exact C.hkk)
    intro u v huv; fin_cases u <;> fin_cases v <;> simp at huv ⊢
    all_goals first
      | exact C.n12 huv | exact C.n12 huv.symm | exact C.n13 huv | exact C.n13 huv.symm
      | exact C.n23 huv | exact C.n23 huv.symm | exact xx huv | exact xx huv.symm
      | exact ax _ C.ha₁.1 huv | exact ax _ C.ha₂.1 huv | exact ax _ C.ha₃.1 huv
      | exact ax' _ C.ha₁.1 huv | exact ax' _ C.ha₂.1 huv | exact ax' _ C.ha₃.1 huv
      | exact ax _ C.ha₁.1 huv.symm | exact ax _ C.ha₂.1 huv.symm | exact ax _ C.ha₃.1 huv.symm
      | exact ax' _ C.ha₁.1 huv.symm | exact ax' _ C.ha₂.1 huv.symm | exact ax' _ C.ha₃.1 huv.symm

noncomputable def emb : Fin 15 → Q.Path := C.hom.mapPath ∘ enW5
noncomputable def ret (p : Q.Path) : Option (Fin 15) := (C.hom.pre p).bind rW5

lemma embeds (k : Type*) [CommRing k] [Fintype Q.Path] : Embeds (Q.pathAlg k) (Q_K3_tgtfork.alg k) C.emb C.ret :=
  (C.hom.embeds k).trans (embeds_W5 k)

/-- every vertex of the subquiver is one of `i, j, k, k'` -/
lemma fv_mem (v : Fin 4) : C.hom.fv v = C.i ∨ C.hom.fv v = C.j ∨ C.hom.fv v = C.k ∨ C.hom.fv v = C.k' := by
  fin_cases v <;> simp [hom]

/-- an arrow joining two vertices of the subquiver is an arrow of it (fullness) -/
lemma arrow_mem (c : Fin Q.na) (u v : Fin 4) (hu : Q.s c = C.hom.fv u) (hv : Q.t c = C.hom.fv v) :
    ∃ a, C.hom.fa a = c := by
  rcases C.full c (by rcases C.fv_mem u with h | h | h | h <;> simp [hu, h])
      (by rcases C.fv_mem v with h | h | h | h <;> simp [hv, h]) with h | h | h | h | h
  exacts [⟨0, by simp [hom, h]⟩, ⟨1, by simp [hom, h]⟩, ⟨2, by simp [hom, h]⟩, ⟨3, by simp [hom, h]⟩,
    ⟨4, by simp [hom, h]⟩]

/-- the remaining Δ condition: off the image, `P Δ(e_l) Pᵀ = 0` for every matrix `P` -/
lemma delta_out [Fintype Q.Path] {k : Type*} [CommRing k] (P : Matrix (Fin 15) (Fin 15) k) (l : Q.Path)
    (hl : C.ret l = none) :
    pushM C.ret P * (Q.pathAlg k).delta (Pi.single l 1) * (pushM C.ret P)ᵀ = 0 := by
  refine conj_delta_eq_zero (C.embeds k) P (fun _ => True) (fun _ _ h => absurd trivial h) l (fun u v _ _ => ?_)
  show ((Q.DZ l (C.emb u) (C.emb v) : ℤ) : k) = 0
  rw [DZ, if_neg, Int.cast_zero]
  rintro ⟨c, hc1, hc2, hc3, hc4⟩
  have tgt : (C.emb u).tgt = C.hom.fv (enW5 u).tgt := C.hom.mapPath_tgt _
  obtain ⟨a, ha⟩ := C.arrow_mem c (enW5 u).tgt (enW5 v).src (by rw [← hc1, tgt]) (by rw [hc2]; rfl)
  obtain ⟨p, hp⟩ := C.hom.mem_range (v := (enW5 u).src) (l := (enW5 u).arrs ++ a :: (enW5 v).arrs)
    (by rw [hc3]; rfl) (by rw [hc4]; simp [emb, Hom.mapPath, ha])
  rw [← hp] at hl
  simp [ret, Hom.pre_map, rW5] at hl

end Fork

open Q_K3_tgtfork in
/-- **Theorem 5 on every acyclic quiver with a full fork:** both fork families, pushed into `k Q`, are folds. -/
theorem fork_general {Q : FinQuiver} [Fintype Q.Path] (C : Fork Q) {k : Type*} [CommRing k] (t : k) :
    (Q.pathAlg k).IsFold (pushM C.ret (evalM forkP t)) ∧ (Q.pathAlg k).IsFold (pushM C.ret (evalM forkMP t)) :=
  ⟨isFold_push (C.embeds k) _ (fork_isFold t) fun l hl => C.delta_out _ l hl,
   isFold_push (C.embeds k) _ (forkM_isFold t) fun l hl => C.delta_out _ l hl⟩

/-- an intertwining in `k Q` restricts to the subquiver when `h` acts there as `S` -/
lemma Fork.restrict {Q : FinQuiver} [Fintype Q.Path] (C : Fork Q) {k : Type*} [CommRing k]
    {h : Matrix Q.Path Q.Path k} {S P P' : Matrix (Fin 15) (Fin 15) k}
    (hh : ∀ v : Fin 15 → k, h *ᵥ pushV C.ret v = pushV C.ret (S *ᵥ v))
    (hc : h * pushM C.ret P = pushM C.ret P' * h) : S * P = P' * S := by
  have he := C.embeds k
  have col : ∀ c : Fin 15, (S * P) *ᵥ Pi.single c 1 = (P' * S) *ᵥ Pi.single c 1 := fun c => by
    have e := congrArg (· *ᵥ pushV C.ret (Pi.single c (1 : k))) hc
    simp only [← Matrix.mulVec_mulVec] at e
    rw [push_mulVec he, hh, hh, push_mulVec he] at e
    rw [← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec]; exact pushV_inj he e
  ext r c
  have := congrFun (col c) r
  simpa [Matrix.mulVec, dotProduct, Pi.single_apply] using this

open Q_K3_tgtfork in
/-- **Theorem 5, rigidity on every acyclic quiver with a full fork.** -/
theorem fork_general_rigid {Q : FinQuiver} [Fintype Q.Path] (C : Fork Q) {k : Type*} [CommRing k]
    {h : Matrix Q.Path Q.Path k} {g00 g01 g02 g10 g11 g12 wa wb wc t t' : k}
    (hh : ∀ v : Fin 15 → k, h *ᵥ pushV C.ret v = pushV C.ret (substK g00 g01 g02 g10 g11 g12 wa wb wc *ᵥ v))
    (hc : h * pushM C.ret (evalM forkP t) = pushM C.ret (evalM forkP t') * h) : t = t' :=
  fork_rigid_GT (C.restrict hh hc)

open Q_K3_tgtfork in
theorem forkM_general_rigid {Q : FinQuiver} [Fintype Q.Path] (C : Fork Q) {k : Type*} [CommRing k]
    {h : Matrix Q.Path Q.Path k} {g00 g01 g02 g10 g11 g12 wa wb wc t t' : k}
    (hh : ∀ v : Fin 15 → k, h *ᵥ pushV C.ret v = pushV C.ret (substK g00 g01 g02 g10 g11 g12 wa wb wc *ᵥ v))
    (hc : h * pushM C.ret (evalM forkMP t) = pushM C.ret (evalM forkMP t') * h) : t = t' :=
  forkM_rigid_GT (C.restrict hh hc)

open Q_K3_tgtfork in
theorem fork_forkM_general_disjoint {Q : FinQuiver} [Fintype Q.Path] (C : Fork Q) {k : Type*} [CommRing k]
    {h : Matrix Q.Path Q.Path k} {g00 g01 g02 g10 g11 g12 wa wb wc t t' : k}
    (hh : ∀ v : Fin 15 → k, h *ᵥ pushV C.ret v = pushV C.ret (substK g00 g01 g02 g10 g11 g12 wa wb wc *ᵥ v))
    (hc : h * pushM C.ret (evalM forkP t) = pushM C.ret (evalM forkMP t') * h) : (1 : k) = 0 :=
  fork_forkM_disjoint (C.restrict hh hc)

end EatTheFibre
