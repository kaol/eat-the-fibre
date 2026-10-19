/-
# 𝔼, the seven-valued logic, derived (v2 §8)

* **L1.** The values are the subsets of {T, F, M}, nonempty in the intended semantics. The connectives are the
  connective-wise (power-set) lift of Bochvar's tables, where M is infectious. A value is designated iff it contains T.
* **L2.** `collapse` sends a value to its classical part: B if it holds T and F, N if it holds neither. It is a
  homomorphism onto the four-valued matrix S_fde (LP on T/B/F, N infectious). It preserves and reflects designation,
  and it has a section. So 𝔼-consequence *is* S_fde-consequence (`E_iff_Sfde`).
* **L3.** `Γ ⊨ φ ⟺ Γ ⊨_LP φ ∧ Var φ ⊆ Var Γ`, for all finite `Γ` and all `φ` (`L3`).
* **L4.** There are no tautologies and no unsatisfiable premise sets (`L4_no_tautology`, `L4_all_satisfiable`).
* **L5.** The coNP-hardness reduction: `γ` is classically unsatisfiable iff `γ ⊨ ⋁_{p ∈ Var γ} (p ∧ ¬p)` (`L5`).
* **Meta-1.** `¬(μ ∨ ν) ∧ (μ ∧ ν)` is designated in exactly 4 of the 49 assignments (`meta1`).
-/
import Mathlib

namespace EatTheFibre.Discordia.Logic

/-! ## Formulas -/

inductive Fml
  | var (n : ℕ)
  | neg (φ : Fml)
  | conj (φ ψ : Fml)
  | disj (φ ψ : Fml)
  deriving DecidableEq

namespace Fml

def vars : Fml → List ℕ
  | var n => [n]
  | neg φ => vars φ
  | conj φ ψ => vars φ ++ vars ψ
  | disj φ ψ => vars φ ++ vars ψ

lemma vars_ne_nil : ∀ φ : Fml, vars φ ≠ []
  | var _ => by simp [vars]
  | neg φ => vars_ne_nil φ
  | conj φ _ => by simp [vars, vars_ne_nil φ]
  | disj φ _ => by simp [vars, vars_ne_nil φ]

/-- evaluation in any matrix given by its three operations -/
def eval {V : Type*} (n : V → V) (c d : V → V → V) (v : ℕ → V) : Fml → V
  | var p => v p
  | neg φ => n (eval n c d v φ)
  | conj φ ψ => c (eval n c d v φ) (eval n c d v ψ)
  | disj φ ψ => d (eval n c d v φ) (eval n c d v ψ)

/-- a homomorphism of matrices commutes with evaluation -/
lemma eval_hom {V W : Type*} {n : V → V} {c d : V → V → V} {n' : W → W} {c' d' : W → W → W} (h : V → W)
    (hn : ∀ a, h (n a) = n' (h a)) (hc : ∀ a b, h (c a b) = c' (h a) (h b)) (hd : ∀ a b, h (d a b) = d' (h a) (h b))
    (v : ℕ → V) : ∀ φ, h (eval n c d v φ) = eval n' c' d' (h ∘ v) φ
  | var _ => rfl
  | neg φ => by simp [eval, hn, eval_hom h hn hc hd v φ]
  | conj φ ψ => by simp [eval, hc, eval_hom h hn hc hd v φ, eval_hom h hn hc hd v ψ]
  | disj φ ψ => by simp [eval, hd, eval_hom h hn hc hd v φ, eval_hom h hn hc hd v ψ]

lemma eval_congr {V : Type*} {n : V → V} {c d : V → V → V} {v w : ℕ → V} :
    ∀ φ, (∀ p ∈ vars φ, v p = w p) → eval n c d v φ = eval n c d w φ
  | var p, h => h p (by simp [vars])
  | neg φ, h => by simp only [eval, eval_congr φ h]
  | conj φ ψ, h => by
      simp only [vars, List.mem_append] at h
      simp only [eval, eval_congr φ fun p hp => h p (.inl hp), eval_congr ψ fun p hp => h p (.inr hp)]
  | disj φ ψ, h => by
      simp only [vars, List.mem_append] at h
      simp only [eval, eval_congr φ fun p hp => h p (.inl hp), eval_congr ψ fun p hp => h p (.inr hp)]

/-- a predicate closed under the operations holds of every formula whose atoms satisfy it -/
lemma eval_closed {V : Type*} {n : V → V} {c d : V → V → V} (S : V → Prop) (hn : ∀ a, S a → S (n a))
    (hc : ∀ a b, S a → S b → S (c a b)) (hd : ∀ a b, S a → S b → S (d a b)) {v : ℕ → V} :
    ∀ φ, (∀ p ∈ vars φ, S (v p)) → S (eval n c d v φ)
  | var p, h => h p (by simp [vars])
  | neg φ, h => hn _ (eval_closed S hn hc hd φ h)
  | conj φ ψ, h => by
      simp only [vars, List.mem_append] at h
      exact hc _ _ (eval_closed S hn hc hd φ fun p hp => h p (.inl hp))
        (eval_closed S hn hc hd ψ fun p hp => h p (.inr hp))
  | disj φ ψ, h => by
      simp only [vars, List.mem_append] at h
      exact hd _ _ (eval_closed S hn hc hd φ fun p hp => h p (.inl hp))
        (eval_closed S hn hc hd ψ fun p hp => h p (.inr hp))

/-- an infectious value `z` (absorbing for all operations) infects every formula containing an atom valued `z` -/
lemma eval_infect {V : Type*} {n : V → V} {c d : V → V → V} (z : V) (hn : n z = z)
    (hc : ∀ a, c z a = z ∧ c a z = z) (hd : ∀ a, d z a = z ∧ d a z = z) {v : ℕ → V} :
    ∀ φ, (∃ p ∈ vars φ, v p = z) → eval n c d v φ = z
  | var p, ⟨q, hq, h⟩ => by simp [vars] at hq; subst hq; exact h
  | neg φ, h => by simp only [eval]; rw [eval_infect z hn hc hd φ h, hn]
  | conj φ ψ, ⟨q, hq, h⟩ => by
      simp only [vars, List.mem_append] at hq
      simp only [eval]; rcases hq with hq | hq
      · rw [eval_infect z hn hc hd φ ⟨q, hq, h⟩, (hc _).1]
      · rw [eval_infect z hn hc hd ψ ⟨q, hq, h⟩, (hc _).2]
  | disj φ ψ, ⟨q, hq, h⟩ => by
      simp only [vars, List.mem_append] at hq
      simp only [eval]; rcases hq with hq | hq
      · rw [eval_infect z hn hc hd φ ⟨q, hq, h⟩, (hd _).1]
      · rw [eval_infect z hn hc hd ψ ⟨q, hq, h⟩, (hd _).2]

end Fml

open Fml

/-! ## The four-valued matrix S_fde -/

inductive V4 | T | B | F | N
  deriving DecidableEq, Repr

instance : Fintype V4 := ⟨{.T, .B, .F, .N}, fun x => by cases x <;> simp⟩

namespace V4

def neg : V4 → V4
  | T => F | F => T | B => B | N => N

/-- LP meet on `F < B < T`, `N` infectious -/
def conj : V4 → V4 → V4
  | N, _ => N | _, N => N
  | F, _ => F | _, F => F
  | B, _ => B | _, B => B
  | T, T => T

def disj (a b : V4) : V4 := neg (conj (neg a) (neg b))

def desig (a : V4) : Prop := a = T ∨ a = B
instance : DecidablePred desig := fun a => inferInstanceAs (Decidable (a = T ∨ a = B))

end V4

abbrev ev4 (v : ℕ → V4) (φ : Fml) : V4 := eval V4.neg V4.conj V4.disj v φ

/-- S_fde consequence -/
def EntSfde (Γ : List Fml) (φ : Fml) : Prop := ∀ v, (∀ ψ ∈ Γ, (ev4 v ψ).desig) → (ev4 v φ).desig
/-- LP consequence: valuations avoiding `N` -/
def EntLP (Γ : List Fml) (φ : Fml) : Prop :=
  ∀ v : ℕ → V4, (∀ p, v p ≠ .N) → (∀ ψ ∈ Γ, (ev4 v ψ).desig) → (ev4 v φ).desig
/-- classical satisfiability -/
def SatCL (γ : Fml) : Prop := ∃ v : ℕ → V4, (∀ p, v p = .T ∨ v p = .F) ∧ ev4 v γ = .T

/-! ## 𝔼: the seven values -/

inductive V3 | T | F | M
  deriving DecidableEq, Repr

instance : Fintype V3 := ⟨{.T, .F, .M}, fun x => by cases x <;> simp⟩

namespace V3
def neg : V3 → V3 | T => F | F => T | M => M
/-- Bochvar: `M` infectious, otherwise classical -/
def conj : V3 → V3 → V3
  | M, _ => M | _, M => M
  | T, T => T | _, _ => F
def disj : V3 → V3 → V3
  | M, _ => M | _, M => M
  | F, F => F | _, _ => T
end V3

/-- the connective-wise lift to subsets -/
def neg7 (S : Finset V3) : Finset V3 := S.image V3.neg
def conj7 (S U : Finset V3) : Finset V3 := (S ×ˢ U).image fun x => V3.conj x.1 x.2
def disj7 (S U : Finset V3) : Finset V3 := (S ×ˢ U).image fun x => V3.disj x.1 x.2

abbrev ev7 (v : ℕ → Finset V3) (φ : Fml) : Finset V3 := eval neg7 conj7 disj7 v φ

/-- 𝔼-consequence: valuations into the seven nonempty subsets -/
def EntE (Γ : List Fml) (φ : Fml) : Prop :=
  ∀ v : ℕ → Finset V3, (∀ p, (v p).Nonempty) → (∀ ψ ∈ Γ, V3.T ∈ ev7 v ψ) → V3.T ∈ ev7 v φ

/-! ## L2: the collapse -/

/-- the classical part of a value -/
def collapse (S : Finset V3) : V4 :=
  if V3.T ∈ S then (if V3.F ∈ S then .B else .T) else (if V3.F ∈ S then .F else .N)

/-- a section of `collapse` with nonempty values -/
def lift : V4 → Finset V3
  | .T => {V3.T} | .F => {V3.F} | .B => {V3.T, V3.F} | .N => {V3.M}

lemma collapse_neg : ∀ S, collapse (neg7 S) = V4.neg (collapse S) := by decide
lemma collapse_conj : ∀ S U, collapse (conj7 S U) = V4.conj (collapse S) (collapse U) := by decide
lemma collapse_disj : ∀ S U, collapse (disj7 S U) = V4.disj (collapse S) (collapse U) := by decide
lemma collapse_desig : ∀ S, V3.T ∈ S ↔ (collapse S).desig := by decide
lemma collapse_lift : ∀ a, collapse (lift a) = a := by decide
lemma lift_nonempty : ∀ a, (lift a).Nonempty := by decide

lemma collapse_ev (v : ℕ → Finset V3) (φ : Fml) : collapse (ev7 v φ) = ev4 (collapse ∘ v) φ :=
  eval_hom collapse collapse_neg collapse_conj collapse_disj v φ

/-- **L2.** 𝔼 and S_fde have the same consequence relation. -/
theorem E_iff_Sfde (Γ : List Fml) (φ : Fml) : EntE Γ φ ↔ EntSfde Γ φ := by
  constructor
  · intro h w hΓ
    have := h (lift ∘ w) (fun p => lift_nonempty _) (fun ψ hψ => by
      rw [collapse_desig, collapse_ev]; convert hΓ ψ hψ; funext p; simp [collapse_lift])
    rw [collapse_desig, collapse_ev] at this; convert this; funext p; simp [collapse_lift]
  · intro h v _ hΓ
    rw [collapse_desig, collapse_ev]
    exact h _ fun ψ hψ => by rw [← collapse_ev, ← collapse_desig]; exact hΓ ψ hψ

/-! ## L3: the right variable-inclusion companion of LP -/

section L3
open V4

lemma N_infects (v : ℕ → V4) (φ : Fml) (h : ∃ p ∈ vars φ, v p = N) : ev4 v φ = N :=
  eval_infect N rfl (fun a => by cases a <;> decide) (fun a => by cases a <;> decide) φ h

lemma B_fixed (v : ℕ → V4) (φ : Fml) (h : ∀ p ∈ vars φ, v p = B) : ev4 v φ = B :=
  eval_closed (· = B) (by rintro _ rfl; rfl) (by rintro _ _ rfl rfl; rfl) (by rintro _ _ rfl rfl; rfl) φ h

lemma desig_ne_N {a : V4} (h : a.desig) : a ≠ N := by rcases h with rfl | rfl <;> decide

/-- `Var φ ⊆ Var Γ` -/
def VarSub (Γ : List Fml) (φ : Fml) : Prop := ∀ p ∈ vars φ, ∃ ψ ∈ Γ, p ∈ vars ψ

/-- **L3.** `Γ ⊨_𝔼 φ ⟺ Γ ⊨_LP φ ∧ Var φ ⊆ Var Γ`. -/
theorem L3 (Γ : List Fml) (φ : Fml) : EntE Γ φ ↔ EntLP Γ φ ∧ VarSub Γ φ := by
  rw [E_iff_Sfde]
  constructor
  · intro h
    refine ⟨fun v _ hΓ => h v hΓ, fun p hp => ?_⟩
    by_contra hno
    push Not at hno
    let v : ℕ → V4 := fun q => if q = p then N else B
    have hΓ : ∀ ψ ∈ Γ, (ev4 v ψ).desig := fun ψ hψ => by
      rw [B_fixed v ψ fun q hq => by simp [v, show q ≠ p from fun e => hno ψ hψ (e ▸ hq)]]; right; rfl
    exact desig_ne_N (h v hΓ) (N_infects v φ ⟨p, hp, by simp [v]⟩)
  · rintro ⟨hLP, hsub⟩ v hΓ
    let v' : ℕ → V4 := fun q => if v q = N then B else v q
    have hv' : ∀ p, v' p ≠ N := fun p => by by_cases h : v p = N <;> simp [v', h]
    -- atoms of designated formulas are not `N`, so `v'` agrees with `v` on them
    have agree : ∀ ψ, (ev4 v ψ).desig → ∀ p ∈ vars ψ, v' p = v p := fun ψ hd p hp => by
      have : v p ≠ N := fun hN => desig_ne_N hd (N_infects v ψ ⟨p, hp, hN⟩)
      simp [v', this]
    have hΓ' : ∀ ψ ∈ Γ, (ev4 v' ψ).desig := fun ψ hψ => by
      rw [show ev4 v' ψ = ev4 v ψ from eval_congr ψ (agree ψ (hΓ ψ hψ))]; exact hΓ ψ hψ
    have := hLP v' hv' hΓ'
    rwa [show ev4 v' φ = ev4 v φ from eval_congr φ fun p hp => by
      obtain ⟨ψ, hψ, hpψ⟩ := hsub p hp; exact agree ψ (hΓ ψ hψ) p hpψ] at this

/-! ## L4 -/

/-- **L4.** No tautologies: the all-`{M}` assignment kills every formula. -/
theorem L4_no_tautology (φ : Fml) : ¬ EntE [] φ := by
  intro h
  have := h (fun _ => {V3.M}) (fun _ => by simp) (by simp)
  rw [collapse_desig, collapse_ev] at this
  have hN : ev4 (collapse ∘ fun _ => ({V3.M} : Finset V3)) φ = N :=
    N_infects _ φ ⟨(vars φ).head (vars_ne_nil φ), List.head_mem _, show collapse {V3.M} = N by decide⟩
  rw [hN] at this; exact absurd this (by decide)

/-- **L4.** No unsatisfiable premise sets: the all-`{T, F}` assignment designates every formula. -/
theorem L4_all_satisfiable (Γ : List Fml) :
    ∃ v : ℕ → Finset V3, (∀ p, (v p).Nonempty) ∧ ∀ ψ ∈ Γ, V3.T ∈ ev7 v ψ :=
  ⟨fun _ => {V3.T, V3.F}, fun _ => by simp, fun ψ _ => by
    rw [collapse_desig, collapse_ev, B_fixed (collapse ∘ fun _ => ({V3.T, V3.F} : Finset V3)) ψ
      fun _ _ => show collapse {V3.T, V3.F} = B by decide]; right; rfl⟩

/-! ## L5: the coNP-hardness reduction -/

/-- `⋁_{p ∈ l} (p ∧ ¬p)` -/
def gluts : List ℕ → Fml
  | [] => conj (var 0) (neg (var 0))
  | [p] => conj (var p) (neg (var p))
  | p :: q :: l => disj (conj (var p) (neg (var p))) (gluts (q :: l))

lemma vars_gluts_sub : ∀ (l : List ℕ), l ≠ [] → ∀ p ∈ vars (gluts l), p ∈ l
  | [], h => absurd rfl h
  | [p], _ => by simp [gluts, vars]
  | p :: q :: l, _ => by
      intro r hr
      simp only [gluts, vars, List.mem_append, List.mem_singleton] at hr
      rcases hr with ((rfl | rfl) | hr)
      · simp
      · simp
      · exact List.mem_cons_of_mem _ (vars_gluts_sub (q :: l) (by simp) r hr)

lemma classical_closed (v : ℕ → V4) (φ : Fml) (h : ∀ p ∈ vars φ, v p = T ∨ v p = F) :
    ev4 v φ = T ∨ ev4 v φ = F :=
  eval_closed (fun a => a = T ∨ a = F) (by rintro _ (rfl | rfl) <;> decide)
    (by rintro _ _ (rfl | rfl) (rfl | rfl) <;> decide) (by rintro _ _ (rfl | rfl) (rfl | rfl) <;> decide) φ h

lemma gluts_classical (v : ℕ → V4) (hv : ∀ p, v p = T ∨ v p = F) : ∀ l, ev4 v (gluts l) = F
  | [] => by rcases hv 0 with h | h <;> simp [gluts, eval, h] <;> decide
  | [p] => by rcases hv p with h | h <;> simp [gluts, eval, h] <;> decide
  | p :: q :: l => by
      rcases hv p with h | h <;> simp only [gluts, eval, h, gluts_classical v hv (q :: l)] <;> decide

lemma aux_B : ∀ a b : V4, a = B → b ≠ N → (V4.disj (V4.conj a (V4.neg a)) b).desig := by decide
lemma aux_rest : ∀ a b : V4, a ≠ N → b.desig → (V4.disj (V4.conj a (V4.neg a)) b).desig := by decide
lemma neg_ne_N : ∀ a : V4, a ≠ N → V4.neg a ≠ N := by decide
lemma conj_ne_N : ∀ a b : V4, a ≠ N → b ≠ N → V4.conj a b ≠ N := by decide
lemma disj_ne_N : ∀ a b : V4, a ≠ N → b ≠ N → V4.disj a b ≠ N := by decide

/-- under a valuation with no `N` on `l` and some `B` on `l`, the glut disjunction is designated -/
lemma gluts_B (v : ℕ → V4) : ∀ l : List ℕ, (∀ p ∈ l, v p ≠ N) → (∃ p ∈ l, v p = B) → (ev4 v (gluts l)).desig
  | [], _, ⟨_, hp, _⟩ => by simp at hp
  | [p], _, ⟨q, hq, hB⟩ => by
      simp at hq; subst hq
      show (V4.conj (v q) (V4.neg (v q))).desig
      rw [hB]; decide
  | p :: q :: l, hN, ⟨r, hr, hB⟩ => by
      show (V4.disj (V4.conj (v p) (V4.neg (v p))) (ev4 v (gluts (q :: l)))).desig
      have hrest : ev4 v (gluts (q :: l)) ≠ N :=
        eval_closed (· ≠ N) neg_ne_N conj_ne_N disj_ne_N _ fun s hs =>
          hN s (List.mem_cons_of_mem _ (vars_gluts_sub _ (by simp) s hs))
      rcases List.mem_cons.mp hr with rfl | hr
      · exact aux_B _ _ hB hrest
      · exact aux_rest _ _ (hN p (by simp))
          (gluts_B v (q :: l) (fun s hs => hN s (List.mem_cons_of_mem _ hs)) ⟨r, hr, hB⟩)

/-- **L5.** `γ` is classically unsatisfiable iff `γ ⊨_𝔼 ⋁_{p ∈ Var γ} (p ∧ ¬p)`. -/
theorem L5 (γ : Fml) : ¬ SatCL γ ↔ EntE [γ] (gluts (vars γ)) := by
  rw [E_iff_Sfde]
  constructor
  · intro hu v hγ
    have hγ' := hγ γ (by simp)
    have hN : ∀ p ∈ vars γ, v p ≠ N := fun p hp hN => desig_ne_N hγ' (N_infects v γ ⟨p, hp, hN⟩)
    by_cases hB : ∃ p ∈ vars γ, v p = B
    · exact gluts_B v _ hN hB
    · push Not at hB
      exfalso; apply hu
      let c : ℕ → V4 := fun q => if v q = F then F else T
      have hc : ∀ p ∈ vars γ, c p = v p := fun p hp => by
        have := hN p hp; have := hB p hp
        cases h : v p <;> simp_all [c]
      refine ⟨c, fun q => by by_cases h : v q = F <;> simp [c, h], ?_⟩
      rw [show ev4 c γ = ev4 v γ from eval_congr γ hc]
      rcases classical_closed v γ (fun p hp => by rw [← hc p hp]; by_cases h : v p = F <;> simp [c, h])
        with h | h
      · exact h
      · rw [h] at hγ'; exact absurd hγ' (by decide)
  · rintro h ⟨c, hc, hT⟩
    have := h c (by simp [hT]; left; rfl)
    rw [gluts_classical c hc] at this; exact absurd this (by decide)

end L3

/-! ## Meta-1 -/

/-- the seven values -/
def vals7 : Finset (Finset V3) := Finset.univ.filter Finset.Nonempty

/-- `¬(μ ∨ ν) ∧ (μ ∧ ν)` -/
def absurdity : Fml := conj (neg (disj (var 0) (var 1))) (conj (var 0) (var 1))

/-- **Meta-1.** Designated in exactly 4 of the 49 assignments: those where `μ` and `ν` both contain T and F. -/
theorem meta1 :
    ((vals7 ×ˢ vals7).filter fun ab => V3.T ∈ ev7 (fun p => if p = 0 then ab.1 else ab.2) absurdity) =
      (vals7 ×ˢ vals7).filter fun ab => V3.T ∈ ab.1 ∧ V3.F ∈ ab.1 ∧ V3.T ∈ ab.2 ∧ V3.F ∈ ab.2 := by
  decide

theorem meta1_card :
    ((vals7 ×ˢ vals7).filter fun ab => V3.T ∈ ev7 (fun p => if p = 0 then ab.1 else ab.2) absurdity).card = 4 ∧
      (vals7 ×ˢ vals7).card = 49 := by
  decide

end EatTheFibre.Discordia.Logic
