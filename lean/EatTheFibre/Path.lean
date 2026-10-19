/-
# Path algebras of finite quivers

A finite quiver has vertices `Fin nv` and arrows `Fin na` with source `s` and target `t`. A path is a start vertex
together with a list of arrows, each starting where the previous one ends (vertices are the paths with no arrows).
The path algebra has these paths as basis:
* `b_p · b_q = b_{pq}` when `p` ends where `q` starts, and `0` otherwise;
* arrow deletion `Δ(b_r) = Σ b_p ⊗ b_q` over the ways of writing `r = p · a · q` with `a` an arrow.

These are the census engine's conventions (`gcensus.build`).

An **acyclicity potential** `h` (`h (s a) < h (t a)` for every arrow) bounds path lengths, so the paths are finite
(`Path.finite`). An **injective quiver morphism** `φ : Q₁ → Q₂` maps paths to paths. The image is closed under products,
contains every prefix and suffix of its paths, and has Δ local to it. So `φ` gives an `Embeds` of path algebras
(`Hom.embeds`), and folds can be pushed along it (`Reduction.lean`).
-/
import EatTheFibre.Reduction

namespace EatTheFibre

structure FinQuiver where
  nv : ℕ
  na : ℕ
  s : Fin na → Fin nv
  t : Fin na → Fin nv

namespace FinQuiver

variable (Q : FinQuiver)

/-- `l` is a walk starting at `x` -/
def Walk : Fin Q.nv → List (Fin Q.na) → Prop
  | _, [] => True
  | x, a :: l => Q.s a = x ∧ Walk (Q.t a) l

instance decWalk : ∀ x l, Decidable (Q.Walk x l)
  | _, [] => isTrue trivial
  | x, a :: l => by unfold Walk; haveI := decWalk (Q.t a) l; infer_instance

/-- where a walk ends -/
def endv : Fin Q.nv → List (Fin Q.na) → Fin Q.nv
  | x, [] => x
  | _, a :: l => endv (Q.t a) l

lemma walk_append : ∀ (x : Fin Q.nv) (l l' : List (Fin Q.na)),
    Q.Walk x (l ++ l') ↔ Q.Walk x l ∧ Q.Walk (Q.endv x l) l'
  | x, [], l' => by simp [Walk, endv]
  | x, a :: l, l' => by simp only [List.cons_append, Walk, endv, walk_append (Q.t a) l l', and_assoc]

lemma endv_append : ∀ (x : Fin Q.nv) (l l' : List (Fin Q.na)), Q.endv x (l ++ l') = Q.endv (Q.endv x l) l'
  | x, [], l' => rfl
  | x, a :: l, l' => endv_append (Q.t a) l l'

structure Path where
  src : Fin Q.nv
  arrs : List (Fin Q.na)
  ok : Q.Walk src arrs
deriving DecidableEq

variable {Q}

lemma Path.ext' {p q : Q.Path} (h1 : p.src = q.src) (h2 : p.arrs = q.arrs) : p = q := by
  cases p; cases q; simp_all

/-- the end vertex of a path -/
def Path.tgt (p : Q.Path) : Fin Q.nv := Q.endv p.src p.arrs

/-- concatenation (when `p` ends where `q` starts) -/
def Path.comp (p q : Q.Path) (h : p.tgt = q.src) : Q.Path :=
  ⟨p.src, p.arrs ++ q.arrs, (Q.walk_append _ _ _).mpr ⟨p.ok, by rw [show Q.endv p.src p.arrs = q.src from h]; exact q.ok⟩⟩

variable (Q)

/-- structure constants of the path algebra: `b_p b_q = b_r` iff `r = p q` -/
def cZ (p q r : Q.Path) : ℤ := if p.tgt = q.src ∧ r.src = p.src ∧ r.arrs = p.arrs ++ q.arrs then 1 else 0

/-- arrow deletion: `Δ b_r ∋ b_p ⊗ b_q` iff `r = p a q` for an arrow `a` -/
def DZ (r p q : Q.Path) : ℤ :=
  if ∃ a, p.tgt = Q.s a ∧ Q.t a = q.src ∧ r.src = p.src ∧ r.arrs = p.arrs ++ a :: q.arrs then 1 else 0

/-- the path algebra, over any commutative ring -/
def pathAlg (k : Type*) [CommRing k] : CoordAlg k Q.Path := ⟨fun p q r => (Q.cZ p q r : k), fun r p q => (Q.DZ r p q : k)⟩

/-! ## Finiteness from an acyclicity potential -/

/-- `h` strictly increases along every arrow -/
def Potential (h : Fin Q.nv → ℕ) : Prop := ∀ a, h (Q.s a) < h (Q.t a)

lemma walk_length {h : Fin Q.nv → ℕ} (hp : Q.Potential h) :
    ∀ (x : Fin Q.nv) (l : List (Fin Q.na)), Q.Walk x l → h x + l.length ≤ h (Q.endv x l)
  | x, [], _ => by simp [endv]
  | x, a :: l, hw => by
    obtain ⟨hs, hw'⟩ := hw
    have := walk_length hp (Q.t a) l hw'
    have ha := hp a; rw [hs] at ha
    simp only [List.length_cons, endv]; omega

lemma Path.finite {h : Fin Q.nv → ℕ} (hp : Q.Potential h) : Finite Q.Path := by
  let M := Finset.univ.sup h
  have hb : ∀ p : Q.Path, p.arrs.length ≤ M := fun p => by
    have := Q.walk_length hp p.src p.arrs p.ok
    have hM : h (Q.endv p.src p.arrs) ≤ M := Finset.le_sup (Finset.mem_univ _)
    omega
  have : Finite {l : List (Fin Q.na) | l.length ≤ M} := (List.finite_length_le _ M).to_subtype
  exact Finite.of_injective (fun p : Q.Path => (p.src, (⟨p.arrs, hb p⟩ : {l : List (Fin Q.na) | l.length ≤ M})))
    (fun p q hpq => by
      simp only [Prod.mk.injEq, Subtype.mk.injEq] at hpq
      exact Path.ext' hpq.1 hpq.2)

/-! ## Quiver morphisms -/

/-- a morphism of quivers, injective on vertices and on arrows -/
structure Hom (Q₁ Q₂ : FinQuiver) where
  fv : Fin Q₁.nv → Fin Q₂.nv
  fa : Fin Q₁.na → Fin Q₂.na
  hs : ∀ a, Q₂.s (fa a) = fv (Q₁.s a)
  ht : ∀ a, Q₂.t (fa a) = fv (Q₁.t a)
  injv : Function.Injective fv
  inja : Function.Injective fa

namespace Hom

variable {Q₁ Q₂ : FinQuiver} (φ : Hom Q₁ Q₂)

lemma walk_map : ∀ (x : Fin Q₁.nv) (l : List (Fin Q₁.na)), Q₂.Walk (φ.fv x) (l.map φ.fa) ↔ Q₁.Walk x l
  | x, [] => by simp [Walk]
  | x, a :: l => by
    simp only [List.map_cons, Walk, φ.hs, φ.ht, walk_map (Q₁.t a) l, φ.injv.eq_iff]

lemma endv_map : ∀ (x : Fin Q₁.nv) (l : List (Fin Q₁.na)), Q₂.endv (φ.fv x) (l.map φ.fa) = φ.fv (Q₁.endv x l)
  | x, [] => rfl
  | x, a :: l => by simp only [List.map_cons, endv, φ.ht, endv_map (Q₁.t a) l]

/-- the image of a path -/
def mapPath (p : Q₁.Path) : Q₂.Path := ⟨φ.fv p.src, p.arrs.map φ.fa, (φ.walk_map _ _).mpr p.ok⟩

lemma mapPath_tgt (p : Q₁.Path) : (φ.mapPath p).tgt = φ.fv p.tgt := φ.endv_map _ _

lemma mapPath_injective : Function.Injective φ.mapPath := fun p q h => by
  simp only [mapPath, Path.mk.injEq] at h
  exact Path.ext' (φ.injv h.1) ((List.map_injective_iff.mpr φ.inja) h.2)

open Classical in
/-- the partial inverse of `mapPath` -/
noncomputable def pre (x : Q₂.Path) : Option Q₁.Path := if h : ∃ p, φ.mapPath p = x then some h.choose else none

lemma pre_map (p : Q₁.Path) : φ.pre (φ.mapPath p) = some p := by
  have h : ∃ q, φ.mapPath q = φ.mapPath p := ⟨p, rfl⟩
  simp only [pre, dif_pos h, Option.some.injEq]; exact φ.mapPath_injective h.choose_spec

lemma map_of_pre {x : Q₂.Path} {p : Q₁.Path} (h : φ.pre x = some p) : φ.mapPath p = x := by
  unfold pre at h; split at h
  · next hx => simp only [Option.some.injEq] at h; rw [← h]; exact hx.choose_spec
  · simp at h

lemma pre_eq_none {x : Q₂.Path} : φ.pre x = none ↔ ∀ p, φ.mapPath p ≠ x := by
  constructor
  · intro h p hp; rw [← hp, φ.pre_map] at h; simp at h
  · intro h; unfold pre; rw [dif_neg (fun ⟨p, hp⟩ => h p hp)]

/-- a path whose arrows are images, starting at an image vertex, is an image -/
lemma mem_range {x : Q₂.Path} {v : Fin Q₁.nv} {l : List (Fin Q₁.na)} (hv : x.src = φ.fv v) (hl : x.arrs = l.map φ.fa) :
    ∃ p, φ.mapPath p = x := by
  have hw : Q₁.Walk v l := (φ.walk_map v l).mp (by rw [← hv, ← hl]; exact x.ok)
  exact ⟨⟨v, l, hw⟩, Path.ext' hv.symm hl.symm⟩

lemma cZ_map (p q r : Q₁.Path) : Q₂.cZ (φ.mapPath p) (φ.mapPath q) (φ.mapPath r) = Q₁.cZ p q r := by
  simp only [cZ]; rw [mapPath_tgt]
  simp only [mapPath, φ.injv.eq_iff, ← List.map_append, (List.map_injective_iff.mpr φ.inja).eq_iff]

lemma DZ_map (r p q : Q₁.Path) : Q₂.DZ (φ.mapPath r) (φ.mapPath p) (φ.mapPath q) = Q₁.DZ r p q := by
  simp only [DZ]; rw [mapPath_tgt]
  simp only [mapPath]
  congr 1
  apply propext; constructor
  · rintro ⟨a, h1, h2, h3, h4⟩
    -- the deleted arrow is an image arrow
    rw [List.map_eq_append_iff] at h4
    obtain ⟨l1, l2, hl, hm1, hm2⟩ := h4
    rw [List.map_eq_cons_iff] at hm2
    obtain ⟨a₀, l3, hl2, ha, hm3⟩ := hm2
    subst ha
    refine ⟨a₀, ?_, ?_, φ.injv h3, ?_⟩
    · rw [φ.hs] at h1; exact φ.injv h1
    · rw [φ.ht] at h2; exact φ.injv h2
    · have inj := List.map_injective_iff.mpr φ.inja
      rw [hl, hl2, inj hm1, inj hm3]
  · rintro ⟨a, h1, h2, h3, h4⟩
    exact ⟨φ.fa a, by rw [h1, φ.hs], by rw [φ.ht, h2], by rw [h3], by rw [h4]; simp⟩

/-- **an injective quiver morphism embeds path algebras** -/
theorem embeds (k : Type*) [CommRing k] [Fintype Q₁.Path] [Fintype Q₂.Path] :
    CoordAlg.Embeds (Q₂.pathAlg k) (Q₁.pathAlg k) φ.mapPath φ.pre where
  re := φ.pre_map
  er _ _ h := φ.map_of_pre h
  c_emb i j m := by simp [pathAlg, φ.cZ_map]
  closed i j m hm := by
    simp only [pathAlg, cZ]
    split_ifs with h
    · exfalso
      obtain ⟨h1, h2, h3⟩ := h
      rw [mapPath_tgt] at h1
      obtain ⟨p, hp⟩ := φ.mem_range (v := i.src) (l := i.arrs ++ j.arrs) (by rw [h2]; rfl)
        (by rw [h3]; simp [mapPath])
      rw [φ.pre_eq_none] at hm; exact hm p hp
    · simp
  ideal x y m hxy hc := by
    simp only [pathAlg, cZ] at hc
    split_ifs at hc with h
    · obtain ⟨h1, h2, h3⟩ := h
      rw [φ.pre_eq_none]; intro r hr
      -- x and y are a prefix and a suffix of an image path, hence images
      have ha : m.arrs = r.arrs.map φ.fa := by rw [← hr]; rfl
      rw [ha, List.map_eq_append_iff] at h3
      obtain ⟨l1, l2, _, hx, hy⟩ := h3
      have hxs : x.src = φ.fv r.src := by rw [← h2, ← hr]; rfl
      obtain ⟨p, hp⟩ := φ.mem_range hxs hx.symm
      have hys : y.src = φ.fv (Q₁.endv r.src l1) := by
        rw [← h1, Path.tgt, hxs, ← hx, φ.endv_map]
      obtain ⟨q, hq⟩ := φ.mem_range hys hy.symm
      rcases hxy with hn | hn
      · rw [φ.pre_eq_none] at hn; exact hn p hp
      · rw [φ.pre_eq_none] at hn; exact hn q hq
    · simp at hc
  D_emb l a b := by simp [pathAlg, φ.DZ_map]
  D_local l x y hxy := by
    simp only [pathAlg, DZ]
    split_ifs with h
    · exfalso
      obtain ⟨a, h1, h2, h3, h4⟩ := h
      have hm : (φ.mapPath l).arrs = l.arrs.map φ.fa := rfl
      rw [hm, List.map_eq_append_iff] at h4
      obtain ⟨l1, l2, hl, hx, hy⟩ := h4
      rw [List.map_eq_cons_iff] at hy
      obtain ⟨a₀, l3, _, ha, hy⟩ := hy
      have hxs : x.src = φ.fv l.src := by rw [← h3]; rfl
      obtain ⟨p, hp⟩ := φ.mem_range hxs hx.symm
      have hys : y.src = φ.fv (Q₁.t a₀) := by rw [← h2, ← ha, φ.ht]
      obtain ⟨q, hq⟩ := φ.mem_range hys hy.symm
      rcases hxy with hn | hn
      · rw [φ.pre_eq_none] at hn; exact hn p hp
      · rw [φ.pre_eq_none] at hn; exact hn q hq
    · simp

end Hom

end FinQuiver

end EatTheFibre
