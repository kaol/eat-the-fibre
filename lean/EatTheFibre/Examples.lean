/-
# The general theorems applied to census quivers

The hypotheses of `thm4_general` and `thm3_general` are met by the census quivers. Finiteness of the paths comes from
the potential `h = Fin.val` (`Path.finite`).
-/
import EatTheFibre.General

open Matrix

namespace EatTheFibre

open FinQuiver CoordAlg Reflect

/-- 1⇉2⇉3 + 1→3 (`two_pairs_short`): arrows a₁ a₂ : 0 → 1, b₁ b₂ : 1 → 2, c : 0 → 2 -/
abbrev Qshort : FinQuiver := ⟨3, 5, ![0, 0, 1, 1, 0], ![1, 1, 2, 2, 2]⟩

lemma Qshort_pot : Qshort.Potential (fun v => v.val) := by unfold FinQuiver.Potential; decide

noncomputable instance : Fintype Qshort.Path :=
  haveI := Path.finite Qshort Qshort_pot; Fintype.ofFinite _

def Cshort : Consecutive Qshort where
  h := fun v => v.val
  pot := Qshort_pot
  u := 0
  v := 1
  w := 2
  a₁ := 0
  a₂ := 1
  b₁ := 2
  b₂ := 3
  ha₁ := by decide
  ha₂ := by decide
  hb₁ := by decide
  hb₂ := by decide
  na := by decide
  nb := by decide

/-- **Theorem 4 on 1⇉2⇉3 + 1→3**, where Theorem 2 fails (`Q_two_pairs_short.thm2_check_fails`). -/
theorem thm4_two_pairs_short {k : Type*} [CommRing k] (t : k) :
    (Qshort.pathAlg k).IsFold (pushM Cshort.ret (evalM Q_two_pairs.thm4P t)) :=
  thm4_general Cshort t

/-- 1⇉2⇇3 + 2→4 (`pairs_in_tail`): a₁ a₂ : 0 → 1, b₁ b₂ : 2 → 1, d : 1 → 3 -/
abbrev Qtail : FinQuiver := ⟨4, 5, ![0, 0, 2, 2, 1], ![1, 1, 1, 1, 3]⟩

/-- a potential: 0 ↦ 0, 1 ↦ 1, 2 ↦ 0, 3 ↦ 2 -/
def hTail : Fin 4 → ℕ := ![0, 1, 0, 2]

lemma Qtail_pot : Qtail.Potential hTail := by unfold FinQuiver.Potential; decide

noncomputable instance : Fintype Qtail.Path :=
  haveI := Path.finite Qtail Qtail_pot; Fintype.ofFinite _

def Ctail : Converging Qtail where
  h := hTail
  pot := Qtail_pot
  u := 0
  v := 1
  w := 2
  a₁ := 0
  a₂ := 1
  b₁ := 2
  b₂ := 3
  ha₁ := by decide
  ha₂ := by decide
  hb₁ := by decide
  hb₂ := by decide
  na := by decide
  nb := by decide
  huw := by decide
  noArrow := by decide

/-- **Theorem 3 on 1⇉2⇇3 + 2→4**: the four families are folds (v2.27 found 18 lines of vertex blocks here;
these are four of them). -/
theorem thm3_pairs_in_tail {k : Type*} [CommRing k] (t : k) :
    (Qtail.pathAlg k).IsFold (pushM Ctail.ret (evalM Q_pairs_in.thm3F1P t)) ∧
    (Qtail.pathAlg k).IsFold (pushM Ctail.ret (evalM Q_pairs_in.thm3F2P t)) ∧
    (Qtail.pathAlg k).IsFold (pushM Ctail.ret (evalM Q_pairs_in.thm3F3P t)) ∧
    (Qtail.pathAlg k).IsFold (pushM Ctail.ret (evalM Q_pairs_in.thm3F4P t)) :=
  thm3_general Ctail t

/-- 4→1⇉2⇇3 (`pairs_in_pre`): an arrow into `u` -/
abbrev Qpre : FinQuiver := ⟨4, 5, ![3, 0, 0, 2, 2], ![0, 1, 1, 1, 1]⟩

def hPre : Fin 4 → ℕ := ![1, 2, 0, 0]

lemma Qpre_pot : Qpre.Potential hPre := by unfold FinQuiver.Potential; decide

noncomputable instance : Fintype Qpre.Path :=
  haveI := Path.finite Qpre Qpre_pot; Fintype.ofFinite _

def Cpre : Converging Qpre where
  h := hPre
  pot := Qpre_pot
  u := 0
  v := 1
  w := 2
  a₁ := 1
  a₂ := 2
  b₁ := 3
  b₂ := 4
  ha₁ := by decide
  ha₂ := by decide
  hb₁ := by decide
  hb₂ := by decide
  na := by decide
  nb := by decide
  huw := by decide
  noArrow := by decide

theorem thm3_pairs_in_pre {k : Type*} [CommRing k] (t : k) :
    (Qpre.pathAlg k).IsFold (pushM Cpre.ret (evalM Q_pairs_in.thm3F1P t)) := (thm3_general Cpre t).1

end EatTheFibre
