/-
# Axiom audit

Every headline result depends only on Lean's standard axioms (`propext`, `Classical.choice`, `Quot.sound`):
no `sorry`, and no `Lean.ofReduceBool` (the kernel checks are `decide +kernel`, not `native_decide`).
The build log prints the lists below.
-/
import EatTheFibre.Examples
import EatTheFibre.Gen.two_pairs_short
import EatTheFibre.Gen.pairs_in_short

open EatTheFibre

#print axioms CoordAlg.isFold_iff_basis
#print axioms CoordAlg.IsFold.conj
#print axioms CoordAlg.vertex_block
#print axioms Reflect.check_sound
#print axioms Q_pairs_in.thm3F1_isFold
#print axioms Q_pairs_in.thm3_rigid_GT
#print axioms Q_two_pairs.thm2_isFold
#print axioms Q_two_pairs.thm4_isFold
#print axioms Q_two_pairs.thm4_rigid_GT
#print axioms Q_two_pairs_short.thm4_isFold
#print axioms Q_two_pairs_short.thm2_check_fails
#print axioms Q_pairs_in_short.thm3F1_check_fails
#print axioms CoordAlg.isFold_of_support
#print axioms CoordAlg.isFold_push
#print axioms CoordAlg.Embeds.trans
#print axioms FinQuiver.Hom.embeds
#print axioms thm4_general
#print axioms thm3_general
#print axioms thm3_general_rigid
#print axioms thm4_general_rigid
#print axioms thm4_two_pairs_short
#print axioms thm3_pairs_in_tail
