/-
# Axiom audit

Every headline result depends only on Lean's standard axioms (`propext`, `Classical.choice`, `Quot.sound`):
no `sorry`, and no `Lean.ofReduceBool` (the kernel checks are `decide +kernel`, not `native_decide`).
The build log prints the lists below.
-/
import EatTheFibre.Examples
import EatTheFibre.Gen.two_pairs_short
import EatTheFibre.Gen.pairs_in_short
import EatTheFibre.Discordia.Bridge
import EatTheFibre.Discordia.Exact
import EatTheFibre.Gen.orbits_thm2
import EatTheFibre.Transport
import EatTheFibre.Gen.two_pairs_fam1
import EatTheFibre.Gen.two_pairs_fam2
import EatTheFibre.Gen.two_pairs_fam3
import EatTheFibre.Gen.two_pairs_fam4
import EatTheFibre.Gen.two_pairs_fam5
import EatTheFibre.Gen.two_pairs_fam6
import EatTheFibre.Gen.two_pairs_fam7
import EatTheFibre.Gen.two_pairs_fam8
import EatTheFibre.Gen.two_pairs_fam9
import EatTheFibre.Gen.two_pairs_fam10
import EatTheFibre.Gen.two_pairs_fam11
import EatTheFibre.Gen.two_pairs_fam12
import EatTheFibre.Gen.two_pairs_fam13
import EatTheFibre.Gen.two_pairs_fam14
import EatTheFibre.Gen.two_pairs_fam15
import EatTheFibre.Gen.two_pairs_fam16
import EatTheFibre.Gen.two_pairs_fam17
import EatTheFibre.Gen.two_pairs_fam18
import EatTheFibre.Gen.two_pairs_fam19
import EatTheFibre.Gen.two_pairs_rigid
import EatTheFibre.Gen.two_pairs_disjoint
import EatTheFibre.Gen.two_pairs_behaviour
import EatTheFibre.Gen.pairs_in_behaviour
import EatTheFibre.AutFilt
import EatTheFibre.Gen.orbits_thm3F2
import EatTheFibre.Gen.orbits_thm3F4
import EatTheFibre.Discordia.Logic
import EatTheFibre.Discordia.Ops

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

-- Discordianism v2, derived (v2.31)
#print axioms Discordia.M4_split
#print axioms Discordia.M4_rotaBaxter
#print axioms Discordia.dL_A1
#print axioms Discordia.dL1_leibniz
#print axioms Discordia.E1_left
#print axioms Discordia.dD_3
#print axioms Discordia.dD1_split
#print axioms Discordia.dD1_hom
#print axioms Discordia.E1star
#print axioms Discordia.M3_isRB
#print axioms Discordia.Completion.baxter
#print axioms Discordia.Completion.line1
#print axioms Discordia.Completion.point0
#print axioms Discordia.Completion.line2
#print axioms Discordia.Completion.pointB
#print axioms Discordia.Completion.line0
#print axioms Discordia.Completion.lineB
#print axioms Discordia.Completion.endo_bimod
#print axioms CoordAlg.isFold_iff_rb
#print axioms CoordAlg.fold_star
#print axioms thm4_isRB
#print axioms thm3_isRB
#print axioms Q_two_pairs.thm4_not_dialgebra
#print axioms Q_two_pairs.thm4_not_multiplicative
#print axioms Q_two_pairs.thm4_not_module
#print axioms Discordia.Logic.E_iff_Sfde
#print axioms Discordia.Logic.L3
#print axioms Discordia.Logic.L4_no_tautology
#print axioms Discordia.Logic.L4_all_satisfiable
#print axioms Discordia.Logic.L5
#print axioms Discordia.Logic.meta1
#print axioms Discordia.Logic.meta1_card
#print axioms Discordia.Ops.M2
#print axioms Discordia.Ops.N5
#print axioms Discordia.Ops.N6_nonassoc
#print axioms Discordia.Ops.N6_no_unit
#print axioms Discordia.Ops.β4
#print axioms Discordia.Ops.B4
#print axioms Discordia.Ops.B4_no_neutral_chaos
#print axioms Discordia.Ops.B1_chaos_uncountable
#print axioms Discordia.Ops.B1_named_null
#print axioms Discordia.Ops.B3_periods
#print axioms Discordia.Ops.hand
#print axioms Discordia.Ops.calendar

-- δ-D2 exactly (v2.32)
#print axioms Discordia.Completion.mult_iff
#print axioms Discordia.Completion.point0_iff
#print axioms Discordia.Completion.line1_iff
#print axioms Discordia.Completion.pointB_iff
#print axioms Discordia.Completion.line2_iff
#print axioms Discordia.Completion.line0_iff
#print axioms Discordia.Completion.lineB_iff

-- the polynomial families are single inner orbits (v2.33)
#print axioms Reflect.orbit_sound
#print axioms Q_two_pairs.thm2_one_orbit
#print axioms Q_pairs_in.thm3F2_one_orbit
#print axioms Q_pairs_in.thm3F4_one_orbit

-- the fibre over a behaviour is a set of coproduct twists (v2.34)
#print axioms CoordAlg.twist_delta
#print axioms CoordAlg.transport

-- all nineteen families of 1⇉2⇉3 (v2.36)
#print axioms Q_two_pairs.fam1_isFold
#print axioms Q_two_pairs.fam1_one_orbit
#print axioms Q_two_pairs.fam2_isFold
#print axioms Q_two_pairs.fam2_one_orbit
#print axioms Q_two_pairs.fam3_isFold
#print axioms Q_two_pairs.fam3_one_orbit
#print axioms Q_two_pairs.fam4_isFold
#print axioms Q_two_pairs.fam4_one_orbit
#print axioms Q_two_pairs.fam5_isFold
#print axioms Q_two_pairs.fam5_one_orbit
#print axioms Q_two_pairs.fam6_isFold
#print axioms Q_two_pairs.fam6_one_orbit
#print axioms Q_two_pairs.fam7_isFold
#print axioms Q_two_pairs.fam7_one_orbit
#print axioms Q_two_pairs.fam8_isFold
#print axioms Q_two_pairs.fam8_one_orbit
#print axioms Q_two_pairs.fam9_isFold
#print axioms Q_two_pairs.fam9_one_orbit
#print axioms Q_two_pairs.fam10_isFold
#print axioms Q_two_pairs.fam10_one_orbit
#print axioms Q_two_pairs.fam11_isFold
#print axioms Q_two_pairs.fam11_one_orbit
#print axioms Q_two_pairs.fam12_isFold
#print axioms Q_two_pairs.fam12_one_orbit
#print axioms Q_two_pairs.fam13_isFold
#print axioms Q_two_pairs.fam13_one_orbit
#print axioms Q_two_pairs.fam14_isFold
#print axioms Q_two_pairs.fam14_one_orbit
#print axioms Q_two_pairs.fam15_isFold
#print axioms Q_two_pairs.fam15_one_orbit
#print axioms Q_two_pairs.fam16_isFold
#print axioms Q_two_pairs.fam16_one_orbit
#print axioms Q_two_pairs.fam17_isFold
#print axioms Q_two_pairs.fam17_one_orbit
#print axioms Q_two_pairs.fam18_isFold
#print axioms Q_two_pairs.fam18_one_orbit
#print axioms Q_two_pairs.fam19_isFold
#print axioms Q_two_pairs.fam19_one_orbit
#print axioms Q_two_pairs.fam1_rigid_GT
#print axioms Q_two_pairs.fam2_rigid_GT
#print axioms Q_two_pairs.fam3_rigid_GT
#print axioms Q_two_pairs.fam4_rigid_GT
#print axioms Q_two_pairs.fam5_rigid_GT
#print axioms Q_two_pairs.fam6_rigid_GT
#print axioms Q_two_pairs.fam7_rigid_GT
#print axioms Q_two_pairs.fam8_rigid_GT
#print axioms Q_two_pairs.fam9_rigid_GT
#print axioms Q_two_pairs.fam10_rigid_GT
#print axioms Q_two_pairs.fam11_rigid_GT
#print axioms Q_two_pairs.fam12_rigid_GT
#print axioms Q_two_pairs.fam13_rigid_GT
#print axioms Q_two_pairs.fam14_rigid_GT
#print axioms Q_two_pairs.fam15_rigid_GT
#print axioms Q_two_pairs.fam16_rigid_GT
#print axioms Q_two_pairs.fam17_rigid_GT
#print axioms Q_two_pairs.fam18_rigid_GT
#print axioms Q_two_pairs.fam19_rigid_GT
#print axioms Q_two_pairs.fam_orbits_injective
#print axioms Q_two_pairs.innM_conj
#print axioms Q_two_pairs.fam_behaviours_distinct
#print axioms Q_pairs_in.innM_conj
#print axioms Q_pairs_in.thm3_behaviours_distinct
#print axioms Q_two_pairs.Mult.filt
#print axioms Q_two_pairs.fam_behaviours_distinct_aut
