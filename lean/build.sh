#!/usr/bin/env bash
# Build the Lean project module by module (each generated file needs ~5 GB while the kernel checks run; lake would
# otherwise run several at once). Fetches Mathlib's prebuilt cache first; nothing is compiled from Mathlib sources.
set -euo pipefail
cd "$(dirname "$0")"
lake exe cache get
for m in Basic Poly Reflect Gen.pairs_in Gen.pairs_in_short Gen.two_pairs Gen.two_pairs_short Invariants Reduction Path General Examples Discordia.Fold Discordia.Completion Discordia.Exact Discordia.Logic Discordia.Ops Discordia.Bridge OrbitReflect Gen.orbits_thm2 Gen.orbits_thm3F2 Gen.orbits_thm3F4 Transport Gen.two_pairs_fam1 Gen.two_pairs_fam2 Gen.two_pairs_fam3 Gen.two_pairs_fam4 Gen.two_pairs_fam5 Gen.two_pairs_fam6 Gen.two_pairs_fam7 Gen.two_pairs_fam8 Gen.two_pairs_fam9 Gen.two_pairs_fam10 Gen.two_pairs_fam11 Gen.two_pairs_fam12 Gen.two_pairs_fam13 Gen.two_pairs_fam14 Gen.two_pairs_fam15 Gen.two_pairs_fam16 Gen.two_pairs_fam17 Gen.two_pairs_fam18 Gen.two_pairs_fam19 Gen.two_pairs_rigid Gen.two_pairs_disjoint Gen.two_pairs_behaviour Gen.pairs_in_behaviour AutFilt AutFiltIn Gen.pairs_in_aut Thm3Behaviours Gen.two_pairs_short_line Gen.two_pairs_short_fam1 Gen.two_pairs_short_fam2 Gen.two_pairs_short_fam3 Gen.two_pairs_short_fam4 Gen.two_pairs_short_fam5 Gen.two_pairs_short_fam6 Gen.two_pairs_short_fam7 Gen.two_pairs_short_fam8 Gen.two_pairs_short_fam9 Gen.two_pairs_short_fam10 Gen.two_pairs_short_fam11 Gen.two_pairs_short_fam12 Gen.two_pairs_short_fam13 Gen.two_pairs_short_fam14 Gen.two_pairs_short_fam15 Gen.two_pairs_short_fam16 Gen.two_pairs_short_fam17 Gen.two_pairs_short_fam18 Audit; do
  echo "== EatTheFibre.$m"; lake build "EatTheFibre.$m"
done
lake build
