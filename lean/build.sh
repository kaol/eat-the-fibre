#!/usr/bin/env bash
# Build the Lean project module by module (each generated file needs ~5 GB while the kernel checks run; lake would
# otherwise run several at once). Fetches Mathlib's prebuilt cache first; nothing is compiled from Mathlib sources.
set -euo pipefail
cd "$(dirname "$0")"
lake exe cache get
for m in Basic Poly Reflect Gen.pairs_in Gen.pairs_in_short Gen.two_pairs Gen.two_pairs_short Invariants Reduction Path General Examples Discordia.Fold Discordia.Completion Discordia.Exact Discordia.Logic Discordia.Ops Discordia.Bridge Audit; do
  echo "== EatTheFibre.$m"; lake build "EatTheFibre.$m"
done
lake build
