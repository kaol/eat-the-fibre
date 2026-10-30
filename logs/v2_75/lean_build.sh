#!/bin/bash
export PATH=$HOME/.elan/bin:$PATH
cd /home/claude/eat-the-fibre/lean
for m in Gen.K3_tgtfork Gen.K3_tgtfork_rigid; do
  echo "== $m $(date +%T)"; lake build "EatTheFibre.$m" 2>&1 | grep -v "^\s*$" | tail -40
done
echo LEANDONE
