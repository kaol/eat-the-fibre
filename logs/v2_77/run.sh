#!/bin/bash
cd /home/claude/eat-the-fibre
for job in "K3_tgtfork_long 2" "K3_tgtfork_long 3" "K3_after_fork 3"; do
  set -- $job; o=logs/v2_77/census_$1_$2
  [ -s $o.json ] || python3 pcensus.py $1 $2 --engine torch --out $o.json --save-reps logs/v2_77/reps_$1_$2.pkl > $o.txt 2>&1
done
echo ALLDONE >> logs/v2_77/census.done
