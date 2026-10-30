#!/bin/bash
cd /home/claude/eat-the-fibre
until [ -s logs/v2_76/census.done ]; do sleep 20; done
for job in "K3_after_pair 2" "K3_after_pair 3" "K3_tgtfork_sq 3"; do
  set -- $job; o=logs/v2_76/census_$1_$2
  [ -s $o.json ] || python3 pcensus.py $1 $2 --engine torch --out $o.json > $o.txt 2>&1
done
echo ALLDONE >> logs/v2_76/census2.done
