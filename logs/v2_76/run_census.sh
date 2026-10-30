#!/bin/bash
cd /home/claude/eat-the-fibre
until [ -s logs/v2_75/census3.done ]; do sleep 20; done
for job in "K3_after_fork 2" "K3_tgtfork_sq 2" "K3_after_fork 3"; do
  set -- $job; o=logs/v2_76/census_$1_$2
  [ -s $o.json ] || python3 pcensus.py $1 $2 --engine torch --out $o.json > $o.txt 2>&1
done
echo ALLDONE >> logs/v2_76/census.done
