#!/bin/bash
cd /home/claude/eat-the-fibre
until [ -s logs/v2_75/census.done ]; do sleep 20; done
for Q in K3_tgtfork_long; do
  for p in 2 3; do
    o=logs/v2_75/census_${Q}_${p}
    [ -s $o.json ] && continue
    python3 pcensus.py $Q $p --engine torch --out $o.json > $o.txt 2>&1
  done
done
echo ALLDONE >> logs/v2_75/census2.done
