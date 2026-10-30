#!/bin/bash
# v2.75 censuses at q = 2, 3 (torch engine), sequential; skips finished ones
cd /home/claude/eat-the-fibre
for Q in K3_tgtpair K3_tgtinout K3_srcout_after pair_tgtfork3 K3_tgtfork3; do
  for p in 2 3; do
    o=logs/v2_75/census_${Q}_${p}
    [ -s $o.json ] && continue
    python3 pcensus.py $Q $p --engine torch --out $o.json > $o.txt 2>&1
  done
done
echo ALLDONE >> logs/v2_75/census.done
