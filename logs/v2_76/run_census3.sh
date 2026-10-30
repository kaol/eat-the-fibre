#!/bin/bash
cd /home/claude/eat-the-fibre
until [ -s logs/v2_76/census2.done ]; do sleep 20; done
o=logs/v2_76/census_K3_after_fork_3
[ -s $o.json ] || python3 pcensus.py K3_after_fork 3 --engine torch --out $o.json > $o.txt 2>&1
echo ALLDONE >> logs/v2_76/census3.done
