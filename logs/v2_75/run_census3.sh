#!/bin/bash
cd /home/claude/eat-the-fibre
until [ -s logs/v2_75/census2.done ]; do sleep 20; done
o=logs/v2_75/census_pair_tgtfork3_3
[ -s $o.json ] || python3 pcensus.py pair_tgtfork3 3 --engine torch --out $o.json > $o.txt 2>&1
echo ALLDONE >> logs/v2_75/census3.done
