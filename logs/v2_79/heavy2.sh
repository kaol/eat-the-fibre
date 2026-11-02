#!/bin/bash
# v2.79 heavy runs for the user's machine (from the repository root, branch claude/great-wright-qp17q6):
#   nohup bash logs/v2_79/heavy2.sh JOB [JOB ...] > logs/v2_79/heavy2.log 2>&1 &
# JOB: pairid | drivers2 | pair5  (or: all = pairid drivers2 pair5). Resumable; text outputs under logs/v2_79/.
#   pairid    the v2.55 pipeline on the identity block 1000,0100,0010,0001 of 1⇛2→3⇉4 (K3_after_pair): components,
#             and for each whether it has finitely many orbits (char 0). Decides whether its census growth
#             2236 → 2240 is a family of orbits or a q = 2 effect. CPU (Singular + python), no GPU.
#   drivers2  the v2.66 driver again on K3_tgtinout and K3_srcout_after with a long pipeline timeout: their blocks
#             0001,0100,0010,0001 resp. 0000,0100,0010,0001 and 0000,0000,0000,00m1 hit the 6000 s limit in both
#             modes. RESUME keeps the finished blocks (needs logs/v2_78/drv/<Q>/ from the v2.78 run).
#   pair5     census of 1⇛2→3⇉4 at q = 5 with representatives, and the block split against q = 2, 3 (GPU or CPU).
set -u
cd "$(dirname "$0")/../.."
O=logs/v2_79; mkdir -p $O/pipe
export ETF_MEM=${ETF_MEM:-1024}
DEV=${CENSUS_DEVICE:+--device $CENSUS_DEVICE}
run () {
  case $1 in
    pairid)
      if [ -e $O/pipe/after_pair_id.done ]; then echo "skip pairid"; return; fi
      echo "== pairid $(date)"
      CERT=$O/pipe/after_pair_id.cert.jsonl python3 -u analysis/v2_55_pipeline.py K3_after_pair 1000,0100,0010,0001 \
        $O/pipe/after_pair_id > $O/pipe/after_pair_id.txt 2>&1 && touch $O/pipe/after_pair_id.done
      echo "   done $(date): $(grep -c 'finitely many orbits: True' $O/pipe/after_pair_id.txt) finite, $(grep -c 'finitely many orbits: False' $O/pipe/after_pair_id.txt) not" ;;
    drivers2)
      for Q in K3_tgtinout K3_srcout_after; do
        echo "== driver $Q $(date)"
        PIPE_TMO=200000 RESUME=1 python3 analysis/v2_66_prove.py $Q logs/v2_78/drv/$Q >> logs/v2_78/drv/$Q.prove.txt 2>&1
        echo "   done $(date): $(grep -h 'blocks [0-9]*:' logs/v2_78/drv/$Q/tally.txt | tail -1)"
      done ;;
    pair5)
      o=$O/census_K3_after_pair_5
      if [ ! -s $o.json ]; then
        echo "== census K3_after_pair q=5 $(date)"
        python3 pcensus.py K3_after_pair 5 --engine torch $DEV --out $o.json --save-reps logs/v2_78/reps_K3_after_pair_5.pkl \
          --checkpoint $O/ckpt_K3_after_pair_5.pkl >> $o.txt 2>&1
      fi
      if [ -s $o.json ] && [ -s logs/v2_78/reps_K3_after_pair_2.pkl ] && [ -s logs/v2_78/reps_K3_after_pair_3.pkl ]; then
        python3 -I analysis/v2_53_block_census.py logs/v2_78/reps_K3_after_pair_{2,3,5}.pkl > $O/K3_after_pair.blockcensus235.txt
        python3 analysis/v2_77_blockdiff.py $O/K3_after_pair.blockcensus235.txt > $O/K3_after_pair.blockdiff235.txt
        cat $O/K3_after_pair.blockdiff235.txt
      fi ;;
    all) for j in pairid drivers2 pair5; do run $j; done ;;
    *) echo "unknown job $1"; exit 2 ;;
  esac
}
for j in "$@"; do run $j; done
echo "ALL DONE $(date)"
