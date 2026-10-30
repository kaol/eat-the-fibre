#!/bin/bash
# v2.78 heavy runs, for a machine with more memory and cores than the session container.
# Usage (from the repository root, on branch claude/great-wright-qp17q6):
#   nohup bash logs/v2_78/heavy.sh JOB [JOB ...] > logs/v2_78/heavy.log 2>&1 &
# JOB is one of: tgtfork5 afterfork afterpair long5 drivers  (or: all).
# Each job is resumable: a finished output (non-empty .json / final driver line) is skipped on a rerun.
# Needs: python3 with numpy, sympy, torch; Singular on PATH (for 'drivers').
# Outputs (all text; the *.pkl representatives are gitignored and stay on the machine):
#   logs/v2_78/census_<Q>_<q>.{json,txt}, logs/v2_78/<Q>.blockcensus.txt, logs/v2_78/<Q>.blockdiff_<a>_<b>.txt,
#   logs/v2_78/drv/<Q>.prove.txt (+ logs/v2_78/drv/<Q>/ without scratch/)
# Tunables: ETF_MEM (MiB per torch working set; default here 8192), OMP/torch threads follow the machine.
set -u
cd "$(dirname "$0")/../.."
export ETF_MEM=${ETF_MEM:-8192}
O=logs/v2_78
census () {   # QUIVER q
  local o=$O/census_$1_$2
  if [ -s $o.json ]; then echo "skip $o"; return; fi
  echo "== census $1 q=$2 $(date)"
  python3 pcensus.py $1 $2 --engine torch --out $o.json --save-reps $O/reps_$1_$2.pkl \
    --checkpoint $O/ckpt_$1_$2.pkl > $o.txt 2>&1
  echo "   done $(date): $(python3 -c "import json; print(json.load(open('$o.json'))['orbits_GxTloc'])")"
}
blocks () {   # QUIVER q1 q2 [q3]
  local Q=$1; shift; local args=""
  for q in "$@"; do args="$args $O/reps_${Q}_$q.pkl"; done
  python3 -I analysis/v2_53_block_census.py $args > $O/$Q.blockcensus.txt
  # pairwise differences (blocks identified mod 2 when q = 2 is involved)
  python3 analysis/v2_77_blockdiff.py $O/$Q.blockcensus.txt > $O/$Q.blockdiff.txt || true
}
driver () {   # QUIVER
  mkdir -p $O/drv
  if [ -s $O/drv/$1.prove.txt ] && tail -1 $O/drv/$1.prove.txt | grep -q '"quiver"'; then echo "skip driver $1"; return; fi
  echo "== driver $1 $(date)"
  RESUME=1 python3 analysis/v2_66_prove.py $1 $O/drv/$1 >> $O/drv/$1.prove.txt 2>&1
  echo "   done $(date)"
}
run () {
  case $1 in
    tgtfork5)  census K3_tgtfork 2; census K3_tgtfork 3; census K3_tgtfork 5; blocks K3_tgtfork 2 3 5 ;;
    afterfork) census K3_after_fork 2; census K3_after_fork 3; blocks K3_after_fork 2 3 ;;
    afterpair) census K3_after_pair 2; census K3_after_pair 3; blocks K3_after_pair 2 3 ;;
    long5)     census K3_tgtfork_long 2; census K3_tgtfork_long 3; census K3_tgtfork_long 5; blocks K3_tgtfork_long 2 3 5 ;;
    drivers)   for Q in K3_mid K3_srcboth K3_tgtinout K3_srcout_after; do driver $Q; done ;;
    all)       for j in afterfork afterpair tgtfork5 drivers long5; do run $j; done ;;
    *) echo "unknown job $1"; exit 2 ;;
  esac
}
for j in "$@"; do run $j; done
echo "ALL DONE $(date)"
