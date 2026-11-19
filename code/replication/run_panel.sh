#!/usr/bin/env bash
# run_panel.sh -- drive the trait-pair panel end to end.
set -u

BASE="${LAVAGATE_BASE:?set LAVAGATE_BASE}"
SUB="$BASE/SUBMISSION"
OUTROOT="$BASE/analysis/out_R/multipair"
PANEL="$SUB/results/panel"
LOGDIR="$BASE/analysis/out_R/multipair/logs"
RE_CSV="$BASE/Section 7_LAVA Analysis/multipair.pair_re.csv"
MAX_PAR=2

mkdir -p "$LOGDIR" "$PANEL"

run_group () {
  local g="$1"; shift
  echo "[$(date +%H:%M:%S)] starting group $g $*"
  Rscript "$SUB/code/replication/R1_multipair.R" "$g" "$@" \
    > "$LOGDIR/R1_$g.log" 2>&1
  echo "[$(date +%H:%M:%S)] group $g finished with status $?"
}

pair_re () {
  python - "$1" <<'PY'
import csv, sys
want = sys.argv[1]
for row in csv.DictReader(open(r"""${LAVAGATE_BASE}/Section 7_LAVA Analysis/multipair.pair_re.csv""")):
    if row["pair"] == want:
        print(row["r_e"]); break
else:
    sys.exit(f"no r_e for {want}")
PY
}

case "${1:-}" in
  pilot)
    for g in P1 P2 P3 P5; do run_group "$g" 200; done
    echo
    echo "Check each $OUTROOT/PANEL_SUMMARY_<G>.csv, then:"
    echo "  rm $OUTROOT/checkpoint_*.rds     # pilots leave checkpoints behind"
    echo "  bash run_panel.sh full"
    ;;

  full)
    if ls "$OUTROOT"/checkpoint_*.rds >/dev/null 2>&1; then
      echo "Checkpoints present in $OUTROOT."
      echo "If these are from a PILOT, delete them first or the full run resumes"
      echo "from locus 201 and stops at 200. If they are from an interrupted full"
      echo "run, keep them. Aborting so the choice is explicit."
      exit 1
    fi
    run_group P2
    for g in P1 P3 P5; do
      while [ "$(jobs -rp | wc -l)" -ge "$MAX_PAR" ]; do wait -n; done
      run_group "$g" &
    done
    wait
    echo "[$(date +%H:%M:%S)] all groups done"
    ;;

  analyse)
    cd "$SUB/code/python"
    python 14_pair_pipeline.py --selftest || {
      echo "SELFTEST FAILED -- stopping. Nothing downstream is trustworthy."; exit 1; }
    while IFS=, read -r pair t1 t2; do
      aud="$OUTROOT/$pair/attrition_audit.csv"
      [ -f "$aud" ] || { echo "skip $pair: no $aud"; continue; }
      re="$(pair_re "$pair")" || exit 1
      python 14_pair_pipeline.py --pair "$pair" --trait1 "$t1" --trait2 "$t2" \
        --r-e "$re" --audit "$aud" --outdir "$PANEL/$pair" --null-sign
    done <<'PAIRS'
P1,SCZ,BIP
P2a,MDD,BIP
P2b,MDD_noUKB,BIP_noUKB
P3,ASD,ADHD
P5,AN,ASD
PAIRS
    cp "$OUTROOT"/PANEL_SUMMARY_*.csv "$PANEL/" 2>/dev/null || true
    python 15_build_table8.py --panel "$PANEL" --group-summaries "$PANEL" \
      --out "$PANEL/table8.md"
    echo
    echo "Now paste table8.md into manuscript.md Table 8, rebuild, and re-audit:"
    echo "  (manuscript build is not part of this deposit)"
    echo "  python $SUB/code/python/16_audit_numbers.py"
    ;;

  *)
    echo "usage: bash run_panel.sh {pilot|full|analyse}"; exit 1;;
esac
