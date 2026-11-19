#!/usr/bin/env bash
# resume_panel.sh -- continue the panel from wherever it stopped.
set -u

BASE="${LAVAGATE_BASE:?set LAVAGATE_BASE}"
SUB="$BASE/SUBMISSION"
ADIR="$BASE/Section 7_LAVA Analysis"
OUTROOT="$BASE/analysis/out_R/multipair"
LOGDIR="$OUTROOT/logs"
mkdir -p "$LOGDIR"

say () { echo "[$(date '+%F %T')] $*"; }

if [ ! -f "$ADIR/multipair.sample.overlap.txt" ] || [ ! -f "$ADIR/multipair.info.txt" ]; then
  say "Panel inputs are missing. Run overnight.sh instead (it builds them)."
  exit 1
fi

done_group () {
  [ -f "$OUTROOT/PANEL_SUMMARY_$1.csv" ] && [ ! -f "$OUTROOT/checkpoint_$1.rds" ]
}

run_group () {
  local g="$1"
  if done_group "$g"; then say "  group $g already complete -- skipping"; return 0; fi
  if [ -f "$OUTROOT/checkpoint_$g.rds" ]; then
    say "  group $g resuming from its checkpoint"
  else
    say "  group $g starting from the beginning"
  fi
  Rscript "$SUB/code/replication/R1_multipair.R" "$g" >> "$LOGDIR/R1_$g.log" 2>&1
  say "  group $g finished"
}

say "RESUMING PANEL"
for g in P2 P1 P3 P5; do
  if done_group "$g"; then say "  $g: complete"; else
    if [ -f "$OUTROOT/checkpoint_$g.rds" ]; then say "  $g: partial (checkpoint present)"
    else say "  $g: not started"; fi
  fi
done

run_group P2

say "  groups P1 and P3 -- concurrently"
run_group P1 & PID1=$!
run_group P3 & PID3=$!
wait $PID1 $PID3

run_group P5

say "PANEL COMPLETE"
for f in "$OUTROOT"/PANEL_SUMMARY_*.csv; do
  [ -f "$f" ] && { echo "--- $f"; cat "$f"; }
done
say "Next: bash $SUB/code/replication/run_panel.sh analyse"
