#!/usr/bin/env bash
# overnight.sh -- rebuild the panel inputs, verify them, then run all four
set -u

BASE="${LAVAGATE_BASE:?set LAVAGATE_BASE}"
SUB="$BASE/SUBMISSION"
ADIR="$BASE/Section 7_LAVA Analysis"
OUTROOT="$BASE/analysis/out_R/multipair"
LOGDIR="$OUTROOT/logs"
mkdir -p "$LOGDIR"

say () { echo "[$(date '+%F %T')] $*"; }

if ls "$OUTROOT"/checkpoint_*.rds >/dev/null 2>&1; then
  say "REFUSING TO RUN. Checkpoints exist in $OUTROOT:"
  ls -la "$OUTROOT"/checkpoint_*.rds
  say "This script would delete them and restart those passes from locus 1."
  say "To continue an interrupted panel:  bash resume_panel.sh"
  say "To genuinely start over:           rm $OUTROOT/checkpoint_*.rds  then rerun"
  exit 1
fi

say "STEP 1/5  re-harmonising the panel (both blocks of the split files)"
Rscript "$SUB/code/R/01b_harmonize_panel.R" > "$LOGDIR/01b_harmonize.log" 2>&1
if ! grep -q "^DONE\.$" "$LOGDIR/01b_harmonize.log"; then
  say "FAILED: harmonise did not reach DONE. See $LOGDIR/01b_harmonize.log"
  tail -20 "$LOGDIR/01b_harmonize.log"; exit 1
fi
grep -E "MALFORMED|wrote    :|N        :" "$LOGDIR/01b_harmonize.log" || true
for T in BIP_noUKB MDD; do
  R=$(grep -A3 "^===== $T =====" "$LOGDIR/01b_harmonize.log" | grep -oP "wrote    : \K[0-9]+" || echo 0)
  say "  $T harmonised rows: $R"
done

say "STEP 2/5  discarding stale munged files for MDD and BIP_noUKB"
rm -f "$ADIR/munged/MDD.sumstats.gz" "$ADIR/munged/BIP_noUKB.sumstats.gz"

say "STEP 3/5  cross-trait LDSC over all eight phenotypes"
Rscript "$SUB/code/R/00b_sample_overlap_ldsc_panel.R" > "$LOGDIR/00b_ldsc.log" 2>&1
if ! grep -q "^DONE\.$" "$LOGDIR/00b_ldsc.log"; then
  say "FAILED: LDSC did not reach DONE. See $LOGDIR/00b_ldsc.log"
  tail -20 "$LOGDIR/00b_ldsc.log"; exit 1
fi
if grep -q "Regression test PASSES" "$LOGDIR/00b_ldsc.log"; then
  say "  LDSC regression test PASSES (ASD x SCZ reproduces the published r_e)"
else
  say "FAILED: LDSC no longer reproduces the published ASD x SCZ overlap. STOPPING."
  grep -A6 "ASD x SCZ r_e" "$LOGDIR/00b_ldsc.log" || true
  exit 1
fi
sed -n '/PER-PAIR SAMPLE OVERLAP/,/^$/p' "$LOGDIR/00b_ldsc.log" || true
cp "$ADIR"/multipair.pair_re.csv "$ADIR"/multipair.sample.overlap.txt \
   "$ADIR"/multipair.info.txt "$ADIR"/multipair.provenance.txt \
   "$SUB/results/ldsc_panel/" 2>/dev/null || true

say "STEP 4/5  piloting R1 on the P2 group (50 loci) to confirm the holes are gone"
rm -f "$OUTROOT"/checkpoint_*.rds
Rscript "$SUB/code/replication/R1_multipair.R" P2 50 > "$LOGDIR/pilot_P2.log" 2>&1
NFAIL=$(grep -c "none of specified SNP IDs are present" "$LOGDIR/pilot_P2.log" || true)
NPROC=$(grep -oP 'processed\s+\K[0-9]+' "$LOGDIR/pilot_P2.log" | head -1 || echo 0)
SNPS=$(grep -oP 'SNPs common to.*: \K[0-9]+' "$LOGDIR/pilot_P2.log" | head -1 || echo 0)
say "  pilot: $SNPS common SNPs, $NPROC/50 loci processed, $NFAIL reference misses"
if [ "${NFAIL:-99}" -gt 5 ]; then
  say "FAILED: still losing loci to missing reference SNPs. NOT launching the panel."
  say "        See $LOGDIR/pilot_P2.log"
  exit 1
fi
if [ "${SNPS:-0}" -lt 3500000 ]; then
  say "FAILED: SNP intersection $SNPS is far below the ~4-5 M expected once both"
  say "        blocks are read. NOT launching the panel. See $LOGDIR/01b_harmonize.log"
  exit 1
fi
say "  pilot looks sane -- launching the full panel"
rm -f "$OUTROOT"/checkpoint_*.rds

say "STEP 5/5  full panel"
say "  group P2 (MDD, BIP, MDD_noUKB, BIP_noUKB) -- alone"
Rscript "$SUB/code/replication/R1_multipair.R" P2 > "$LOGDIR/R1_P2.log" 2>&1
say "  group P2 finished (exit $?)"

say "  groups P1 and P3 -- concurrently"
Rscript "$SUB/code/replication/R1_multipair.R" P1 > "$LOGDIR/R1_P1.log" 2>&1 &
PID1=$!
Rscript "$SUB/code/replication/R1_multipair.R" P3 > "$LOGDIR/R1_P3.log" 2>&1 &
PID3=$!
wait $PID1; say "  group P1 finished (exit $?)"
wait $PID3; say "  group P3 finished (exit $?)"

say "  group P5 -- last"
Rscript "$SUB/code/replication/R1_multipair.R" P5 > "$LOGDIR/R1_P5.log" 2>&1
say "  group P5 finished (exit $?)"

say "PANEL COMPLETE"
for f in "$OUTROOT"/PANEL_SUMMARY_*.csv; do
  [ -f "$f" ] && { echo "--- $f"; cat "$f"; }
done
say "Next: bash $SUB/code/replication/run_panel.sh analyse"
