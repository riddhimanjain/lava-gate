# Replication scripts — run these on the machine holding the LD reference

Neither R script can run in a cloud container: they need
`Section 2_LD Reference/ukb_ref_consolidated/` (~15 GB) and LAVA's C++ loader. The Python
steps run anywhere.

| Script | Answers the reviewer question | Cost |
|---|---|---|
| `R2_tierB_extended.R` | "Your headline slope rests on 18 survivor loci." | ~1.5 h, 120 candidate loci |
| `R1_multipair.R` | "Is this LAVA-specific, or ASD-specific?" | **~12–15 h per phenotype group, four groups** |
| `../python/14_pair_pipeline.py` | the applied protocol, per pair | minutes (hours with `--null-sign`) |
| `../python/15_build_table8.py` | assembles Table 8 | seconds |
| `V0_exactness_check.R` | "Is your fast reimplementation really LAVA?" | ~1 min, no LD reference needed |
| `V1_gate_attribution.R` | **"Is the attenuation the filter, or just `cap.estimates`?"** | ~10 min, no LD reference needed |
| `../python/V2_gate_null_baseline.py` | "Is your 5% null the right null?" | seconds |
| `../python/V3_null_calibration.py` | "Does the weaker trait carry any signal at all?" | seconds |
| `../python/V4_closedform_vs_tierB.py` | "Does the theory predict the measurement?" | ~2 min |
| `../python/V5_null_recalibration.py` | "Is the chi-square reference calibrated?" | seconds |

---

## Order

1. **`R2_tierB_extended.R` first.** Cheapest, and it tests the paper's own headline number.
   If the extended slope differs materially from 0.5512, that number changes in the
   manuscript before anything else is done.
2. **Download the panel** — `../../data/DATA_SOURCES_MULTIPAIR.md` has the figshare DOIs.
   Note the download host: `https://ndownloader.figshare.com/files/<id>` with a browser
   User-Agent. `Start-BitsTransfer` returns HTTP 403 against that host; `curl` works.
3. **`../R/01b_harmonize_panel.R`** — converts all six releases to LAVA input and writes
   `multipair.info.txt`. Auto-detects Ricopili `daner` and PGC sumstats-VCF layouts, and
   takes the effect allele from a per-file verified mapping rather than by convention.
4. **`../R/00b_sample_overlap_ldsc_panel.R`** — munges the new traits and runs cross-trait
   LDSC over all eight, writing the single `multipair.sample.overlap.txt` and the per-pair
   `multipair.pair_re.csv`. It regression-tests the ASD × SCZ cell against the published
   0.02756 and warns loudly if the panel does not reproduce it.
5. **`R1_multipair.R <GROUP>`** — pilot each group first (see below), then run in full.
6. **`../python/14_pair_pipeline.py`** per pair, with that pair's `r_e` from step 4.
   **Run `--selftest` before trusting any of it.**
7. **`../python/15_build_table8.py`** to assemble Table 8.
8. **Verification, in order.** `V0_exactness_check.R` **must pass before `V1`** — it exits
   non-zero otherwise, because V1's whole argument depends on its direct formulas being
   LAVA's. Neither needs the LD reference: both work from the `(K, theta1, theta2, r_e)`
   tuples Tier B already wrote, which is sufficient because rho.hat and the gate statistic
   are invariant to rescaling Sigma and Omega together. Then `V2`-`V5` in numeric order.
9. **`../python/17_figure6_validation.py`** draws Figure 6 from V1 and V4.
10. **`../python/16_audit_numbers.py`** last. It checks all 269 numeric claims in both
    directions and exits non-zero on any failure.

---

## One pass per group — NOT one pass for the panel

An earlier design ran a single genome pass over all eight phenotypes, on the grounds that
`process.locus()` returns `omega` and `sigma` for every phenotype at once. That reasoning
checked the two mechanisms that make `K` depend on the phenotype set — the
`max.K = floor(max.prop.K * min(N))` cap and the `K/N > 0.1` rule — and correctly found
that neither binds here. Both are still asserted at startup, and the script stops if you
add a GWAS small enough to break them (roughly *N* < 6,100).

**A third mechanism binds, and it is decisive.** `process.input()` reduces every phenotype
to the SNPs common to *all* phenotypes in the object and the reference. Measured:

| SNP set | Count |
|---|---|
| ASD × SCZ + reference — what the published pair used | 5,774,835 |
| All eight phenotypes + reference | **2,306,339** |

The first figure is exactly what LAVA printed for the published run, so this is the real
quantity. *K* is close to proportional to a block's SNP count and the univariate
noncentrality is *K*θ, so a combined pass would change the attrition cascade it exists to
measure, and the completed pair would not reproduce its own numbers.

So: one pass per group, where a group is the smallest phenotype set that must share a SNP
set.

| Group | Phenotypes in the pass | Pairs scored |
|---|---|---|
| `P1` | SCZ, BIP | P1 |
| `P2` | MDD, BIP, MDD_noUKB, BIP_noUKB | P2a **and** P2b |
| `P3` | ASD, ADHD | P3 |
| `P5` | AN, ASD | P5 |

**P2a and P2b must share a pass.** The two PGC releases have genuinely different SNP
coverage, so separate passes would add a second difference between the arms and destroy
the control that makes P2 an experiment rather than an observation. Sharing one pass forces
an identical SNP set on both arms, at the cost of a smaller SNP set — a cost confined to P2.

Groups are independent processes and can run concurrently. `process.locus()` is dominated
by single-threaded C++ work in `load_ld.cpp` (block eigen-reconstruction), not by I/O, so
concurrency scales with cores — but each process holds its own input object, so **RAM is
the limit**, not CPU. On an 8 GB machine, run at most two groups at once, and start the
four-phenotype `P2` group on its own.

Adding the LD reference folder to the antivirus exclusion list materially reduces
`process.locus()` time.

---

## Pilot before committing a group

```bash
Rscript R1_multipair.R P1 200
```

The second argument caps the number of loci. Check the printed SNP intersection and the
`PANEL_SUMMARY_P1.csv` cascade look sane, then run the full pass:

```bash
Rscript R1_multipair.R P1
```

Checkpoints land every 25 loci in `out_R/multipair/checkpoint_<GROUP>.rds` and the script
auto-resumes, so interruption costs at most 25 loci. A pilot leaves its checkpoint in
place deliberately; delete it before the full run of that group.

---

## The Channel 2 trap

`process.locus()` drops phenotypes individually on a negative variance estimate, so in the
`P2` group a locus can be analysable for P2a and not for P2b. `R1_multipair.R` records the
surviving phenotype set per locus (`channel2_survivors_<GROUP>.csv`) and scores each pair
independently against it.

Do not shortcut this. Code that assumes both members of a pair are present will score the
filter on whichever trait survived and emit a clean-looking result for a comparison that
was never made — the exact failure the paper documents.

---

## Generic per-pair analysis, and its regression test

Scripts `09`–`12` implement the applied protocol with ASD × SCZ baked in (`p_ASD`,
`p_SCZ`, `R_E = 0.02756`, `assert len(rep) == 78`). They are the record of the completed
pair and are **not** edited. `14_pair_pipeline.py` is the same protocol parameterised.

```bash
python 14_pair_pipeline.py --selftest
```

reruns the completed pair through the generic code and checks it reproduces the published
P4 numbers to floating point, and that the frozen diagnostic thresholds independently
return the published "not applicable" verdict. **If the selftest fails, nothing downstream
is trustworthy.** Then, per pair:

```bash
python 14_pair_pipeline.py --pair P1 --trait1 SCZ --trait2 BIP --r-e <from multipair.pair_re.csv> --audit <out_R>/multipair/P1/attrition_audit.csv --outdir <results>/panel/P1
```

Add `--null-sign` for the filtered sign test; it is much slower because each replicate
needs its own bootstrap interval.

Applicability thresholds are frozen in `../../docs/DIAGNOSTIC_THRESHOLDS_frozen.md`, set on
the completed pair *before* any other pair was computed, with a per-pair expectation
recorded for each. Do not tune them to an outcome.

---

## Placeholders

The manuscript marks every dependent point with `[PENDING R1]` or `[PENDING R2]`, each
naming the script and output file its number will come from. `15_build_table8.py` writes
`PENDING` for any cell whose source file is absent, in code rather than by discipline.
Search for `PENDING` before submitting — nothing there is invented.
