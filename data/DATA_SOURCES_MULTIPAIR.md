# Datasets for the trait-pair panel

All are public, case–control (binary), European-ancestry, and downloadable from the PGC
results page. **Binary matters**: LAVA gates binary traits with a χ² test and continuous
traits with an F-test (`analysis_functions.R:64`), and the paper's closed form is derived
for the χ² case. Quantitative phenotypes — including the PGC's quantitative anxiety
release — are therefore out of scope and are not in this panel.

Download index: <https://pgc.unc.edu/for-researchers/download-results/>

| Trait | Study | figshare DOI | Cases | Controls | Verified |
|---|---|---|---|---|---|
| SCZ | Trubetskoy 2022, *Nature* | `10.6084/m9.figshare.19426775` | 53,386 | 77,258 | ✓ already in project |
| ASD | Grove 2019, *Nat Genet* | `10.6084/m9.figshare.14671989` | 18,381 | 27,969 | ✓ already in project |
| BIP | Mullins 2021, *Nat Genet* | `10.6084/m9.figshare.14102594` | 41,917 | 371,549 | ✓ from paper |
| BIP (no UKBB) | Mullins 2021, alt release | `10.6084/m9.figshare.22564402` | read from file | read from file | ✗ not published separately |
| MDD | Wray/Howard 2018, *Nat Genet* | `10.6084/m9.figshare.14672085` | read from file | read from file | ✗ |
| MDD (rmUKBB) | Wray 2018, alt release | `10.6084/m9.figshare.21655784` | read from file | read from file | ✗ |
| ADHD | Demontis 2023, *Nat Genet* | `10.6084/m9.figshare.22564390` | 38,691 | 186,843 | ✓ from paper |
| AN | Watson 2019, *Nat Genet* | `10.6084/m9.figshare.14671980` | 16,992 | 55,525 | ✓ from paper |

Cells marked "read from file" are **not** filled in from memory. The PGC does not publish
separate sample counts for the without-UKBB releases, and `R1_multipair.R` reads N from
the summary statistics themselves and writes it into the output. Do not fill these in by
hand.

---

## The panel, and why each pair is in it

The panel spans the two dimensions along which the bias varies — the power of the weaker
trait, and the sample overlap — not a disease area.

| # | Pair | Regime probed | What it should show |
|---|---|---|---|
| **P1** | SCZ × BIP | Both well powered, moderate overlap | The case where the correction **works**: diagnostic passes, corrected magnitude reportable. Without this the paper never demonstrates its own method succeeding. |
| **P2a** | MDD × BIP, **full releases** | High overlap (both include UK Biobank) | Large manufactured intercept |
| **P2b** | MDD × BIP, **both without UKBB** | Low overlap, *same two traits* | Intercept should collapse |
| **P3** | ASD × ADHD | High overlap from shared iPSYCH cohort, moderate power | High *r*~e~ arising naturally in psychiatric data, not from a contrived biobank pair |
| **P4** | ASD × SCZ | Extreme asymmetric power, low overlap | The diagnostic **fails** — already complete |
| **P5** | AN × ASD | Both weakly powered | The floor case: how little of the genome is testable when neither trait is strong |

### P2 is the centrepiece and is worth explaining

P2a and P2b are the **same two disorders** — same true local genetic architecture, same
LD reference, same blocks. The only thing that differs is whether UK Biobank participants
are in both GWAS. So the difference between P2a and P2b isolates the effect of sample
overlap **experimentally**, holding biology fixed.

The paper currently predicts the manufactured intercept is ≈ 0.68 × *r*~e~ from theory
and confirms it at one value of *r*~e~ (0.0276). P2 tests that prediction the way a
methods paper should: by changing *r*~e~ and nothing else, and checking the intercept
moves as predicted. If it does, the mechanism is established rather than merely fitted.

This is a controlled comparison, not a survey, and it is the strongest argument in the
revised paper. It exists only because the PGC happens to ship both releases.

---

## Compute: one pass, not six

The obvious approach is one `process.input()` / `process.locus()` run per pair — six runs
at roughly 13 hours each, about 78 hours, because LAVA's C++ loader reopens the ~1 GB
`.bcor` file and rebuilds the chromosome index on every `process.locus()` call.

**This is avoidable here**, and the source audit is what shows it. `process.locus()`
accepts an input object containing *any* number of phenotypes and returns `omega` and
`sigma` for all of them at once, so one pass can serve every pair — provided the two
mechanisms that make K depend on the phenotype set do not bind:

1. `max.K = floor(max.prop.K * min(loc$N))` (`input_processing.R:125–128`) uses the
   **minimum N across all phenotypes in the input object**. Smallest N in this panel is
   ASD at 46,350, giving `max.K` = 34,762. The largest real *K* in the partition is 607.
   **Never binds.**
2. The `K/N > 0.1` rule (`input_processing.R:161–163`) nulls `h2.obs`. At ASD's N,
   607/46,350 = 0.013. **Never fires.**

Both are checked and asserted at the top of `R1_multipair.R` rather than assumed. If you
later add a GWAS with N below about 6,100, condition 2 starts firing and the combined pass
is no longer safe — the script will stop and tell you.

So: **one genome pass, ~13–15 hours, all pairs.**

One thing the combined pass makes unavoidable, which is a feature: Channel 2 drops
phenotypes individually, so a given locus may be analysable for some pairs and not others.
The script records the surviving phenotype set per locus and scores each pair
independently against it. Code that assumed both members of a pair were present would
mis-score the filter here — exactly the trap the paper documents.

---

## Order of work

1. **`R2_tierB_extended.R`** first, unchanged. Cheapest, and it tests the paper's own
   headline slope of 0.551.
2. Download the six datasets above. ASD and SCZ you already have.
3. **`R1_multipair.R`** — one combined pass.
4. Feed each pair's `attrition_audit.csv` to `08_correction.py` and `09_null_sign.py`,
   updating `R_E` to that pair's cross-trait LDSC intercept (from
   `00_sample_overlap_ldsc.R`).
5. Fill the `PENDING` markers in `manuscript.md`.

## Sources

- PGC download index: <https://pgc.unc.edu/for-researchers/download-results/>
- Mullins *et al.* 2021, *Nat Genet* **53**:817–29: <https://www.nature.com/articles/s41588-021-00857-4>
- Demontis *et al.* 2023, *Nat Genet* **55**:198–208: <https://www.nature.com/articles/s41588-022-01285-8>
- Watson *et al.* 2019, *Nat Genet* **51**:1207–14: <https://pmc.ncbi.nlm.nih.gov/articles/PMC6779477/>
- MDD without-UKBB release notes: <https://psychiatric-genomics-consortium.github.io/mdd-rmUKBB/>
