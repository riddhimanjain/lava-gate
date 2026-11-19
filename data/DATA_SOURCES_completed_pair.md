# Data sources

Every dataset used by this study, with its landing page and citation. All public. GRCh37 throughout.

## Primary GWAS

- **ASD** — Grove et al. 2019, *Nat Genet* 51:431–444. iPSYCH–PGC, 18,381 cases / 27,969 controls
  (N = 46,350). figshare: https://figshare.com/articles/dataset/asd2019/14671989
  On disk: `Section 1_Primary GWAS/Autism Spectrum (Grove et al. 2019)/iPSYCH-PGC_ASD_Nov2017`
- **SCZ** — Trubetskoy et al. 2022, *Nature* 604:502–508. PGC3 European, 53,386 / 77,258
  (N = 130,644). figshare: https://figshare.com/articles/dataset/scz2022/19426775
  On disk: `Section 1_Primary GWAS/Schizophrenia (Trubetskoy et al. 2022)/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv`

Case/control counts are recorded in `Section 7_LAVA Analysis/input.info.txt` and are taken from the
publications, not derived here.

## LD reference and genome partition

- **LAVA UK Biobank European LD reference (v1.1)**, N ≈ 99,339 — Werme et al. 2022,
  *Nat Genet* 54:274–282. https://github.com/josefin-werme/LAVA/blob/main/REFERENCE.md
  On disk: `Section 2_LD Reference/ukb_ref_consolidated/lava-ukb-v1.1*`
- **Genome partition, 2,495 approximately independent blocks** — shipped with LAVA.
  Scripts use the formatted locus file
  `LAVA/support_data/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile`, **not** the raw
  `Section 2_LD Reference/lava-partitioning/LAVA_s2500_m25_f1_w200.blocks` (which lacks the `LOC`
  column `read.loci()` requires — this was a real bug, see `TASKS_EXECUTION_LOG.md` Task 2).
- **LD scores / weights** (1000G Phase 3 EUR, HapMap3 no-MHC) — Zenodo record 10515792. Used by
  cross-trait LDSC for the sample-overlap estimate. Lives outside this folder, in the sibling
  `DSM..Neurodivergence` project; the path is set by `ldbase` in `analysis/R/00_sample_overlap_ldsc.R`.
- **HapMap3 SNP list** — via the GENESIS package data object, used at the munge step.

## Software

- **LAVA v0.1.5** — Werme et al. 2022 (as above). Vendored at `LAVA/`; this is the source every
  code-level claim in the paper is read from.
- **GenomicSEM** — Grotzinger et al. 2019, *Nat Hum Behav* 3:513.
  https://github.com/GenomicSEM/GenomicSEM — supplies the R-native cross-trait LDSC used for the
  sample-overlap estimate.
- **LDSC method** — Bulik-Sullivan et al. 2015, *Nat Genet* 47:291 and 47:1236.
- R 4.6.1, Python 3.13.3 (numpy / scipy / pandas / matplotlib).

## Ancestry constraint

European ancestry only, forced jointly by the LD reference panel and the ASD GWAS. Stated as a
limitation in the paper rather than left implicit.
