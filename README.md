# 161 TCM Herbs — Untargeted Metabolomics Pipeline

This repository contains the source analysis code for **CHA (Chinese Herbs Metabolic
Perturbation Atlas)**; for more information, visit <https://cha.qfxulab.com>.

Six-batch untargeted LC-MS/MS metabolomics study of 161 TCM herbs. Config-driven R
pipeline: raw mzML → feature extraction (XCMS) → QC & preprocessing → PCA/OPLS-DA →
differential analysis (limma, FDR) → multi-database annotation (MS2 libraries + SIRIUS)
→ KEGG pathway enrichment → cross-batch integration (feature alignment + ComBat/limma
batch correction).

## Code layout (`code/this_project/`)

| File | Purpose |
|---|---|
| `metabolomics_pipeline.R` | Main pipeline (steps 1–10); `--qualitative-only` mode |
| `run_all_batches.R` | One-shot runner for the 6 batches; triggers integration; applies the QC quality gate |
| `make_config.R` | Generate a batch config YAML from a metadata CSV |
| `cross_batch_integration.R` | Cross-batch feature alignment (Union-Find, m/z ±12 ppm / RT ±20 s) + ComBat/limma/median-centering correction + meta-analysis |
| `oplsda_functions.R` | OPLS-DA |
| `metaboanalyst_advanced.R` | MSEA, pathway topology, biomarker (RF/ROC), chemical class enrichment |
| `performance_utils.R` | Parallel backend, caching, fast I/O helpers |
| `export_for_web.R` / `refresh_web_export.R` | JSON web export (stable annotation contract) |
| `cck8_analysis/` | CCK8 viability × metabolomics downstream module (own README) |
| `concentration_analysis/` | Dose-response downstream module (own README) |
| `herb_category_analysis/` | TCM category × metabolomics downstream module (own README) |

## Inputs

- Raw mzML: `161herbs_raw_data/Batch{1..6}(...)/`
- Per-batch metadata: `code/metadata/mzml_file_list_Batch{N}_convert.csv`
  (columns: `batch_id,file_name,mode,MS_level,drug_id,replicate_id,drug_concentration,id_match`)
- Drug↔batch map: `metadata6.xlsx`
- Internal standards: `code/metadata/internal_standards.csv`
- CCK8 viability screen: `大批量cck8实验_v3.xlsx`

Sample naming: `{drug}_{replicate}_{H|L}.mzML`; controls `CT{n}_{n}`; QC `QC_{n}`; MS2 `QCMSMS_{n}`.

## Run all batches

```bash
cd code/this_project

# 0. Sanity check: paths, metadata, raw mzML all resolve. Fast.
Rscript run_all_batches.R --validate-only

# 1. Run all 6 batches + integration (~30 h wall clock on 20 cores)
Rscript run_all_batches.R

# Variants
Rscript run_all_batches.R --batches 1,3,5        # subset
Rscript run_all_batches.R --resume-from 3        # skip batches < 3
Rscript run_all_batches.R --integration-only     # integration only, no pipelines
Rscript run_all_batches.R --dry-run              # print commands, execute nothing
Rscript run_all_batches.R --qualitative-only 4   # batch 4, annotation-only
```

Output: `<project root>/results/Batch{1..6}/detail/web_export/` plus per-step dirs
(`03_PCA/`, `04_Differential/`, `06_Annotation/`, `07_KEGG/`, ...).

## Single batch run

```bash
Rscript code/this_project/make_config.R \
  --metadata "code/metadata/mzml_file_list_Batch1_convert.csv" \
  --raw-dir "161herbs_raw_data/Batch1(1-16)" \
  --output config_batch1.yaml \
  --batch-id 1 --project-name Batch1 --polarity negative --split-concentration

Rscript code/this_project/metabolomics_pipeline.R --config config_batch1.yaml
```

## Downstream modules

Each reads the existing pipeline + integration outputs (no re-processing) and has its
own README with step tables and assumptions.

```bash
# CCK8 × metabolomics (~40 s): toxicity classification, viability↔perturbation linkage
cd cck8_analysis && Rscript run_cck8_analysis.R --config cck8_config.yaml

# Concentration-response (~90 s): dose-gradient phenotyping, DRI
cd concentration_analysis && Rscript run_concentration_analysis.R --config concentration_config.yaml

# Herb category × metabolomics (~30–40 min): category association, prediction, signatures
cd herb_category_analysis && Rscript run_herb_category.R --config herb_category_config.yaml
```

All three support `--steps 0,1` to run a subset of steps.

## Notes for running

- **Batch 4** is a degraded acquisition — it must run with `--qualitative-only`
  (annotation + QC-diagnostic PCA only) and is excluded from integration automatically
  by the QC gate. The batch needs re-acquisition; there is no code fix.
- **Peak-table caching**: feature extraction is skipped whenever
  `Peak_table_for_cleaning.csv` exists in the batch's `peak_table/` dir — delete that
  directory to force a re-extract.
- **Module resolution**: `metabolomics_pipeline.R` loads its helpers
  (`oplsda_functions.R`, `performance_utils.R`, `metaboanalyst_advanced.R`,
  `export_for_web.R`) from the working dir's `code/` or the script's own directory —
  run from anywhere.
- **Injection order** for SVR/LOESS drift correction comes from each mzML's
  `startTimeStamp`, not column order.
- **SIRIUS v6+ requires login**: run `sirius login` once (academic token cached in
  `~/.sirius-6.3/.rtoken`).

## Key analysis settings (already in the batch configs)

| Setting | Value | Why |
|---|---|---|
| Treatment grouping | `--split-concentration` → `{drug}_{High\|Low}` | Keeps the dose contrast |
| Significance metric | BH FDR (`adj_p_value`) | Thousands of features at n=3/group |
| `peakwidth` | `c(10, 60)` s | Matches real chromatography (10 min RP gradient) |
| Detection floor | ≥2 per group, `any()` over groups | Keeps herb-specific features |
| `annotation.inhouse_enabled` | `false` | In-house DB's mz/Formula columns are inconsistent |
| `sirius.enabled` | `true` (`HMDB,CHEBI`) | Structure/formula/CANOPUS annotation |
| Parallelism | `performance.cores: 0` (auto) | OPLS-DA permutations are the dominant cost |

## Environment

R 4.5.3 (`/opt/R/4.5.3`). Core: tidymass 2.0.10 (massprocesser / masscleaner /
massdataset), xcms 4.8.0, limma, ropls, openxlsx, tidyverse, yaml, optparse.
Optional: MetaboAnalystR 4.2.0, CAMERA, metfRag, future/future.apply/furrr, data.table,
qs/fst, digest. External: SIRIUS CLI 6.3 (login required), Java ≤17 for metfRag.
