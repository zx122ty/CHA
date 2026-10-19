# 161 TCM Herbs — Untargeted Metabolomics Pipeline

This repository contains the source analysis code for **CHA (Chinese Herbs Metabolic
Perturbation Atlas)**; for more information, visit <https://cha.qfxulab.com>.

Six-batch untargeted LC-MS/MS metabolomics study of 161 TCM herbs (HILIC, ESI−,
Agilent 6545B Q-TOF). Config-driven R pipeline: raw mzML → feature extraction (XCMS) →
QC & preprocessing → PCA/OPLS-DA → differential analysis (limma, FDR) → annotation
(MS2 spectral libraries + SIRIUS + endogenous accurate-mass candidates, MSI levels) →
pathway enrichment (mummichog, KEGG mouse) → cross-batch integration (feature alignment
+ ComBat/limma batch correction).

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
Rscript run_all_batches.R --reuse-oplsda --qualitative-only 4
                                                 # redo annotation onward (~20-30 min/batch),
                                                 # keep the saved OPLS-DA models
Rscript run_all_batches.R --reuse-oplsda --reuse-pathways --qualitative-only 4
                                                 # annotation-only change: also keep the
                                                 # mummichog results (Batch5 Step 7 ~80 min)
                                                 # (not after a change to Step 2: OPLS-DA,
                                                 # pathways and features must be recomputed)
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
- **No blank filter**: no metadata CSV lists blank runs, so Step 2a never runs and
  background ions stay in the feature tables (`blank_filter_applied: false` in each
  manifest). Blank runs exist as `KB*.d` in Batches 3–6 (`其他文件（调谐、KB、WASH）`) and
  as `KB*.mzML` in Batches 5–6; adding them requires re-extracting those batches.
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
| `peakwidth` | `c(10, 60)` s | Matches real chromatography (10 min HILIC gradient) |
| Detection floor | ≥2 per group, `any()` over groups | Keeps herb-specific features |
| `annotation.inhouse_enabled` | `false` | In-house DB's mz/Formula columns are inconsistent |
| `sirius.enabled` | `true` (`HMDB,CHEBI`) | Structure/formula/CANOPUS annotation |
| `sirius.min_structure_confidence` | `0.5` | A CSI:FingerID hit names a feature only above this COSMIC confidence, and only if the structure is in `sirius.structure_db` |
| `annotation.column` | `hilic` | BEH Amide; selects metid's acetate (not fluoride) adduct table for MS2 matching |
| `annotation.ms1_naming_sources` | `HMDB`, `KEGG` | Endogenous naming set (see "Annotation" below) |
| `annotation.ms1_adducts` | `(M-H)-`, `(M+CH3COO)-` | Ammonium acetate buffer; other adducts matched mostly by chance |
| `kegg.organism` | `mmu` | Mouse KEGG pathways for mummichog |
| Parallelism | `performance.cores: 0` (auto) | OPLS-DA permutations are the dominant cost |

## Split peaks: one ion in two elution waves (Step 2-pre)

On this LC system each phospholipid class elutes at its own RT (PE/PC ~60–66 s, PI
~64–72 s), and then all classes elute again in a second wave at one fixed RT (~97 s;
~88 s in Batch1). Both waves are present in every sample, so one ion became two features
24–31 s apart (gap fixed per batch), with intensities correlated at r ≈ 0.8. Taurocholate
and lyso-lipids do the same. `merge_split_peaks()` sums such pairs into one feature that
keeps the earlier, class-specific RT, so the same lipid still aligns across batches. The
gap is estimated per batch from strongly linked same-m/z pairs (r ≥ 0.8, co-detected in
≥ 80% of samples). A pair is merged when its gap is within ± 6 s of that, r ≥ 0.7 and
co-detection ≥ 50%; other same-m/z pairs (isomers) are kept. This merges ~150–190 pairs
per batch. The merges are listed in `02_Preprocessing/split_peak_merges.csv`, and merged
features carry `split_members` in `features.json` (config block `split_peaks`). Before
the merge, the duplicated taurocholate made "Primary bile acid biosynthesis" significant
in 7 of 30 Batch1 comparisons.

## Annotation (Step 6) and what the levels mean

Each feature gets one MSI level (`06_Annotation/feature_annotation.csv`, one row per
feature; `all_annotated_metabolites.csv` holds up to 5 ranked candidates per feature):

| Level | Evidence | Shown as |
|---|---|---|
| 2 — putatively identified | MS2 spectral-library match (HMDB/GNPS/MassBank/MoNA), or SIRIUS CSI:FingerID with COSMIC confidence ≥ 0.5 on a structure from HMDB/ChEBI | the compound name |
| 3 — putatively annotated | accurate mass only, against the endogenous naming set; or a SIRIUS/CANOPUS class only | `Uridine? (M243T122_NEG)` / `M763T88_NEG [Glycerophospholipids]` |
| 4 — unknown | nothing | the feature ID |

- Only Level 1–2 features carry `compound_name`, `hmdb_id`, `kegg_id` in tables,
  figures and the web export; Level 3 candidates are in `putative_*` columns.
- **SIRIUS structures** count only if they are in the searched databases
  (`sirius.structure_db`). SIRIUS 6 "expansive search" returns a PubChem structure when
  HMDB/ChEBI score low — about half of all top hits here — and COSMIC confidence stays
  high for some of them (0.82–0.85: a bromo-glucosinolate, a decyloxy-pyridine amide, a
  lipopeptide). They are dropped on import; the feature keeps the SIRIUS formula and
  CANOPUS class. SIRIUS assumes [M-H]- (acetate adducts are not in its adduct list), so
  formulas of [M+CH3COO]- ions such as PCs are wrong.
- **Endogenous naming set** for accurate-mass candidates: HMDB entries with status
  quantified/detected *and* a curated cellular location (8,000), plus the compounds of
  the mouse KEGG metabolic pathways (mummichog library). ChEBI, FooDB, PubChem and
  HMDB expected/predicted entries are not searched: against all ~490k compounds, m/z
  values shifted by a few Da ("decoys") were annotated as often as real features.
  `06_Annotation/annotation_qc.txt` reports the decoy match rate for every batch — it is
  the expected share of chance hits among Level 3 candidates (Batch1: ~40%).
- **m/z recalibration**: before annotation, 29 lock-mass ions (abundant cell metabolites
  and lipids, TFA) are located in the batch and a low-degree m/z-error curve is fitted
  through the ones that agree. If their median error is ≥ 5 ppm, feature and MS2 m/z are
  corrected (`06_Annotation/mass_calibration.txt`, original m/z kept as `mz_raw`).
  Batch5 drifts from −199 ppm (m/z < 180) to +8 ppm, Batch6 ~−9 ppm, Batch2 ~+6 ppm;
  Batches 1/3/4 are within 3 ppm and left unchanged. Features outside the lock-mass range
  whose correction exceeds 20 ppm (Batch5, m/z < ~180) are flagged
  `mz_calibrated = FALSE` and get no accurate-mass candidates.
- **SIRIUS** (Step 6b) is told the ion types: `[M-H]-` and acetate `[M+C2H4O2-H]-`.
  With MS2-only input it otherwise assumes `[M-H]-`. Element constraints also apply to
  bottom-up formula candidates (otherwise C/H/N/O/F/Si/B formulas win for phospholipids),
  and the PubChem expansive search is off. All settings go into the first SIRIUS call: a
  `config` on a later call invalidates the stored results, and the structure search then
  returns nothing. A settings stamp (`sirius_project/pipeline_sirius_settings.txt`)
  forces a recompute when the settings change. Even so, MS2 alone often cannot
  determine PC acetate formulas: PC(34:1) gets the right ion formula but as a formate
  adduct, and for PC(34:2) the correct formula ranks 9th. The formula shown for an
  identified compound is therefore the library's/HMDB's (`identified_formula`), and
  SIRIUS's own is in `sirius_formula`.
- **MS2 libraries without spectra for the acquisition polarity are skipped** (GNPS NIH
  Clinical Collection has no negative-mode spectra; metid used to fall back to MS1
  matching against this drug library).
- `tidymass2/database/MS1_Database/HMDB/hmdb_ms1_database.rda` was rebuilt with
  monoisotopic masses on 2026-10-09: massdatabase 1.0.14 `convert_hmdb2metid()` stores
  the *average* molecular weight in `mz` (median +676 ppm). The original is kept as
  `hmdb_ms1_database.avgmass_backup_20261009.rda`; the pipeline refuses an HMDB file whose
  masses disagree with its formulas.

## Pathway enrichment (Step 7)

mummichog v2 (MetaboAnalystR) on all features against the mouse KEGG pathways — it does
not need compound IDs. Significant features use the Step 4 p column and cutoff (FDR <
0.05); the reported p is mummichog's Fisher exact test on empirical compounds, q is BH
over the tested pathways (the gamma and permutation p-values are kept as extra columns;
neither was calibrated on this data).
Results: `07_KEGG/mummichog_<cmp>_pathways.csv` (all tested pathways) and
`kegg_pathway_<cmp>.xlsx` (q < 0.05; read by the web export). The library is downloaded
once into `cache/mummichog/` — metaboanalyst.ca serves qs2 files, so the R package
`qs2` must be installed.

## Advanced analysis (Step 10, `10_MetaboAnalystR_Advanced/`)

- Pathway ORA on identified compounds: significant vs all Level 1–2 compounds with a
  KEGG ID in the mouse KEGG pathways (Fisher, BH). It needs ≥ 3 significant identified
  compounds, so most comparisons are skipped (Batch1 has 6 identified compounds in the
  library).
- Biomarker analysis: Random Forest on the significant features plus a ROC curve for
  each of the top 20, labelled with display names.
- Chemical-class enrichment: one class per feature (the identified compound's, the
  CANOPUS prediction, or the class shared by all of the feature's accurate-mass
  candidates), Fisher + BH over classes.
- MetaboAnalystR is no longer used here. Its MSEA/topology took the Level 3 accurate-mass
  candidates as compound IDs, the 4.2 API was called with files where it expects vectors,
  and its reference libraries ship in qs2 format, which 4.2 cannot read. Pathway topology
  is not computed.

## Web export

`metabolites.json` lists annotated features (Level 1–3); `features.json` lists every
feature in the same record format and is what cross-batch integration and the
downstream modules align on.

## Environment

R 4.5.3 (`/opt/R/4.5.3`). Core: tidymass 2.0.10 (massprocesser / masscleaner /
massdataset), xcms 4.8.0, limma, ropls, openxlsx, tidyverse, yaml, optparse.
Optional: MetaboAnalystR 4.2.0, CAMERA, metfRag, future/future.apply/furrr, data.table,
qs/fst, digest. External: SIRIUS CLI 6.3 (login required), Java ≤17 for metfRag.
