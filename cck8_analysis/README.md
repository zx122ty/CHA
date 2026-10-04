# CCK8 × Metabolomics Downstream Analysis

Integrates the 161-herb CCK8 cell-viability screen with the six-batch untargeted
LC-MS/MS metabolomics results. CCK8 provides the **phenotypic outcome** (viability);
the metabolomics pipeline provides the **mechanistic profile** (intrahepatocellular
metabolic perturbation). Together they answer: *which herbs are toxic, and do toxic
herbs perturb the metabolic profile differently?*

- Design doc & assumptions: `prompt/下游分析/CCK8_analysis_plan.md`
- Results: `results/cck8_analysis/` (data / figures / reports / web_export)
- Website integration: see § Website below.

## Inputs

| Input | Path | Notes |
|---|---|---|
| CCK8 results | `大批量cck8实验_v3.xlsx` | 161 herbs × 2 doses (low/high) × 2 replicates; values assumed control-normalized (=1) |
| Drug metadata | `code/metadata/Browse_Table_v2.csv` | join key `match_id` (161/161 match); read with `readr` (base R chokes on a quote in this file) |
| Per-batch file lists | `code/metadata/mzml_file_list_Batch{N}_convert.csv` | drug→batch map; `drug_id == id_match` for all 5 valid batches |
| Batch web exports | `results/Batch{N}/detail/web_export/*.json` | volcano / expression / pathway JSONs (batch 4 excluded) |
| Cross-batch integration | `results/cross_batch_integration/web_export/integrated_logFC_matrix.csv` | 16,521 feature groups × 266 herb-dose comparisons |

**Assumptions** (all in `cck8_config.yaml`, change there — never in code):

- CCK8 doses are the same as the metabolomics experiment: **0.2 mg/mL (Low) / 1.0 mg/mL (High)**,
  rat primary hepatocytes, 24 h (unconfirmed with the lab — Q1/Q3 of the plan doc).
- Toxicity threshold: viability **< 0.90** = toxic; **> 1.10** = proliferative
  (sensitivity at 0.85/0.90/0.95 is reported, not just the default).
- Batch 4 herbs (28) keep their CCK8 values but are excluded from all MS-linked steps
  (`no_quant_ms = TRUE`).

## Running

```bash
cd code/this_project/cck8_analysis
Rscript run_cck8_analysis.R --config cck8_config.yaml                 # all steps, ~40 s
Rscript run_cck8_analysis.R --config cck8_config.yaml --steps 0,1     # CCK8-only (no MS data needed)
```

No packages beyond what the main pipeline uses: `readxl dplyr readr ggplot2 jsonlite yaml`.

| Step | Script | What it does | Main outputs |
|---|---|---|---|
| 0a | `load_cck8.R` | Read xlsx, join metadata, derive per-herb means / replicate discordance / dose delta + QC flags | `data/cck8_normalized.csv` |
| 0b | `qc_cck8.R` | QC report (replicate discordance > 0.2, extreme values < 0.4 / > 1.25; no plate layout available) | `reports/qc_report.md` |
| 1 | `classify_toxicity.R` | Threshold classification + sensitivity table | `data/toxicity_class.csv`, `reports/toxicity_summary.md`, `figures/toxicity_classification.*` |
| 2 | `link_metabolomics.R` | Herb-level: Spearman(viability vs perturbation magnitude) per dose, bootstrap CIs; metrics = n_significant, median \|logFC\| of sig features, max −log10 p, global perturbation distance (mean \|z_treat − z_ctrl\| across features) | `data/{metabolomics_metrics,herb_linkage,linkage_correlations}.csv`, `reports/linkage_report.md`, 2 scatter figures |
| 3 | `toxicity_features.R` | **Core**: per integrated feature group — (A) Spearman vs high-dose viability & dose delta with BH-FDR; (B) toxic-vs-neutral Wilcoxon contrast; (C) overlap of A/B hits with per-herb significant features (deduplicated, with chance baseline). Group IDs are re-derived from the same `metabolites.json` inputs + parameters as `cross_batch_integration.R` (`mz_ppm: 12`, `rt_tol: 20`) so they match the stored matrix exactly | `data/toxicity_features.csv` (16,521 rows), `reports/toxicity_features_report.md`, 2 volcano figures |
| 4 | `pathway_toxicity.R` | Fisher test per pathway: toxic vs neutral herbs with ≥1 hit; exploratory only — stored enrichment is sparse (102/266 comparisons, disease-pathway names) | `data/pathway_*.csv`, `reports/pathway_toxicity_report.md`, figure |
| 6 | `export_cck8_web.R` | Stable JSON contract for the website: per-herb rows + class counts + top-50 toxicity-associated features + per-herb linkage metrics | `web_export/cck8.json` |

`cck8_utils.R` holds shared helpers (config loading, Spearman + bootstrap CI, feature
alignment, plotting). `features_of()` / `align_features()` are adapted from
`../cross_batch_integration.R` — **keep them in sync** if the integration parameters change.

## Key results (2026-09-27 run, default thresholds)

- Classification: **toxic_high 29 / toxic_low 12 / proliferative 13 / neutral 107**;
  8 of the 41 toxic herbs sit in batch 4 (CCK8-only).
- Herb-level linkage (n = 133): viability correlates with `median_abs_logFC_sig`
  (High: ρ = 0.239, p = 0.006; Low: ρ = 0.190, p = 0.029) — toxic herbs do **not** show
  stronger global perturbation; no link to n_significant or perturbation distance.
- Feature level: 4,539 groups FDR-significant vs high-dose viability (|ρ| ≥ 0.3),
  164 vs dose delta, 16 in the toxic-vs-neutral contrast. Robustness: 36.6% of the
  toxic herbs' significant feature groups are also flagged by A/B vs a 27.5% chance baseline.
- Pathways: nothing FDR-significant (sparse stored enrichment).

## Caveats

- **n = 2 replicates per dose, no control-well raw values** → threshold-based screening only;
  no per-herb p-values. Re-run CCK8 (≥3 reps, 4–5 concentrations) for any herb before
  mechanistic claims.
- Control group differs per batch (CT1/CT2/CT3/CT6); all parsers use `CT[0-9]+`.
- OPLS-DA exports are all `available: false`, so R²Y/Q² are not used in Step 2.
- 15,320 of 16,521 integrated matrix rows carry annotations; the rest are single-batch
  features without mz/RT (not part of any aligned group).
- Step 5 (elastic-net ranking) from the plan doc was not implemented — with p ≫ n it
  adds little over the correlation screen.

## Website integration

`export_cck8_web.R` writes `results/cck8_analysis/web_export/cck8.json`. On the machine
that runs the site:

```bash
cd website/backend
python manage.py migrate metabo    # applies metabo/migrations/0010_cck8result.py
python manage.py ingest_cck8       # loads cck8.json into CCK8Result (table replaced; idempotent)
```

Backend endpoints (`website/backend/metabo/api/cck8.py`):

- `GET /api/v1/cck8` — all rows; optional `?toxicity_class=toxic_high,toxic_low&batch_id=N`
- `GET /api/v1/cck8/herbs/{match_id}` — single herb (detail page)
- `GET /api/v1/cck8/scatter` — 266 herb-dose points for the integrate page

Frontend (`website/frontend/cha_frontend/`): `/integrate` gains a
"CCK8 Toxicity × Metabolic Perturbation" scatter panel; Browse gains a CCK8 toxicity
filter facet + per-card badge; the herb detail page gains a "CCK8 Cell Viability" section.
