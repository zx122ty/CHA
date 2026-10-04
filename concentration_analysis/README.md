# Concentration-Response Analysis (Concentration-Response)

Exploits the **intra-drug dose gradient** (High − Low) that the CCK8 × metabolomics module
(`cck8_analysis/`) did not use. CCK8 provides a two-point phenotypic slope; the integrated
logFC matrix provides a two-point metabolic gradient per feature group. Together they answer:
*how does each drug's response change with concentration — linear, threshold, or insensitive —
and which metabolites/pathways track that dose dependence?*

- Design document & assumptions: `prompt/下游分析/concentration_response_plan.md` (user-confirmed 2026-09-27)
- Results: `results/concentration_analysis/` (data / figures / reports / web_export)
- Website contract: `web_export/concentration.json` (same flat style as `cck8.json`; a future `ingest_concentration` migration can load it)

## Inputs

| Input | Path | Notes |
|---|---|---|
| CCK8 normalized values + toxicity class | `results/cck8_analysis/data/{cck8_normalized.csv,toxicity_class.csv}` | run the cck8 module first (its steps 0–1) |
| Integrated logFC matrix | `results/cross_batch_integration/web_export/integrated_logFC_matrix.csv` | 16,521 groups × 266 comparisons; batches 1/2/3/5/6 |
| Group→member annotation | `results/Batch{N}/detail/web_export/metabolites.json` | mz±12 ppm / RT±20 s window (must match the integration) |
| KEGG compound→pathway map | `cache/kegg_link_pathway_cpd.tsv` | auto-downloaded once from KEGGREST (~2 MB); pathway names in `cache/kegg_list_pathway.txt` |

**Assumptions** (all in `concentration_config.yaml`, change there — never in code):
doses 0.2/1 mg/mL (A1, unconfirmed with the lab), single CCK8 plate (A2), proliferative herbs
counted as active (A3), batch-4 herbs phenotype-only (A4), per-herb classes descriptive at n=2
(A5), paired NA masking on the integrated logFC (A6).

## Running

```bash
cd code/this_project/concentration_analysis
Rscript run_concentration_analysis.R --config concentration_config.yaml            # all steps, ~3 min (bootstrap)
Rscript run_concentration_analysis.R --config concentration_config.yaml --steps 0,1   # CCK8-only (phenotype side)
```

| Step | Script | Output |
|---|---|---|
| 0 | `load_inputs.R` | per-herb z_viab / effects / QC flags → `data/conc_inputs.csv`, `reports/qc_report.md` |
| 1 | `phenotype_dose_response.R` | Q1 response-shape classes (linear / threshold_high / saturating_low / weak_inhibitory / proliferative / insensitive) + TCM enrichment |
| 2 | `metabolic_dose_gradient.R` | per-herb × group gradient matrix, robust z, per-herb magnitude → `data/metabolic_gradient.csv`, `gradient_z.rds` |
| 3 | `feature_viability_linkage.R` | Q2 pooled Spearman (group gradient vs z_viab) + BH-FDR + direction consistency + overlap & rank-based robustness |
| 4 | `discordant_herbs.R` | Q3 two-axis quadrants (active_both / active_phenotype_only / active_metabolism_only / inactive) + TCM enrichment |
| 5 | `pathway_dose_response.R` | Q4 KEGG pathway-level gradient association + toxic_high-vs-neutral contrast (exploratory, ~13% KEGG coverage) |
| 6 | `export_conc_web.R` | Q5 DRI + `web_export/concentration.json` |

## Why this design suits two concentrations

At n=2 per dose each herb yields exactly one slope (1 degree of freedom): a simple difference
with SE from replicate variance (`fit_dose_slope()` in `conc_utils.R`). No Hill curves, no
per-herb p-values — classification is descriptive with QC flags, and all formal inference is
cross-herb (pooled Spearman + BH-FDR, Fisher enrichment, bootstrap CIs). Metabolic gradients
(High logFC − Low logFC) cancel herb-specific baseline shifts, so the remaining signal is the
dose-dependent component itself. The code is parameterized on `doses:` in the config: adding a
third concentration extends every step (slopes switch to OLS-on-log-dose automatically).

## Extensibility notes

- New dose level → extend `doses:` list + provide the new CCK8/logFC columns upstream; no code changes.
- KEGG caches are fetched once and reused offline thereafter (`cache/`).
- If the integration parameters (mz_ppm / rt_tol) change, re-run — group IDs must match.
