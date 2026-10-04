#!/usr/bin/env Rscript
###############################################################################
# export_for_web.R — Web-Ready JSON Exporter for herbMetabo Website
#
# Reads the pipeline's internal objects (diff_results, app3, object2, etc.)
# after all steps complete and exports ECharts-compatible JSON files into a
# 'web_export/' subdirectory under the output root.
#
# Each batch run produces its own web_export/ directory. The Django ingestion
# command (ingest_r_data.py) reads these JSON files to populate PostgreSQL.
#
# Output files (all in web_export/):
#   manifest.json          — Metadata: batch info, treatments, file index
#   volcano_{cmp}.json     — Per-comparison volcano plot data
#   pca_scores.json        — PCA scores for all samples in this batch
#   pathway_enrichment_{cmp}.json  — KEGG enrichment per comparison
#   msea_{cmp}.json        — MetaboAnalystR MSEA results (if available)
#   pathway_topology_{cmp}.json    — MetaboAnalystR topology results (if available)
#   biomarker_{cmp}.json   — MetaboAnalystR biomarker results (if available)
#   chem_class_{cmp}.json  — MetaboAnalystR chemical class enrichment (if available)
#   metabolites.json       — All annotated metabolites with IDs
#   expression_{cmp}.json  — Expression matrix for heatmap/boxplot
#
# Usage:
#   source("export_for_web.R")
#   export_all(diff_results, app3, object2, cfg, msea_results = NULL, ...)
#
# Dependencies: jsonlite, dplyr, tidyr (all part of tidyverse)
###############################################################################

suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
  library(tidyr)
})

# ══════════════════════════════════════════════════════════════════════════════
# Helper: safe JSON write (pretty-print, create dir)
# ══════════════════════════════════════════════════════════════════════════════

#' Write a list to a pretty-printed JSON file, creating dirs as needed
write_json_export <- function(data, file_path) {
  dir.create(dirname(file_path), showWarnings = FALSE, recursive = TRUE)
  # na = "null" is required, not cosmetic: the default serializes NA as the
  # STRING "NA", so a downstream `is True` / truthiness check on a missing
  # metric would read "NA" as present-and-valid rather than absent.
  write_json(data, path = file_path, pretty = TRUE, auto_unbox = TRUE,
             na = "null", null = "null")
  cat(sprintf("  => Exported: %s\n", file_path))
  invisible(TRUE)
}

# ══════════════════════════════════════════════════════════════════════════════
# 1. Manifest — batch metadata + treatment index
# ══════════════════════════════════════════════════════════════════════════════

#' Build manifest.json — batch metadata and treatment list
#'
#' @param cfg Pipeline config list
#' @param diff_results Output of run_differential()
#' @param comparisons List of comparison objects from generate_comparisons()
#' @param qc_quality Output of compute_qc_quality()/assess_quantification_quality()
#' @param qualitative_only TRUE when the run skipped intensity-dependent steps
#' @return List ready for JSON export
build_manifest <- function(cfg, diff_results, comparisons,
                           qc_quality = NULL, qualitative_only = FALSE) {
  # Batch info from config
  batch_name <- cfg$project$name %||% "unknown"
  batch_id <- cfg$project$batch_id %||% NA_integer_
  output_dir <- cfg$project$output_dir

  # Build treatment list from diff_results
  treatments <- list()
  for (cmp_name in names(diff_results)) {
    dr <- diff_results[[cmp_name]]
    if (is.null(dr)) next

    cmp <- dr$comparison
    treat_name <- cmp$name %||% cmp_name
    treat_label <- cmp$label %||% treat_name

    # Extract drug_id from treatment group name (e.g., "10" from "10_vs_CT1")
    treat_group <- strsplit(treat_name, "_vs_")[[1]][1]
    control_group <- strsplit(treat_name, "_vs_")[[1]][2] %||% "CT1"

    # Determine concentration from the treatment group name
    # R pipeline names groups like "10_High" or "10_Low" depending on config
    # We store the raw comparison info and let the python parser figure it out
    n_sig <- if (!is.null(dr$sig)) nrow(dr$sig) else 0
    n_up <- dr$n_up %||% 0
    n_down <- dr$n_down %||% 0

    treatments[[treat_name]] <- list(
      comparison_name   = treat_name,
      display_label     = treat_label,
      treatment_group   = treat_group,
      control_group     = control_group,
      n_significant     = n_sig,
      n_up              = n_up,
      n_down            = n_down,
      has_volcano       = TRUE,
      has_pathway       = TRUE,
      has_msea          = FALSE,  # updated by caller if MetaboAnalystR ran
      has_expression    = TRUE
    )
  }

  # ── QC quality gate ─────────────────────────────────────────────────────
  # Cross-batch integration reads quantification_valid to decide whether this
  # batch may enter the combined quantitative analysis. A batch whose QC
  # reproducibility is far worse than its peers can still contribute m/z-based
  # annotation, but must not contribute intensities. Missing metrics are
  # treated as invalid-but-explicit rather than silently valid.
  valid <- isTRUE(qc_quality$quantification_valid)
  qc_block <- list(
    quantification_valid = valid,
    qualitative_only     = isTRUE(qualitative_only),
    exclude_reason       = if (valid) NULL else if (!is.null(qc_quality$exclude_reason) &&
                                                   !is.na(qc_quality$exclude_reason)) {
      qc_quality$exclude_reason
    } else if (isTRUE(qualitative_only)) {
      "run in qualitative-only mode"
    } else {
      "QC quality metrics unavailable"
    },
    n_qc_samples         = qc_quality$n_qc_samples %||% NA_integer_,
    n_features           = qc_quality$n_features %||% NA_integer_,
    qc_rsd_median        = qc_quality$qc_rsd_median %||% NA_real_,
    qc_rsd_pass_rate     = qc_quality$qc_rsd_pass_rate %||% NA_real_,
    qc_qc_correlation    = qc_quality$qc_qc_correlation %||% NA_real_
  )

  list(
    # Marks a batch that produced annotation but no quantitative comparison.
    # Consumers must not read an absent treatments list as "no data found".
    qualitative_only = isTRUE(qualitative_only),
    batch = list(
      id          = batch_id,
      name        = batch_name,
      output_dir  = output_dir,
      export_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ),
    config = list(
      polarity            = cfg$feature_extraction$polarity,
      ref_group           = cfg$differential$reference_group,
      p_value_cutoff      = cfg$differential$p_value_cutoff,
      fc_threshold        = cfg$differential$fc_threshold,
      significance_metric = cfg$differential$significance_metric %||% "p_value"
    ),
    qc_quality = qc_block,
    treatments = treatments,
    file_index = list(
      manifest             = "manifest.json",
      metabolites          = "metabolites.json",
      pca_scores           = "pca_scores.json",
      volcano_pattern      = "volcano_{cmp}.json",
      pathway_pattern      = "pathway_enrichment_{cmp}.json",
      expression_pattern   = "expression_{cmp}.json",
      msea_pattern         = "msea_{cmp}.json",
      topology_pattern     = "pathway_topology_{cmp}.json",
      biomarker_pattern    = "biomarker_{cmp}.json",
      chem_class_pattern   = "chem_class_{cmp}.json"
    )
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# 2. Volcano Plot JSON — {name, logFC, negLog10P, group, significance}
# ══════════════════════════════════════════════════════════════════════════════

#' Export volcano plot data for a single comparison
#'
#' Returns an ECharts-friendly array of points:
#'   {name, logFC, negLog10P, p_value, fdr, group, significance, vip_score, compound_name, hmdb_id, kegg_id}
#'
#' @param cmp_name Comparison name (e.g., "10_vs_CT1")
#' @param dr Differential result list from run_differential()
#' @param app3 Annotation table (to add compound names)
#' @param cfg Config list (for thresholds)
#' @return List with structure for JSON
build_volcano_json <- function(cmp_name, dr, app3, cfg) {
  if (is.null(dr) || is.null(dr$all)) {
    return(list(comparison = cmp_name, points = list(), n_up = 0, n_down = 0, n_not_sig = 0))
  }

  diff_all <- dr$all
  alpha <- cfg$differential$p_value_cutoff
  logfc_cutoff <- log2(cfg$differential$fc_threshold)

  # Determine significance metric
  sig_metric <- cfg$differential$significance_metric %||% "p_value"
  p_col <- if (identical(sig_metric, "adj_p_value")) "adj.P.Val" else "P.Value"

  # Build annotation lookup
  annot_lookup <- list()
  if (!is.null(app3) && nrow(app3) > 0) {
    for (i in seq_len(nrow(app3))) {
      vid <- as.character(app3$variable_id[i])
      annot_lookup[[vid]] <- list(
        compound_name = if (!is.null(app3$Compound.name[i]) && !is.na(app3$Compound.name[i]) && app3$Compound.name[i] != "") app3$Compound.name[i] else NA,
        hmdb_id = if (!is.null(app3$HMDB.ID[i]) && !is.na(app3$HMDB.ID[i]) && app3$HMDB.ID[i] != "") app3$HMDB.ID[i] else NA,
        kegg_id = if (!is.null(app3$KEGG.ID[i]) && !is.na(app3$KEGG.ID[i]) && app3$KEGG.ID[i] != "") app3$KEGG.ID[i] else NA
      )
    }
  }

  # Build points array
  points <- list()
  for (i in seq_len(nrow(diff_all))) {
    row <- diff_all[i, ]
    vid <- rownames(diff_all)[i]

    p_val <- if (!is.null(row[[p_col]])) as.numeric(row[[p_col]]) else NA
    logfc <- if (!is.null(row$logFC)) as.numeric(row$logFC) else NA
    neg_log10p <- if (!is.na(p_val) && p_val > 0) -log10(p_val) else 0

    # Significance classification
    is_sig <- !is.na(p_val) && !is.na(logfc) && p_val < alpha && abs(logfc) > logfc_cutoff
    sig_group <- if (is_sig) {
      if (logfc > logfc_cutoff) "Up" else "Down"
    } else {
      "Not"
    }

    # VIP score (if available from OPLS-DA)
    vip <- if (!is.null(dr$vip_scores) && vid %in% names(dr$vip_scores)) dr$vip_scores[[vid]] else NA

    # Annotation
    ann <- annot_lookup[[vid]] %||% list(compound_name = NA, hmdb_id = NA, kegg_id = NA)

    point <- list(
      name = vid,
      logFC = round(logfc, 4),
      negLog10P = round(neg_log10p, 4),
      p_value = if (!is.na(p_val)) format(p_val, scientific = TRUE) else NA,
      fdr = if (!is.null(row$adj.P.Val)) format(as.numeric(row$adj.P.Val), scientific = TRUE) else NA,
      group = sig_group,
      significance = if (is_sig) "significant" else "not_significant",
      vip_score = if (!is.na(vip)) round(vip, 4) else NA,
      compound_name = ann$compound_name,
      hmdb_id = ann$hmdb_id,
      kegg_id = ann$kegg_id
    )
    points[[i]] <- point
  }

  n_up <- sum(sapply(points, function(p) p$group == "Up"))
  n_down <- sum(sapply(points, function(p) p$group == "Down"))
  n_not <- sum(sapply(points, function(p) p$group == "Not"))

  list(
    comparison = cmp_name,
    label = dr$comparison$label %||% cmp_name,
    thresholds = list(
      logFC_cutoff = logfc_cutoff,
      p_value_cutoff = alpha,
      significance_metric = sig_metric
    ),
    summary = list(
      total = length(points),
      n_up = n_up,
      n_down = n_down,
      n_not_significant = n_not,
      n_significant = n_up + n_down
    ),
    points = points
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. PCA JSON — {sample_id, PC1, PC2, PC3, group, batch}
# ══════════════════════════════════════════════════════════════════════════════

#' Export PCA scores as ECharts scatter plot data
#'
#' @param object2 mass_dataset object (after preprocessing)
#' @param sample_info data.frame with sample_id and class
#' @param cfg Config list
#' @return List with scores array and variance explained
build_pca_json <- function(object2, sample_info, cfg) {
  # Run PCA (same as run_pca() in pipeline)
  group <- sample_info$class[match(colnames(object2@expression_data),
                                    sample_info$sample_id)]
  pca_input <- t(object2@expression_data)
  pca_input <- pca_input[, apply(pca_input, 2, var, na.rm = TRUE) > 0, drop = FALSE]

  if (ncol(pca_input) < 2) {
    return(list(scores = list(), variance_explained = list(), n_components = 0))
  }

  pca1 <- prcomp(pca_input, scale. = TRUE)
  summ1 <- summary(pca1)
  scores_df <- as.data.frame(pca1$x)

  # Build scores array for ECharts
  scores <- list()
  for (i in seq_len(nrow(scores_df))) {
    sid <- rownames(scores_df)[i]
    scores[[i]] <- list(
      sample_id = sid,
      PC1 = round(scores_df$PC1[i], 4),
      PC2 = round(scores_df$PC2[i], 4),
      PC3 = round(scores_df$PC3[i], 4),
      group = as.character(group[i]),
      batch = cfg$project$batch_id %||% NA
    )
  }

  # Variance explained
  var_exp <- list(
    PC1 = round(summ1$importance[2, 1] * 100, 2),
    PC2 = round(summ1$importance[2, 2] * 100, 2),
    PC3 = round(summ1$importance[2, 3] * 100, 2),
    PC4 = round(summ1$importance[2, 4] * 100, 2),
    PC5 = round(summ1$importance[2, 5] * 100, 2)
  )

  list(
    scores = scores,
    variance_explained = var_exp,
    n_components = min(5, ncol(scores_df)),
    n_samples = nrow(scores_df),
    n_features = ncol(pca_input)
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# 3b. Internal Standard QC JSON — per-batch instrument stability
# ══════════════════════════════════════════════════════════════════════════════

#' Build internal_standard_qc.json — IS peak areas, RSD% and detection rates
#'
#' @param is_qc_stats Output of extract_and_remove_internal_standards() (or NULL)
#' @param cfg Validated config list
build_internal_standard_qc_json <- function(is_qc_stats, cfg) {
  if (is.null(is_qc_stats)) return(NULL)

  rounds <- function(x) if (is.null(x) || all(is.na(x))) x else round(x, 4)
  standards <- lapply(unname(is_qc_stats), function(s) {
    list(
      name               = s$name,
      variable_id        = s$variable_id,
      mz                 = s$mz,
      rt_expected_sec    = rounds(s$rt_expected_sec),
      rt_observed_sec    = rounds(s$rt_observed_sec),
      mean_area_qc       = rounds(s$mean_area_qc),
      rsd_pct_qc         = rounds(s$rsd_pct_qc),
      mean_area_all      = rounds(s$mean_area_all),
      rsd_pct_all        = rounds(s$rsd_pct_all),
      n_samples_detected = s$n_samples_detected,
      n_samples_total    = s$n_samples_total,
      detection_rate     = rounds(s$detection_rate)
    )
  })

  list(
    batch_id           = cfg$project$batch_id %||% NA_integer_,
    batch_name         = cfg$project$name %||% "unknown",
    n_internal_standards = length(standards),
    internal_standards = standards
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. Pathway Enrichment JSON — bubble chart data
# ══════════════════════════════════════════════════════════════════════════════

#' Export KEGG pathway enrichment results for a single comparison
#'
#' Returns ECharts bubble chart data:
#'   {pathway_name, p_value, neg_log10_p, mapped_count, background_count, impact_score}
#'
#' @param cmp_name Comparison name
#' @param kegg_dir Path to KEGG output directory (07_KEGG/)
#' @return List with pathway array
build_pathway_json <- function(cmp_name, kegg_dir) {
  # Try to find the KEGG enrichment result file
  # Pattern: kegg_pathway_{cmp_name}.xlsx or similar
  xlsx_path <- file.path(kegg_dir, paste0("kegg_pathway_", cmp_name, ".xlsx"))

  if (!file.exists(xlsx_path)) {
    # Try CSV fallback
    csv_path <- file.path(kegg_dir, paste0("kegg_pathway_", cmp_name, ".csv"))
    if (file.exists(csv_path)) {
      xlsx_path <- csv_path
    }
  }

  if (!file.exists(xlsx_path)) {
    # Try the annotated CSV as fallback
    alt_path <- file.path(kegg_dir, paste0("annotated_", cmp_name, ".csv"))
    if (file.exists(alt_path)) {
      # annotated file has different structure — return empty
      return(list(comparison = cmp_name, pathways = list(), n_pathways = 0))
    }
    return(list(comparison = cmp_name, pathways = list(), n_pathways = 0))
  }

  # Read the enrichment results
  if (grepl("\\.xlsx$", xlsx_path)) {
    tryCatch({
      pathway_df <- openxlsx::read.xlsx(xlsx_path)
    }, error = function(e) {
      pathway_df <- NULL
    })
  } else {
    pathway_df <- tryCatch({
      read.csv(xlsx_path, stringsAsFactors = FALSE)
    }, error = function(e) NULL)
  }

  if (is.null(pathway_df) || nrow(pathway_df) == 0) {
    return(list(comparison = cmp_name, pathways = list(), n_pathways = 0))
  }

  # Standardize column names (HMDB KEGG enrichment output)
  # Expected columns: pathway_name, p_value, mapped_number, background_number, ...
  p_name_col <- if ("pathway_name" %in% colnames(pathway_df)) "pathway_name" else
                if ("Pathway" %in% colnames(pathway_df)) "Pathway" else NULL
  p_val_col <- if ("p_value" %in% colnames(pathway_df)) "p_value" else
               if ("P.Value" %in% colnames(pathway_df)) "P.Value" else
               if ("p.value" %in% colnames(pathway_df)) "p.value" else NULL
  mapped_col <- if ("mapped_number" %in% colnames(pathway_df)) "mapped_number" else
                if ("Mapped" %in% colnames(pathway_df)) "Mapped" else
                if ("mapped" %in% colnames(pathway_df)) "mapped" else NULL
  bg_col <- if ("background_number" %in% colnames(pathway_df)) "background_number" else
            if ("Background" %in% colnames(pathway_df)) "Background" else NULL
  impact_col <- if ("impact_score" %in% colnames(pathway_df)) "impact_score" else
                if ("Impact" %in% colnames(pathway_df)) "Impact" else NULL

  if (is.null(p_name_col) || is.null(p_val_col)) {
    return(list(comparison = cmp_name, pathways = list(), n_pathways = 0))
  }

  pathways <- list()
  for (i in seq_len(nrow(pathway_df))) {
    p_val <- as.numeric(pathway_df[[p_val_col]][i])
    neg_log10 <- if (!is.na(p_val) && p_val > 0) -log10(p_val) else 0
    mapped <- if (!is.null(mapped_col)) as.numeric(pathway_df[[mapped_col]][i]) %||% NA else NA
    bg <- if (!is.null(bg_col)) as.numeric(pathway_df[[bg_col]][i]) %||% NA else NA
    impact <- if (!is.null(impact_col)) as.numeric(pathway_df[[impact_col]][i]) %||% NA else NA

    pathways[[i]] <- list(
      pathway_name = as.character(pathway_df[[p_name_col]][i]),
      p_value = format(p_val, scientific = TRUE),
      neg_log10_p = round(neg_log10, 4),
      mapped_count = if (!is.na(mapped)) mapped else 0,
      background_count = if (!is.na(bg)) bg else 0,
      impact_score = if (!is.na(impact)) round(impact, 4) else NA
    )
  }

  # Sort by p-value ascending
  p_vals <- sapply(pathways, function(p) as.numeric(p$p_value))
  ord <- order(p_vals, na.last = TRUE)
  pathways <- pathways[ord]

  list(
    comparison = cmp_name,
    # The tidymass step runs enrich_hmdb() against the HMDB/SMPDB pathway
    # collection — the "KEGG" folder/file naming is a legacy label only.
    database = "HMDB/SMPDB",
    n_pathways = length(pathways),
    pathways = pathways
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# 5. Metabolites JSON — master metabolite reference table
# ══════════════════════════════════════════════════════════════════════════════

#' Export all annotated metabolites as a reference table
#'
#' Each metabolite gets a unique key (variable_id) and carries
#' HMDB ID, KEGG ID, compound name, formula, m/z, RT, super class.
#'
#' Build the metabolites.json payload
#'
#' Emits the standard annotation contract: every record carries an identifier
#' (`inchikey`), the adduct it was matched as, its MSI level, the engine that
#' produced it (`source`), and a confidence `score` — so downstream consumers do
#' not have to know which annotation tool ran.
#'
#' @param app3 Annotation table from run_annotation() (optionally enriched by
#'   run_sirius_analysis() / run_metfrag_analysis())
#' @param object2 mass_dataset (for expression data)
#' @return List of metabolite records
build_metabolites_json <- function(app3, object2) {
  if (is.null(app3) || nrow(app3) == 0) {
    return(list(metabolites = list(), n_metabolites = 0))
  }

  # Deduplicate by variable_id (keep first annotation for each feature)
  app3_dedup <- app3 %>%
    group_by(variable_id) %>%
    slice(1) %>%
    ungroup()

  # `%||%` is defined in metabolomics_pipeline.R; standalone callers may not
  # have it in scope.
  if (!exists("%||%", mode = "function")) {
    `%||%` <<- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a
  }

  # Safe scalar getter: returns NA when the column is absent, NULL or blank.
  pick <- function(row, col) {
    if (is.null(col) || !col %in% colnames(row)) return(NA)
    v <- row[[col]]
    if (is.null(v) || length(v) == 0 || is.na(v[1]) || identical(as.character(v[1]), "")) return(NA)
    v[1]
  }

  metabolites <- list()
  for (i in seq_len(nrow(app3_dedup))) {
    row <- app3_dedup[i, ]
    vid <- as.character(row$variable_id)

    # Prefer a SIRIUS CSI:FingerID structure hit when present (it is what the
    # identification report promotes to Level 2); fall back to the library hit.
    sirius_name    <- pick(row, "sirius_name")
    library_name   <- pick(row, "library_compound_name")
    library_name   <- if (is.na(library_name)) pick(row, "Compound.name") else library_name
    compound_name  <- if (!is.na(sirius_name)) sirius_name else library_name

    # Formula: SIRIUS predicts it even when no MS1 library hit carries one.
    formula <- pick(row, "formula")
    if (is.na(formula)) formula <- pick(row, "sirius_formula")
    if (is.na(formula)) formula <- pick(row, "Formula")

    # Confidence score: prefer the spectral score (SS) for MS2 hits, else the
    # composite score, else the SIRIUS formula score.
    score <- pick(row, "SS")
    if (is.na(score)) score <- pick(row, "Total.score")
    if (is.na(score)) score <- pick(row, "sirius_formula_score")

    source <- pick(row, "compound_source")
    if (is.na(source)) {
      # `compound_source` is only set by generate_identification_report(); when
      # the export runs off the pre-report table, derive provenance from what
      # actually supplied the name — otherwise a SIRIUS-named compound would be
      # attributed to whichever MS1 library happened to match its m/z.
      source <- if (!is.na(sirius_name)) "SIRIUS_CSI:FingerID" else pick(row, "database_source")
    }

    metabolites[[i]] <- list(
      variable_id = vid,
      compound_name = compound_name,
      # ── standard annotation contract ──────────────────────────────────────
      inchikey = pick(row, "sirius_InChIkey"),
      smiles = pick(row, "sirius_smiles"),
      formula = formula,
      adduct = pick(row, "Adduct"),
      msi_level = pick(row, "confidence_level"),
      source = source,
      score = if (is.na(score)) NA else round(as.numeric(score), 4),
      # ── provenance ────────────────────────────────────────────────────────
      library_compound_name = library_name,
      sirius_name = sirius_name,
      canopus_class = pick(row, "canopus_cf_class"),
      canopus_superclass = pick(row, "canopus_cf_superclass"),
      hmdb_id = pick(row, "HMDB.ID"),
      kegg_id = pick(row, "KEGG.ID"),
      mz = if (!is.null(row$mz)) round(as.numeric(row$mz), 6) else NA,
      rt = if (!is.null(row$rt)) round(as.numeric(row$rt), 2) else NA,
      database_source = pick(row, "database_source"),
      match_type = pick(row, "match_type"),
      confidence_level = pick(row, "confidence_level")
    )
  }

  list(
    metabolites = metabolites,
    n_metabolites = length(metabolites),
    n_unique_compounds = length(unique(
      sapply(metabolites, function(m) m$compound_name)
    ))
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# 6. Expression Matrix JSON — for heatmap / boxplot
# ══════════════════════════════════════════════════════════════════════════════

#' Export expression matrix for a comparison
#'
#' Returns the normalized expression values for significant metabolites
#' across all samples, ready for ECharts heatmap.
#'
#' @param cmp_name Comparison name
#' @param dr Differential result
#' @param object2 mass_dataset
#' @param app3 Annotation table
#' @return List with expression data
build_expression_json <- function(cmp_name, dr, object2, app3) {
  if (is.null(dr) || is.null(dr$sig) || nrow(dr$sig) == 0) {
    return(list(comparison = cmp_name, metabolites = list(), samples = list(), data = list()))
  }

  # Get significant metabolite IDs
  sig_ids <- rownames(dr$sig)

  # Expression data for these metabolites
  expr <- object2@expression_data
  sig_expr <- expr[rownames(expr) %in% sig_ids, , drop = FALSE]

  if (nrow(sig_expr) == 0) {
    return(list(comparison = cmp_name, metabolites = list(), samples = list(), data = list()))
  }

  # Sample info
  samples <- colnames(sig_expr)

  # Metabolite names (from annotation)
  name_map <- list()
  if (!is.null(app3) && nrow(app3) > 0) {
    for (i in seq_len(nrow(app3))) {
      vid <- as.character(app3$variable_id[i])
      if (vid %in% names(name_map)) next
      cname <- if (!is.null(app3$Compound.name[i]) && !is.na(app3$Compound.name[i]) && app3$Compound.name[i] != "") app3$Compound.name[i] else vid
      name_map[[vid]] <- cname
    }
  }

  # Build metabolite list
  metabolites <- list()
  for (vid in rownames(sig_expr)) {
    metabolites[[length(metabolites) + 1]] <- list(
      variable_id = vid,
      compound_name = name_map[[vid]] %||% vid
    )
  }

  # Build data matrix (list of rows)
  data_rows <- list()
  for (i in seq_len(nrow(sig_expr))) {
    vid <- rownames(sig_expr)[i]
    vals <- as.numeric(sig_expr[i, ])
    names(vals) <- colnames(sig_expr)

    row <- list(
      variable_id = vid,
      compound_name = name_map[[vid]] %||% vid,
      values = as.list(setNames(round(vals, 4), colnames(sig_expr)))
    )
    data_rows[[i]] <- row
  }

  list(
    comparison = cmp_name,
    n_metabolites = nrow(sig_expr),
    n_samples = length(samples),
    metabolites = metabolites,
    samples = samples,
    data = data_rows
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# 7. MetaboAnalystR Results Export (if available)
# ══════════════════════════════════════════════════════════════════════════════

#' Export MetaboAnalystR MSEA results
#'
#' @param msea_results Output from run_msea() in metaboanalyst_advanced.R
#' @param cmp_name Comparison name
#' @return List for JSON
build_msea_json <- function(msea_results, cmp_name) {
  if (is.null(msea_results) || !is.list(msea_results) || length(msea_results) == 0) {
    return(list(comparison = cmp_name, sets = list(), n_sets = 0))
  }

  # msea_results is typically a list with entries per comparison
  cmp_result <- msea_results[[cmp_name]]
  if (is.null(cmp_result)) {
    return(list(comparison = cmp_name, sets = list(), n_sets = 0))
  }

  # Extract the result table (varies by MSEA method)
  result_table <- cmp_result$result
  if (is.null(result_table) || nrow(result_table) == 0) {
    return(list(comparison = cmp_name, sets = list(), n_sets = 0))
  }

  sets <- list()
  for (i in seq_len(nrow(result_table))) {
    sets[[i]] <- list(
      set_name = as.character(result_table$set_name[i] %||% result_table$ pathway[i] %||% paste0("Set_", i)),
      p_value = format(as.numeric(result_table$p_value[i] %||% 1), scientific = TRUE),
      neg_log10_p = if (!is.na(as.numeric(result_table$p_value[i] %||% 1)) && as.numeric(result_table$p_value[i] %||% 1) > 0)
                      round(-log10(as.numeric(result_table$p_value[i])), 4) else 0,
      hit_count = as.integer(result_table$hit_count[i] %||% result_table$matched[i] %||% 0),
      total_count = as.integer(result_table$total_count[i] %||% result_table$size[i] %||% 0),
      fdr = if (!is.null(result_table$fdr[i])) format(as.numeric(result_table$fdr[i]), scientific = TRUE) else NA
    )
  }

  list(
    comparison = cmp_name,
    method = cmp_result$method %||% "globaltest",
    n_sets = length(sets),
    sets = sets
  )
}

#' Export MetaboAnalystR Pathway Topology results
#'
#' @param topology_results Output from run_pathway_topology()
#' @param cmp_name Comparison name
#' @return List for JSON
build_topology_json <- function(topology_results, cmp_name) {
  if (is.null(topology_results) || !is.list(topology_results) || length(topology_results) == 0) {
    return(list(comparison = cmp_name, pathways = list(), n_pathways = 0))
  }

  cmp_result <- topology_results[[cmp_name]]
  if (is.null(cmp_result)) {
    return(list(comparison = cmp_name, pathways = list(), n_pathways = 0))
  }

  result_table <- cmp_result$result
  if (is.null(result_table) || nrow(result_table) == 0) {
    return(list(comparison = cmp_name, pathways = list(), n_pathways = 0))
  }

  pathways <- list()
  for (i in seq_len(nrow(result_table))) {
    pathways[[i]] <- list(
      pathway_name = as.character(result_table$pathway_name[i] %||% paste0("Pathway_", i)),
      p_value = format(as.numeric(result_table$p_value[i] %||% 1), scientific = TRUE),
      neg_log10_p = if (!is.na(as.numeric(result_table$p_value[i] %||% 1)) && as.numeric(result_table$p_value[i] %||% 1) > 0)
                      round(-log10(as.numeric(result_table$p_value[i])), 4) else 0,
      impact_score = round(as.numeric(result_table$impact[i] %||% result_table$impact_score[i] %||% 0), 4),
      hit_count = as.integer(result_table$hit_count[i] %||% result_table$matched[i] %||% 0),
      total_count = as.integer(result_table$total_count[i] %||% result_table$size[i] %||% 0)
    )
  }

  list(
    comparison = cmp_name,
    method = cmp_result$metric %||% "rbc",
    n_pathways = length(pathways),
    pathways = pathways
  )
}

#' Export MetaboAnalystR Biomarker Analysis results
#'
#' @param biomarker_results Output from run_biomarker_analysis()
#' @param cmp_name Comparison name
#' @return List for JSON
build_biomarker_json <- function(biomarker_results, cmp_name) {
  if (is.null(biomarker_results) || !is.list(biomarker_results) || length(biomarker_results) == 0) {
    return(list(comparison = cmp_name, features = list(), n_features = 0, roc_auc = NA))
  }

  cmp_result <- biomarker_results[[cmp_name]]
  if (is.null(cmp_result)) {
    return(list(comparison = cmp_name, features = list(), n_features = 0, roc_auc = NA))
  }

  features <- list()
  if (!is.null(cmp_result$feature_importance) && nrow(cmp_result$feature_importance) > 0) {
    fi <- cmp_result$feature_importance
    for (i in seq_len(nrow(fi))) {
      features[[i]] <- list(
        variable_id = rownames(fi)[i] %||% as.character(fi$variable_id[i] %||% i),
        importance = round(as.numeric(fi$importance[i] %||% fi$MeanDecreaseAccuracy[i] %||% 0), 4)
      )
    }
  }

  list(
    comparison = cmp_name,
    method = "Random Forest",
    n_features = length(features),
    roc_auc = if (!is.null(cmp_result$roc_auc)) round(as.numeric(cmp_result$roc_auc), 4) else NA,
    features = features
  )
}

#' Export MetaboAnalystR Chemical Class Enrichment results
#'
#' @param chem_class_results Output from run_chem_class_enrichment()
#' @param cmp_name Comparison name
#' @return List for JSON
build_chem_class_json <- function(chem_class_results, cmp_name) {
  if (is.null(chem_class_results) || !is.list(chem_class_results) || length(chem_class_results) == 0) {
    return(list(comparison = cmp_name, classes = list(), n_classes = 0))
  }

  cmp_result <- chem_class_results[[cmp_name]]
  if (is.null(cmp_result)) {
    return(list(comparison = cmp_name, classes = list(), n_classes = 0))
  }

  result_table <- cmp_result$result
  if (is.null(result_table) || nrow(result_table) == 0) {
    return(list(comparison = cmp_name, classes = list(), n_classes = 0))
  }

  classes <- list()
  for (i in seq_len(nrow(result_table))) {
    classes[[i]] <- list(
      class_name = as.character(result_table$class_name[i] %||% result_table$super_class[i] %||% paste0("Class_", i)),
      p_value = format(as.numeric(result_table$p_value[i] %||% 1), scientific = TRUE),
      neg_log10_p = if (!is.na(as.numeric(result_table$p_value[i] %||% 1)) && as.numeric(result_table$p_value[i] %||% 1) > 0)
                      round(-log10(as.numeric(result_table$p_value[i])), 4) else 0,
      hit_count = as.integer(result_table$hit_count[i] %||% result_table$matched[i] %||% 0),
      total_count = as.integer(result_table$total_count[i] %||% result_table$size[i] %||% 0),
      fdr = if (!is.null(result_table$fdr[i])) format(as.numeric(result_table$fdr[i]), scientific = TRUE) else NA
    )
  }

  list(
    comparison = cmp_name,
    taxon_level = cmp_result$taxon_level %||% "super_class",
    n_classes = length(classes),
    classes = classes
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# 8. OPLS-DA Results (if available)
# ══════════════════════════════════════════════════════════════════════════════

#' Export OPLS-DA scores and VIP values
#'
#' @param oplsda_dir Path to OPLS-DA output directory
#' @param cmp_name Comparison name
#' @return List for JSON
build_oplsda_json <- function(oplsda_dir, cmp_name) {
  # Try to find OPLS-DA score files
  score_file <- file.path(oplsda_dir, paste0(cmp_name, "_oplsda_scores.csv"))
  vip_file <- file.path(oplsda_dir, paste0(cmp_name, "_vip_scores.csv"))

  if (!file.exists(score_file)) {
    return(list(comparison = cmp_name, available = FALSE))
  }

  scores <- list()
  if (file.exists(score_file)) {
    score_df <- read.csv(score_file, stringsAsFactors = FALSE)
    for (i in seq_len(nrow(score_df))) {
      scores[[i]] <- list(
        sample_id = as.character(score_df$sample_id[i] %||% rownames(score_df)[i]),
        score_1 = round(as.numeric(score_df$p1[i] %||% score_df$t1[i] %||% 0), 4),
        score_2 = round(as.numeric(score_df$o1[i] %||% score_df$to1[i] %||% 0), 4),
        group = as.character(score_df$group[i] %||% "Unknown")
      )
    }
  }

  vip <- list()
  if (file.exists(vip_file)) {
    vip_df <- read.csv(vip_file, stringsAsFactors = FALSE)
    for (i in seq_len(nrow(vip_df))) {
      vip[[i]] <- list(
        variable_id = as.character(vip_df$variable_id[i] %||% rownames(vip_df)[i]),
        vip_score = round(as.numeric(vip_df$vip[i] %||% vip_df$VIP[i] %||% 0), 4)
      )
    }
  }

  list(
    comparison = cmp_name,
    available = TRUE,
    n_samples = length(scores),
    n_vip_features = length(vip),
    scores = scores,
    vip = vip
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# 9. Summary Statistics JSON
# ══════════════════════════════════════════════════════════════════════════════

#' Export per-comparison summary statistics
#'
#' @param diff_results Output of run_differential()
#' @param app3 Annotation table
#' @return List of summary stats
build_summary_json <- function(diff_results, app3) {
  comparisons <- list()
  for (cmp_name in names(diff_results)) {
    dr <- diff_results[[cmp_name]]
    if (is.null(dr)) next

    n_total <- if (!is.null(dr$all)) nrow(dr$all) else 0
    n_sig <- if (!is.null(dr$sig)) nrow(dr$sig) else 0
    n_up <- dr$n_up %||% 0
    n_down <- dr$n_down %||% 0

    # Count how many significant metabolites have annotations
    n_annotated <- 0
    if (!is.null(dr$sig) && nrow(dr$sig) > 0 && !is.null(app3) && nrow(app3) > 0) {
      sig_ids <- rownames(dr$sig)
      ann_ids <- unique(app3$variable_id)
      n_annotated <- sum(sig_ids %in% ann_ids)
    }

    comparisons[[cmp_name]] <- list(
      comparison_name = cmp_name,
      label = dr$comparison$label %||% cmp_name,
      n_total_features = n_total,
      n_significant = n_sig,
      n_up = n_up,
      n_down = n_down,
      n_annotated = n_annotated,
      n_ctrl_samples = length(dr$ctrl_samples %||% c()),
      n_treat_samples = length(dr$treat_samples %||% c())
    )
  }

  list(
    n_comparisons = length(comparisons),
    comparisons = comparisons
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# 10. Master Export Function
# ══════════════════════════════════════════════════════════════════════════════

#' Export all web-ready JSON files
#'
#' Call this function at the end of the pipeline's main() function.
#' It reads the pipeline's internal objects and writes ECharts-compatible
#' JSON files into <output_dir>/web_export/.
#'
#' @param diff_results Output of run_differential() — list of comparison results
#' @param app3 Annotation table from run_annotation() — data.frame
#' @param object2 mass_dataset object after preprocessing
#' @param cfg Validated config list from load_and_validate_config()
#' @param sample_info data.frame with sample_id and class (from preprocessing)
#' @param comparisons List of comparison objects (from generate_comparisons())
#' @param msea_results Optional — MetaboAnalystR MSEA results
#' @param topology_results Optional — MetaboAnalystR pathway topology results
#' @param biomarker_results Optional — MetaboAnalystR biomarker results
#' @param chem_class_results Optional — MetaboAnalystR chemical class enrichment
#' @param kegg_dir Optional — path to KEGG output directory (default: <output_dir>/07_KEGG)
#' @param oplsda_dir Optional — path to OPLS-DA output directory
#' @return Invisible list of all exported data
export_all <- function(diff_results, app3, object2, cfg,
                       sample_info = NULL,
                       comparisons = NULL,
                       msea_results = NULL,
                       topology_results = NULL,
                       biomarker_results = NULL,
                       chem_class_results = NULL,
                       kegg_dir = NULL,
                       oplsda_dir = NULL,
                       is_qc_stats = NULL,
                       qc_quality = NULL,
                       qualitative_only = FALSE) {
  cat("\n========================================\n")
  cat("     Web Export: Generating JSON files\n")
  cat("========================================\n")

  out_root <- cfg$project$output_dir
  web_dir <- file.path(out_root, "web_export")
  dir.create(web_dir, showWarnings = FALSE, recursive = TRUE)
  cat(sprintf("=> Export directory: %s\n", web_dir))

  # ── Use defaults if not provided ────────────────────────────────────────
  if (is.null(sample_info) && exists("object2", inherits = FALSE)) {
    sample_info <- object2@sample_info
  }
  if (is.null(comparisons) && exists("diff_results", inherits = FALSE)) {
    comparisons <- lapply(diff_results, function(dr) dr$comparison)
  }
  if (is.null(kegg_dir)) {
    kegg_dir <- file.path(out_root, "07_KEGG")
  }
  if (is.null(oplsda_dir)) {
    oplsda_dir <- file.path(out_root, "03.5_OPLSDA")
  }

  # ── 1. Manifest ─────────────────────────────────────────────────────────
  cat("\n[1/7] Building manifest.json...\n")
  manifest <- build_manifest(cfg, diff_results, comparisons,
                             qc_quality = qc_quality,
                             qualitative_only = qualitative_only)
  write_json_export(manifest, file.path(web_dir, "manifest.json"))
  if (!isTRUE(qc_quality$quantification_valid)) {
    cat(sprintf("  ⚠️  manifest marks quantification_valid = false (%s)\n",
                qc_quality$exclude_reason %||% "no QC metrics"))
    cat("      This batch will be excluded from cross-batch integration.\n")
  }

  # ── 2. Metabolites reference ────────────────────────────────────────────
  cat("[2/7] Building metabolites.json...\n")
  metabolites <- build_metabolites_json(app3, object2)
  write_json_export(metabolites, file.path(web_dir, "metabolites.json"))

  # ── 3. PCA scores ───────────────────────────────────────────────────────
  cat("[3/7] Building pca_scores.json...\n")
  pca_data <- build_pca_json(object2, sample_info, cfg)
  write_json_export(pca_data, file.path(web_dir, "pca_scores.json"))

  # ── 3b. Internal standard QC (optional) ─────────────────────────────────
  if (!is.null(is_qc_stats)) {
    cat("[3b/7] Building internal_standard_qc.json...\n")
    is_qc_json <- build_internal_standard_qc_json(is_qc_stats, cfg)
    write_json_export(is_qc_json, file.path(web_dir, "internal_standard_qc.json"))
  }

  # ── 4. Summary statistics ───────────────────────────────────────────────
  cat("[4/7] Building summary.json...\n")
  summary_data <- build_summary_json(diff_results, app3)
  write_json_export(summary_data, file.path(web_dir, "summary.json"))

  # ── 5. Per-comparison files ─────────────────────────────────────────────
  cat("[5/7] Building per-comparison files...\n")

  # Update manifest with presence of MetaboAnalystR results
  has_msea <- !is.null(msea_results)
  has_topology <- !is.null(topology_results)
  has_biomarker <- !is.null(biomarker_results)
  has_chem_class <- !is.null(chem_class_results)

  for (cmp_name in names(diff_results)) {
    dr <- diff_results[[cmp_name]]
    if (is.null(dr)) next

    cat(sprintf("   Processing: %s\n", cmp_name))

    # Volcano
    volcano <- build_volcano_json(cmp_name, dr, app3, cfg)
    write_json_export(volcano, file.path(web_dir, paste0("volcano_", cmp_name, ".json")))

    # Pathway enrichment (KEGG)
    pathway <- build_pathway_json(cmp_name, kegg_dir)
    write_json_export(pathway, file.path(web_dir, paste0("pathway_enrichment_", cmp_name, ".json")))

    # Expression matrix
    expr <- build_expression_json(cmp_name, dr, object2, app3)
    write_json_export(expr, file.path(web_dir, paste0("expression_", cmp_name, ".json")))

    # OPLS-DA (if available)
    if (!is.null(oplsda_dir) && dir.exists(oplsda_dir)) {
      oplsda <- build_oplsda_json(oplsda_dir, cmp_name)
      write_json_export(oplsda, file.path(web_dir, paste0("oplsda_", cmp_name, ".json")))
    }

    # Update manifest treatment entry with available analyses
    if (cmp_name %in% names(manifest$treatments)) {
      manifest$treatments[[cmp_name]]$has_msea <- has_msea && !is.null(msea_results[[cmp_name]])
      manifest$treatments[[cmp_name]]$has_pathway_topology <- has_topology && !is.null(topology_results[[cmp_name]])
      manifest$treatments[[cmp_name]]$has_biomarker <- has_biomarker && !is.null(biomarker_results[[cmp_name]])
      manifest$treatments[[cmp_name]]$has_chem_class <- has_chem_class && !is.null(chem_class_results[[cmp_name]])
    }
  }

  # ── 6. MetaboAnalystR Advanced (per-comparison) ─────────────────────────
  if (has_msea || has_topology || has_biomarker || has_chem_class) {
    cat("[6/7] Building MetaboAnalystR files...\n")

    for (cmp_name in names(diff_results)) {
      # MSEA
      if (has_msea) {
        msea <- build_msea_json(msea_results, cmp_name)
        write_json_export(msea, file.path(web_dir, paste0("msea_", cmp_name, ".json")))
      }

      # Pathway Topology
      if (has_topology) {
        topology <- build_topology_json(topology_results, cmp_name)
        write_json_export(topology, file.path(web_dir, paste0("pathway_topology_", cmp_name, ".json")))
      }

      # Biomarker
      if (has_biomarker) {
        biomarker <- build_biomarker_json(biomarker_results, cmp_name)
        write_json_export(biomarker, file.path(web_dir, paste0("biomarker_", cmp_name, ".json")))
      }

      # Chemical Class
      if (has_chem_class) {
        chem_class <- build_chem_class_json(chem_class_results, cmp_name)
        write_json_export(chem_class, file.path(web_dir, paste0("chem_class_", cmp_name, ".json")))
      }
    }
  }

  # ── 7. Re-write manifest with updated flags ─────────────────────────────
  cat("[7/7] Finalizing manifest.json...\n")
  write_json_export(manifest, file.path(web_dir, "manifest.json"))

  # ── Summary ─────────────────────────────────────────────────────────────
  n_metabolites <- metabolites$n_metabolites
  n_comparisons <- length(diff_results)
  n_files <- length(list.files(web_dir, pattern = "\\.json$"))

  cat(sprintf("\n=> Web export complete!\n"))
  cat(sprintf("   Directory: %s\n", web_dir))
  cat(sprintf("   JSON files: %d\n", n_files))
  cat(sprintf("   Comparisons: %d\n", n_comparisons))
  cat(sprintf("   Metabolites: %d\n", n_metabolites))
  cat(sprintf("   PCA samples: %d\n", pca_data$n_samples %||% 0))
  cat(sprintf("   MetaboAnalystR: MSEA=%s, Topology=%s, Biomarker=%s, ChemClass=%s\n",
              has_msea, has_topology, has_biomarker, has_chem_class))
  cat("========================================\n")

  invisible(list(
    manifest     = manifest,
    metabolites  = metabolites,
    pca          = pca_data,
    summary      = summary_data,
    web_dir      = web_dir
  ))
}