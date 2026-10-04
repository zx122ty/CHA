#!/usr/bin/env Rscript
###############################################################################
# metabolomics_pipeline.R — Config-driven Untargeted Metabolomics Pipeline
#
# Refactored from FuZi_QinPi_synergy_analysis.R to be fully configurable via
# an external YAML configuration file. Supports any study design with MS1/MS2
# separation, limma differential analysis, OPLS-DA, annotation, KEGG enrichment,
# and overlap analysis.
#
# Enhanced with advanced compound identification modules:
#   - SIRIUS (molecular formula prediction + CSI:FingerID structure search)
#   - MetFrag (in silico fragmentation validation)
#   - CAMERA (adduct/isotope grouping)
#   - MSI confidence level scoring (Level 1-4)
#   - Extended MS2 spectral libraries (GNPS, MoNA, etc.)
#   - MS1-only file isolation for clean feature extraction
#
# Usage:
#   Rscript metabolomics_pipeline.R --config /path/to/config.yaml
#   Rscript metabolomics_pipeline.R --config /path/to/config.yaml --validate-only
#
# Dependencies:
#   Required:  tidymass, tidyverse, yaml, optparse, limma, ropls, openxlsx,
#              ggplot2, ggrepel, ggsci, ggpubr, pheatmap, VennDiagram, UpSetR
#   Optional:  CAMERA, xcms (Bioconductor, for adduct grouping)
#              metfRag (for in silico fragmentation)
#   External:  SIRIUS CLI (Java, for formula/structure prediction)
###############################################################################

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(tidymass)
  library(tidyverse)
  library(openxlsx)
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(ggsci)
  library(ggpubr)
  library(pheatmap)
  library(VennDiagram)
  library(ropls)
})

# ══════════════════════════════════════════════════════════════════════════════
# Utility functions
# ══════════════════════════════════════════════════════════════════════════════

#' Null-coalescing operator
`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

#' Directory of this script, captured before any setwd().
#'
#' load_and_validate_config() chdirs to cfg$project$working_dir early, so a
#' later normalizePath() of the relative --file argument would resolve against
#' the wrong directory whenever working_dir != the script's directory, and
#' every optional module would be silently skipped.
PIPELINE_SCRIPT_DIR <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE))
  } else {
    getwd()
  }
})

#' Compute QC-derived quality metrics for a batch
#'
#' QC samples are repeated injections of one pooled reference, so they measure
#' technical reproducibility without confounding by biology. These metrics are
#' recorded in the batch manifest and gate cross-batch integration: a batch
#' whose QC-RSD median far exceeds the others cannot support quantitative
#' comparison, even when its m/z-based annotation is still usable.
#'
#' @param expression_data Features x samples matrix (pre-filter peak table)
#' @param sample_classes Character vector of class labels. May be named by
#'   sample_id or supplied in the same order as the matrix columns.
#' @param qc_name Class label identifying MS1 QC samples
#' @param rsd_threshold RSD % above which a feature is considered unreliable
#' @return List of metrics, or NULL when fewer than 3 QC samples are present
compute_qc_quality <- function(expression_data, sample_classes, qc_name,
                               rsd_threshold = 30) {
  # Callers may pass a bare factor/vector (as dplyr pipelines do), in which case
  # names() is empty. Fall back to positional alignment with the matrix columns
  # rather than silently returning NULL and losing the whole quality gate.
  if (is.null(names(sample_classes)) &&
      length(sample_classes) == ncol(expression_data)) {
    names(sample_classes) <- colnames(expression_data)
  }
  qc_samples <- names(sample_classes)[
    !is.na(sample_classes) & sample_classes == qc_name]
  qc_samples <- intersect(qc_samples, colnames(expression_data))
  if (length(qc_samples) < 3) return(NULL)

  m <- as.matrix(expression_data[, qc_samples, drop = FALSE])
  # Match the pipeline's own RSD definition (measured values only, >= 3 points)
  rsd <- apply(m, 1, function(x) {
    x <- x[!is.na(x) & x > 0]
    if (length(x) < 3) return(NA_real_)
    sd(x) / mean(x) * 100
  })
  rsd <- rsd[is.finite(rsd)]
  if (length(rsd) == 0) return(NULL)

  # Correlation needs a complete matrix; NA pairs would otherwise collapse the
  # estimate. Minimum-value imputation matches what the pipeline does downstream.
  # Features that are entirely NA across QC cannot be imputed (min() of an empty
  # vector is Inf, which would poison the whole correlation matrix) so they are
  # dropped first — they carry no reproducibility information by definition.
  mi <- m[rowSums(!is.na(m)) > 0, , drop = FALSE]
  for (i in seq_len(nrow(mi))) {
    v <- mi[i, ]
    if (anyNA(v)) {
      mn <- suppressWarnings(min(v, na.rm = TRUE))
      if (!is.finite(mn)) next
      mi[i, is.na(v)] <- mn
    }
  }
  cm <- suppressWarnings(cor(mi, use = "pairwise.complete.obs"))
  # A QC sample that is constant across all features (or entirely NA) yields an
  # all-NA column in the matrix; drop such columns before taking the median so a
  # single degenerate sample cannot erase the whole estimate.
  cm <- cm[!apply(cm, 1, function(r) all(is.na(r))),
           !apply(cm, 2, function(r) all(is.na(r))), drop = FALSE]
  qc_corr <- if (nrow(cm) >= 2) {
    suppressWarnings(median(cm[upper.tri(cm)], na.rm = TRUE))
  } else NA_real_
  if (!is.finite(qc_corr)) qc_corr <- NA_real_

  list(
    n_qc_samples      = length(qc_samples),
    n_features        = nrow(m),
    qc_rsd_median     = round(median(rsd), 1),
    qc_rsd_pass_rate  = round(mean(rsd <= rsd_threshold), 3),
    qc_qc_correlation = round(qc_corr, 3)
  )
}

#' Decide whether a batch's quantification is trustworthy enough to integrate
#'
#' Thresholds are deliberately loose: they catch a batch that is an order of
#' magnitude worse than its peers, not one that is merely average.
#'
#' @param qc_quality Output of compute_qc_quality()
#' @param max_rsd_median Reject when median QC-RSD exceeds this (%)
#' @param min_qc_corr Reject when QC-QC correlation falls below this
#' @return List(valid = logical, reasons = character vector)
assess_quantification_quality <- function(qc_quality,
                                          max_rsd_median = 40,
                                          min_qc_corr = 0.85) {
  if (is.null(qc_quality)) {
    return(list(valid = FALSE, reasons = "no QC samples available"))
  }
  reasons <- character(0)
  if (is.finite(qc_quality$qc_rsd_median) &&
      qc_quality$qc_rsd_median > max_rsd_median) {
    reasons <- c(reasons, sprintf("QC-RSD median %.1f%% exceeds %.0f%%",
                                  qc_quality$qc_rsd_median, max_rsd_median))
  }
  if (is.finite(qc_quality$qc_qc_correlation) &&
      qc_quality$qc_qc_correlation < min_qc_corr) {
    reasons <- c(reasons, sprintf("QC-QC correlation %.3f below %.2f",
                                  qc_quality$qc_qc_correlation, min_qc_corr))
  }
  list(valid = length(reasons) == 0, reasons = reasons)
}

#' Assign MSI confidence level based on available evidence
#'
#' Follows the Metabolomics Standards Initiative (MSI) guidelines:
#'   Level 1 — Confirmed:      MS2 + RT match against in-house database
#'   Level 2 — Putatively identified: MS2 spectral match against external database
#'   Level 3 — Putatively annotated:  MS1 accurate mass match only
#'   Level 4 — Unknown:        no match in any database
#'
#' @param has_ms2_match Logical, whether the feature has an MS2 spectral match
#' @param has_ms1_match Logical, whether the feature has an MS1 mass match
#' @param has_rt_match  Logical, whether RT was matched (reserved for future use)
#' @param is_inhouse    Logical, whether the match is from an in-house database
#' @param is_sirius     Logical, whether SIRIUS provided formula/structure evidence
#' @param is_metfrag    Logical, whether MetFrag validated the fragmentation
#' @return Character string with MSI confidence level
assign_confidence_level <- function(has_ms2_match, has_ms1_match,
                                     has_rt_match = FALSE, is_inhouse = FALSE,
                                     is_sirius = FALSE, is_metfrag = FALSE) {
  if (has_ms2_match && has_rt_match && is_inhouse) {
    return("Level 1 — Confirmed")
  } else if (is_metfrag && has_ms2_match) {
    return("Level 1 — Confirmed")
  } else if (has_ms2_match) {
    return("Level 2 — Putatively identified")
  } else if (is_sirius) {
    return("Level 2 — Putatively identified")
  } else if (has_ms1_match) {
    return("Level 3 — Putatively annotated")
  } else {
    return("Level 4 — Unknown")
  }
}

#' Save a ggplot object as both PDF and PNG
#'
#' @param filename Output filename without extension
#' @param plot ggplot object (defaults to last_plot())
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Resolution for PNG output
save_plot <- function(filename, plot = NULL, width = 8, height = 7, dpi = 300) {
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  ggsave(paste0(filename, ".pdf"), plot = plot, width = width, height = height)
  ggsave(paste0(filename, ".png"), plot = plot, width = width, height = height, dpi = dpi)
  cat(sprintf("  => Saved: %s (.pdf/.png)\n", filename))
}

# ══════════════════════════════════════════════════════════════════════════════
# 0. CLI & Config Loading
# ══════════════════════════════════════════════════════════════════════════════

#' Parse command-line arguments
parse_cli <- function() {
  option_list <- list(
    make_option(c("-c", "--config"), type = "character", default = NULL,
                help = "Path to YAML configuration file", metavar = "FILE"),
    make_option(c("--validate-only"), action = "store_true", default = FALSE,
                help = "Validate config and paths only, do not run pipeline"),
    make_option(c("--qualitative-only"), dest = "qualitative_only",
                action = "store_true", default = FALSE,
                help = paste("Skip every intensity-dependent step and keep only",
                             "analyses valid on a degraded batch (annotation +",
                             "QC diagnostic PCA). For batches whose QC-RSD makes",
                             "quantification unreliable."))
  )
  opt <- parse_args(OptionParser(
    option_list = option_list,
    usage = paste("Rscript metabolomics_pipeline.R --config <FILE>",
                  "[--validate-only] [--qualitative-only]")
  ))
  if (is.null(opt$config)) {
    stop("--config is required. Use --help for usage.")
  }
  opt
}

#' Load and validate the YAML configuration file
#'
#' Checks that all required top-level keys are present, paths exist,
#' and threshold values are in valid ranges. Stops with clear messages on failure.
#'
#' @param config_path Path to the YAML config file
#' @return A validated configuration list
load_and_validate_config <- function(config_path) {
  if (!file.exists(config_path)) {
    stop(sprintf("Config file not found: %s", config_path))
  }
  cat(sprintf("===== Loading config: %s =====\n", normalizePath(config_path)))

  cfg <- read_yaml(config_path)

  # ── Required top-level keys ──────────────────────────────────────────────
  required_keys <- c("project", "data", "metadata_columns", "sample_roles",
                     "feature_extraction", "filtering", "normalization",
                     "differential", "annotation", "visualization")
  missing_keys <- setdiff(required_keys, names(cfg))
  if (length(missing_keys) > 0) {
    stop(sprintf("Config missing required keys: %s",
                 paste(missing_keys, collapse = ", ")))
  }

  # ── Set working directory ────────────────────────────────────────────────
  wd <- cfg$project$working_dir
  if (is.null(wd) || !dir.exists(wd)) {
    stop(sprintf("Working directory not found: %s", wd))
  }
  setwd(wd)
  cat(sprintf("=> Working directory: %s\n", wd))
  cat(sprintf("=> Project: %s\n", cfg$project$name))

  # ── Validate file paths ──────────────────────────────────────────────────
  metadata_path <- cfg$data$metadata_file
  if (is.null(metadata_path) || !file.exists(metadata_path)) {
    stop(sprintf("Metadata file not found: %s", metadata_path))
  }
  cat(sprintf("=> Metadata: %s (%d rows)\n", metadata_path,
              nrow(read.csv(metadata_path, stringsAsFactors = FALSE))))

  raw_dir <- cfg$data$raw_mzml_dir
  if (is.null(raw_dir) || !dir.exists(raw_dir)) {
    stop(sprintf("Raw mzML directory not found: %s", raw_dir))
  }
  mzml_count <- length(list.files(raw_dir, pattern = "\\.mzML$",
                                   ignore.case = TRUE, recursive = TRUE))
  cat(sprintf("=> Raw mzML dir: %s (%d .mzML files)\n", raw_dir, mzml_count))
  if (mzml_count == 0) {
    stop(sprintf("No .mzML files found in: %s", raw_dir))
  }

  # ── Validate MS level column mapping ─────────────────────────────────────
  ms_col <- cfg$metadata_columns$ms_level
  if (!is.null(ms_col) && ms_col != "") {
    cat(sprintf("=> MS level column: '%s' (MS1/MS2 separation enabled)\n", ms_col))
  } else {
    cat("=> MS level column not specified; treating all samples as MS1\n")
  }

  # ── Validate thresholds ──────────────────────────────────────────────────
  f <- cfg$filtering
  if (is.null(f$blank_fold_change) || f$blank_fold_change < 0) {
    stop("filtering.blank_fold_change must be >= 0")
  }
  if (is.null(f$rsd_threshold) || f$rsd_threshold < 0 || f$rsd_threshold > 100) {
    stop("filtering.rsd_threshold must be between 0 and 100")
  }
  cat(sprintf("=> Polarity: %s\n", cfg$feature_extraction$polarity))
  cat(sprintf("=> Reference group: %s\n", cfg$differential$reference_group))
  cat(sprintf("=> Config validation passed.\n"))

  cfg
}

# ══════════════════════════════════════════════════════════════════════════════
# 1. Metadata Loading & Sample Info Construction
# ══════════════════════════════════════════════════════════════════════════════

#' Load metadata CSV using column-name mapping from config
#'
#' Reads the CSV, maps columns via metadata_columns config, separates MS1/MS2.
#'
#' @param cfg Validated config list
#' @return List with components: ms1 (data.frame), ms2 (data.frame), all (data.frame)
load_metadata <- function(cfg) {
  cat("\n===== Loading metadata =====\n")

  meta <- read.csv(cfg$data$metadata_file, stringsAsFactors = FALSE)

  col_map <- cfg$metadata_columns
  required_cols <- c("filename", "group_id")
  for (rc in required_cols) {
    cn <- col_map[[rc]]
    if (is.null(cn) || !(cn %in% colnames(meta))) {
      stop(sprintf("Required metadata column '%s' (mapped as '%s') not found in CSV. Available columns: %s",
                   rc, cn, paste(colnames(meta), collapse = ", ")))
    }
  }

  # Standardize column names internally
  result <- list()
  result$all <- meta

  # MS1/MS2 separation
  ms_col <- col_map$ms_level
  if (!is.null(ms_col) && ms_col %in% colnames(meta)) {
    ms1 <- meta[meta[[ms_col]] == "MS", , drop = FALSE]
    ms2 <- meta[meta[[ms_col]] == "MSMS", , drop = FALSE]
    cat(sprintf("   MS1 samples: %d, MS2 samples: %d\n", nrow(ms1), nrow(ms2)))
  } else {
    ms1 <- meta
    ms2 <- meta[0, , drop = FALSE]
    cat(sprintf("   MS1 samples: %d (MS_level column not present; treating all as MS1)\n", nrow(ms1)))
  }

  result$ms1 <- ms1
  result$ms2 <- ms2
  result$col_map <- col_map
  result
}

#' Build sample_info data.frame for tidymass from metadata
#'
#' Assigns sample class labels from the group_id column. MS2 files get
#' special class tags (e.g., "QC-MSMS") to ensure they are excluded from
#' MS1-only steps (RSD filtering, differential analysis).
#'
#' @param meta_list Output of load_metadata()
#' @param cfg Validated config list
#' @return data.frame with columns: sample_id, class, injection.order
build_sample_info_all <- function(meta_list, cfg) {
  cat("\n===== Building sample info =====\n")

  col_map <- meta_list$col_map
  fn_col  <- col_map$filename
  grp_col <- col_map$group_id
  ms_col  <- col_map$ms_level

  # Combine MS1 and MS2
  meta <- meta_list$all

  # Strip .mzML extension for sample_id
  sample_ids <- gsub("\\.mzML$", "", meta[[fn_col]], ignore.case = TRUE)

  # Base class from group_id
  classes <- as.character(meta[[grp_col]])

  # Optional: split treatment groups by a second column (e.g. dose High/Low).
  # Only non-special samples get the suffix, so CT1/QC/Blank keep their names.
  suffix_col <- col_map$group_suffix
  if (!is.null(suffix_col) && suffix_col %in% colnames(meta)) {
    special <- c(cfg$sample_roles$control, cfg$sample_roles$qc_ms1, cfg$sample_roles$blank)
    suffix_vals <- as.character(meta[[suffix_col]])
    has_suffix <- !(classes %in% special) & !is.na(suffix_vals) & suffix_vals != ""
    classes[has_suffix] <- paste0(classes[has_suffix], "_", suffix_vals[has_suffix])
    cat(sprintf("   Group suffix from '%s': %d treatment classes\n",
                suffix_col, length(unique(classes[has_suffix]))))
  }

  # Append -MSMS suffix for MS2 samples
  if (!is.null(ms_col) && ms_col %in% colnames(meta)) {
    is_ms2 <- meta[[ms_col]] == "MSMS"
    classes[is_ms2] <- paste0(classes[is_ms2], "-MSMS")
  }

  sample_info <- data.frame(
    sample_id       = sample_ids,
    class           = classes,
    injection.order = seq_len(nrow(meta)),
    stringsAsFactors = FALSE
  )

  cat("   Sample types:\n")
  print(table(sample_info$class))

  sample_info
}

#' Generate differential comparisons from config
#'
#' If auto_generate_comparisons is true, creates "treat_vs_ref" for every
#' unique treatment group (excluding control, QC, and Blank).
#' Merges with any explicit extra_comparisons from the config.
#'
#' @param cfg Validated config list
#' @param all_classes Character vector of all unique sample classes
#' @return List of comparison objects, each with name, ctrl, treat
generate_comparisons <- function(cfg, all_classes) {
  cat("\n===== Generating comparisons =====\n")

  ref     <- cfg$differential$reference_group
  control <- cfg$sample_roles$control
  qc_ms1  <- cfg$sample_roles$qc_ms1
  blank   <- cfg$sample_roles$blank

  # Identify treatment groups: all non-special classes
  special <- c(control, qc_ms1, blank)
  # Also exclude any class ending in -MSMS
  treatment_groups <- setdiff(all_classes, special)
  treatment_groups <- treatment_groups[!grepl("-MSMS$", treatment_groups)]

  comparisons <- list()

  if (isTRUE(cfg$differential$auto_generate_comparisons)) {
    for (trt in sort(treatment_groups)) {
      cmp_name <- sprintf("%s_vs_%s", trt, ref)
      comparisons[[cmp_name]] <- list(
        name  = cmp_name,
        ctrl  = ref,
        treat = trt
      )
    }
    cat(sprintf("   Auto-generated %d comparisons (%s_vs_%s)\n",
                length(comparisons), "<treatment>", ref))
  }

  # Add explicit extra comparisons
  extras <- cfg$differential$extra_comparisons
  if (is.list(extras) && length(extras) > 0) {
    for (ex in extras) {
      if (is.list(ex) && !is.null(ex$name)) {
        comparisons[[ex$name]] <- ex
        cat(sprintf("   Extra comparison: %s (%s vs %s)\n", ex$name, ex$treat, ex$ctrl))
      }
    }
  }

  # Build display labels (use index-based loop to modify list IN PLACE)
  labels <- cfg$visualization$group_labels
  for (i in seq_along(comparisons)) {
    t_label <- if (!is.null(labels[[comparisons[[i]]$treat]])) labels[[comparisons[[i]]$treat]] else comparisons[[i]]$treat
    c_label <- if (!is.null(labels[[comparisons[[i]]$ctrl]]))  labels[[comparisons[[i]]$ctrl]]  else comparisons[[i]]$ctrl
    comparisons[[i]]$label <- sprintf("%s vs %s", t_label, c_label)
  }

  cat(sprintf("   Total comparisons: %d\n", length(comparisons)))
  # Convert to simple list-of-lists (not named for easier iteration)
  unname(comparisons)
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. Step 1: Feature Extraction (with MS1-only isolation)
# ══════════════════════════════════════════════════════════════════════════════

#' Run massprocesser-based feature extraction (XCMS peak picking)
#'
#' Uses symlink-based MS1-only file isolation to prevent MS/MS DDA files
#' from polluting the peak table. Writes sample_info.csv (MS1 only) to
#' the raw data directory. Skips processing if the peak table already exists.
#'
#' @param cfg Validated config list
#' @param meta_list Output of load_metadata()
run_feature_extraction <- function(cfg, meta_list) {
  cat("\n===== Step 1: Feature extraction (MS1-only peak picking) =====\n")

  out_dir <- cfg$data$peak_table_dir
  peak_file <- file.path(out_dir, "Peak_table_for_cleaning.csv")

  if (file.exists(peak_file)) {
    cat(sprintf("=> Found existing peak table: %s\n=> Skipping feature extraction.\n",
                peak_file))
    return(invisible(TRUE))
  }

  cat("=> Starting feature extraction from mzML files (MS1 only)...\n")

  # Write sample_info.csv with MS1 samples only — critical MS1/MS2 gate
  meta_ms1 <- meta_list$ms1
  col_map  <- meta_list$col_map
  fn_col   <- col_map$filename
  grp_col  <- col_map$group_id
  suffix_col <- col_map$group_suffix

  classes_ms1 <- as.character(meta_ms1[[grp_col]])
  if (!is.null(suffix_col) && suffix_col %in% colnames(meta_ms1)) {
    special <- c(cfg$sample_roles$control, cfg$sample_roles$qc_ms1, cfg$sample_roles$blank)
    suffix_vals <- as.character(meta_ms1[[suffix_col]])
    has_suffix <- !(classes_ms1 %in% special) & !is.na(suffix_vals) & suffix_vals != ""
    classes_ms1[has_suffix] <- paste0(classes_ms1[has_suffix], "_", suffix_vals[has_suffix])
  }

  sample_info_ms1 <- data.frame(
    sample_id       = gsub("\\.mzML$", "", meta_ms1[[fn_col]], ignore.case = TRUE),
    class           = classes_ms1,
    injection.order = seq_len(nrow(meta_ms1)),
    stringsAsFactors = FALSE
  )

  raw_dir <- cfg$data$raw_mzml_dir

  # ── MS1-only file isolation via symlinks into group-specific subfolders ─────────
  # massprocesser::process_data(path=".") groups samples by parent folder name:
  # unlist(lapply(str_split(f.in, "/"), function(x) x[length(x)-1])).
  # If symlinks are placed directly in tmp_ms1_only, all samples become group ".",
  # and min_fraction = 0.5 requires a feature to appear in >=50% of the entire cohort!
  # Fix: organize symlinks into subfolders by their biological class (e.g. 10_High, CT1, QC).
  # This allows xcms::PeakDensityParam to check minFraction WITHIN EACH GROUP (e.g. 2 of 3 replicates).
  ms1_filenames <- paste0(sample_info_ms1$sample_id, ".mzML")
  ms1_existing  <- ms1_filenames[file.exists(file.path(raw_dir, ms1_filenames))]
  cat(sprintf("=> MS1 .mzML files found: %d / %d named in sample_info\n",
              length(ms1_existing), length(ms1_filenames)))

  tmp_ms1_dir <- file.path(raw_dir, "tmp_ms1_only")
  if (dir.exists(tmp_ms1_dir)) unlink(tmp_ms1_dir, recursive = TRUE)
  dir.create(tmp_ms1_dir, showWarnings = FALSE, recursive = TRUE)

  # Create subdirectories per class and link mzML files inside them
  for (i in seq_along(ms1_existing)) {
    fn <- ms1_existing[i]
    sid <- gsub("\\.mzML$", "", fn, ignore.case = TRUE)
    cls <- sample_info_ms1$class[sample_info_ms1$sample_id == sid][1]
    if (is.na(cls) || cls == "") cls <- "Subject"
    grp_dir <- file.path(tmp_ms1_dir, cls)
    if (!dir.exists(grp_dir)) dir.create(grp_dir, showWarnings = FALSE, recursive = TRUE)
    file.symlink(file.path(raw_dir, fn), file.path(grp_dir, fn))
  }
  cat(sprintf("=> Created temp MS1-only dir: %s with %d class subdirectories (%d symlinks)\n",
              tmp_ms1_dir, length(unique(sample_info_ms1$class)), length(ms1_existing)))

  # Write sample_info.csv into the temp dir as well
  write.csv(sample_info_ms1,
            file.path(tmp_ms1_dir, "sample_info.csv"),
            row.names = FALSE)

  # Run massprocesser on the MS1-only temp directory
  orig_wd <- getwd()
  setwd(tmp_ms1_dir)
  cat(sprintf("=> Changed working directory to: %s\n", tmp_ms1_dir))

  fe <- cfg$feature_extraction
  tryCatch({
    massprocesser::process_data(
      path                     = ".",
      polarity                 = fe$polarity,
      ppm                      = fe$ppm,
      peakwidth                = unlist(fe$peakwidth),
      snthresh                 = fe$snthresh,
      noise                    = fe$noise,
      threads                  = fe$threads,
      output_tic               = isTRUE(fe$output_tic),
      output_bpc               = isTRUE(fe$output_bpc),
      output_rt_correction_plot = isTRUE(fe$output_rt_correction_plot),
      min_fraction             = fe$min_fraction,
      fill_peaks               = isTRUE(fe$fill_peaks)
    )
  }, error = function(e) {
    # If massprocesser::process_data dies during peak table writing due to dplyr::select column mismatch,
    # recover cleanly from xdata3 intermediate data
    cat(sprintf("=> massprocesser encountered: %s. Generating peak table directly from xdata3...\n", conditionMessage(e)))
    xdata3_path <- file.path(tmp_ms1_dir, "Result/intermediate_data/xdata3")
    if (file.exists(xdata3_path)) {
      load(xdata3_path)
      values <- xcms::featureValues(xdata3, value = "into")
      definition <- xcms::featureDefinitions(object = xdata3)
      definition <- definition[, -ncol(definition)]
      definition <- definition[names(definition) != "peakidx"]
      definition <- definition@listData %>% do.call(cbind, .) %>% as.data.frame()
      peak_name <- xcms::groupnames(xdata3)
      peak_name <- paste(peak_name, ifelse(fe$polarity == "positive", "POS", "NEG"), sep = "_")
      colnames(values) <- stringr::str_replace(string = colnames(values), pattern = "\\.mz[X|x]{0,1}[M|m][L|l]", replacement = "")

      sample_group_cols <- setdiff(colnames(definition), c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax", "npeaks"))
      peak_table_for_cleaning <- definition %>%
        dplyr::select(-dplyr::any_of(c("mzmin", "mzmax", "rtmin", "rtmax", "npeaks", sample_group_cols))) %>%
        dplyr::rename(mz = mzmed, rt = rtmed) %>%
        data.frame(variable_id = peak_name, ., values, stringsAsFactors = FALSE, check.names = FALSE)

      peak_table <- data.frame(peak.name = peak_name, definition, values, stringsAsFactors = FALSE, check.names = FALSE)
      readr::write_csv(peak_table, file = file.path(tmp_ms1_dir, "Result/Peak_table.csv"))
      readr::write_csv(peak_table_for_cleaning, file = file.path(tmp_ms1_dir, "Result/Peak_table_for_cleaning.csv"))
      cat("=> Successfully generated Peak_table.csv and Peak_table_for_cleaning.csv\n")
    } else {
      stop(e)
    }
  })

  setwd(orig_wd)

  # Move Result/ from temp MS1-only dir to configured peak_table_dir
  result_in_raw <- file.path(tmp_ms1_dir, "Result")
  target_dir <- cfg$data$peak_table_dir
  if (dir.exists(result_in_raw) &&
      normalizePath(target_dir) != normalizePath(result_in_raw)) {
    dir.create(target_dir, showWarnings = FALSE, recursive = TRUE)
    result_files <- list.files(result_in_raw, full.names = TRUE)
    file.copy(result_files, target_dir, recursive = TRUE, overwrite = TRUE)
    unlink(result_in_raw, recursive = TRUE)
    cat(sprintf("=> Moved Result/ to: %s\n", target_dir))
  }

  # Clean up temp MS1-only directory
  unlink(tmp_ms1_dir, recursive = TRUE)
  cat("=> Cleaned up temp MS1-only directory\n")

  cat("=> Feature extraction complete\n")
  invisible(TRUE)
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. Step 2: Preprocessing (Blank filter → Missing value → QC-RSD → Impute → SVR)
# ══════════════════════════════════════════════════════════════════════════════

#' Extract internal-standard QC stats and remove IS features from the peak table.
#'
#' Runs BEFORE any filtering: the IS mix (isotope-labeled amino acids, 2 uM) is
#' present in blanks too, so a future blank filter would treat IS features as
#' background and drop them before QC stats could be computed.
#'
#' @param expression_data features x samples matrix (raw peak areas)
#' @param variable_info data.frame with variable_id, mz, rt (rt in seconds)
#' @param sample_info_all data.frame with sample_id and class
#' @param cfg Validated config list
#' @return list(expression_data, variable_info, is_qc_stats) — is_qc_stats is
#'   NULL when the block is disabled or nothing matched
extract_and_remove_internal_standards <- function(expression_data, variable_info,
                                                  sample_info_all, cfg) {
  is_cfg <- cfg$internal_standards
  if (is.null(is_cfg) || !isTRUE(is_cfg$enabled)) {
    return(list(expression_data = expression_data, variable_info = variable_info,
                is_qc_stats = NULL))
  }

  if (is.null(is_cfg$file) || !file.exists(is_cfg$file)) {
    cat(sprintf("   [IS] Reference file not found: %s — skipping\n",
                is_cfg$file %||% "<NULL>"))
    return(list(expression_data = expression_data, variable_info = variable_info,
                is_qc_stats = NULL))
  }

  cat("\n----- Step 2-pre: Internal standards -----\n")
  is_ref <- read.csv(is_cfg$file, stringsAsFactors = FALSE)
  qc_samples <- sample_info_all$sample_id[sample_info_all$class == cfg$sample_roles$qc_ms1]
  all_samples <- sample_info_all$sample_id[!grepl("-MSMS$", sample_info_all$class)]

  mz_tol_ppm <- is_cfg$mz_tol_ppm %||% 15
  rt_tol_sec <- is_cfg$rt_tol_sec %||% 20

  matched_ids <- character(0)
  qc_stats <- list()

  for (i in seq_len(nrow(is_ref))) {
    ref <- is_ref[i, ]
    mz_tol <- ref$mz * mz_tol_ppm / 1e6
    hits <- which(abs(variable_info$mz - ref$mz) <= mz_tol &
                  abs(variable_info$rt - ref$rt_sec) <= rt_tol_sec)

    if (length(hits) == 0) {
      cat(sprintf("   [IS] %s: no match (mz=%.4f, rt=%.1fs) — skipped\n",
                  ref$name, ref$mz, ref$rt_sec))
      next
    }
    if (length(hits) > 1) {
      dist <- abs(variable_info$mz[hits] - ref$mz) / ref$mz +
              abs(variable_info$rt[hits] - ref$rt_sec) / rt_tol_sec
      hits <- hits[order(dist)]
      cat(sprintf("   [IS] %s: %d candidates matched, using nearest (%s)\n",
                  ref$name, length(hits), variable_info$variable_id[hits[1]]))
    }

    vid <- variable_info$variable_id[hits[1]]
    row_vals <- as.numeric(expression_data[hits[1], ])
    names(row_vals) <- colnames(expression_data)

    qc_vals  <- row_vals[intersect(qc_samples, names(row_vals))]
    all_vals <- row_vals[intersect(all_samples, names(row_vals))]
    calc_rsd <- function(x) {
      x <- x[!is.na(x) & x > 0]
      if (length(x) < 3) return(NA_real_)
      sd(x) / mean(x) * 100
    }
    n_detected <- sum(!is.na(all_vals) & all_vals > 0)

    qc_stats[[ref$name]] <- list(
      name               = ref$name,
      variable_id        = vid,
      mz                 = ref$mz,
      rt_expected_sec    = ref$rt_sec,
      rt_observed_sec    = variable_info$rt[hits[1]],
      mean_area_qc       = if (length(qc_vals)) mean(qc_vals, na.rm = TRUE) else NA_real_,
      rsd_pct_qc         = calc_rsd(qc_vals),
      mean_area_all      = if (length(all_vals)) mean(all_vals, na.rm = TRUE) else NA_real_,
      rsd_pct_all        = calc_rsd(all_vals),
      n_samples_detected = n_detected,
      n_samples_total    = length(all_vals),
      detection_rate     = if (length(all_vals)) n_detected / length(all_vals) else NA_real_
    )
    matched_ids <- c(matched_ids, vid)
  }

  cat(sprintf("   [IS] Matched %d/%d internal standards; removing from feature table\n",
              length(matched_ids), nrow(is_ref)))

  keep <- !(variable_info$variable_id %in% matched_ids)
  list(
    expression_data = expression_data[keep, , drop = FALSE],
    variable_info   = variable_info[keep, , drop = FALSE],
    is_qc_stats     = if (length(qc_stats) > 0) qc_stats else NULL
  )
}

#' Read the acquisition start timestamp from an mzML file header.
#'
#' `<run startTimeStamp="...">` sits ~18 KB into the file, so a bounded read is
#' enough — no need to parse the whole (often hundreds of MB) file.
read_mzml_start_time <- function(path, n_bytes = 131072) {
  if (!file.exists(path)) return(NA_real_)
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  raw_bytes <- readBin(con, "raw", n = n_bytes)
  txt <- rawToChar(raw_bytes[raw_bytes != as.raw(0)])
  hit <- regmatches(txt, regexpr('startTimeStamp="[^"]*"', txt))
  if (length(hit) == 0) return(NA_real_)
  ts <- sub('startTimeStamp="', "", sub('"$', "", hit))
  as.numeric(as.POSIXct(ts, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC"))
}

#' Derive the true acquisition order for the samples in the peak table.
#'
#' The peak-table column order is the alphabetically sorted file listing, NOT
#' the acquisition order: it places every QC sample at one end of the run, so
#' the SVR/LOESS QC-drift model would be fitted on a narrow x-range and
#' extrapolate onto all biological samples. Read the real order from the mzML
#' headers; fall back to the given order (with a warning) if unavailable.
get_injection_order <- function(sample_ids, raw_dir) {
  ts <- vapply(sample_ids, function(sid) {
    tryCatch(read_mzml_start_time(file.path(raw_dir, paste0(sid, ".mzML"))),
             error = function(e) NA_real_)
  }, numeric(1))

  if (all(is.na(ts))) {
    cat("   ⚠️  No mzML acquisition timestamps found; using peak-table order for drift correction\n")
    return(seq_along(sample_ids))
  }
  n_missing <- sum(is.na(ts))
  if (n_missing > 0) {
    cat(sprintf("   ⚠️  %d sample(s) have no acquisition timestamp; ordered last\n", n_missing))
  }
  # rank() maps each sample to its position in the run (order() would return
  # the permutation instead, mis-assigning every sample's injection order).
  rank(ts, na.last = TRUE)
}

#' Full preprocessing pipeline: load peak table, filter, impute, correct, subset
#'
#' Returns a list with: object2 (mass_dataset), sample_info (core groups only),
#' variable_info, expression_data (matrix), core_groups (vector).
#'
#' @param cfg Validated config list
#' @param meta_list Output of load_metadata()
#' @param sample_info_all Output of build_sample_info_all()
run_preprocessing <- function(cfg, meta_list, sample_info_all,
                             qualitative_only = FALSE) {
  cat("\n===== Step 2: Data preprocessing =====\n")

  # ── Load peak table ──────────────────────────────────────────────────────
  peak_dir <- cfg$data$peak_table_dir
  peak_file <- ifelse(
    file.exists(file.path(peak_dir, "Peak_table_for_cleaning.csv")),
    file.path(peak_dir, "Peak_table_for_cleaning.csv"),
    file.path(peak_dir, "peak_table_for_cleaning.csv"))
  if (!file.exists(peak_file)) {
    stop(sprintf("Peak table not found: %s", peak_file))
  }

  raw_data <- read.csv(peak_file, row.names = 1, header = TRUE, check.names = FALSE)
  cat(sprintf("=> Peak table: %d features x %d samples\n", nrow(raw_data), ncol(raw_data)))

  expression_data <- raw_data[, -1:-2, drop = FALSE]
  variable_info <- raw_data[, 1:2, drop = FALSE] %>%
    mutate(variable_id = rownames(raw_data))
  variable_info <- variable_info[, c(3, 1, 2)]
  rownames(variable_info) <- NULL

  # ── 2-pre: Internal standards (before any filter — see function docs) ────
  is_result <- extract_and_remove_internal_standards(
    expression_data, variable_info, sample_info_all, cfg)
  expression_data <- is_result$expression_data
  variable_info   <- is_result$variable_info
  is_qc_stats     <- is_result$is_qc_stats

  # ── Determine group membership per sample ────────────────────────────────
  col_map <- meta_list$col_map
  fn_col  <- col_map$filename

  match_group <- function(sid, meta, fn_col, grp_col) {
    idx <- which(meta[[fn_col]] == paste0(sid, ".mzML"))
    if (length(idx) > 0) return(meta[[grp_col]][idx[1]])
    idx <- which(gsub("\\.mzML$", "", meta[[fn_col]]) == sid)
    if (length(idx) > 0) return(meta[[grp_col]][idx[1]])
    NA_character_
  }

  sample_ids_in_data <- colnames(expression_data)
  classes <- sapply(sample_ids_in_data, function(sid) {
    match_group(sid, meta_list$all, fn_col, col_map$group_id)
  })

  # A peak-table column that cannot be traced back to the metadata would
  # otherwise be silently assigned a fabricated class (e.g. "QC_" for "QC_1"),
  # enter the SVR as a subject and spawn a bogus comparison. Fail loudly instead.
  unmatched <- sample_ids_in_data[is.na(classes)]
  if (length(unmatched) > 0) {
    stop(sprintf(
      "Peak-table samples not found in the metadata (%d): %s%s\n  Check that the peak table was built from this batch's mzML files.",
      length(unmatched),
      paste(head(unmatched, 10), collapse = ", "),
      if (length(unmatched) > 10) ", ..." else ""))
  }

  # Mirror build_sample_info_all()'s optional dose split (10 -> 10_High / 10_Low).
  # This class vector is what generate_comparisons() sees, so without it the
  # High and Low replicates would silently be pooled again here.
  suffix_col <- col_map$group_suffix
  if (!is.null(suffix_col) && suffix_col %in% colnames(meta_list$all)) {
    special <- c(cfg$sample_roles$control, cfg$sample_roles$qc_ms1, cfg$sample_roles$blank)
    meta_fn <- gsub("\\.mzML$", "", meta_list$all[[fn_col]], ignore.case = TRUE)
    suffix_vals <- as.character(meta_list$all[[suffix_col]])[match(sample_ids_in_data, meta_fn)]
    has_suffix <- !(classes %in% special) & !is.na(suffix_vals) & suffix_vals != ""
    classes[has_suffix] <- paste0(classes[has_suffix], "_", suffix_vals[has_suffix])
    cat(sprintf("   Group suffix from '%s': %d treatment classes\n",
                suffix_col, length(unique(classes[has_suffix]))))
  }

  # Tag MS2 samples
  ms_col <- col_map$ms_level
  if (!is.null(ms_col) && ms_col %in% colnames(meta_list$all)) {
    ms2_filenames <- meta_list$all[[fn_col]][meta_list$all[[ms_col]] == "MSMS"]
    ms2_ids <- gsub("\\.mzML$", "", ms2_filenames, ignore.case = TRUE)
    is_ms2 <- sample_ids_in_data %in% ms2_ids
    if (any(is_ms2)) {
      classes[is_ms2] <- paste0(classes[is_ms2], "-MSMS")
    }
  }

  sample_info_all <- data.frame(
    sample_id = sample_ids_in_data,
    class     = classes,
    stringsAsFactors = FALSE
  )
  cat("=> Sample types:\n")
  print(table(sample_info_all$class))

  f <- cfg$filtering

  # ── QC quality metrics (pre-filter) ──────────────────────────────────────
  # Computed on the raw peak table, before any filtering, so the numbers
  # describe the batch as acquired rather than the state of the surviving
  # subset. Written into the manifest to gate cross-batch integration.
  qc_quality <- compute_qc_quality(
    expression_data, sample_info_all$class, cfg$sample_roles$qc_ms1,
    rsd_threshold = f$rsd_threshold)
  if (!is.null(qc_quality)) {
    cat(sprintf("\n----- QC quality (pre-filter) -----\n"))
    cat(sprintf("   QC samples: %d, features: %d\n",
                qc_quality$n_qc_samples, qc_quality$n_features))
    cat(sprintf("   QC-RSD median: %.1f%%   QC-RSD pass rate (<=%.0f%%): %.1f%%\n",
                qc_quality$qc_rsd_median, f$rsd_threshold,
                100 * qc_quality$qc_rsd_pass_rate))
    cat(sprintf("   QC-QC correlation (median): %.3f\n", qc_quality$qc_qc_correlation))
    qa <- assess_quantification_quality(qc_quality)
    qc_quality$quantification_valid <- qa$valid
    qc_quality$exclude_reason <- if (qa$valid) NA_character_ else
      paste(qa$reasons, collapse = "; ")
    if (!qa$valid) {
      cat(sprintf("   ⚠️  Quantification NOT valid for cross-batch comparison: %s\n",
                  qc_quality$exclude_reason))
    }
  }

  # ── 2a: Blank background filter ──────────────────────────────────────────
  cat("\n----- Step 2a: Blank background filter -----\n")
  blank_name <- cfg$sample_roles$blank
  blank_samples <- sample_info_all$sample_id[
    grepl(paste0("^", blank_name), sample_info_all$class) &
      !grepl("-MSMS$", sample_info_all$class)]
  bio_samples <- sample_info_all$sample_id[
    !(sample_info_all$class %in% c(blank_name,
                                   paste0(blank_name, "-MSMS"),
                                   cfg$sample_roles$qc_ms1,
                                   paste0(cfg$sample_roles$qc_ms1, "-MSMS")))]

  cat(sprintf("   Blank samples: %d, Biological samples: %d\n",
              length(blank_samples), length(bio_samples)))

  if (length(blank_samples) >= 2) {
    # Per-feature blank level uses the MEDIAN across blanks (not the mean): it is
    # robust to unstable early/column-conditioning blank injections, so a few
    # noisy blanks no longer distort the background estimate.
    # Biological side stays a MEAN on purpose — a median across all samples would
    # be ~0 for a feature present in only a few of the many drug groups and would
    # wrongly drop it as background.
    blank_med  <- apply(as.matrix(expression_data[, blank_samples, drop = FALSE]),
                        1, median, na.rm = TRUE)
    bio_mean   <- rowMeans(as.matrix(expression_data[, bio_samples, drop = FALSE]),
                           na.rm = TRUE)
    blank_pass <- bio_mean >= f$blank_fold_change * blank_med
    blank_pass[is.na(blank_pass) | blank_med == 0] <- TRUE
    cat(sprintf("   Retained: %d/%d (removed %d background features)\n",
                sum(blank_pass), length(blank_pass), sum(!blank_pass)))
    expression_data <- expression_data[blank_pass, , drop = FALSE]
    variable_info   <- variable_info[blank_pass, , drop = FALSE]
  } else {
    cat("   Insufficient blank samples, skipping blank filter\n")
  }

  # ── 2b: Missing value filter ─────────────────────────────────────────────
  cat("\n----- Step 2b: Missing value filter -----\n")
  bio_class <- sample_info_all$class[match(bio_samples, sample_info_all$sample_id)]
  names(bio_class) <- bio_samples

  grp_levels <- unique(bio_class)
  grp_n <- vapply(grp_levels, function(grp) sum(bio_class == grp), integer(1))
  # per-feature missing COUNT per group (features x groups, integer-valued)
  miss_n <- vapply(grp_levels, function(grp) {
    grp_samps <- names(bio_class)[bio_class == grp]
    if (length(grp_samps) == 0) return(rep(0, nrow(expression_data)))
    apply(expression_data[, grp_samps, drop = FALSE], 1,
          function(x) sum(is.na(x)))
  }, numeric(nrow(expression_data)))
  miss_n <- as.matrix(miss_n)

  # (a) relative cap: for herb metabolomics, herb components are group-specific!
  # A feature is retained if it is detected in AT LEAST ONE biological group
  # with missing fraction <= max_missing_per_group AND detected count >= min_det.
  # (Alternatively, if detected in CT or QC samples).
  min_det <- f$min_detected_per_group %||% 2
  detected <- grp_n[col(miss_n)] - miss_n
  missing_frac <- miss_n / grp_n[col(miss_n)]

  # Pass if AT LEAST ONE biological group meets both the fraction and detection floor
  group_pass_matrix <- (missing_frac <= f$max_missing_per_group) & (detected >= min_det)
  det_pass <- apply(group_pass_matrix, 1, any, na.rm = TRUE)

  if (qualitative_only) {
    # In qualitative mode, keep if detected in >= 2 samples anywhere in cohort
    tot_det <- rowSums(!is.na(expression_data[, bio_samples, drop = FALSE]))
    cohort_pass <- tot_det >= min_det
    cat(sprintf("   [qualitative-only] per-group floor disabled; cohort floor (>=%d samples) applied\n",
                min_det))
    det_pass <- cohort_pass
  }

  mv_pass <- det_pass
  cat(sprintf("   Retained: %d/%d (removed %d features not meeting >= %d detections in any group)\n",
              sum(mv_pass, na.rm = TRUE), length(mv_pass),
              sum(!mv_pass, na.rm = TRUE), min_det))
  expression_data <- expression_data[mv_pass, , drop = FALSE]
  variable_info   <- variable_info[mv_pass, , drop = FALSE]

  # ── 2c: QC-RSD filter (MS1 only!) ────────────────────────────────────────
  cat("\n----- Step 2c: QC-RSD filter -----\n")
  qc_name <- cfg$sample_roles$qc_ms1
  qc_samples <- sample_info_all$sample_id[
    sample_info_all$class == qc_name]  # exact match - MS1 QC only, not QC-MSMS
  cat(sprintf("   QC samples (MS1 only): %d\n", length(qc_samples)))

  if (qualitative_only) {
    cat("   [qualitative-only] QC-RSD filter SKIPPED — it exists only to guard\n")
    cat("                      quantitative comparison, which this run does not do.\n")
    cat(sprintf("                      (QC-RSD median here is %.1f%%; see manifest qc_quality)\n",
                if (!is.null(qc_quality)) qc_quality$qc_rsd_median else NA_real_))
  } else if (length(qc_samples) >= 3) {
    qc_data <- as.matrix(expression_data[, qc_samples, drop = FALSE])
    calc_rsd <- function(x) {
      x <- x[!is.na(x) & x > 0]
      if (length(x) < 3) return(NA)
      sd(x) / mean(x) * 100
    }
    qc_rsd <- apply(qc_data, 1, calc_rsd)
    rsd_pass <- is.na(qc_rsd) | qc_rsd <= f$rsd_threshold
    cat(sprintf("   Retained: %d/%d (removed %d RSD > %.0f%% features)\n",
                sum(rsd_pass, na.rm = TRUE), length(rsd_pass),
                sum(!rsd_pass, na.rm = TRUE), f$rsd_threshold))
    expression_data <- expression_data[rsd_pass, , drop = FALSE]
    variable_info   <- variable_info[rsd_pass, , drop = FALSE]
  } else {
    cat("   Insufficient QC samples, skipping RSD filter\n")
  }

  # ── 2d: Imputation + low-intensity filter + SVR correction ──────────────
  cat("\n----- Step 2d: Imputation + normalization -----\n")
  norm_cfg <- cfg$normalization

  # Every path below ends in normalize_data(). If filtering has already removed
  # every feature, the median fallback dies deep inside apply()/median() with
  # "missing value where TRUE/FALSE needed" — median() of a zero-row column is
  # NA, and masscleaner's `if (median_x == 0)` then rejects the NA. That message
  # hides the real cause, an unusable peak table, so fail here instead.
  assert_features_remain <- function(object, stage) {
    if (nrow(object@expression_data) > 0) return(invisible(NULL))
    if (qualitative_only) {
      # Downstream steps here are annotation and a diagnostic PCA, neither of
      # which needs a normalized matrix. Surface the condition and continue.
      cat(sprintf("   ⚠️  No features remain after %s — continuing in qualitative-only mode\n",
                  stage))
      return(invisible(NULL))
    }
    stop(sprintf(paste0(
      "No features remain after %s for batch '%s' — cannot normalize.\n",
      "  Feature extraction most likely produced a degraded peak table; check ",
      "the\n  log above for peak-picking/alignment errors before re-running.\n",
      "  Delete %s to force a fresh feature extraction.\n",
      "  If this batch is only usable qualitatively, re-run with --qualitative-only."),
      stage, cfg$project$name %||% "unknown", cfg$data$peak_table_dir))
  }

  sample_info_svr <- sample_info_all
  # injection.order must reflect real acquisition order: the peak table's
  # column order is alphabetical and would cluster all QC at the end of the run.
  sample_info_svr$injection.order <- get_injection_order(
    sample_info_svr$sample_id, cfg$data$raw_mzml_dir)
  qc_pos <- sample_info_svr$injection.order[sample_info_svr$class == qc_name]
  cat(sprintf("   Injection order: %d samples, QC positions %s\n",
              nrow(sample_info_svr),
              if (length(qc_pos) > 0)
                sprintf("%d-%d (median %.0f)", min(qc_pos), max(qc_pos), median(qc_pos))
              else "n/a"))
  sample_info_svr$class[!sample_info_svr$class %in%
    c(blank_name, paste0(blank_name, "-MSMS"),
      qc_name, paste0(qc_name, "-MSMS"))] <- "Subject"
  sample_info_svr$class[sample_info_svr$class == qc_name] <- "QC"

  object <- create_mass_dataset(
    expression_data = expression_data,
    sample_info     = sample_info_svr,
    variable_info   = variable_info
  )

  mv_filled <- impute_mv(object = object, method = "minimum")
  cat("   Missing value imputation complete (minimum method)\n")

  isOK <- apply(mv_filled@expression_data, 1, function(row) {
    !any(row < f$min_intensity, na.rm = TRUE)
  })
  cat(sprintf("   Low-signal filter (< %g): retained %d/%d features\n",
              f$min_intensity, sum(isOK), length(isOK)))
  mv_filled@expression_data <- mv_filled@expression_data[isOK, , drop = FALSE]
  mv_filled@variable_info   <- mv_filled@variable_info[isOK, , drop = FALSE]
  assert_features_remain(mv_filled, "low-signal filtering")

  # QC drift correction
  qc_count <- sum(sample_info_all$class == qc_name)
  if (qc_count >= norm_cfg$qc_min_samples) {
    cat(sprintf("   Running QC drift correction (%d QC samples)...\n", qc_count))

    # Remove non-finite features
    data_ok <- apply(mv_filled@expression_data, 1, function(r) all(is.finite(r)))
    if (sum(data_ok) < nrow(mv_filled@expression_data)) {
      cat(sprintf("   Removing %d non-finite features\n", sum(!data_ok)))
      mv_filled@expression_data <- mv_filled@expression_data[data_ok, , drop = FALSE]
      mv_filled@variable_info   <- mv_filled@variable_info[data_ok, , drop = FALSE]
    }

    # For SVR drift correction, only features with non-zero variance across QC can be corrected.
    # Features with zero variance in QC (e.g. absent in QC, imputed with constant) are kept intact
    # rather than discarded from the biological dataset.
    qc_in_obj <- mv_filled@sample_info$sample_id[mv_filled@sample_info$class == "QC"]
    qc_data2  <- mv_filled@expression_data[, qc_in_obj, drop = FALSE]
    qc_var    <- apply(qc_data2, 1, function(x) var(as.numeric(x)))
    zero_var  <- qc_var == 0 | is.na(qc_var)

    if (nrow(mv_filled@expression_data) > 10) {
      if (any(zero_var)) {
        cat(sprintf("   Notice: %d features have zero variance in QC (absent or constant in QC).\n", sum(zero_var)))
        cat(sprintf("           Running SVR on %d QC-active features; uncorrected features preserved.\n", sum(!zero_var)))

        # Split into QC-active and QC-inactive
        obj_active <- mv_filled[!zero_var, ]
        obj_inactive <- mv_filled[zero_var, ]

        obj_active_norm <- tryCatch({
          normalize_data(object = obj_active, method = "svr",
                         optimization = TRUE, threads = norm_cfg$svr_threads)
        }, error = function(e1) {
          cat("   SVR on active features failed, falling back to LOESS...\n")
          tryCatch({
            normalize_data(object = obj_active, method = "loess",
                           optimization = TRUE, threads = norm_cfg$svr_threads)
          }, error = function(e2) {
            normalize_data(obj_active, method = norm_cfg$fallback_method)
          })
        })

        # Combine active normalized + inactive back together
        combined_expr <- rbind(obj_active_norm@expression_data, obj_inactive@expression_data)
        combined_var  <- rbind(obj_active_norm@variable_info, obj_inactive@variable_info)
        # reorder to match mv_filled
        ord <- match(mv_filled@variable_info$variable_id, combined_var$variable_id)
        object2 <- mv_filled
        object2@expression_data <- combined_expr[ord, , drop = FALSE]
        object2@variable_info   <- combined_var[ord, , drop = FALSE]
      } else {
        object2 <- tryCatch({
          normalize_data(object = mv_filled, method = "svr",
                         optimization = TRUE, threads = norm_cfg$svr_threads)
        }, error = function(e1) {
          cat("   SVR failed, falling back to LOESS...\n")
          tryCatch({
            normalize_data(object = mv_filled, method = "loess",
                           optimization = TRUE, threads = norm_cfg$svr_threads)
          }, error = function(e2) {
            cat("   LOESS also failed, falling back to median normalization\n")
            normalize_data(mv_filled, method = norm_cfg$fallback_method)
          })
        })
      }
    } else {
      cat("   Insufficient features, using median normalization\n")
      object2 <- normalize_data(mv_filled, method = norm_cfg$fallback_method)
    }
    cat("   Normalization complete\n")
  } else {
    cat(sprintf("   Insufficient QC samples (%d < %d), using median normalization\n",
                qc_count, norm_cfg$qc_min_samples))
    assert_features_remain(mv_filled, "low-signal filtering")
    object2 <- normalize_data(mv_filled, method = norm_cfg$fallback_method)
  }

  # Restore original class labels
  orig_class <- sample_info_all$class
  names(orig_class) <- sample_info_all$sample_id
  object2@sample_info$class <- orig_class[object2@sample_info$sample_id]

  # ── 2e: Extract core groups (non-QC, non-Blank, non-MS2) ────────────────
  cat("\n----- Step 2e: Extract core groups -----\n")
  core_classes <- setdiff(unique(sample_info_all$class),
                          c(blank_name, paste0(blank_name, "-MSMS"),
                            qc_name, paste0(qc_name, "-MSMS")))
  core_samples <- sample_info_all$sample_id[sample_info_all$class %in% core_classes]
  cat(sprintf("   Core samples: %d (%s)\n", length(core_samples),
              paste(core_classes, collapse = ", ")))

  sample_info <- object2@sample_info
  sample_info <- sample_info[sample_info$sample_id %in% core_samples, , drop = FALSE]
  rownames(sample_info) <- NULL

  expression_data <- object2@expression_data[, core_samples, drop = FALSE]
  variable_info   <- object2@variable_info

  cat(sprintf("   Final: %d features x %d samples\n", nrow(expression_data), ncol(expression_data)))
  print(table(sample_info$class))

  object2 <- create_mass_dataset(
    expression_data = expression_data,
    sample_info     = sample_info,
    variable_info   = variable_info
  )

  # ── 2f: Load MS2 spectra (from MSMS mzML files) ─────────────────────────
  cat("\n----- Step 2f: Load MS2 spectra -----\n")
  ms2_loaded <- FALSE

  meta_ms2 <- meta_list$ms2

  # Exclude MS2 BLANK runs from the annotation spectra source. A blank MS2 scan
  # that happens to match a feature's m/z + RT would attach a noise spectrum and
  # cause misannotation, so only QC/pooled/sample MS2 should feed mutate_ms2().
  if (nrow(meta_ms2) >= 1) {
    blank_name <- cfg$sample_roles$blank
    grp_col_ms2 <- meta_list$col_map$group_id
    is_blank_ms2 <- meta_ms2[[grp_col_ms2]] == blank_name
    if (any(is_blank_ms2)) {
      cat(sprintf("   Excluding %d MS2 blank run(s) from annotation spectra\n",
                  sum(is_blank_ms2)))
      meta_ms2 <- meta_ms2[!is_blank_ms2, , drop = FALSE]
    }
  }

  if (nrow(meta_ms2) >= 1) {
    ms2_files <- meta_ms2[[meta_list$col_map$filename]]
    ms2_paths <- file.path(cfg$data$raw_mzml_dir, ms2_files)
    ms2_exists <- ms2_paths[file.exists(ms2_paths)]

    if (length(ms2_exists) >= 1) {
      tmp_ms2_dir <- file.path(tempdir(), "ms2_spectra")
      dir.create(tmp_ms2_dir, showWarnings = FALSE, recursive = TRUE)
      file.copy(ms2_exists, tmp_ms2_dir, overwrite = TRUE)

      ms2_loaded <- tryCatch({
        object2 <- mutate_ms2(
          object = object2,
          column = cfg$annotation$column,
          polarity = cfg$feature_extraction$polarity,
          ms1.ms2.match.mz.tol = cfg$annotation$ms1_match_ppm,
          ms1.ms2.match.rt.tol = cfg$annotation$rt_tol_inhouse,
          path = tmp_ms2_dir
        )
        cat("   MS2 spectra loaded successfully\n")
        TRUE
      }, error = function(e) {
        cat(sprintf("   MS2 spectra loading failed: %s\n", conditionMessage(e)))
        FALSE
      })
      unlink(tmp_ms2_dir, recursive = TRUE)
    } else {
      cat("   No MS2 mzML files found, skipping\n")
    }
  } else {
    cat("   No MS2 samples in metadata, skipping\n")
  }

  list(
    object2         = object2,
    sample_info     = sample_info,
    variable_info   = variable_info,
    expression_data = expression_data,
    core_classes    = core_classes,
    ms2_loaded      = ms2_loaded,
    is_qc_stats     = is_qc_stats,
    qc_quality      = qc_quality
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# 5. Step 3: PCA
# ══════════════════════════════════════════════════════════════════════════════

run_pca <- function(object2, sample_info, cfg) {
  cat("\n===== Step 3: PCA =====\n")

  out_dir <- file.path(cfg$project$output_dir, "03_PCA")
  group <- sample_info$class[match(colnames(object2@expression_data),
                                    sample_info$sample_id)]
  pca_input <- t(object2@expression_data)
  pca_input <- pca_input[, apply(pca_input, 2, var, na.rm = TRUE) > 0, drop = FALSE]
  cat(sprintf("   PCA input: %d features\n", ncol(pca_input)))

  pca1 <- prcomp(pca_input, scale. = TRUE)
  df1 <- as.data.frame(pca1$x)
  summ1 <- summary(pca1)
  xlab1 <- paste0("PC1(", round(summ1$importance[2, 1] * 100, 2), "%)")
  ylab1 <- paste0("PC2(", round(summ1$importance[2, 2] * 100, 2), "%)")

  df1$group <- group

  # Build color map from config + auto-assign for unlisted groups
  viz_colors <- cfg$visualization$group_colors
  all_groups <- unique(group)
  auto_palette <- c("#D20A13", "#088247", "#FFD121", "#7E6148FF", "#5BC0EB",
                    "#F39B7FFF", "#BC3C29FF", "#0072B5FF", "#E18727FF",
                    "#20854EFF", "#7876B1FF", "#6F99ADFF", "#FFDC91FF",
                    "#EE4C97FF", "#00A087FF", "#8491B4FF", "#CC3333",
                    "#8C510A", "#01665E", "#5C88DA", "#A6CEE3", "#B2DF8A",
                    "#FB9A99", "#FF7F00", "#CAB2D6", "#33A02C", "#E31A1C",
                    "#1F78B4", "#B15928", "#FDBF6F", "#6A3D9A")
  plot_colors <- sapply(all_groups, function(g) {
    if (!is.null(viz_colors[[g]])) viz_colors[[g]]
    else auto_palette[(which(all_groups == g) - 1) %% length(auto_palette) + 1]
  })
  names(plot_colors) <- all_groups

  viz_labels <- cfg$visualization$group_labels

  p <- ggplot(data = df1, aes(x = PC1, y = PC2, color = group)) +
    stat_ellipse(aes(fill = group), type = "norm", geom = "polygon",
                 alpha = 0.15, color = NA, level = 0.8) +
    scale_color_manual(values = plot_colors,
                       labels = viz_labels[names(plot_colors)]) +
    scale_fill_manual(values = plot_colors) +
    theme_bw() +
    geom_point(size = 5) +
    geom_point(aes(x = PC1, y = PC2, color = group),
               shape = 21, color = "black", size = 5) +
    labs(x = xlab1, y = ylab1, color = "Group", title = "PCA Scores Plot") +
    guides(fill = "none") +
    labs(title = "") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 15),
      axis.title.x = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      axis.title.y = element_text(size = 16),
      axis.text.y = element_text(size = 14),
      panel.grid = element_blank(),
      aspect.ratio = 1,
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 13),
      plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), 'cm')
    )

  save_plot(file.path(out_dir, "PCA_scores"), plot = p, width = 10, height = 8)
  invisible(list(pca = pca1, scores = df1))
}

# ══════════════════════════════════════════════════════════════════════════════
# 6. Step 4: Differential Analysis (limma)
# ══════════════════════════════════════════════════════════════════════════════

run_differential <- function(object2, sample_info, comparisons, cfg) {
  cat("\n===== Step 4: Differential analysis (limma) =====\n")

  out_dir <- file.path(cfg$project$output_dir, "04_Differential")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  diff_cfg <- cfg$differential
  alpha <- diff_cfg$p_value_cutoff
  fc_cut <- diff_cfg$fc_threshold
  logfc_cutoff <- log2(fc_cut)

  # Significance metric: "adj_p_value" (FDR, topTable's adj.P.Val) or "p_value"
  # (raw P.Value, default). Controls both the sig-feature filter and the volcano
  # plot's dashed threshold / point coloring so they stay consistent.
  sig_metric <- if (!is.null(diff_cfg$significance_metric)) diff_cfg$significance_metric else "p_value"
  p_col <- if (identical(sig_metric, "adj_p_value")) "adj.P.Val" else "P.Value"
  p_lab <- if (identical(sig_metric, "adj_p_value")) "adj. p-value (FDR)" else "p-value"
  cat(sprintf("   Significance: |log2FC| >= %.3f & %s < %.3g\n", logfc_cutoff, p_lab, alpha))

  Nomalize_data <- object2@expression_data
  log_data <- log2(Nomalize_data)
  core_group_vec <- sample_info$class[match(colnames(Nomalize_data),
                                             sample_info$sample_id)]
  all_samples <- colnames(Nomalize_data)

  diff_results <- list()
  viz_labels <- cfg$visualization$group_labels

  for (cmp in comparisons) {
    cmp_name  <- cmp$name
    ctrl_grp  <- cmp$ctrl
    treat_grp <- cmp$treat
    label     <- cmp$label

    cat(sprintf("   Analyzing: %s\n", label))

    ctrl_samples  <- all_samples[core_group_vec == ctrl_grp]
    treat_samples <- all_samples[core_group_vec == treat_grp]

    if (length(ctrl_samples) < 2 || length(treat_samples) < 2) {
      cat(sprintf("   ⚠️  Insufficient samples (ctrl=%d, treat=%d), skipping\n",
                  length(ctrl_samples), length(treat_samples)))
      next
    }

    ctrl_data  <- log_data[, ctrl_samples, drop = FALSE]
    treat_data <- log_data[, treat_samples, drop = FALSE]
    ana_data   <- cbind(ctrl_data, treat_data)
    ctl_num    <- ncol(ctrl_data)
    treat_num  <- ncol(treat_data)

    Type <- c(rep("CTL", ctl_num), rep("TREAT", treat_num))
    design <- model.matrix(~0 + factor(Type))
    colnames(design) <- c("CTL", "TREAT")

    fit <- lmFit(ana_data, design = design)
    cont.matrix <- makeContrasts(TREAT - CTL, levels = design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2)
    Diff <- topTable(fit2, adjust = diff_cfg$p_adjust_method, number = nrow(ana_data))

    condition <- abs(Diff$logFC) >= logfc_cutoff & Diff[[p_col]] < alpha
    Diff_sig <- Diff[condition, , drop = FALSE]

    diff_results[[cmp_name]] <- list(
      all           = Diff,
      sig           = Diff_sig,
      ctrl_samples  = ctrl_samples,
      treat_samples = treat_samples,
      n_up          = sum(Diff_sig$logFC > 0, na.rm = TRUE),
      n_down        = sum(Diff_sig$logFC < 0, na.rm = TRUE),
      comparison    = cmp
    )

    cat(sprintf("     Significant: %d (up: %d, down: %d)\n",
                nrow(Diff_sig),
                diff_results[[cmp_name]]$n_up,
                diff_results[[cmp_name]]$n_down))

    write.csv(Diff_sig,
              file = file.path(out_dir, paste0(cmp_name, "_differential_peaks.csv")),
              row.names = TRUE)

    # Volcano plot — y-axis, coloring and threshold line all use the chosen
    # significance metric (p_col) so the plot matches the sig-feature table.
    Diff$.plot_p <- Diff[[p_col]]
    Significant <- ifelse(
      (Diff$.plot_p < alpha & abs(Diff$logFC) > logfc_cutoff),
      ifelse(Diff$logFC > logfc_cutoff, "Up", "Down"),
      "Not"
    )

    p <- ggplot(Diff, aes(x = logFC, y = -log10(.plot_p))) +
      geom_point(alpha = 0.4, size = 3.5, aes(color = Significant)) +
      ylab(sprintf("-log10(%s)", p_lab)) +
      scale_color_manual(values = c("blue4", "grey", "red3")) +
      geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff),
                 lty = 4, col = "black", lwd = 0.8) +
      geom_hline(yintercept = -log10(alpha),
                 lty = 4, col = "black", lwd = 0.8) +
      labs(title = label) +
      theme_bw() +
      theme(aspect.ratio = 1, panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 14))

    save_plot(file.path(out_dir, paste0(cmp_name, "_volcano")),
              plot = p, width = 8, height = 7, dpi = cfg$visualization$dpi)
  }

  cat(sprintf("   Completed %d comparisons\n", length(diff_results)))
  invisible(diff_results)
}

# ══════════════════════════════════════════════════════════════════════════════
# 7. Step 5: Overlap Analysis (generic, replaces synergy step)
# ══════════════════════════════════════════════════════════════════════════════

run_overlap_analysis <- function(diff_results, cfg) {
  cat("\n===== Step 5: Overlap analysis =====\n")

  if (!isTRUE(cfg$overlap_analysis$enabled)) {
    cat("   Overlap analysis disabled in config, skipping\n")
    return(invisible(NULL))
  }

  out_dir <- file.path(cfg$project$output_dir, "05_Overlap")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  oa_cfg <- cfg$overlap_analysis

  # Build set list: for each comparison, collect significant variable_ids
  sig_sets <- lapply(diff_results, function(dr) {
    if (!is.null(dr$sig) && nrow(dr$sig) > 0) rownames(dr$sig) else character(0)
  })

  # Use display labels for set names
  set_labels <- sapply(diff_results, function(dr) {
    if (!is.null(dr$comparison$label)) dr$comparison$label
    else if (!is.null(dr$comparison$name)) dr$comparison$name
    else "unknown"
  })

  # Filter to comparisons with at least min_comparisons_per_metabolite hits.
  # n_sig must be filtered alongside sig_sets/set_labels: the UpSet branch below
  # indexes the filtered list with ranks taken from n_sig, so leaving n_sig at
  # its original length selects the wrong comparisons (or NULLs) whenever any
  # comparison is dropped here.
  n_sig <- sapply(sig_sets, length)
  keep <- n_sig >= oa_cfg$min_comparisons_per_metabolite
  sig_sets <- sig_sets[keep]
  set_labels <- set_labels[keep]
  n_sig <- n_sig[keep]

  if (length(sig_sets) < 2) {
    cat("   < 2 comparisons with significant hits, skipping overlap analysis\n")
    return(invisible(NULL))
  }

  cat(sprintf("   %d comparisons with >= %d significant metabolites\n",
              length(sig_sets), oa_cfg$min_comparisons_per_metabolite))

  n_sets <- length(sig_sets)
  venn_max <- oa_cfg$venn_max_sets

  if (n_sets <= venn_max) {
    # Use VennDiagram (up to 5 sets)
    cat("   Drawing Venn diagram...\n")
    venn_colors <- c("#D20A13", "#088247", "#FFD121", "#0072B5FF", "#7E6148FF")
    for (fmt in c("pdf", "png")) {
      fpath <- file.path(out_dir, paste0("Venn_overlap.", fmt))
      if (fmt == "pdf") {
        pdf(fpath, width = 8, height = 7)
      } else {
        png(fpath, width = 8, height = 7, units = "in", res = cfg$visualization$dpi)
      }
      venn_plot <- venn.diagram(
        x = sig_sets,
        category.names = set_labels,
        filename = NULL, output = TRUE,
        col = venn_colors[1:n_sets],
        fill = alpha(venn_colors[1:n_sets], 0.3),
        cat.cex = 1.0, cex = 1.2,
        cat.default.pos = "outer"
      )
      grid.draw(venn_plot)
      dev.off()
    }
    cat("   Venn diagram saved\n")
  } else {
    # Use UpSetR for larger sets
    cat(sprintf("   %d comparisons > Venn max (%d), using UpSet plot...\n",
                n_sets, venn_max))

    # Select top comparisons by n_significant for UpSet
    upset_max <- oa_cfg$upset_max_sets
    if (n_sets > upset_max) {
      top_idx <- order(n_sig, decreasing = TRUE)[1:upset_max]
      sig_sets <- sig_sets[top_idx]
      set_labels <- set_labels[top_idx]
    }

    if (requireNamespace("UpSetR", quietly = TRUE)) {
      # Build binary membership matrix
      all_metabs <- unique(unlist(sig_sets))
      upset_mat <- as.data.frame(sapply(sig_sets, function(s) {
        as.integer(all_metabs %in% s)
      }))
      rownames(upset_mat) <- all_metabs
      colnames(upset_mat) <- set_labels

      for (fmt in c("pdf", "png")) {
        fpath <- file.path(out_dir, paste0("upset_overlap.", fmt))
        if (fmt == "pdf") {
          pdf(fpath, width = 14, height = 8)
        } else {
          png(fpath, width = 14, height = 8, units = "in", res = cfg$visualization$dpi)
        }
        UpSetR::upset(upset_mat,
                      nsets = min(ncol(upset_mat), upset_max),
                      order.by = "freq",
                      main.bar.color = "#58CDD9",
                      sets.bar.color = "#D20A13")
        dev.off()
      }
      cat("   UpSet plot saved\n")
    } else {
      cat("   UpSetR package not available, saving membership table only\n")
    }

    # Save metabolite membership table
    membership <- data.frame(
      variable_id = all_metabs,
      upset_mat,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    write.csv(membership, file.path(out_dir, "metabolite_overlap_membership.csv"),
              row.names = FALSE)
    cat("   Membership table saved\n")
  }

  invisible(NULL)
}

# ══════════════════════════════════════════════════════════════════════════════
# 8. Step 6: Enhanced Annotation (multi-database MS1/MS2 + MSI confidence)
# ══════════════════════════════════════════════════════════════════════════════

#' Run comprehensive metabolite annotation against MS1 and MS2 databases
#'
#' Loads all available databases (MS1 accurate-mass and MS2 spectral libraries),
#' annotates features, merges by priority, and assigns MSI confidence levels.
#' Extended with GNPS sub-libraries (IH, NIH) and MoNA for broader coverage.
#'
#' @param object2 mass_dataset object with MS2 spectra loaded
#' @param cfg Validated config list
#' @return data.frame of merged annotations with confidence levels, or NULL
run_annotation <- function(object2, cfg) {
  cat("\n===== Step 6: Metabolite annotation (enhanced) =====\n")

  out_dir <- file.path(cfg$project$output_dir, "06_Annotation")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ann_cfg <- cfg$annotation

  # ── Build comprehensive database source list ─────────────────────────────
  db_sources <- list()
  # The in-house DB is opt-in: its mz column does not match its own formulas,
  # so when enabled it is loaded but must not be given merge priority.
  if (isTRUE(ann_cfg$inhouse_enabled)) {
    db_sources <- c(db_sources, list(
      list(name = "inhouse",  label = "In-house",  path = ann_cfg$inhouse_db_path, obj_name = "inhouse_Metabolite.database", type = "MS1")
    ))
  } else {
    cat("   In-house DB disabled (annotation$inhouse_enabled = FALSE)\n")
  }
  db_sources <- c(db_sources, list(
    # MS1 accurate-mass databases
    list(name = "HMDB",       label = "HMDB",        path = file.path(ann_cfg$ms1_db_dir, "HMDB/hmdb_ms1_database.rda"),    obj_name = "hmdb_ms1", type = "MS1"),
    list(name = "KEGG",       label = "KEGG",        path = file.path(ann_cfg$ms1_db_dir, "KEGG/kegg_ms1_database.rda"),    obj_name = "kegg_ms1", type = "MS1"),
    list(name = "ChEBI",      label = "ChEBI",       path = file.path(ann_cfg$ms1_db_dir, "ChEBI/chebi_ms1_database.rda"),  obj_name = "chebi_ms1", type = "MS1"),
    list(name = "FooDB",      label = "FooDB",       path = file.path(ann_cfg$ms1_db_dir, "FooDB/foodb_ms1_database.rda"),  obj_name = "foodb_ms1", type = "MS1"),
    list(name = "PubChem",    label = "PubChem",     path = file.path(ann_cfg$ms1_db_dir, "PubChem/pubchem_ms1_database.rda"), obj_name = "pubchem_ms1", type = "MS1"),
    # MS2 spectral libraries
    list(name = "MS2_HMDB",     label = "HMDB MS2",     path = file.path(ann_cfg$ms2_db_dir, "hmdb_ms2_merged.rda"),       obj_name = "hmdb_ms2", type = "MS2"),
    list(name = "MS2_GNPS",     label = "GNPS MS2",     path = file.path(ann_cfg$ms2_db_dir, "gnps_ms2_merged.rda"),       obj_name = "gnps_ms2", type = "MS2"),
    list(name = "MS2_MassBank", label = "MassBank MS2", path = file.path(ann_cfg$ms2_db_dir, "massbank_ms2_merged.rda"),   obj_name = "massbank_ms2", type = "MS2"),
    list(name = "MS2_MoNA",     label = "MoNA MS2",     path = file.path(ann_cfg$ms2_db_dir, "mona_ms2_merged.rda"),       obj_name = "mona_ms2", type = "MS2"),
    # GNPS sub-libraries (specialised collections)
    list(name = "MS2_GNPS_IH",  label = "GNPS IH MS2",  path = file.path(ann_cfg$ms2_db_dir, "gnps_iobanhc_ms2_merged.rda"),  obj_name = "gnps_iobanhc_ms2", type = "MS2"),
    list(name = "MS2_GNPS_NIH", label = "GNPS NIH MS2", path = file.path(ann_cfg$ms2_db_dir, "gnps_nihclinicalcollection1_ms2_merged.rda"), obj_name = "gnps_nihclinicalcollection1_ms2", type = "MS2")
  ))

  databases <- list()
  for (src in db_sources) {
    if (file.exists(src$path)) {
      tryCatch({
        load(src$path)
        db_obj <- get(src$obj_name)
        if (inherits(db_obj, "databaseClass")) {
          if ("Lab.ID" %in% colnames(db_obj@spectra.info)) {
            db_obj@spectra.info$Lab.ID <- as.character(db_obj@spectra.info$Lab.ID)
          }
          # Ensure standard annotation columns exist
          std_cols <- c("Compound.name", "CAS.ID", "HMDB.ID", "KEGG.ID", "Formula", "RT")
          for (sc in setdiff(std_cols, colnames(db_obj@spectra.info))) {
            db_obj@spectra.info[[sc]] <- NA
          }
          databases[[src$name]] <- db_obj
          cat(sprintf("   Loaded: %s | %d compounds | type: %s\n", src$label,
                      nrow(databases[[src$name]]@spectra.info), src$type))
        }
      }, error = function(e) {
        cat(sprintf("   Load failed: %s - %s\n", src$label, conditionMessage(e)))
      })
    } else {
      # Only warn for known databases; new ones may legitimately not exist
      cat(sprintf("   DB not found: %s (%s), skipping\n", src$label, src$path))
    }
  }
  cat(sprintf("   Total databases loaded: %d\n", length(databases)))

  if (length(databases) == 0) {
    cat("   No databases available, skipping annotation\n")
    return(invisible(NULL))
  }

  # ── Annotate per database ────────────────────────────────────────────────
  all_annotations <- list()
  polarity <- cfg$feature_extraction$polarity

  for (db_name in names(databases)) {
    cat(sprintf("   Annotating with %s...\n", db_name))

    tryCatch({
      rt_tol <- ifelse(db_name == "inhouse", ann_cfg$rt_tol_inhouse, ann_cfg$rt_tol_external)
      has_ms2 <- length(databases[[db_name]]@spectra.data) > 0 &&
        any(sapply(databases[[db_name]]@spectra.data, function(x) length(x) > 0))

      if (has_ms2) {
        cat("     MS2 spectra detected, using MS1+MS2 joint matching\n")
        # Use configurable scoring weights (with defaults)
        ms1_wt <- ann_cfg$ms1_match_weight %||% 0.25
        rt_wt  <- ann_cfg$rt_match_weight  %||% 0.25
        ms2_wt <- ann_cfg$ms2_match_weight %||% 0.5
        object_annotated <- annotate_metabolites_mass_dataset(
          object = object2,
          ms1.match.ppm = ann_cfg$ms1_match_ppm,
          ms2.match.ppm = ann_cfg$ms2_match_ppm,
          mz.ppm.thr = 400,
          ms2.match.tol = ann_cfg$ms2_match_tol,
          rt.match.tol = rt_tol,
          polarity = polarity,
          ms1.match.weight = ms1_wt,
          rt.match.weight = rt_wt,
          ms2.match.weight = ms2_wt,
          total.score.tol = ann_cfg$total_score_tol %||% 0.5,
          candidate.num = ann_cfg$candidate_num %||% 3,
          database = databases[[db_name]]
        )
      } else {
        object_annotated <- annotate_metabolites_mass_dataset(
          object = object2,
          ms1.match.ppm = ann_cfg$ms1_match_ppm,
          rt.match.tol = rt_tol,
          polarity = polarity,
          database = databases[[db_name]]
        )
      }

      ann_table <- extract_annotation_table(object_annotated)
      if (!is.null(ann_table) && nrow(ann_table) > 0) {
        all_annotations[[db_name]] <- ann_table
        cat(sprintf("     Matched %d metabolites\n", nrow(ann_table)))
      } else {
        cat("     No matches\n")
      }
    }, error = function(e) {
      cat(sprintf("     Annotation error: %s\n", conditionMessage(e)))
    })
  }

  # ── Merge annotations by priority ────────────────────────────────────────
  if (length(all_annotations) == 0) {
    cat("   No annotations from any database\n")
    return(invisible(NULL))
  }

  db_priority_order <- ann_cfg$db_priority
  db_order <- db_priority_order[db_priority_order %in% names(all_annotations)]
  # Append databases not in priority list (don't silently drop evidence)
  db_remaining <- setdiff(names(all_annotations), db_order)
  if (length(db_remaining) > 0) {
    cat(sprintf("   Appending %d non-priority DB(s): %s\n",
                length(db_remaining), paste(db_remaining, collapse = ", ")))
    db_order <- c(db_order, db_remaining)
  }
  cat(sprintf("   Merging %d databases, priority: %s\n",
              length(db_order), paste(db_order, collapse = " > ")))

  # MS2 evidence requires an actual spectral match (SS). Some MS2 libraries hold
  # no spectra for this polarity (e.g. gnps_nihclinicalcollection1 has none for
  # negative mode); metid then silently falls back to MS1-only identification,
  # and labelling that as an MS2 hit would promote every such row to MSI Level 2.
  ms2_match_type <- function(df, db_name) {
    if (!grepl("^MS2_", db_name) || !"SS" %in% colnames(df)) {
      return(rep("MS1", nrow(df)))
    }
    ifelse(!is.na(df$SS) & df$SS > 0, "MS2", "MS1")
  }

  # First priority database
  merged_ann <- all_annotations[[db_order[1]]]
  merged_ann$database_source <- db_order[1]
  merged_ann$match_type <- ms2_match_type(merged_ann, db_order[1])

  if (length(db_order) > 1) {
    for (i in 2:length(db_order)) {
      next_db <- all_annotations[[db_order[i]]]
      if (is.null(next_db) || nrow(next_db) == 0) next
      next_db$database_source <- db_order[i]
      next_db$match_type <- ms2_match_type(next_db, db_order[i])
      new_ann <- next_db[!next_db$variable_id %in% merged_ann$variable_id, , drop = FALSE]
      if (nrow(new_ann) > 0) {
        merged_ann <- bind_rows(merged_ann, new_ann)
      }
    }
  }

  # ── Assign MSI confidence levels ─────────────────────────────────────────
  merged_ann <- merged_ann %>%
    rowwise() %>%
    mutate(
      confidence_level = assign_confidence_level(
        has_ms2_match = (.data$match_type == "MS2"),
        has_ms1_match = (.data$match_type == "MS1"),
        has_rt_match  = FALSE,
        is_inhouse    = .data$database_source == "inhouse",
        is_sirius     = FALSE,
        is_metfrag    = FALSE
      )
    ) %>%
    ungroup()

  n_level1 <- sum(merged_ann$confidence_level == "Level 1 — Confirmed", na.rm = TRUE)
  n_level2 <- sum(merged_ann$confidence_level == "Level 2 — Putatively identified", na.rm = TRUE)
  n_level3 <- sum(merged_ann$confidence_level == "Level 3 — Putatively annotated", na.rm = TRUE)
  cat(sprintf("   MSI Confidence breakdown:\n"))
  cat(sprintf("     Level 1 — Confirmed:             %d\n", n_level1))
  cat(sprintf("     Level 2 — Putatively identified: %d\n", n_level2))
  cat(sprintf("     Level 3 — Putatively annotated:  %d\n", n_level3))

  # ── Add m/z and RT info ──────────────────────────────────────────────────
  if (nrow(merged_ann) > 0) {
    vi <- object2@variable_info %>%
      dplyr::select(variable_id, mz, rt)
    merged_ann <- merged_ann %>%
      left_join(vi, by = "variable_id") %>%
      dplyr::select(variable_id, mz, rt, database_source, match_type,
                    confidence_level, Compound.name, everything())
  }

  # Merge mz/rt from peak table (for mzmin/mzmax/rtmin/rtmax)
  peak_file2 <- file.path(cfg$data$peak_table_dir, "Peak_table.csv")
  if (!file.exists(peak_file2)) peak_file2 <- file.path(cfg$data$peak_table_dir, "peak_table.csv")
  if (file.exists(peak_file2)) {
    mzrt <- read.csv(peak_file2, stringsAsFactors = FALSE)[, 1:7]
    colnames(mzrt) <- c("variable_id", "mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax")
    # Avoid duplicate mz/rt columns by dropping them from mzrt before join
    mzrt_clean <- mzrt %>% dplyr::select(variable_id, mzmin, mzmax, rtmin, rtmax)
    app3 <- left_join(merged_ann, mzrt_clean, by = "variable_id")
    # Reorder: put core columns first
    app3 <- app3 %>%
      dplyr::select(variable_id, mz, rt, mzmin, mzmax, rtmin, rtmax,
                    database_source, match_type, confidence_level, everything())
  } else {
    app3 <- merged_ann
  }

  n_unique <- length(unique(app3$Compound.name[!is.na(app3$Compound.name) &
                                                app3$Compound.name != ""]))
  cat(sprintf("   Merged: %d rows, %d unique compound names\n", nrow(app3), n_unique))

  write.csv(app3, file.path(out_dir, "all_annotated_metabolites.csv"), row.names = FALSE)
  cat(sprintf("   Saved: %s\n", file.path(out_dir, "all_annotated_metabolites.csv")))

  # Per-database stats
  db_stats <- table(app3$database_source)
  for (nm in names(db_stats)) {
    cat(sprintf("     %s: %d\n", nm, db_stats[nm]))
  }

  # Save per-confidence-level CSVs
  for (lvl in unique(app3$confidence_level)) {
    lvl_short <- gsub(" .*$", "", gsub("Level ", "L", lvl))
    lvl_data <- app3[app3$confidence_level == lvl, , drop = FALSE]
    if (nrow(lvl_data) > 0) {
      write.csv(lvl_data,
                file.path(out_dir, paste0("annotations_", lvl_short, ".csv")),
                row.names = FALSE)
    }
  }
  # Save unmatched features
  all_var_ids <- rownames(object2@expression_data)
  matched_ids <- unique(app3$variable_id)
  unmatched <- setdiff(all_var_ids, matched_ids)
  if (length(unmatched) > 0) {
    unmatched_df <- data.frame(variable_id = unmatched, stringsAsFactors = FALSE)
    write.csv(unmatched_df, file.path(out_dir, "unmatched_features.csv"), row.names = FALSE)
    cat(sprintf("   Unmatched features: %d -> unmatched_features.csv\n", length(unmatched)))
  }

  invisible(app3)
}

# ══════════════════════════════════════════════════════════════════════════════
# 8a. Step 6a: CAMERA Adduct & Isotope Grouping (optional)
# ══════════════════════════════════════════════════════════════════════════════

#' Run CAMERA for adduct and isotope annotation
#'
#' CAMERA groups features based on RT proximity and peak-shape correlation,
#' then annotates isotope peaks and adducts. This helps determine the true
#' parent ion mass for each compound. Requires xcms and CAMERA Bioconductor
#' packages.
#'
#' @param object mass_dataset object
#' @param cfg Validated config list
#' @return Invisible list with object (updated) and camera_an (xsAnnotate object)
run_camera_annotation <- function(object, cfg) {
  cat("\n===== Step 6a: CAMERA feature grouping (optional) =====\n")

  if (!isTRUE(cfg$camera$enabled)) {
    cat("   CAMERA disabled in config, skipping\n")
    return(invisible(NULL))
  }

  # Check if CAMERA and xcms are available
  if (!requireNamespace("CAMERA", quietly = TRUE)) {
    cat("   CAMERA package not installed. Skipping.\n")
    cat("   Install with: BiocManager::install('CAMERA')\n")
    return(invisible(NULL))
  }
  if (!requireNamespace("xcms", quietly = TRUE)) {
    cat("   xcms package not installed. Skipping CAMERA.\n")
    cat("   Install with: BiocManager::install('xcms')\n")
    return(invisible(NULL))
  }

  cat("   Loading CAMERA...\n")
  suppressPackageStartupMessages(library(CAMERA))
  suppressPackageStartupMessages(library(xcms))

  # Try to load xcmsSet from massprocesser Result directory
  peak_dir <- cfg$data$peak_table_dir
  xcms_set_path <- file.path(peak_dir, "xcms_set.rda")
  xcms_object_path <- file.path(peak_dir, "xcms_object.rda")
  xset <- NULL

  for (fpath in c(xcms_set_path, xcms_object_path)) {
    if (file.exists(fpath)) {
      tryCatch({
        load(fpath)
        cat(sprintf("   Loaded xcmsSet from: %s\n", fpath))
        if (exists("xset") && class(xset) == "xcmsSet") {
          # already have xset
        } else if (exists("xcms_object") && class(xcms_object) == "xcmsSet") {
          xset <- xcms_object
        }
        break
      }, error = function(e) {
        cat(sprintf("   Failed to load xcmsSet: %s\n", conditionMessage(e)))
      })
    }
  }

  if (is.null(xset)) {
    cat("   xcmsSet not found in Result directory.\n")
    cat("   Attempting to create one from raw data...\n")

    raw_dir <- cfg$data$raw_mzml_dir
    ms1_files <- list.files(raw_dir, pattern = "\\.mzML$", full.names = TRUE)
    # Exclude MS/MS files
    ms1_files <- ms1_files[!grepl("MSMS", basename(ms1_files), ignore.case = TRUE)]
    cat(sprintf("   Using %d MS1-only .mzML files for xcmsSet\n", length(ms1_files)))
    fe <- cfg$feature_extraction
    tryCatch({
      xset <- xcms::xcmsSet(
        files = ms1_files,
        method = "centWave",
        ppm = fe$ppm,
        peakwidth = unlist(fe$peakwidth),
        snthresh = fe$snthresh,
        noise = fe$noise,
        polarity = fe$polarity,
        nSlaves = fe$threads
      )
      xset <- xcms::group(xset)
      xset <- xcms::retcor(xset)
      xset <- xcms::group(xset)
      cat(sprintf("   xcmsSet created: %d features\n", nrow(xset@peaks)))
    }, error = function(e) {
      cat(sprintf("   xcmsSet creation failed: %s\n", conditionMessage(e)))
      xset <<- NULL
    })
  }

  if (is.null(xset)) {
    cat("   No xcmsSet available, skipping CAMERA\n")
    return(invisible(NULL))
  }

  # Run CAMERA xsAnnotate
  cat("   Running CAMERA xsAnnotate...\n")
  camera_cfg <- cfg$camera
  an <- tryCatch({
    xsAnnotate(xset, polarity = camera_cfg$polarity)
  }, error = function(e) {
    cat(sprintf("   xsAnnotate failed: %s\n", conditionMessage(e)))
    return(NULL)
  })

  if (is.null(an)) {
    cat("   CAMERA xsAnnotate failed, skipping\n")
    return(invisible(NULL))
  }

  # Group features by FWHM
  an <- groupFWHM(an, perfwhm = camera_cfg$perfwhm)

  # Correlate EIC
  an <- findIsotopes(an, ppm = camera_cfg$ppm)
  an <- groupCorr(an, cor_eic_th = camera_cfg$cor_eic_th,
                  graphMethod = camera_cfg$graphMethod,
                  pval = camera_cfg$pval,
                  calcCiS = camera_cfg$calcCiS,
                  calcCaS = camera_cfg$calcCaS)
  an <- findAdducts(an, polarity = camera_cfg$polarity)

  # Extract results and map back to mass_dataset
  peaklist <- getPeaklist(an)

  if (!is.null(peaklist) && nrow(peaklist) > 0) {
    camera_info <- peaklist %>%
      mutate(
        camera_isotopes = as.character(.data$isotopes),
        camera_adduct   = as.character(.data$adduct),
        camera_pcgroup  = .data$pcgroup
      ) %>%
      dplyr::select(mz, rt, camera_isotopes, camera_adduct, camera_pcgroup)

    # Merge with variable_info by m/z and RT proximity
    camera_info_renamed <- camera_info %>%
      dplyr::rename(camera_mz = mz, camera_rt = rt)

    merged <- object@variable_info %>%
      dplyr::cross_join(camera_info_renamed) %>%
      mutate(
        mz_diff = abs(mz - camera_mz),
        rt_diff = abs(rt - camera_rt)
      ) %>%
      filter(mz_diff < 0.01, rt_diff < 30) %>%
      group_by(variable_id) %>%
      slice_min(order_by = mz_diff + rt_diff / 100, n = 1) %>%
      ungroup() %>%
      dplyr::select(variable_id, camera_isotopes, camera_adduct, camera_pcgroup)

    object@variable_info <- object@variable_info %>%
      left_join(merged, by = "variable_id")

    n_annotated <- sum(!is.na(object@variable_info$camera_adduct))
    n_isotope <- sum(!is.na(object@variable_info$camera_isotopes) &
                       object@variable_info$camera_isotopes != "")
    cat(sprintf("   CAMERA: %d features with adduct annotation, %d with isotopes\n",
                n_annotated, n_isotope))
  } else {
    cat("   CAMERA returned no peaklist\n")
  }

  cat("   CAMERA annotation complete\n")
  invisible(list(object = object, camera_an = an))
}

# ══════════════════════════════════════════════════════════════════════════════
# 8b. Step 6b: SIRIUS Molecular Formula & Structure Prediction (optional)
# ══════════════════════════════════════════════════════════════════════════════

#' List the CSI:FingerID structure databases actually installed in the SIRIUS
#' workspace.
#'
#' SIRIUS 6.3 bundles only a small default set; the natural-product collections
#' (PLANTCYC, KNAPSACK, GNPS, ...) must be downloaded first. Requesting an
#' uninstalled database does not error — it simply matches nothing — so callers
#' need this list to avoid silently degrading a structure search.
#'
#' @param sirius_path Path to the SIRIUS executable
#' @return Character vector of installed database names; character(0) if the
#'   list could not be determined (callers should then pass the request through)
sirius_available_databases <- function(sirius_path) {
  out <- tryCatch(
    suppressWarnings(system2(sirius_path, args = c("custom-db", "show"),
                             stdout = TRUE, stderr = TRUE)),
    error = function(e) character(0))
  if (length(out) == 0) return(character(0))
  # No custom databases installed — SIRIUS prints this and then its citation
  # block, which must not be mistaken for database names.
  if (any(grepl("No Custom database", out))) return(character(0))
  # Otherwise `custom-db show` prints a table whose rows start with the
  # database name. Accept only lines whose first whitespace-delimited token is
  # an all-caps identifier (e.g. PLANTCYC, KNAPSACK, GNPS); this rejects the
  # log preamble and the trailing citation prose.
  keep <- grepl("^\\s*[A-Z][A-Z0-9_]*\\s+\\S", out)
  if (!any(keep)) return(character(0))
  db <- trimws(sub("\\s+.*$", "", out[keep]))
  unique(db[nzchar(db)])
}

#' Export MS2 spectra to MGF format and optionally run SIRIUS
#'
#' SIRIUS is a Java CLI tool for molecular formula prediction and structure
#' identification from MS/MS data. This step:
#'  1. Exports MS2 spectra to MGF for SIRIUS
#'  2. If SIRIUS is available, launches it to predict formula, structure, and fingerprint
#'  3. Imports the results back into the annotation table
#'
#' Requires SIRIUS CLI (v6.x) installed separately.
#' Download from: https://bio.informatik.uni-jena.de/software/sirius/
#'
#' @param object mass_dataset object with MS2 spectra loaded
#' @param app3 Current annotation data.frame
#' @param cfg Validated config list
#' @return Updated annotation data.frame with SIRIUS columns (if available)
run_sirius_analysis <- function(object, app3, cfg) {
  cat("\n===== Step 6b: SIRIUS molecular formula / structure prediction (optional) =====\n")

  if (!isTRUE(cfg$sirius$enabled)) {
    cat("   SIRIUS disabled in config, skipping\n")
    return(app3)
  }

  sirius_path <- cfg$sirius$path
  if (is.null(sirius_path) || !file.exists(sirius_path)) {
    cat(sprintf("   SIRIUS not found at: %s\n", sirius_path))
    cat("   Set sirius.path in config or install SIRIUS from https://bio.informatik.uni-jena.de/software/sirius/\n")
    cat("   Skipping SIRIUS analysis\n")
    return(app3)
  }

  cat(sprintf("   SIRIUS found at: %s\n", sirius_path))

  # Export MS2 spectra to MGF
  mgf_dir <- file.path(cfg$project$output_dir, "sirius_input")
  dir.create(mgf_dir, showWarnings = FALSE, recursive = TRUE)
  mgf_file <- file.path(mgf_dir, "ms2_spectra.mgf")

  n_exported <- 0
  spectra_data <- object@ms2_data

  if (is.null(spectra_data) || length(spectra_data) == 0) {
    cat("   No MS2 spectra in object, cannot run SIRIUS\n")
    return(app3)
  }

  # Write MGF file
  con <- file(mgf_file, "w")
  on.exit(close(con), add = TRUE)

  for (i in seq_along(spectra_data)) {
    md <- spectra_data[[i]]
    if (is.null(md) || length(md@ms2_spectra) == 0) next

    for (j in seq_along(md@ms2_spectra)) {
      spec <- md@ms2_spectra[[j]]
      if (is.null(spec) || (!is.data.frame(spec) && !is.matrix(spec)) ||
          nrow(as.data.frame(spec)) == 0) next

      var_id <- md@variable_id[j]
      precursor_mz <- md@ms2_mz[j]
      precursor_rt <- md@ms2_rt[j]

      if (is.null(var_id) || is.na(var_id) || is.null(precursor_mz) || is.na(precursor_mz)) next

      cat("BEGIN IONS\n", file = con)
      cat(sprintf("TITLE=%s_scan%d\n", var_id, j), file = con)
      cat(sprintf("PEPMASS=%.6f\n", precursor_mz), file = con)
      cat(sprintf("RTINSECONDS=%.2f\n", if (is.null(precursor_rt) || is.na(precursor_rt)) 0 else precursor_rt), file = con)
      cat(sprintf("CHARGE=%s\n", ifelse(cfg$feature_extraction$polarity == "positive", "1+", "1-")), file = con)

      # Write fragment ions
      ms2_df <- as.data.frame(spec)
      for (k in seq_len(nrow(ms2_df))) {
        cat(sprintf("%.6f %.4f\n", ms2_df$mz[k], ms2_df$intensity[k]), file = con)
      }

      cat("END IONS\n", file = con)
      n_exported <- n_exported + 1
    }
  }
  cat(sprintf("   Exported %d MS2 spectra to: %s\n", n_exported, mgf_file))

  if (n_exported == 0) {
    cat("   No MS2 spectra to export\n")
    return(app3)
  }

  # Run SIRIUS analysis
  project_space <- cfg$sirius$project_space
  # SIRIUS 6.x expects -o to be a .sirius FILE path (a Nitrite project), not a
  # directory — passing a directory aborts with "is a directory, must be a file".
  sirius_output <- ifelse(grepl("\\.sirius$", project_space),
                          project_space,
                          file.path(project_space, "project.sirius"))
  dir.create(dirname(sirius_output), showWarnings = FALSE, recursive = TRUE)

  # Summary TSVs are NOT written into the project automatically in 6.3; they
  # must be exported with the `summaries` subcommand into a separate directory.
  summary_dir <- file.path(dirname(sirius_output), "summaries")

  # Check for existing SIRIUS project to enable caching. A pending structure
  # recompute (sirius.recompute_structure) bypasses the cache, because the
  # structure step must run again to pick up a changed database list.
  force_structure <- isTRUE(cfg$sirius$recompute_structure)
  if (file.exists(sirius_output) && !force_structure &&
      length(list.files(summary_dir, pattern = "\\.tsv$")) > 0) {
    cat(sprintf("   Existing SIRIUS project + summaries found: %s\n=> Skipping SIRIUS computation (cached).\n",
                sirius_output))
  } else {
    if (force_structure && file.exists(sirius_output)) {
      cat("   sirius.recompute_structure = TRUE — re-running the tool chain\n")
      cat("   (SIRIUS reuses formula/fingerprint/CANOPUS results and only redoes\n")
      cat("    the structure search against the new database list).\n")
    }
    # SIRIUS CLI --threads: global option (before the sub-command), read
    # from performance.sirius_threads (fallback 4).
    sirius_threads <- as.character(max(1L, as.integer(cfg$performance$sirius_threads %||% 4)))
    cat(sprintf("   SIRIUS threads: %s\n", sirius_threads))

    # SIRIUS 6.3 tool chain, in dependency order:
    #   formula      — fragmentation trees + isotope pattern → molecular formula
    #   fingerprints — ML-predicted molecular fingerprints (needed by structures)
    #   canopus      — compound class prediction (ClassyFire / NPC)
    #   structures   — CSI:FingerID database search (needs formula + fingerprints)
    # Each step is a separate invocation against the same project; SIRIUS
    # preserves prior results, so re-running is incremental.
    #
    # NOTE: the subcommand names are plural (`structures`, `fingerprints`); the
    # singular `structure`/`fingerprint` forms do not exist and abort the CLI.
    sirius_run <- function(tool_args, label) {
      cat(sprintf("   Launching SIRIUS: %s...\n", label))
      # SIRIUS 6.3.7 can finish writing the project and then never exit (observed:
      # structure search wrote project.sirius in 4 min, then sat idle for 7 h
      # with the JVM alive). Without a timeout that hang blocks the whole batch.
      # `timeout` sends SIGTERM after the budget, which is safe here because the
      # project is already persisted — the next invocation resumes from it.
      timeout_sec <- as.integer(cfg$sirius$timeout_sec %||% 3600)
      rc <- system2("timeout",
                    args = c(as.character(timeout_sec), sirius_path,
                             "--log", "WARNING",
                             "--threads", sirius_threads,
                             "-o", sirius_output,
                             tool_args),
                    stdout = TRUE, stderr = TRUE)
      status <- attr(rc, "status")
      if (!is.null(status) && status == 124) {
        cat(sprintf("   WARNING: SIRIUS %s hit the %ds timeout and was terminated.\n",
                    label, timeout_sec))
        cat("            Results written so far are kept; raise sirius.timeout_sec if needed.\n")
      } else if (!is.null(status) && status != 0) {
        cat(sprintf("   WARNING: SIRIUS %s returned exit code %d\n", label, status))
      }
      invisible(status %||% 0L)
    }

    # `-i` (the input MGF) is only needed on the first call that imports the
    # spectra; later tool calls reuse the project's stored features.
    cat(sprintf("   SIRIUS input: %s\n   SIRIUS output: %s\n", mgf_file, sirius_output))
    sirius_run(c("-i", mgf_file, "formula", "-p", "qtof",
                 "--ppm-max", as.character(cfg$sirius$ppm_max %||% 15)),
               "formula prediction")
    sirius_run(c("fingerprints"), "fingerprint prediction")
    sirius_run(c("canopus"), "compound class (CANOPUS) prediction")
    # CSI:FingerID structure database selection. This study is 161 TCM herbs, so
    # the natural-product collections (PLANTCYC, KNAPSACK, GNPS) would be
    # valuable — HMDB/CHEBI are human-centric and hold few plant secondary
    # metabolites.
    #
    # CAVEAT: those collections are NOT bundled with SIRIUS 6.3. Requesting one
    # that is not installed makes `structures` find zero compounds for every
    # feature, and combined with --recompute it can hang for hours with no
    # output. So the requested list is intersected with what the workspace
    # actually has, and a missing database is reported rather than silently
    # degrading the search.
    structure_db <- cfg$sirius$structure_db %||% "HMDB,CHEBI"
    requested <- trimws(strsplit(structure_db, ",")[[1]])
    available <- sirius_available_databases(sirius_path)
    if (length(available) > 0) {
      missing_db <- setdiff(requested, available)
      if (length(missing_db) > 0) {
        cat(sprintf("   NOTE: SIRIUS database(s) not installed, skipping: %s\n",
                    paste(missing_db, collapse = ", ")))
        cat("         Install with: sirius custom-db-downloader  (then `custom-db add`)\n")
      }
      keep <- intersect(requested, available)
      if (length(keep) == 0) {
        cat("   No requested structure database is installed; falling back to SIRIUS default (BIO)\n")
        keep <- "BIO"
      }
      structure_db <- paste(keep, collapse = ",")
    }
    cat(sprintf("   Structure search databases: %s\n", structure_db))

    # SIRIUS preserves already-computed results and skips those instances, so a
    # changed `-d` list silently has NO effect unless the search is forced. Set
    # sirius.recompute_structure=true once after switching databases.
    structure_args <- c("structures", "-d", structure_db)
    if (isTRUE(cfg$sirius$recompute_structure)) {
      structure_args <- c("--recompute", structure_args)
      cat("   (structure search forced to recompute)\n")
    }
    sirius_run(structure_args,
               sprintf("structure database search [%s]", structure_db))

    cat(sprintf("   SIRIUS analysis complete. Project saved to: %s\n", sirius_output))
  }

  # Export summary TSVs (formula/structure/canopus identifications).
  # Re-export when the structure search was just redone, otherwise the stale
  # TSVs from the previous database list would be merged instead.
  dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)
  if (force_structure || length(list.files(summary_dir, pattern = "\\.tsv$")) == 0) {
    if (force_structure) {
      unlink(list.files(summary_dir, pattern = "\\.tsv$", full.names = TRUE))
    }
    cat("   Exporting SIRIUS summary tables...\n")
    # `summaries` has also been seen to hang after writing its TSVs; bound it too.
    sum_timeout <- as.integer(cfg$sirius$timeout_sec %||% 3600)
    rc_sum <- system2("timeout",
                      args = c(as.character(sum_timeout), sirius_path,
                               "--log", "WARNING",
                               "-o", sirius_output,
                               "summaries", "--output", summary_dir),
                      stdout = TRUE, stderr = TRUE)
    status <- attr(rc_sum, "status")
    if (!is.null(status) && status == 124) {
      cat(sprintf("   WARNING: SIRIUS summaries hit the %ds timeout and was terminated.\n",
                  sum_timeout))
    } else if (!is.null(status) && status != 0) {
      cat(sprintf("   WARNING: SIRIUS summaries returned exit code %d\n", status))
    }
  }

  # Import SIRIUS results
  summary_files <- list.files(summary_dir, pattern = "\\.tsv$", full.names = TRUE)
  if (length(summary_files) > 0) {
    cat(sprintf("   Found %d SIRIUS summary file(s)\n", length(summary_files)))

    # SIRIUS 6.3 identifies each compound by `mappingFeatureId`
    # ("<internalId>_UNKNOWN_<MGF TITLE>"), not by a `title` column. The MGF
    # TITLE is "<variable_id>_scan<k>", so strip the leading internal id and
    # the trailing scan suffix to recover variable_id.
    sirius_var_id <- function(x) {
      x <- sub("^[0-9]+_UNKNOWN_", "", x)   # drop "<id>_UNKNOWN_" prefix
      sub("_scan[0-9]+$", "", x)            # drop "_scan<k>" suffix
    }

    for (sf in summary_files) {
      sirius_results <- tryCatch({
        read.delim(sf, stringsAsFactors = FALSE, check.names = FALSE)
      }, error = function(e) NULL)

      if (is.null(sirius_results) || nrow(sirius_results) == 0) next
      cat(sprintf("   SIRIUS results: %d rows from %s\n", nrow(sirius_results), basename(sf)))

      id_col <- intersect(c("mappingFeatureId", "compoundId"), colnames(sirius_results))[1]
      if (is.na(id_col)) {
        cat(sprintf("   No feature-id column in %s, skipping\n", basename(sf)))
        next
      }
      sirius_results$variable_id <- sirius_var_id(as.character(sirius_results[[id_col]]))

      if (is.null(app3) || nrow(app3) == 0) next

      merge_kind <- if (grepl("^structure_identifications", basename(sf))) {
        "structure"
      } else if (grepl("^formula_identifications", basename(sf))) {
        "formula"
      } else if (grepl("^canopus_structure_summary", basename(sf))) {
        "canopus"
      } else {
        NA_character_
      }
      if (is.na(merge_kind)) next

      rank_col <- switch(merge_kind,
                         structure = "structurePerIdRank",
                         formula   = "formulaRank",
                         canopus   = "formulaRank")
      if (rank_col %in% colnames(sirius_results)) {
        sirius_results <- sirius_results[order(sirius_results$variable_id,
                                               sirius_results[[rank_col]]), , drop = FALSE]
        sirius_results <- sirius_results[!duplicated(sirius_results$variable_id), , drop = FALSE]
      }

      sirius_merge <- switch(merge_kind,
        structure = sirius_results %>%
          transmute(
            variable_id,
            sirius_name      = .data$name,
            sirius_formula   = .data$molecularFormula,
            sirius_InChIkey  = .data$InChIkey2D,
            sirius_InChI     = .data$InChI,
            sirius_smiles    = .data$smiles,
            sirius_csi_score = .data[["CSI:FingerIDScore"]],
            sirius_confidence = .data$ConfidenceScoreExact
          ),
        formula = sirius_results %>%
          transmute(
            variable_id,
            sirius_formula       = .data$molecularFormula,
            sirius_formula_score = .data$SiriusScoreNormalized,
            sirius_zodiac_score  = .data$ZodiacScore
          ),
        # CANOPUS class prediction — valuable for TCM chemistry (flavonoids,
        # alkaloids, terpenoids, ...) even when no structure hit exists.
        canopus = sirius_results %>%
          transmute(
            variable_id,
            canopus_npc_pathway      = .data[["NPC#pathway"]],
            canopus_npc_superclass   = .data[["NPC#superclass"]],
            canopus_npc_class        = .data[["NPC#class"]],
            canopus_cf_superclass    = .data[["ClassyFire#superclass"]],
            canopus_cf_class         = .data[["ClassyFire#class"]],
            canopus_cf_most_specific = .data[["ClassyFire#most specific class"]]
          )
      )

      # Drop any columns already merged from an earlier summary file so the
      # left_join cannot create .x/.y duplicates (summary files are iterated in
      # filesystem order, and structure/formula/canopus columns overlap).
      dup_cols <- intersect(setdiff(colnames(sirius_merge), "variable_id"), colnames(app3))
      if (length(dup_cols) > 0) app3 <- app3[, setdiff(colnames(app3), dup_cols), drop = FALSE]

      app3 <- app3 %>% left_join(sirius_merge, by = "variable_id")

      if (merge_kind == "structure") {
        n_struct <- sum(!is.na(app3$sirius_name))
        cat(sprintf("   SIRIUS structure hits merged: %d features with a CSI:FingerID candidate\n",
                    n_struct))
        if (n_struct > 0) {
          # Upgrade Level 3/4 features that now have SIRIUS structure evidence.
          # A CSI:FingerID structure hit is MSI Level 2 (putatively identified);
          # a formula-only hit stays Level 3.
          app3 <- app3 %>%
            mutate(
              confidence_level = ifelse(
                !is.na(.data$sirius_name) &
                  .data$confidence_level %in% c("Level 3 — Putatively annotated",
                                                "Level 4 — Unknown"),
                "Level 2 — Putatively identified",
                .data$confidence_level
              )
            )
          cat("   Upgraded confidence levels for SIRIUS-validated features\n")
        }
      } else if (merge_kind == "formula") {
        cat(sprintf("   SIRIUS formula hits merged: %d features with a predicted formula\n",
                    sum(!is.na(app3$sirius_formula))))
      } else {
        cat(sprintf("   CANOPUS class predictions merged: %d features\n",
                    sum(!is.na(app3$canopus_cf_superclass))))
      }
    }
  } else {
    cat("   No SIRIUS summary files found. Results may need manual review.\n")
  }

  app3
}

# ══════════════════════════════════════════════════════════════════════════════
# 8c. Step 6c: MetFrag In Silico Fragmentation Validation (optional)
# ══════════════════════════════════════════════════════════════════════════════

#' Run MetFrag for in silico fragmentation validation
#'
#' MetFrag is a Java tool for in silico fragmentation search. It scores
#' candidate molecules by matching their theoretical fragments to observed
#' MS2 spectra against PubChem.
#'
#' Note: MetFrag via the metfRag R package may have Java/CDK version
#' incompatibilities (especially with Java ≥ 17). If fragmentation scoring
#' returns 0 results, consider using MetFrag standalone CLI instead.
#'
#' @param app3 Current annotation data.frame
#' @param object mass_dataset object with MS2 spectra
#' @param cfg Validated config list
#' @return Updated annotation data.frame with MetFrag columns (if available)
run_metfrag_analysis <- function(app3, object, cfg) {
  cat("\n===== Step 6c: MetFrag in silico fragmentation validation (optional) =====\n")

  if (!isTRUE(cfg$metfrag$enabled)) {
    cat("   MetFrag disabled in config, skipping\n")
    return(app3)
  }

  if (!requireNamespace("metfRag", quietly = TRUE)) {
    cat("   metfRag package not installed. Skipping.\n")
    cat("   Install with: install.packages('metfRag', repos='https://cloud.r-project.org')\n")
    return(app3)
  }
  suppressPackageStartupMessages(library(metfRag))

  if (is.null(app3) || nrow(app3) == 0) {
    cat("   No identifications to validate, skipping\n")
    return(app3)
  }

  # Check if we have MS2 spectra
  spectra_data <- object@ms2_data
  if (is.null(spectra_data) || length(spectra_data) == 0) {
    cat("   No MS2 spectra in object, cannot run MetFrag\n")
    return(app3)
  }

  cat(sprintf("   Validating %d identified features with MetFrag...\n", nrow(app3)))

  ann_cfg <- cfg$annotation
  metfrag_cfg <- cfg$metfrag
  polarity <- cfg$feature_extraction$polarity

  # Track which features get MetFrag scores
  app3$metfrag_score <- NA_real_
  app3$metfrag_num_explained_peaks <- NA_integer_
  app3$metfrag_fragmenter_score <- NA_real_

  # Process features that have both annotation and MS2
  all_ms2_features <- c()
  for (md in spectra_data) {
    if (isS4(md) && inherits(md, "ms2_data") && length(md@variable_id) > 0) {
      all_ms2_features <- c(all_ms2_features, md@variable_id)
    }
  }

  n_processed <- 0
  n_scored <- 0

  for (i in seq_len(nrow(app3))) {
    var_id <- app3$variable_id[i]
    if (!var_id %in% all_ms2_features) next

    # Find which ms2_data object contains this feature and get its spectrum
    ms2_spec <- NULL
    for (md in spectra_data) {
      if (isS4(md) && inherits(md, "ms2_data")) {
        idx <- which(md@variable_id == var_id)
        if (length(idx) > 0) {
          ms2_spec <- md@ms2_spectra[[idx[1]]]
          break
        }
      }
    }

    if (is.null(ms2_spec) || (!is.data.frame(ms2_spec) && !is.matrix(ms2_spec))) next

    # Get precursor m/z
    precursor_mz <- app3$mz[i]
    if (is.na(precursor_mz) || is.null(precursor_mz)) next

    # Get compound name
    compound_name <- app3$Compound.name[i]
    if (is.na(compound_name) || compound_name == "") next

    formula <- ""
    if ("Formula" %in% colnames(app3)) {
      formula <- as.character(app3$Formula[i])
      if (is.na(formula)) formula <- ""
    }

    ms2_df <- as.data.frame(ms2_spec)
    mzs <- ms2_df$mz
    ints <- ms2_df$intensity

    n_processed <- n_processed + 1

    # Build MetFrag settings
    s <- create.settings.sample()
    s[["MetFragDatabaseType"]] <- "PubChem"
    s[["DatabaseSearchRelativeMassDeviation"]] <- metfrag_cfg$ppm_tol %||% ann_cfg$ms1_match_ppm
    s[["FragmentPeakMatchAbsoluteMassDeviation"]] <- metfrag_cfg$mzabs %||% 0.005
    s[["FragmentPeakMatchRelativeMassDeviation"]] <- metfrag_cfg$ppm_tol %||% 15
    s[["NeutralPrecursorMass"]] <- precursor_mz
    s[["SampleName"]] <- compound_name

    # Set precursor compound IDs from formula if available
    if (!is.na(formula) && formula != "") {
      s[["NeutralPrecursorMolecularFormula"]] <- formula
    }

    # Peak list as matrix
    peak_mat <- cbind(mz = mzs, intensity = ints)
    s[["PeakList"]] <- peak_mat

    # Run MetFrag
    result <- tryCatch({
      run.metfrag(s)
    }, error = function(e) {
      if (n_processed <= 5) {
        cat(sprintf("     MetFrag error for %s: %s\n", compound_name, conditionMessage(e)))
      }
      return(NULL)
    })

    if (!is.null(result) && nrow(result) > 0) {
      best_row <- result[1, , drop = FALSE]

      if ("FragmenterScore" %in% colnames(best_row)) {
        app3$metfrag_fragmenter_score[i] <- as.numeric(best_row$FragmenterScore[1])
        n_scored <- n_scored + 1
      }
      if ("NumberExplainedPeaks" %in% colnames(best_row)) {
        app3$metfrag_num_explained_peaks[i] <- as.integer(best_row$NumberExplainedPeaks[1])
      }
      if ("Score" %in% colnames(best_row)) {
        app3$metfrag_score[i] <- as.numeric(best_row$Score[1])
      }

      if (n_processed %% 10 == 0) {
        cat(sprintf("   MetFrag: %d processed, %d scored...\n", n_processed, n_scored))
      }
    }
  }

  cat(sprintf("   MetFrag complete: %d features processed, %d scored\n", n_processed, n_scored))
  if (n_scored == 0 && n_processed > 0) {
    cat("   Note: MetFrag fragment scoring may have failed due to Java/CDK version incompatibility.\n")
    cat("   The 'metfRag' jar was compiled against CDK with JNI-InChI that is incompatible with Java 21.\n")
    cat("   PubChem candidate retrieval worked, but fragmentation scoring was skipped.\n")
    cat("   To fix: install a compatible Java version (<=17) or use MetFrag standalone CLI.\n")
  }

  # Update confidence levels for MetFrag-validated features
  if (n_scored > 0) {
    app3 <- app3 %>%
      rowwise() %>%
      mutate(
        confidence_level = ifelse(
          !is.na(metfrag_fragmenter_score) && metfrag_fragmenter_score > 0,
          "Level 1 — Confirmed",
          confidence_level
        )
      ) %>%
      ungroup()
    cat(sprintf("   Upgraded %d features to Level 1 based on MetFrag validation\n",
                sum(!is.na(app3$metfrag_fragmenter_score) & app3$metfrag_fragmenter_score > 0)))
  }

  app3
}

# ══════════════════════════════════════════════════════════════════════════════
# 8d. Step 6d: Identification Report Generation
# ══════════════════════════════════════════════════════════════════════════════

#' Generate comprehensive identification report
#'
#' Creates a detailed identification report with:
#'  - Multi-sheet Excel workbook (Summary, All, per-Level sheets)
#'  - Per-confidence-level CSV files
#'  - Text summary with identification statistics
#'
#' @param app3 Merged annotation data.frame with confidence levels
#' @param object mass_dataset object
#' @param cfg Validated config list
#' @return Invisible NULL
generate_identification_report <- function(app3, object, cfg) {
  cat("\n===== Step 6d: Generating identification report =====\n")

  out_dir <- file.path(cfg$project$output_dir, "06_Annotation")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (is.null(app3) || nrow(app3) == 0) {
    cat("   No identifications to report\n")
    return(invisible(NULL))
  }

  # ── Summary statistics ────────────────────────────────────────────────────
  # Counts are per FEATURE, not per candidate row: app3 may carry several
  # candidate annotations for the same variable_id, and row counts (4309 in
  # an earlier run for 1816 features) massively overstate identification.
  total_features <- nrow(object@expression_data)
  app3_feat <- app3 %>% distinct(variable_id, .keep_all = TRUE)
  n_candidates <- nrow(app3)
  n_identified <- nrow(app3_feat)
  n_level1 <- sum(app3_feat$confidence_level == "Level 1 — Confirmed", na.rm = TRUE)
  n_level2 <- sum(app3_feat$confidence_level == "Level 2 — Putatively identified", na.rm = TRUE)
  n_level3 <- sum(app3_feat$confidence_level == "Level 3 — Putatively annotated", na.rm = TRUE)
  n_level4 <- sum(app3_feat$confidence_level == "Level 4 — Unknown", na.rm = TRUE)
  n_unmatched <- total_features - n_identified

  unique_compounds <- unique(app3$Compound.name[!is.na(app3$Compound.name) &
                                                  app3$Compound.name != ""])

  # Per-database stats
  db_stats <- table(app3$database_source, useNA = "ifany")

  # SIRIUS stats (if available; per feature). The base annotation table (app3)
  # only keeps the MS1/MS2-database hit, so a feature that SIRIUS identified but
  # no MS1 library matched would be invisible in the report — and the exported
  # Compound.name would still show the (weaker) library hit. Prefer the SIRIUS
  # structure hit where one exists, and record its provenance.
  has_sirius <- "sirius_name" %in% colnames(app3_feat)
  if (has_sirius) {
    n_sirius_structure <- sum(!is.na(app3_feat$sirius_name))
    n_sirius_formula   <- sum(!is.na(app3_feat$sirius_formula))
    app3_feat$library_compound_name <- app3_feat$Compound.name
    use_sirius <- !is.na(app3_feat$sirius_name) & nzchar(app3_feat$sirius_name)
    app3_feat$Compound.name[use_sirius] <- app3_feat$sirius_name[use_sirius]
    app3_feat$compound_source <- ifelse(use_sirius, "SIRIUS_CSI:FingerID", "MS_library")
    app3_feat$formula[use_sirius] <- app3_feat$sirius_formula[use_sirius]
    # Mirror the same preference into the full candidate table so the CSV/Excel
    # and the web export agree with the summary counts.
    if ("sirius_name" %in% colnames(app3)) {
      use_sirius_all <- !is.na(app3$sirius_name) & nzchar(app3$sirius_name)
      app3$library_compound_name <- app3$Compound.name
      app3$Compound.name[use_sirius_all] <- app3$sirius_name[use_sirius_all]
      app3$compound_source <- ifelse(use_sirius_all, "SIRIUS_CSI:FingerID", "MS_library")
      if (!"formula" %in% colnames(app3)) app3$formula <- NA_character_
      app3$formula[use_sirius_all] <- app3$sirius_formula[use_sirius_all]
    }
  } else {
    n_sirius_structure <- 0
    n_sirius_formula <- 0
  }

  # SIRIUS stats (if available; per feature)
  n_sirius <- if ("sirius_formula" %in% colnames(app3_feat)) {
    sum(!is.na(app3_feat$sirius_formula))
  } else 0

  # MetFrag stats (if available; per feature)
  n_metfrag <- if ("metfrag_fragmenter_score" %in% colnames(app3_feat)) {
    sum(!is.na(app3_feat$metfrag_fragmenter_score) &
          app3_feat$metfrag_fragmenter_score > 0)
  } else 0

  # ── Write text summary ────────────────────────────────────────────────────
  summary_lines <- c(
    "============================================",
    sprintf("  Compound Identification Report — %s", cfg$project$name),
    sprintf("  Generated: %s", Sys.time()),
    "============================================",
    "",
    sprintf("Total MS1 features: %d", total_features),
    sprintf("Features with identification: %d", n_identified),
    sprintf("Candidate annotation rows: %d", n_candidates),
    sprintf("Unmatched features: %d", n_unmatched),
    sprintf("Unique compound names: %d", length(unique_compounds)),
    "",
    "--- MSI Confidence Level Breakdown ---",
    sprintf("  Level 1 — Confirmed (MS2 + RT + in-house):  %d", n_level1),
    sprintf("  Level 2 — Putatively identified (MS2/SIRIUS): %d", n_level2),
    sprintf("  Level 3 — Putatively annotated (MS1 only):  %d", n_level3),
    sprintf("  Level 4 — Unknown:                           %d", n_level4),
    "",
    sprintf("  SIRIUS-validated features: %d", n_sirius),
    sprintf("    of which CSI:FingerID structure hits: %d", n_sirius_structure),
    sprintf("    of which formula-only:                 %d", n_sirius_formula),
    sprintf("  MetFrag-validated features: %d", n_metfrag),
    "",
    "--- Per-Database Matches ---",
    paste(capture.output(print(db_stats)), collapse = "\n"),
    ""
  )

  writeLines(summary_lines, file.path(out_dir, "identification_summary.txt"))
  cat(sprintf("   Summary saved: %s\n", file.path(out_dir, "identification_summary.txt")))

  # ── Build multi-sheet Excel report ────────────────────────────────────────
  wb <- createWorkbook()

  # Summary sheet
  addWorksheet(wb, "Summary")
  summary_df <- data.frame(
    Metric = c("Total MS1 features", "Features identified",
               "Candidate annotation rows", "Unmatched features",
               "Unique compound names",
               "Level 1 — Confirmed", "Level 2 — Putatively identified",
               "Level 3 — Putatively annotated", "Level 4 — Unknown",
               "SIRIUS-validated (formula)", "  of which CSI:FingerID structure",
               "MetFrag-validated"),
    Value = c(total_features, n_identified, n_candidates, n_unmatched,
              length(unique_compounds),
              n_level1, n_level2, n_level3, n_level4,
              n_sirius, n_sirius_structure, n_metfrag),
    stringsAsFactors = FALSE
  )
  writeData(wb, "Summary", summary_df)

  # All identifications sheet
  addWorksheet(wb, "All_Identifications")
  writeData(wb, "All_Identifications", app3)

  # Per-level sheets
  for (lvl in unique(app3$confidence_level)) {
    lvl_short <- gsub(" .*$", "", gsub("Level ", "L", lvl))
    lvl_data <- app3[app3$confidence_level == lvl, , drop = FALSE]
    if (nrow(lvl_data) > 0) {
      addWorksheet(wb, lvl_short)
      writeData(wb, lvl_short, lvl_data)
    }
  }

  # Per-database sheets
  for (db_src in unique(app3$database_source)) {
    db_data <- app3[app3$database_source == db_src, , drop = FALSE]
    if (nrow(db_data) > 0) {
      sheet_name <- substr(db_src, 1, 31)  # Excel sheet name limit
      addWorksheet(wb, sheet_name)
      writeData(wb, sheet_name, db_data)
    }
  }

  saveWorkbook(wb, file.path(out_dir, "Compound_identification_report.xlsx"), overwrite = TRUE)
  cat(sprintf("   Excel report saved: %s\n", file.path(out_dir, "Compound_identification_report.xlsx")))

  # Print summary to console
  for (line in summary_lines) {
    if (nchar(line) > 0) cat(sprintf("   %s\n", line))
  }

  invisible(NULL)
}

# ══════════════════════════════════════════════════════════════════════════════
# 9. Step 6.5: Annotated volcano plots (top 5 labeled metabolites)
# ══════════════════════════════════════════════════════════════════════════════

run_annotated_volcano <- function(diff_results, app3, cfg) {
  cat("\n===== Step 6.5: Annotated volcano plots =====\n")

  if (is.null(app3) || nrow(app3) == 0) {
    cat("   No annotations available, skipping\n")
    return(invisible(NULL))
  }

  out_dir <- file.path(cfg$project$output_dir, "04_Differential")
  diff_cfg <- cfg$differential
  alpha <- diff_cfg$p_value_cutoff
  logfc_cutoff <- log2(diff_cfg$fc_threshold)

  # Mirror the significance metric used in run_differential() so labeled
  # volcanoes color/threshold on the same p-column as the sig-feature tables.
  sig_metric <- if (!is.null(diff_cfg$significance_metric)) diff_cfg$significance_metric else "p_value"
  p_col <- if (identical(sig_metric, "adj_p_value")) "adj.P.Val" else "P.Value"
  p_lab <- if (identical(sig_metric, "adj_p_value")) "adj. p-value (FDR)" else "p-value"

  # one annotation per feature (first candidate) — multi-candidate rows in
  # app3 otherwise duplicate volcano points via the left_join below
  annot_lookup <- app3 %>%
    dplyr::select(variable_id, Compound.name) %>%
    filter(!is.na(Compound.name) & Compound.name != "") %>%
    distinct(variable_id, .keep_all = TRUE)

  for (cmp_name in names(diff_results)) {
    dr <- diff_results[[cmp_name]]
    if (is.null(dr)) next

    diff_all <- dr$all
    label <- dr$comparison$label
    cat(sprintf("   Annotating: %s\n", label))

    diff_all$variable_id <- rownames(diff_all)
    diff_all <- diff_all %>%
      left_join(annot_lookup, by = "variable_id") %>%
      mutate(label = ifelse(is.na(Compound.name) | Compound.name == "",
                            variable_id, Compound.name))

    diff_all$.plot_p <- diff_all[[p_col]]
    diff_all$Significant <- ifelse(
      diff_all$.plot_p < alpha & abs(diff_all$logFC) > logfc_cutoff,
      ifelse(diff_all$logFC > logfc_cutoff, "Up", "Down"), "Not")

    sig_with_name <- diff_all %>%
      filter(Significant != "Not") %>%
      arrange(.plot_p) %>%
      head(20) %>%
      filter(!is.na(Compound.name) & Compound.name != "") %>%
      head(5)

    p <- ggplot(diff_all, aes(x = logFC, y = -log10(.plot_p))) +
      geom_point(alpha = 0.4, size = 3.5, aes(color = Significant)) +
      ylab(sprintf("-log10(%s)", p_lab)) +
      scale_color_manual(values = c("blue4", "grey", "red3")) +
      geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff),
                 lty = 4, col = "black", lwd = 0.8) +
      geom_hline(yintercept = -log10(alpha),
                 lty = 4, col = "black", lwd = 0.8) +
      labs(title = label) +
      theme_bw() +
      theme(aspect.ratio = 1, panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 14))

    if (nrow(sig_with_name) > 0) {
      p <- p + geom_text_repel(
        data = sig_with_name, aes(label = label),
        size = 3.5, max.overlaps = 10,
        box.padding = 0.5, point.padding = 0.3,
        force = 2, segment.color = "grey50")
    }

    save_plot(file.path(out_dir, paste0(cmp_name, "_volcano_labeled")),
              plot = p, width = 10, height = 8, dpi = cfg$visualization$dpi)
  }

  invisible(NULL)
}

# ══════════════════════════════════════════════════════════════════════════════
# 10. Step 7: KEGG Pathway Enrichment
# ══════════════════════════════════════════════════════════════════════════════

run_kegg_enrichment <- function(app3, diff_results, cfg) {
  cat("\n===== Step 7: Pathway enrichment (HMDB/SMPDB; legacy dir name 07_KEGG) =====\n")

  if (!exists("hmdb_pathway")) {
    cat("   hmdb_pathway not loaded, skipping KEGG enrichment\n")
    return(invisible(NULL))
  }
  if (is.null(app3) || nrow(app3) == 0) {
    cat("   No annotations, skipping\n")
    return(invisible(NULL))
  }

  out_dir <- file.path(cfg$project$output_dir, "07_KEGG")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  for (cmp_name in names(diff_results)) {
    dr <- diff_results[[cmp_name]]
    if (is.null(dr) || is.null(dr$sig) || nrow(dr$sig) == 0) next

    label <- dr$comparison$label
    cat(sprintf("   Enrichment: %s\n", label))

    diff_peaks <- dr$sig
    diff_peaks$variable_id <- rownames(diff_peaks)

    # one annotation per feature before joining, so a multi-candidate feature
    # contributes exactly one row (same first candidate as the volcano plots)
    annotated <- merge(app3 %>% distinct(variable_id, .keep_all = TRUE),
                       diff_peaks, by = "variable_id") %>%
      distinct(Compound.name, .keep_all = TRUE)
    if (nrow(annotated) == 0) {
      cat("     No annotated matches, skipping\n")
      next
    }

    write.csv(annotated, file.path(out_dir, paste0("annotated_", cmp_name, ".csv")),
              row.names = FALSE)

    HMDB_ids <- unique(annotated$HMDB.ID)
    HMDB_ids <- HMDB_ids[!is.na(HMDB_ids) & HMDB_ids != ""]

    if (length(HMDB_ids) > 0) {
      result <- tryCatch(
        enrich_hmdb(
          query_id = HMDB_ids, query_type = "compound",
          id_type = "HMDB", pathway_database = hmdb_pathway,
          only_primary_pathway = TRUE, p_cutoff = 0.99,
          p_adjust_method = "BH"
        ),
        error = function(e) NULL
      )

      # enrich_hmdb() returns NULL when no pathway is enriched for this query
      # (and may error on degenerate input), so guard before slot access.
      if (is.null(result)) {
        cat("     No enriched pathways for this comparison, skipping\n")
        next
      }

      HMDB_pathway <- result@result
      if (is.null(HMDB_pathway) || nrow(HMDB_pathway) == 0) {
        cat("     No enriched pathways for this comparison, skipping\n")
        next
      }
      # Filter on BH-adjusted q-values. The tidymass enrich_hmdb result carries
      # no adjusted column (verified 2026-09-10 batch1 run), so compute BH over
      # this query's pathways ourselves; use the result's own column if it ever
      # appears. Cutoff: kegg.q_value_cutoff (default 0.05).
      path_alpha <- cfg$kegg$q_value_cutoff %||% 0.05
      adj_col <- NULL
      for (cc in c("q_value", "p_adjust", "fdr", "p_fdr")) {
        if (cc %in% colnames(HMDB_pathway)) { adj_col <- cc; break }
      }
      if (is.null(adj_col)) {
        HMDB_pathway$q_value <- p.adjust(HMDB_pathway$p_value, method = "BH")
        adj_col <- "q_value"
        cat("     (adjusted-p computed in-pipeline: BH over query pathways)\n")
      }
      HMDB_pathway <- HMDB_pathway[HMDB_pathway[[adj_col]] < path_alpha, , drop = FALSE]
      HMDB_pathway <- arrange(HMDB_pathway, desc(mapped_number))

      if (nrow(HMDB_pathway) > 0) {
        write.xlsx(HMDB_pathway, file.path(out_dir, paste0("kegg_pathway_", cmp_name, ".xlsx")))

        HMDB_plot <- if (nrow(HMDB_pathway) > 10) HMDB_pathway[1:10, ] else HMDB_pathway
        HMDB_plot$pathway_name <- factor(HMDB_plot$pathway_name,
                                         levels = rev(unique(HMDB_plot$pathway_name)))

        p <- ggplot(HMDB_plot, aes(x = mapped_number, y = pathway_name)) +
          geom_point(aes(size = p_value, color = mapped_number)) +
          scale_color_gradientn(colours = c("#f7ca64", "#46bac2", "#7e62a3")) +
          labs(color = "Mapped number", size = "p-value",
               x = "Count Number", title = label) +
          theme_bw() +
          theme(axis.text.y = element_text(size = rel(1.5)),
                axis.title.x = element_text(size = rel(1.5)),
                axis.title.y = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 14)) +
          scale_size(range = c(5, 10))

        save_plot(file.path(out_dir, paste0("kegg_pathway_", cmp_name)),
                  plot = p, width = 13, height = 8, dpi = cfg$visualization$dpi)
        cat(sprintf("     %d significant pathways\n", nrow(HMDB_pathway)))
      } else {
        cat("     No significant pathways\n")
      }
    } else {
      cat("     No HMDB IDs for enrichment\n")
    }
  }

  invisible(NULL)
}

# ══════════════════════════════════════════════════════════════════════════════
# 11. Step 8: Heatmap
# ══════════════════════════════════════════════════════════════════════════════

run_heatmap <- function(object2, sample_info, app3, diff_results, cfg) {
  cat("\n===== Step 8: Heatmap =====\n")

  if (is.null(app3) || nrow(app3) == 0) {
    cat("   No annotations, skipping heatmaps\n")
    return(invisible(NULL))
  }

  out_dir <- file.path(cfg$project$output_dir, "08_Heatmap")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  viz_cfg <- cfg$visualization

  all_samples <- colnames(object2@expression_data)
  expr_for_merge <- object2@expression_data
  expr_for_merge$variable_id <- rownames(expr_for_merge)

  for (cmp_name in names(diff_results)) {
    dr <- diff_results[[cmp_name]]
    if (is.null(dr) || is.null(dr$sig) || nrow(dr$sig) == 0) next

    diff_peaks <- dr$sig
    diff_peaks$variable_id <- rownames(diff_peaks)

    # one annotation per feature before joining (see KEGG step for rationale)
    annotated <- merge(app3 %>% distinct(variable_id, .keep_all = TRUE),
                       diff_peaks, by = "variable_id") %>%
      distinct(Compound.name, .keep_all = TRUE)
    if (nrow(annotated) == 0) next

    annotated <- merge(annotated, expr_for_merge, by = "variable_id")
    annotated <- annotated %>% arrange(desc(logFC))

    top_n <- viz_cfg$heatmap_top_n
    top_up <- annotated %>% filter(logFC > 0) %>% head(top_n)
    top_down <- annotated %>% filter(logFC < 0) %>% arrange(logFC) %>% head(top_n)
    annotated_top <- bind_rows(top_up, top_down)

    if (nrow(annotated_top) < 2) {
      cat(sprintf("   Heatmap: %s - insufficient metabolites, skipping\n", dr$comparison$label))
      next
    }

    core_cols <- all_samples[all_samples %in% colnames(annotated_top)]
    if (length(core_cols) < 2) next

    heat_data <- annotated_top[, core_cols, drop = FALSE]
    rownames(heat_data) <- annotated_top$Compound.name
    heat_data <- heat_data[rowSums(is.na(heat_data)) < ncol(heat_data), , drop = FALSE]
    if (nrow(heat_data) < 2) next

    annotation_col <- data.frame(
      Group = factor(sample_info$class[match(core_cols, sample_info$sample_id)]),
      row.names = core_cols
    )

    # Build annotation colors from config
    all_classes <- unique(sample_info$class[match(core_cols, sample_info$sample_id)])
    viz_colors <- cfg$visualization$group_colors
    ann_colors <- list(Group = setNames(
      sapply(all_classes, function(g) {
        if (!is.null(viz_colors[[g]])) viz_colors[[g]] else "#999999"
      }),
      all_classes
    ))

    cell_w <- 14; cell_h <- 12
    n_rows <- nrow(heat_data); n_cols <- ncol(heat_data)
    hm_height <- n_rows * cell_h / 72 + 2.5
    hm_width  <- max(n_cols * cell_w / 72 + 4.5, 9)

    tryCatch({
      for (fmt in c("pdf", "png")) {
        fpath <- file.path(out_dir, paste0("heatmap_", cmp_name, ".", fmt))
        if (fmt == "pdf") {
          pdf(fpath, width = hm_width, height = hm_height)
        } else {
          png(fpath, width = hm_width, height = hm_height,
              units = "in", res = viz_cfg$dpi)
        }
        pheatmap(heat_data,
                 cluster_cols = TRUE, cluster_rows = TRUE,
                 show_colnames = TRUE, scale = "row",
                 annotation_col = annotation_col,
                 annotation_colors = ann_colors,
                 fontsize = 8, fontsize_row = 8, fontsize_col = 8,
                 color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                 border_color = "grey", cellwidth = cell_w, cellheight = cell_h,
                 main = dr$comparison$label)
        dev.off()
      }
      cat(sprintf("   Heatmap: %s (%d metabolites)\n", dr$comparison$label, nrow(heat_data)))
    }, error = function(e) {
      cat(sprintf("   Heatmap error: %s - %s\n", dr$comparison$label, conditionMessage(e)))
    })
  }

  invisible(NULL)
}

# ══════════════════════════════════════════════════════════════════════════════
# 12. Step 9: Boxplot (key metabolite expression)
# ══════════════════════════════════════════════════════════════════════════════

run_boxplot <- function(object2, sample_info, app3, diff_results, cfg) {
  cat("\n===== Step 9: Boxplot =====\n")

  if (is.null(app3) || nrow(app3) == 0) {
    cat("   No annotations, skipping boxplots\n")
    return(invisible(NULL))
  }

  out_dir <- file.path(cfg$project$output_dir, "09_Boxplot")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  viz_cfg <- cfg$visualization

  all_samples <- colnames(object2@expression_data)
  core_group_vec <- sample_info$class[match(all_samples, sample_info$sample_id)]

  expre <- object2@expression_data
  expre$variable_id <- rownames(expre)
  app3_expr <- left_join(app3, expre, by = "variable_id")

  # Collect top-N metabolites from all comparisons
  metabolite_set <- character(0)
  for (cmp_name in names(diff_results)) {
    dr <- diff_results[[cmp_name]]
    if (is.null(dr) || is.null(dr$sig) || nrow(dr$sig) == 0) next
    top_sig <- dr$sig %>%
      arrange(P.Value) %>%
      head(5)
    metabolite_set <- c(metabolite_set, rownames(top_sig))
  }
  metabolite_set <- unique(metabolite_set)

  # Deduplicate by compound name
  if (length(metabolite_set) > 0) {
    name_map <- app3 %>%
      filter(!is.na(Compound.name) & Compound.name != "") %>%
      dplyr::select(variable_id, Compound.name) %>%
      distinct(Compound.name, .keep_all = TRUE)

    has_name <- metabolite_set %in% name_map$variable_id
    named_vars <- metabolite_set[has_name]
    unnamed_vars <- metabolite_set[!has_name]

    named_dedup <- name_map %>%
      filter(variable_id %in% named_vars) %>%
      distinct(Compound.name, .keep_all = TRUE) %>%
      pull(variable_id)

    plot_metabolites <- unique(c(named_dedup, unnamed_vars))
    plot_metabolites <- head(plot_metabolites, viz_cfg$boxplot_top_n)
    cat(sprintf("   Drawing %d metabolite boxplots\n", length(plot_metabolites)))

    # Color map
    viz_colors <- cfg$visualization$group_colors
    all_groups <- unique(core_group_vec)
    auto_palette <- c("#D20A13", "#088247", "#FFD121", "#7E6148FF", "#5BC0EB",
                      "#F39B7FFF", "#BC3C29FF", "#0072B5FF", "#E18727FF")
    mycol <- sapply(all_groups, function(g) {
      if (!is.null(viz_colors[[g]])) viz_colors[[g]]
      else auto_palette[(which(all_groups == g) - 1) %% length(auto_palette) + 1]
    })
    names(mycol) <- all_groups

    viz_labels <- cfg$visualization$group_labels

    for (vid in plot_metabolites) {
      meta_row <- app3_expr[app3_expr$variable_id == vid, ]
      if (nrow(meta_row) == 0) next

      compound_name <- ifelse(
        !is.na(meta_row$Compound.name[1]) && meta_row$Compound.name[1] != "",
        meta_row$Compound.name[1], vid)

      box_values <- as.numeric(meta_row[1, all_samples])
      box_df <- data.frame(
        value = box_values,
        group = core_group_vec,
        stringsAsFactors = FALSE
      )
      box_df <- box_df[!is.na(box_df$value), ]
      if (nrow(box_df) < 4 || length(unique(box_df$group)) < 2) next

      p <- ggplot(data = box_df, aes(x = group, y = value, color = group)) +
        geom_boxplot(aes(fill = group), color = NA, outlier.color = NA, alpha = 1) +
        stat_summary(fun = mean,
                     fun.min = function(x) mean(x) - sd(x),
                     fun.max = function(x) mean(x) + sd(x),
                     geom = "errorbar", width = 0.6, size = 1.5, alpha = 1) +
        stat_summary(fun = median, geom = "errorbar",
                     aes(ymax = after_stat(y), ymin = after_stat(y)),
                     width = 0.8, color = "white", size = 1) +
        scale_fill_manual(values = mycol, labels = viz_labels[names(mycol)]) +
        scale_color_manual(values = mycol, labels = viz_labels[names(mycol)]) +
        labs(title = compound_name, y = "Normalized intensity", x = "") +
        theme_bw() +
        theme(panel.grid = element_blank(),
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14),
              aspect.ratio = 1,
              plot.title = element_text(hjust = 0.5, size = 12))

      save_plot(file.path(out_dir, paste0("boxplot_", gsub("[^A-Za-z0-9]", "_", compound_name))),
                plot = p, width = 8, height = 7, dpi = viz_cfg$dpi)
    }
    cat("   Boxplots complete\n")
  } else {
    cat("   No significant metabolites to plot\n")
  }

  invisible(NULL)
}

# ══════════════════════════════════════════════════════════════════════════════
# 13. Main Orchestrator
# ══════════════════════════════════════════════════════════════════════════════

#' Resolve a helper-module path (e.g. "export_for_web.R").
#'
#' Modules are loaded with paths relative to the working directory
#' (<working_dir>/code/<file>), which is how test_run/ is laid out. When the
#' pipeline runs from code/this_project (the documented run_all_batches.R
#' entry point) there is no code/ subdir there, so fall back to the directory
#' of the running script — otherwise every optional module (OPLS-DA,
#' performance, MetaboAnalystR, web export) is silently skipped.
resolve_module_path <- function(filename) {
  candidates <- file.path("code", filename)
  if (exists("PIPELINE_SCRIPT_DIR", inherits = TRUE)) {
    candidates <- c(candidates, file.path(PIPELINE_SCRIPT_DIR, filename))
  }
  for (candidate in candidates) {
    if (file.exists(candidate)) return(candidate)
  }
  candidates[1]
}

#' Main pipeline entry point
#'
#' Orchestrates all analysis steps. Sources oplsda_functions.R for
#' OPLS-DA analysis (Step 3.5) and annotation update (Step 6.5).
#'
#' @param cfg Validated configuration list from load_and_validate_config()
main <- function(cfg, qualitative_only = FALSE) {
  cat("\n========================================\n")
  cat(sprintf("  Metabolomics Pipeline: %s\n", cfg$project$name))
  if (qualitative_only) {
    cat("  MODE: qualitative-only — intensity-dependent steps are SKIPPED\n")
  }
  cat("========================================\n")

  if (qualitative_only) {
    cat("\n⚠️  --qualitative-only is active for this batch.\n")
    cat("    No quantitative result produced here can be compared between groups:\n")
    cat("      * QC-RSD filter     SKIPPED (its only purpose is guarding quantitation)\n")
    cat("      * group-detection   SKIPPED (per-group counts need reproducible peaks)\n")
    cat("      * differential/OPLS-DA/KEGG/heatmap/boxplot/MetaboAnalystR  SKIPPED\n")
    cat("      * annotation and the QC diagnostic PCA  KEPT\n")
    cat("    The manifest will mark quantification_valid = false, which excludes\n")
    cat("    this batch from cross-batch integration automatically.\n\n")
  }

  # ── Load helper source ─────────────────────────────────────────────────
  oplsda_source <- resolve_module_path("oplsda_functions.R")
  if (file.exists(oplsda_source)) {
    source(oplsda_source)
    cat(sprintf("=> Sourced OPLS-DA functions: %s\n", oplsda_source))
  } else {
    cat(sprintf("⚠️  OPLS-DA functions not found at %s, Step 3.5 will be skipped\n",
                oplsda_source))
  }

  # ── Initialize parallel backend (performance module) ───────────────────
  perf_source <- resolve_module_path("performance_utils.R")
  if (file.exists(perf_source)) {
    source(perf_source)
    init_parallel_backend(cfg)
    # Install fast I/O wrappers if enabled
    install_fast_io(cfg)
  } else {
    cat(sprintf("⚠️  Performance module not found at %s, running sequentially\n",
                perf_source))
  }

  # ── Create output directory tree ────────────────────────────────────────
  out_root <- cfg$project$output_dir
  for (sub in c("03_PCA", "04_Differential", "05_Overlap",
                "06_Annotation", "07_KEGG", "08_Heatmap", "09_Boxplot",
                "03.5_OPLSDA", "06.5_OPLSDA_annotated",
                "sirius_input", "sirius_project",
                "10_MetaboAnalystR_Advanced")) {
    dir.create(file.path(out_root, sub), showWarnings = FALSE, recursive = TRUE)
  }
  cat(sprintf("=> Output directory: %s\n", out_root))

  # ── Load metadata ──────────────────────────────────────────────────────
  meta_list <- load_metadata(cfg)
  sample_info_all <- build_sample_info_all(meta_list, cfg)

  # ── Step 1: Feature extraction (MS1-only isolation) ─────────────────────
  run_feature_extraction(cfg, meta_list)

  # ── Step 2: Preprocessing ──────────────────────────────────────────────
  prep <- run_preprocessing(cfg, meta_list, sample_info_all,
                            qualitative_only = qualitative_only)
  object2     <- prep$object2
  sample_info <- prep$sample_info
  qc_quality  <- prep$qc_quality

  # ── Generate comparisons ───────────────────────────────────────────────
  core_classes <- prep$core_classes
  comparisons <- generate_comparisons(cfg, core_classes)

  # ── Step 3: PCA ────────────────────────────────────────────────────────
  run_pca(object2, sample_info, cfg)

  # Persist the normalized (SVR/LOESS/median-corrected) matrix for the core
  # samples. cross_batch_integration.R discovers normalized_expression_matrix
  # .csv under each batch's output tree and uses it for ComBat / median batch
  # correction and the integrated PCA; without this file that stage silently
  # skips. Values are on the original intensity scale (cross_batch auto-log2s).
  expr_out <- file.path(out_root, "03_PCA", "normalized_expression_matrix.csv")
  dir.create(dirname(expr_out), showWarnings = FALSE, recursive = TRUE)
  write.csv(as.data.frame(object2@expression_data), expr_out, row.names = TRUE)
  cat(sprintf("=> Saved normalized expression matrix: %s\n", basename(expr_out)))

  # ── Step 3.5: OPLS-DA ─────────────────────────────────────────────────
  if (qualitative_only) {
    cat("\n===== Step 3.5: OPLS-DA SKIPPED (qualitative-only) =====\n")
  } else if (isTRUE(cfg$oplsda$enabled) && exists("run_oplsda_analysis")) {
    # Filter to subset if specified
    oplsda_comps <- comparisons
    subset_comp <- cfg$oplsda$comparison_subset
    if (!is.null(subset_comp) && length(subset_comp) > 0) {
      subset_names <- as.character(unlist(subset_comp))
      oplsda_comps <- comparisons[sapply(comparisons, function(c) c$name %in% subset_names)]
      cat(sprintf("\n   OPLS-DA subset: %d comparisons\n", length(oplsda_comps)))
    }

    if (length(oplsda_comps) > 0) {
      comparison_labels <- setNames(
        vapply(comparisons, function(c) if (is.null(c$label) || is.na(c$label)) c$name else c$label, character(1)),
        vapply(comparisons, function(c) c$name, character(1))
      )
      oplsda_metrics <- run_oplsda_analysis(
        expression_data   = object2@expression_data,
        sample_info       = sample_info,
        comparisons       = oplsda_comps,
        comparison_labels = comparison_labels,
        output_dir        = file.path(cfg$project$output_dir, "03.5_OPLSDA"),
        n_perm            = cfg$oplsda$n_permutations,
        vip_threshold     = cfg$oplsda$vip_threshold
      )
    }
  } else {
    cat("\n===== Step 3.5: OPLS-DA disabled, skipping =====\n")
  }

  # ── Step 4: Differential analysis ─────────────────────────────────────
  if (qualitative_only) {
    # limma compares group means; with QC-RSD median in the tens of percent the
    # fold changes are dominated by injection variability, so any p-value here
    # would be spurious. Skipping is the point of this mode — not a limitation.
    cat("\n===== Step 4: Differential analysis SKIPPED (qualitative-only) =====\n")
    diff_results <- list()
  } else {
    diff_results <- run_differential(object2, sample_info, comparisons, cfg)
  }

  # ── Step 5: Overlap analysis ──────────────────────────────────────────
  if (qualitative_only) {
    cat("===== Step 5: Overlap analysis SKIPPED (qualitative-only) =====\n")
  } else {
    run_overlap_analysis(diff_results, cfg)
  }

  # ════════════════════════════════════════════════════════════════════════
  # Step 6: Enhanced Annotation & Compound Identification
  # ════════════════════════════════════════════════════════════════════════

  # ── Step 6a: CAMERA adduct/isotope grouping (optional) ─────────────────
  run_camera_annotation(object2, cfg)

  # ── Step 6: Multi-database annotation + MSI confidence ─────────────────
  # Kept in qualitative-only mode: matching is on the m/z axis, which depends
  # on mass calibration rather than on signal intensity, so a sensitivity-
  # degraded batch still yields meaningful putative annotations.
  app3 <- run_annotation(object2, cfg)

  # ── Step 6b: SIRIUS formula/structure prediction (optional) ────────────
  app3 <- run_sirius_analysis(object2, app3, cfg)

  # ── Step 6c: MetFrag in silico validation (optional) ───────────────────
  app3 <- run_metfrag_analysis(app3, object2, cfg)

  # ── Step 6d: Identification report ─────────────────────────────────────
  generate_identification_report(app3, object2, cfg)

  # ── Step 6.5: Annotated volcano + OPLS-DA update ──────────────────────
  if (qualitative_only || length(diff_results) == 0) {
    cat("===== Step 6.5: Annotated volcano SKIPPED (no differential results) =====\n")
  } else {
    run_annotated_volcano(diff_results, app3, cfg)
  }

  if (!qualitative_only && !is.null(app3) && nrow(app3) > 0 &&
      exists("update_oplsda_with_annotation") &&
      isTRUE(cfg$oplsda$enabled)) {
    update_oplsda_with_annotation(
      oplsda_dir    = file.path(out_root, "03.5_OPLSDA"),
      output_dir    = file.path(out_root, "06.5_OPLSDA_annotated"),
      app3          = app3,
      vip_threshold = cfg$oplsda$vip_threshold
    )
  }

  # ── Step 7: KEGG enrichment ───────────────────────────────────────────
  if (qualitative_only || length(diff_results) == 0) {
    cat("===== Step 7: KEGG enrichment SKIPPED (no differential results) =====\n")
  } else {
    run_kegg_enrichment(app3, diff_results, cfg)
  }

  # ── Step 10: MetaboAnalystR Advanced Analysis ─────────────────────────
  if (qualitative_only) {
    cat("===== Step 10: MetaboAnalystR SKIPPED (qualitative-only) =====\n")
  } else {
    # Source the MetaboAnalystR module (separate file for modularity)
    ma_source <- resolve_module_path("metaboanalyst_advanced.R")
    if (file.exists(ma_source)) {
      source(ma_source)
      run_metaboanalyst_advanced(object2, sample_info, diff_results, app3, cfg)
    } else {
      cat(sprintf("⚠️  MetaboAnalystR module not found at %s, skipping Step 10\n",
                  ma_source))
    }
  }

  # ── Step 8: Heatmap ──────────────────────────────────────────────────
  if (qualitative_only || length(diff_results) == 0) {
    cat("===== Step 8: Heatmap SKIPPED (no differential results) =====\n")
  } else {
    run_heatmap(object2, sample_info, app3, diff_results, cfg)
  }

  # ── Step 9: Boxplot ──────────────────────────────────────────────────
  if (qualitative_only || length(diff_results) == 0) {
    cat("===== Step 9: Boxplot SKIPPED (no differential results) =====\n")
  } else {
    run_boxplot(object2, sample_info, app3, diff_results, cfg)
  }

  # ── Step 10b: Web Export (JSON for website) ──────────────────────────
  web_export_source <- resolve_module_path("export_for_web.R")
  if (file.exists(web_export_source)) {
    source(web_export_source)
    # Collect MetaboAnalystR results if available
    msea_results <- NULL
    topology_results <- NULL
    biomarker_results <- NULL
    chem_class_results <- NULL
    if (exists("msea_results_global")) msea_results <- msea_results_global
    if (exists("topology_results_global")) topology_results <- topology_results_global
    if (exists("biomarker_results_global")) biomarker_results <- biomarker_results_global
    if (exists("chem_class_results_global")) chem_class_results <- chem_class_results_global

    export_all(
      diff_results       = diff_results,
      app3               = app3,
      object2            = object2,
      cfg                = cfg,
      sample_info        = sample_info,
      comparisons        = comparisons,
      msea_results       = msea_results,
      topology_results   = topology_results,
      biomarker_results  = biomarker_results,
      chem_class_results = chem_class_results,
      is_qc_stats        = prep$is_qc_stats,
      qc_quality         = qc_quality,
      qualitative_only   = qualitative_only
    )
  } else {
    cat(sprintf("⚠️  Web export module not found at %s, skipping\n",
                web_export_source))
  }

  # ── Summary report ────────────────────────────────────────────────────
  cat("\n========================================\n")
  cat("     Pipeline complete!\n")
  cat("========================================\n")
  cat(sprintf("Output directory: %s\n", out_root))
  cat(sprintf("Comparisons: %d\n", length(diff_results)))
  cat("Outputs:\n")
  cat("  ├── 03_PCA/                       - PCA scores plot\n")
  cat("  ├── 03.5_OPLSDA/                  - OPLS-DA score/S-plot/VIP/permutation\n")
  cat("  ├── 04_Differential/              - Volcano + differential peak tables\n")
  cat("  ├── 05_Overlap/                   - Venn/UpSet overlap analysis\n")
  cat("  ├── 06_Annotation/                - Metabolite annotation + ID report\n")
  cat("  │   ├── all_annotated_metabolites.csv\n")
  cat("  │   ├── annotations_L[1-3].csv    - Per-confidence-level annotations\n")
  cat("  │   ├── unmatched_features.csv\n")
  cat("  │   ├── identification_summary.txt\n")
  cat("  │   └── Compound_identification_report.xlsx\n")
  cat("  ├── 06.5_OPLSDA_annotated/        - Annotated OPLS-DA S-plots\n")
  cat("  ├── 07_KEGG/                      - Pathway enrichment (HMDB/SMPDB)\n")
  cat("  ├── 10_MetaboAnalystR_Advanced/   - MSEA, topology, biomarker, chem class\n")
  cat("  │   ├── MSEA/                     - Metabolite Set Enrichment (GlobalTest)\n")
  cat("  │   ├── Pathway_Topology/         - KEGG topology (RBC/Impact)\n")
  cat("  │   ├── Biomarker_Analysis/       - RF feature selection + ROC curves\n")
  cat("  │   └── Chemical_Class_Enrichment/- HMDB taxonomy enrichment\n")
  cat("  ├── 08_Heatmap/                   - Differential metabolite heatmaps\n")
  cat("  ├── 09_Boxplot/                   - Boxplots for key metabolites\n")
  cat("  ├── sirius_input/                 - MS2 spectra exported for SIRIUS\n")
  cat("  ├── sirius_project/               - SIRIUS results\n")
  cat("  └── web_export/                   - JSON export for herbMetabo website\n")
  cat("========================================\n")

  # ── Teardown parallel backend ──────────────────────────────────────────
  if (exists("teardown_parallel_backend")) {
    teardown_parallel_backend(cfg)
  }

  invisible(list(
    diff_results = diff_results,
    app3         = app3,
    object2      = object2,
    qc_quality   = qc_quality
  ))
}

# ══════════════════════════════════════════════════════════════════════════════
# 14. Entry Point
# ══════════════════════════════════════════════════════════════════════════════

if (sys.nframe() == 0) {
  opt <- parse_cli()
  cfg <- load_and_validate_config(opt$config)

  if (isTRUE(opt$`validate-only`)) {
    cat("\n===== Validation-only mode: config looks good! =====\n")
    cat("Run without --validate-only to execute the full pipeline.\n")
    quit(save = "no", status = 0)
  }

  main(cfg, qualitative_only = isTRUE(opt$qualitative_only))
}
