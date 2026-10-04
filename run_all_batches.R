#!/usr/bin/env Rscript
###############################################################################
# run_all_batches.R — Master Orchestrator for 6-Batch herbMetabo Pipeline
#
# Runs the full metabolomics pipeline for all 6 batches sequentially,
# then optionally runs cross-batch integration analysis.
#
# Each batch:
#   1. Generates config.yaml from metadata CSV (make_config.R)
#   2. Runs the full pipeline (metabolomics_pipeline.R)
#   3. Generates web_export/ JSON files (export_for_web.R, called inside pipeline)
#
# After all batches complete:
#   4. Runs cross-batch integration analysis (batch correction + merged results)
#
# Usage:
#   Rscript run_all_batches.R                          # Run all 6 batches
#   Rscript run_all_batches.R --batches 1,3,5          # Run specific batches
#   Rscript run_all_batches.R --skip-pipeline           # Skip pipeline, only integration
#   Rscript run_all_batches.R --integration-only        # Only cross-batch integration
#   Rscript run_all_batches.R --dry-run                 # Show what would be done
#   Rscript run_all_batches.R --validate-only           # Validate configs only
#
# Dependencies:
#   - All pipeline dependencies (see metabolomics_pipeline.R)
#   - jsonlite (for integration manifest)
#   - future / furrr (optional, for parallel integration)
###############################################################################

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(jsonlite)
})

# ══════════════════════════════════════════════════════════════════════════════
# 0. Configuration — Batch Directory Mapping
# ══════════════════════════════════════════════════════════════════════════════

#' Batch directory configuration
#'
#' Maps batch_id → raw data directory name, metadata file, and other
#' per-batch settings. Update these paths if directory structure changes.
BATCH_CONFIG <- list(
  `1` = list(
    raw_dir   = "161herbs_raw_data/Batch1(1-16)",
    metadata  = "code/metadata/mzml_file_list_Batch1_convert.csv",
    polarity  = "negative",
    name      = "Batch1"
  ),
  `2` = list(
    raw_dir   = "161herbs_raw_data/Batch2(17-41)",
    metadata  = "code/metadata/mzml_file_list_Batch2_convert.csv",
    polarity  = "negative",
    name      = "Batch2"
  ),
  `3` = list(
    raw_dir   = "161herbs_raw_data/Batch3(42-71)",
    metadata  = "code/metadata/mzml_file_list_Batch3_convert.csv",
    polarity  = "negative",
    name      = "Batch3"
  ),
  `4` = list(
    raw_dir   = "161herbs_raw_data/Batch4(73-101)",
    metadata  = "code/metadata/mzml_file_list_Batch4_convert.csv",
    polarity  = "negative",
    name      = "Batch4"
  ),
  `5` = list(
    raw_dir   = "161herbs_raw_data/Batch5(102-140)",
    metadata  = "code/metadata/mzml_file_list_Batch5_convert.csv",
    polarity  = "negative",
    name      = "Batch5"
  ),
  `6` = list(
    raw_dir   = "161herbs_raw_data/Batch6(5,140-164)",
    metadata  = "code/metadata/mzml_file_list_Batch6_convert.csv",
    polarity  = "negative",
    name      = "Batch6"
  )
)

# Paths relative to the project root
# Detect project root: script is at code/this_project/run_all_batches.R
# Project root is ../../ (two levels up from the script dir)
script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) getwd()
)
if (is.null(script_dir) || script_dir == ".") {
  script_dir <- getwd()
}
# Go up two levels: code/this_project/ -> code/ -> project root
PROJECT_ROOT <- normalizePath(file.path(script_dir, "..", ".."))
# Verify we found the right root
if (!file.exists(file.path(PROJECT_ROOT, "161herbs_raw_data"))) {
  # Try one more level up
  PROJECT_ROOT <- normalizePath(file.path(PROJECT_ROOT, ".."))
}
cat(sprintf("  Project root: %s\n", PROJECT_ROOT))

# Output root for all batch results
OUTPUT_ROOT <- file.path(PROJECT_ROOT, "results")
INTEGRATION_DIR <- file.path(OUTPUT_ROOT, "cross_batch_integration")

# ══════════════════════════════════════════════════════════════════════════════
# 1. CLI Parsing
# ══════════════════════════════════════════════════════════════════════════════

parse_cli <- function() {
  option_list <- list(
    make_option(c("--batches"), type = "character", default = NULL,
                help = "Comma-separated batch IDs to run, e.g. '1,3,5'. Default: all 6"),
    make_option(c("--skip-pipeline"), action = "store_true", default = FALSE,
                help = "Skip the pipeline steps (useful if already run)"),
    make_option(c("--integration-only"), action = "store_true", default = FALSE,
                help = "Only run cross-batch integration, skip all pipelines"),
    make_option(c("--dry-run"), dest = "dry_run", action = "store_true", default = FALSE,
                help = "Print what would be done without executing"),
    make_option(c("--validate-only"), action = "store_true", default = FALSE,
                help = "Validate configs and paths only, do not run"),
    make_option(c("--resume-from"), dest = "resume_from", type = "integer", default = NULL,
                help = "Resume from a specific batch ID (skip earlier batches)",
                metavar = "INT"),
    make_option(c("--parallel"), action = "store_true", default = FALSE,
                help = "Run batches in parallel (use with caution — memory heavy)"),
    make_option(c("--log-dir"), dest = "log_dir", type = "character", default = NULL,
                help = "Directory for log files. Default: <OUTPUT_ROOT>/logs",
                metavar = "DIR"),
    make_option(c("--qualitative-only"), dest = "qualitative_only",
                type = "character", default = NULL,
                help = paste("Comma-separated batch IDs to run in qualitative-only mode",
                             "(annotation + QC diagnostic PCA; no quantitative steps).",
                             "Use for batches whose QC reproducibility makes",
                             "quantification unreliable, e.g. '4'."),
                metavar = "IDS"),
    make_option(c("--oplsda-permutations"), dest = "oplsda_permutations",
                type = "integer", default = NULL,
                help = "OPLS-DA permutation count passed to make_config.R (default 1000; use 100 for a quick pilot batch)",
                metavar = "INT")
  )

  opt <- parse_args(OptionParser(
    option_list = option_list,
    usage = "Rscript run_all_batches.R [options]",
    description = "Master orchestrator for 6-batch herbMetabo metabolomics pipeline."
  ))

  # Parse --batches
  if (!is.null(opt$batches)) {
    opt$batch_ids <- as.integer(strsplit(opt$batches, ",")[[1]])
  } else {
    opt$batch_ids <- 1:6
  }

  # Parse --qualitative-only batch list
  if (!is.null(opt$qualitative_only)) {
    opt$qualitative_only_ids <- as.integer(strsplit(opt$qualitative_only, ",")[[1]])
  } else {
    opt$qualitative_only_ids <- integer(0)
  }

  # Validate batch IDs
  invalid <- setdiff(opt$batch_ids, 1:6)
  if (length(invalid) > 0) {
    stop(sprintf("Invalid batch IDs: %s. Valid range: 1-6",
                 paste(invalid, collapse = ", ")))
  }
  invalid_q <- setdiff(opt$qualitative_only_ids, 1:6)
  if (length(invalid_q) > 0) {
    stop(sprintf("Invalid --qualitative-only batch IDs: %s. Valid range: 1-6",
                 paste(invalid_q, collapse = ", ")))
  }

  opt
}

# ══════════════════════════════════════════════════════════════════════════════
# 2. Helper Functions
# ══════════════════════════════════════════════════════════════════════════════

#' Log a message with timestamp
log_msg <- function(..., level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] [%s] ", timestamp, level), ..., "\n")
}

#' Null-coalescing. This orchestrator runs standalone, so it cannot rely on the
#' definition in metabolomics_pipeline.R.
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a

#' Run a system command with logging
run_cmd <- function(cmd, args, log_file = NULL, dry_run = FALSE) {
  # Ensure dry_run is a logical scalar
  dry_run <- isTRUE(dry_run)

  # Build a proper shell command string with double-quoted args
  # Using double quotes (") handles special chars like () in paths
  quoted_args <- paste(sprintf('"%s"', gsub('"', '\\"', args)), collapse = " ")
  cmd_str <- paste(cmd, quoted_args)

  if (dry_run) {
    log_msg("[DRY RUN] Would execute:", cmd_str)
    return(TRUE)
  }

  log_msg("Running:", cmd_str)

  if (!is.null(log_file)) {
    dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
    exit_code <- system(paste(cmd_str, ">", log_file, "2>&1"))
  } else {
    exit_code <- system(cmd_str)
  }

  if (exit_code != 0) {
    log_msg(sprintf("Command failed with exit code %d: %s", exit_code, cmd_str),
            level = "ERROR")
    if (!is.null(log_file)) {
      log_msg(sprintf("See log: %s", log_file), level = "ERROR")
    }
    return(FALSE)
  }

  log_msg("Completed successfully")
  TRUE
}

#' Resolve the web_export directory for a batch output root.
#'
#' The pipeline writes <output_dir>/detail/web_export (make_config sets
#' output_dir = <finalOutput>/detail). Older/hand-made layouts may put
#' web_export directly under the batch root, so accept both.
resolve_web_export_dir <- function(output_dir) {
  candidates <- c(file.path(output_dir, "detail", "web_export"),
                  file.path(output_dir, "web_export"))
  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "manifest.json"))) return(candidate)
  }
  candidates[1]
}

#' Check if a batch has already been processed (has web_export/)
check_batch_done <- function(batch_id, output_dir) {
  web_dir <- resolve_web_export_dir(output_dir)
  manifest_path <- file.path(web_dir, "manifest.json")
  done <- file.exists(manifest_path)
  if (done) {
    n_json <- length(list.files(web_dir, pattern = "\\.json$"))
    log_msg(sprintf("Batch %d already processed: %s (%d JSON files)",
                    batch_id, web_dir, n_json))
  }
  done
}

#' Get the output directory for a batch
get_batch_output_dir <- function(batch_id, cfg) {
  # The pipeline's make_config.R generates output_dir based on project name
  # Default pattern: <metadata_dir>/Result_<project_name>
  # We override to put all results under OUTPUT_ROOT
  batch_name <- cfg$name
  file.path(OUTPUT_ROOT, batch_name)
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. Run Single Batch
# ══════════════════════════════════════════════════════════════════════════════

#' Run the full pipeline for a single batch
#'
#' @param batch_id Batch ID (1-6)
#' @param cfg Batch configuration from BATCH_CONFIG
#' @param opts CLI options
#' @return TRUE if successful, FALSE otherwise
run_batch <- function(batch_id, cfg, opts) {
  log_msg(sprintf("===== Starting Batch %d: %s =====", batch_id, cfg$name))

  # Resolve absolute paths
  metadata_abs <- file.path(PROJECT_ROOT, cfg$metadata)
  raw_dir_abs  <- file.path(PROJECT_ROOT, cfg$raw_dir)
  output_dir   <- get_batch_output_dir(batch_id, cfg)
  config_file  <- file.path(PROJECT_ROOT, sprintf("config_batch%d.yaml", batch_id))
  log_dir      <- opts$log_dir %||% file.path(OUTPUT_ROOT, "logs")
  log_file     <- file.path(log_dir, sprintf("batch%d.log", batch_id))

  # ── Validate paths ──────────────────────────────────────────────────────
  if (!file.exists(metadata_abs)) {
    log_msg(sprintf("Metadata not found: %s", metadata_abs), level = "ERROR")
    return(FALSE)
  }
  if (!dir.exists(raw_dir_abs)) {
    log_msg(sprintf("Raw data dir not found: %s", raw_dir_abs), level = "ERROR")
    return(FALSE)
  }

  log_msg(sprintf("  Metadata:  %s", metadata_abs))
  log_msg(sprintf("  Raw data:  %s", raw_dir_abs))
  log_msg(sprintf("  Output:    %s", output_dir))
  log_msg(sprintf("  Config:    %s", config_file))
  log_msg(sprintf("  Log:       %s", log_file))

  # ── Step 1: Generate config ─────────────────────────────────────────────
  log_msg("[Step 1] Generating config.yaml...")
  config_args <- c(
    "make_config.R",
    "--metadata", metadata_abs,
    "--raw-dir", raw_dir_abs,
    "--output", config_file,
    "--batch-id", as.character(batch_id),
    "--project-name", cfg$name,
    "--polarity", cfg$polarity,
    "--finalOutput", output_dir,
    # One treatment per drug × dose (10_High_vs_CT1 / 10_Low_vs_CT1) instead of
    # pooling the 3 High + 3 Low replicates into a single 6-sample group.
    "--split-concentration"
  )
  if (!is.null(opts$oplsda_permutations)) {
    config_args <- c(config_args,
                     "--oplsda-permutations", as.character(opts$oplsda_permutations))
  }

  success <- run_cmd("Rscript", config_args,
                     log_file = gsub("\\.log$", "_config.log", log_file),
                     dry_run = isTRUE(opts$dry_run))
  if (!success) {
    log_msg(sprintf("Batch %d: config generation failed", batch_id), level = "ERROR")
    return(FALSE)
  }

  # ── Step 2: Validate config ─────────────────────────────────────────────
  log_msg("[Step 2] Validating config...")
  validate_args <- c(
    "metabolomics_pipeline.R",
    "--config", config_file,
    "--validate-only"
  )
  success <- run_cmd("Rscript", validate_args,
                     log_file = gsub("\\.log$", "_validate.log", log_file),
                     dry_run = isTRUE(opts$dry_run))
  if (!success) {
    log_msg(sprintf("Batch %d: config validation failed", batch_id), level = "ERROR")
    return(FALSE)
  }

  if (isTRUE(opts$`validate-only`)) {
    log_msg(sprintf("Batch %d: validation passed", batch_id))
    return(TRUE)
  }

  # ── Step 3: Run pipeline ────────────────────────────────────────────────
  log_msg("[Step 3] Running metabolomics pipeline...")
  pipeline_args <- c(
    "metabolomics_pipeline.R",
    "--config", config_file
  )
  # Batches listed in --qualitative-only skip their intensity-dependent steps.
  # Their manifest then records quantification_valid = false, which excludes
  # them from cross-batch integration automatically.
  if (batch_id %in% (opts$qualitative_only_ids %||% integer(0))) {
    pipeline_args <- c(pipeline_args, "--qualitative-only")
    log_msg(sprintf("  Batch %d: qualitative-only mode (no quantitative steps)", batch_id))
  }
  success <- run_cmd("Rscript", pipeline_args,
                     log_file = log_file,
                     dry_run = isTRUE(opts$dry_run))
  if (!success) {
    log_msg(sprintf("Batch %d: pipeline failed", batch_id), level = "ERROR")
    return(FALSE)
  }

  # ── Step 4: Verify web export ───────────────────────────────────────────
  log_msg("[Step 4] Verifying web export...")
  web_dir <- resolve_web_export_dir(output_dir)
  manifest_path <- file.path(web_dir, "manifest.json")
  if (file.exists(manifest_path)) {
    n_json <- length(list.files(web_dir, pattern = "\\.json$"))
    log_msg(sprintf("Batch %d: web export OK (%s, %d JSON files)",
                    batch_id, web_dir, n_json))
  } else {
    log_msg(sprintf("Batch %d: web_export/manifest.json not found under %s!",
                    batch_id, output_dir), level = "WARNING")
    log_msg("  The pipeline may have failed silently or export_for_web.R is missing.")
    log_msg("  Check the batch log for details.")
    return(FALSE)
  }

  log_msg(sprintf("===== Batch %d complete =====", batch_id))
  TRUE
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. Cross-Batch Integration
# ══════════════════════════════════════════════════════════════════════════════

#' Run cross-batch integration analysis
#'
#' After all 6 batches complete, this function:
#'   1. Collects all web_export/manifest.json files
#'   2. Builds a cross-batch integration manifest
#'   3. Calls the R integration script (if available)
#'   4. Generates web_export JSON for the integrated results
#'
#' @param opts CLI options
#' @return TRUE if successful
run_integration <- function(opts) {
  log_msg("\n==========================================")
  log_msg("  Cross-Batch Integration Analysis")
  log_msg("==========================================")

  # ── Collect completed batches ───────────────────────────────────────────
  completed_batches <- list()
  excluded_batches  <- list()
  for (bid in opts$batch_ids) {
    cfg <- BATCH_CONFIG[[as.character(bid)]]
    output_dir <- get_batch_output_dir(bid, cfg)
    web_dir <- resolve_web_export_dir(output_dir)
    manifest_path <- file.path(web_dir, "manifest.json")

    if (file.exists(manifest_path)) {
      manifest <- tryCatch({
        read_yaml(manifest_path)
      }, error = function(e) NULL)

      if (!is.null(manifest)) {
        # Mirror the quantitative gate applied inside cross_batch_integration.R.
        # Without this the integration manifest would advertise every batch that
        # merely produced a manifest, including ones excluded from the maths —
        # and the website ingests this file.
        qq <- manifest$qc_quality
        if (!is.null(qq) && !isTRUE(qq$quantification_valid)) {
          reason <- qq$exclude_reason %||% "quantification_valid = false"
          excluded_batches[[as.character(bid)]] <- list(
            batch_id = bid,
            name = cfg$name,
            web_dir = web_dir,
            reason = as.character(reason),
            qc_rsd_median = qq$qc_rsd_median %||% NA_real_,
            n_features = qq$n_features %||% NA_integer_
          )
          log_msg(sprintf("  Batch %d (%s): EXCLUDED from quantitative integration — %s",
                          bid, cfg$name, reason), level = "WARNING")
          next
        }

        completed_batches[[as.character(bid)]] <- list(
          batch_id = bid,
          name = cfg$name,
          output_dir = output_dir,
          web_dir = web_dir,
          manifest = manifest,
          n_treatments = length(manifest$treatments %||% list()),
          n_metabolites = 0  # will be populated
        )
        log_msg(sprintf("  Batch %d (%s): %d treatments",
                        bid, cfg$name, completed_batches[[as.character(bid)]]$n_treatments))
      }
    } else {
      log_msg(sprintf("  Batch %d: not found (no web_export), skipping", bid),
              level = "WARNING")
    }
  }

  if (length(completed_batches) < 2) {
    log_msg("Need at least 2 completed batches for integration, skipping",
            level = "WARNING")
    return(FALSE)
  }

  log_msg(sprintf("  Total batches for integration: %d", length(completed_batches)))

  # ── Create integration output directory ─────────────────────────────────
  dir.create(INTEGRATION_DIR, showWarnings = FALSE, recursive = TRUE)
  web_export_dir <- file.path(INTEGRATION_DIR, "web_export")
  dir.create(web_export_dir, showWarnings = FALSE, recursive = TRUE)

  # ── Build integration manifest ──────────────────────────────────────────
  integration_manifest <- list(
    integration = list(
      name = "Cross-Batch Integration",
      date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      n_batches = length(completed_batches),
      batches = names(completed_batches),
      methods_available = list(
        combat = TRUE,    # ComBat batch correction
        limma = TRUE,     # Limma with batch as covariate
        meta_analysis = TRUE  # Fisher/Stouffer meta-analysis
      )
    ),
    included_batches = lapply(completed_batches, function(b) list(
      batch_id = b$batch_id,
      name = b$name,
      n_treatments = b$n_treatments,
      web_dir = b$web_dir
    )),
    # Batches that produced annotation but failed the quantitative gate. Listed
    # explicitly rather than omitted, so a downstream reader can tell "excluded
    # for cause" apart from "never ran".
    excluded_batches = if (length(excluded_batches) == 0) list() else
      lapply(excluded_batches, function(b) list(
        batch_id = b$batch_id,
        name = b$name,
        web_dir = b$web_dir,
        reason = b$reason,
        qc_rsd_median = b$qc_rsd_median,
        n_features = b$n_features
      ))
  )

  write_json(integration_manifest,
             file.path(web_export_dir, "integration_manifest.json"),
             pretty = TRUE, auto_unbox = TRUE)
  log_msg(sprintf("  Integration manifest: %s",
                  file.path(web_export_dir, "integration_manifest.json")))

  # ── Check if there's an integration R script to run ─────────────────────
  integration_script <- file.path(
    script_dir, "cross_batch_integration.R"
  )

  if (file.exists(integration_script)) {
    log_msg("Running cross_batch_integration.R...")
    if (!opts$dry_run) {
      integration_args <- c(
        integration_script,
        "--input-dir", web_export_dir,
        "--output-dir", INTEGRATION_DIR,
        "--batches", paste(names(completed_batches), collapse = ",")
      )
      success <- run_cmd("Rscript", integration_args,
                         log_file = file.path(INTEGRATION_DIR, "integration.log"),
                         dry_run = isTRUE(opts$dry_run))
      if (success) {
        log_msg("Cross-batch integration complete")
      } else {
        log_msg("Cross-batch integration script failed", level = "ERROR")
        log_msg("  You can run it manually later.")
      }
    }
  } else {
    log_msg("Integration script not found at:", integration_script)
    log_msg("  Creating placeholder. Implement cross_batch_integration.R later.")
    log_msg("  The integration manifest is ready for manual processing.")

    # Create a placeholder for the integration script
    if (!opts$dry_run) {
      integration_placeholder <- file.path(
        dirname(integration_script), "cross_batch_integration.R"
      )
      if (!file.exists(integration_placeholder)) {
        writeLines(
          c(
            "#!/usr/bin/env Rscript",
            "# Cross-batch integration analysis — to be implemented",
            "#",
            "# This script should:",
            "#   1. Read all batch web_export/ files",
            "#   2. Perform batch effect correction (ComBat, limma)",
            "#   3. Merge metabolite-level results across batches",
            "#   4. Generate integrated web_export JSON files",
            "#",
            "# Usage:",
            "#   Rscript cross_batch_integration.R \\",
            "#     --input-dir <web_export_dir> \\",
            "#     --output-dir <integration_dir> \\",
            "#     --batches 1,2,3,4,5,6",
            "",
            "suppressPackageStartupMessages(library(optparse))",
            "suppressPackageStartupMessages(library(jsonlite))",
            "",
            "option_list <- list(",
            "  make_option('--input-dir', type='character', help='Input web_export directory'),",
            "  make_option('--output-dir', type='character', help='Output directory'),",
            "  make_option('--batches', type='character', help='Comma-separated batch IDs')",
            ")",
            "opt <- parse_args(OptionParser(option_list = option_list))",
            "",
            "cat('Cross-batch integration analysis placeholder\\n')",
            "cat(sprintf('  Input: %s\\n', opt$`input-dir`))",
            "cat(sprintf('  Output: %s\\n', opt$`output-dir`))",
            "cat(sprintf('  Batches: %s\\n', opt$batches))",
            "cat('\\nTODO: Implement batch correction and integration\\n')",
            "",
            "# Save a placeholder result",
            "write_json(",
            "  list(status = 'placeholder', message = 'Integration not yet implemented'),",
            "  file.path(opt$`output-dir`, 'web_export', 'integration_results.json'),",
            "  pretty = TRUE, auto_unbox = TRUE",
            ")",
            ""
          ),
          integration_placeholder
        )
        log_msg(sprintf("  Created placeholder: %s", integration_placeholder))
      }
    }
  }

  log_msg("Cross-batch integration phase complete")
  TRUE
}

# ══════════════════════════════════════════════════════════════════════════════
# 5. Main Entry Point
# ══════════════════════════════════════════════════════════════════════════════

main <- function() {
  opts <- parse_cli()

  cat("\n")
  cat("========================================\n")
  cat("  herbMetabo — 6-Batch Pipeline Orchestrator\n")
  cat("========================================\n")
  cat(sprintf("  Project root:   %s\n", PROJECT_ROOT))
  cat(sprintf("  Output root:    %s\n", OUTPUT_ROOT))
  cat(sprintf("  Batches:        %s\n", paste(opts$batch_ids, collapse = ", ")))
  cat(sprintf("  Dry run:        %s\n", opts$dry_run))
  cat(sprintf("  Validate only:  %s\n", opts$`validate-only`))
  cat(sprintf("  Parallel:       %s\n", opts$parallel))
  cat("========================================\n\n")

  # Create output directories
  dir.create(OUTPUT_ROOT, showWarnings = FALSE, recursive = TRUE)
  log_dir <- opts$log_dir %||% file.path(OUTPUT_ROOT, "logs")
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

  # ── Phase 1: Run pipelines ─────────────────────────────────────────────
  if (isFALSE(opts$`integration-only`)) {
    results <- list()

    if (opts$parallel && length(opts$batch_ids) > 1) {
      # ── Parallel execution ──────────────────────────────────────────────
      log_msg("Running batches in parallel (this may use a lot of memory)...")

      if (!requireNamespace("future", quietly = TRUE) ||
          !requireNamespace("furrr", quietly = TRUE)) {
        log_msg("future/furrr packages not available. Falling back to sequential.",
                level = "WARNING")
        opts$parallel <- FALSE
      } else {
        library(future)
        library(furrr)

        # Use multisession plan (separate R processes)
        n_cores <- min(length(opts$batch_ids), future::availableCores() - 1)
        plan(multisession, workers = n_cores)
        log_msg(sprintf("  Using %d parallel workers", n_cores))

        # Resume support
        batch_ids_to_run <- opts$batch_ids
        if (!is.null(opts$resume_from)) {
          batch_ids_to_run <- batch_ids_to_run[batch_ids_to_run >= opts$resume_from]
        }

        results <- furrr::future_map(setNames(batch_ids_to_run, batch_ids_to_run), function(bid) {
          cfg <- BATCH_CONFIG[[as.character(bid)]]
          run_batch(bid, cfg, opts)
        }, .options = furrr_options(seed = TRUE))

        plan(sequential)
      }
    }

    if (!opts$parallel) {
      # ── Sequential execution (default) ──────────────────────────────────
      for (bid in opts$batch_ids) {
        # Resume support: skip batches before resume_from
        if (!is.null(opts$resume_from) && bid < opts$resume_from) {
          log_msg(sprintf("Skipping Batch %d (resume from %d)", bid, opts$resume_from))
          next
        }

        # Skip if already processed (resume mode)
        if (!is.null(opts$resume_from)) {
          cfg <- BATCH_CONFIG[[as.character(bid)]]
          output_dir <- get_batch_output_dir(bid, cfg)
          if (check_batch_done(bid, output_dir)) {
            log_msg(sprintf("Batch %d already done, skipping", bid))
            next
          }
        }

        cfg <- BATCH_CONFIG[[as.character(bid)]]
        results[[as.character(bid)]] <- run_batch(bid, cfg, opts)
      }
    }

    # ── Report results ────────────────────────────────────────────────────
    cat("\n========================================\n")
    cat("  Pipeline Results Summary\n")
    cat("========================================\n")
    n_success <- sum(unlist(results), na.rm = TRUE)
    n_total <- length(results)
    for (bid in names(results)) {
      status <- if (isTRUE(results[[bid]])) "✓" else "✗"
      cfg <- BATCH_CONFIG[[bid]]
      cat(sprintf("  Batch %s (%s): %s\n", bid, cfg$name, status))
    }
    cat(sprintf("\n  %d / %d batches succeeded\n", n_success, n_total))
    cat("========================================\n")

  } else {
    log_msg("Skipping pipeline (--integration-only mode)")
  }

  # ── Phase 2: Cross-batch integration ───────────────────────────────────
  if (isTRUE(opts$`integration-only`) && isFALSE(opts$dry_run) && isFALSE(opts$`validate-only`)) {
    run_integration(opts)
  } else if (isFALSE(opts$`validate-only`) && isFALSE(opts$`skip-pipeline`) && isFALSE(opts$dry_run)) {
    run_integration(opts)
  } else if (isTRUE(opts$`validate-only`)) {
    log_msg("Skipping integration (--validate-only mode)")
  } else if (isTRUE(opts$`skip-pipeline`)) {
    # If only integration was requested, still run it
    if (isTRUE(opts$`integration-only`)) {
      run_integration(opts)
    } else {
      log_msg("Skipping integration (--skip-pipeline mode)")
    }
  }

  # ── Final summary ──────────────────────────────────────────────────────
  cat("\n========================================\n")
  cat("  All Done!\n")
  cat("========================================\n")
  cat(sprintf("  Results:     %s\n", OUTPUT_ROOT))
  cat(sprintf("  Logs:        %s\n", log_dir))
  cat(sprintf("  Integration: %s\n", INTEGRATION_DIR))
  cat("\n  Next steps:\n")
  cat("  1. Start PostgreSQL + Django:\n")
  cat("      cd website && docker compose up -d db\n")
  cat("      cd website/backend && python manage.py migrate\n")
  cat("      python manage.py seed_metadata\n")
  cat("  2. Import data into Django:\n")
  cat(sprintf("      python manage.py ingest_r_data --base-dir %s\n", OUTPUT_ROOT))
  cat("  3. Start the API server:\n")
  cat("      python manage.py runserver\n")
  cat("========================================\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# 6. Entry Point
# ══════════════════════════════════════════════════════════════════════════════

if (sys.nframe() == 0) {
  main()
}