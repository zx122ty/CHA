#!/usr/bin/env Rscript
###############################################################################
# performance_utils.R — High-Performance Optimizations for Metabolomics Pipeline
#
# Provides:
#   1. Parallel backend management (future / future.apply / furrr)
#   2. High-performance I/O wrappers (data.table::fread/fwrite, qs serialization)
#   3. File-hash caching system (MD5/SHA256 keyed cache)
#   4. SQLite cache for network/API calls (MetaboAnalystR KEGG)
#   5. Memory management utilities
#   6. Refactored bottleneck functions:
#      - run_metfrag_analysis_parallel (parallelized over features)
#      - run_boxplot_parallel (parallelized over metabolites)
#      - run_annotated_volcano_parallel (parallelized over comparisons)
#      - run_differential_parallel (parallelized over comparisons)
#      - run_sirius_analysis_cached (with hash-based caching)
#      - run_preprocessing_fast (with matrix ops + fread/fwrite)
#
# All functions fall back to sequential execution if parallelization fails.
#
# Usage: source("code/performance_utils.R")  # sourced by main()
# Dependencies: future, future.apply, furrr, data.table, qs, RSQLite (optional)
###############################################################################

suppressPackageStartupMessages({
  library(future)
  library(future.apply)
  library(furrr)
  library(data.table)
})

#' List-safe null-coalescing operator
#'
#' Unlike the main script's `%||%`, this version handles lists correctly:
#' `is.na()` on a list returns a vector; we only call `is.na()` on scalars.
`%||%` <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.list(a) || length(a) > 1) return(a)
  if (is.na(a)) return(b)
  a
}

# ══════════════════════════════════════════════════════════════════════════════
# 1. Parallel Backend Management
# ══════════════════════════════════════════════════════════════════════════════

#' Initialize parallel backend based on configuration
#'
#' Sets up the future plan (multisession, multicore, or sequential).
#' Logs the number of workers and backend type. Falls back to sequential
#' if the requested backend fails.
#'
#' @param cfg Config list (performance sub-block)
#' @return Invisible TRUE/FALSE indicating whether parallel was enabled
init_parallel_backend <- function(cfg) {
  perf <- cfg$performance %||% list()
  n_cores <- perf$cores %||% 4
  if (n_cores == 0) n_cores <- max(1, parallel::detectCores() - 1)
  backend <- perf$parallel_backend %||% "multisession"
  timeout <- perf$parallel_setup_timeout %||% 30

  cat("\n===== Parallel Backend Initialization =====\n")
  cat(sprintf("   Requested: %s with %d workers\n", backend, n_cores))

  if (identical(backend, "sequential")) {
    plan("sequential")
    cat("   Sequential mode (config override)\n")
    return(invisible(FALSE))
  }

  # Set timeout for future resolution
  options(future.wait.timeout = timeout)

  parallel_ok <- tryCatch({
    switch(backend,
      multisession = {
        plan("multisession", workers = n_cores)
        TRUE
      },
      multicore = {
        plan("multicore", workers = n_cores)
        TRUE
      },
      cluster = {
        cl <- parallel::makeCluster(n_cores, timeout = timeout)
        plan("cluster", workers = cl)
        TRUE
      },
      {
        plan("multisession", workers = n_cores)
        TRUE
      }
    )
  }, error = function(e) {
    cat(sprintf("   WARNING: Parallel backend '%s' failed: %s\n", backend, conditionMessage(e)))
    cat("   Falling back to sequential execution\n")
    plan("sequential")
    FALSE
  })

  if (parallel_ok) {
    actual_workers <- nbrOfWorkers()
    cat(sprintf("   Active backend: %s with %d worker(s)\n",
                class(plan())[1], actual_workers))
  } else {
    plan("sequential")
    cat("   Sequential fallback active\n")
  }

  invisible(parallel_ok)
}

#' Shut down parallel backend and clean up workers
#'
#' @param cfg Optional config list (for logging level)
teardown_parallel_backend <- function(cfg = NULL) {
  cat("\n===== Tearing down parallel backend =====\n")
  tryCatch({
    future::plan("sequential")
    cat("   Parallel backend shut down, reset to sequential\n")
  }, error = function(e) {
    cat(sprintf("   Cleanup warning: %s\n", conditionMessage(e)))
    # Force reset
    tryCatch(future::plan("sequential"), error = function(e2) NULL)
  })
  invisible(TRUE)
}


# ══════════════════════════════════════════════════════════════════════════════
# 2. High-Performance I/O Wrappers
# ══════════════════════════════════════════════════════════════════════════════

#' Fast read CSV — uses fread or falls back to read.csv
#'
#' @param file Path to CSV file
#' @param ... Additional arguments passed to fread or read.csv
#' @return data.frame/data.table
fast_read <- function(file, ...) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    as.data.frame(data.table::fread(file, ...))
  } else {
    read.csv(file, stringsAsFactors = FALSE, ...)
  }
}

#' Fast write CSV — uses fwrite or falls back to write.csv
#'
#' @param x Data to write
#' @param file Output path
#' @param ... Additional arguments passed to fwrite or write.csv
fast_write <- function(x, file, ...) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(x, file, ...)
  } else {
    write.csv(x, file, row.names = FALSE, ...)
  }
  invisible(file)
}

#' Fast serialization of intermediate R objects
#'
#' Uses qs (fastest), fst (fast), or base saveRDS (fallback).
#'
#' @param object R object to save
#' @param file Path to output file (extension .qs, .fst, or .rds)
#' @param cfg Config list (performance sub-block for preset)
fast_save_object <- function(object, file, cfg = NULL) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  engine <- if (!is.null(cfg)) cfg$performance$serialize_engine %||% "qs" else "qs"

  saved <- switch(engine,
    qs = {
      if (!requireNamespace("qs", quietly = TRUE)) {
        cat("   qs not available, falling back to saveRDS\n")
        saveRDS(object, file)
        TRUE
      } else {
        preset <- if (!is.null(cfg)) cfg$performance$qs_preset %||% "high" else "high"
        qs::qsave(object, file, preset = preset)
        TRUE
      }
    },
    fst = {
      if (!requireNamespace("fst", quietly = TRUE)) {
        cat("   fst not available, falling back to saveRDS\n")
        saveRDS(object, file)
        TRUE
      } else {
        fst::write_fst(object, file)
        TRUE
      }
    },
    {
      saveRDS(object, file)
      TRUE
    }
  )
  cat(sprintf("   Saved: %s (%s engine)\n", file, engine))
  invisible(saved)
}

#' Fast deserialization of intermediate R objects
#'
#' @param file Path to input file
#' @param cfg Config list (performance sub-block)
#' @return Deserialized R object
fast_load_object <- function(file, cfg = NULL) {
  if (!file.exists(file)) {
    stop(sprintf("Cache file not found: %s", file))
  }
  engine <- if (!is.null(cfg)) cfg$performance$serialize_engine %||% "qs" else "qs"

  result <- switch(engine,
    qs = {
      if (!requireNamespace("qs", quietly = TRUE)) {
        readRDS(file)
      } else {
        qs::qread(file)
      }
    },
    fst = {
      if (!requireNamespace("fst", quietly = TRUE)) {
        readRDS(file)
      } else {
        fst::read_fst(file)
      }
    },
    readRDS(file)
  )
  cat(sprintf("   Loaded: %s (%s engine)\n", file, engine))
  result
}


# ══════════════════════════════════════════════════════════════════════════════
# 3. File-Hash Caching System
# ══════════════════════════════════════════════════════════════════════════════

#' Compute a hash for a file or a set of parameters
#'
#' @param input Path to file, or a list/character vector of parameters
#' @param algo Hash algorithm ("md5" or "sha256")
#' @return Character hash string
compute_hash <- function(input, algo = "md5") {
  if (is.character(input) && length(input) == 1 && file.exists(input)) {
    # Hash a file
    if (algo == "sha256") {
      digest::digest(input, algo = "sha256", file = TRUE)
    } else {
      digest::digest(input, algo = "md5", file = TRUE)
    }
  } else {
    # Hash parameters
    serialized <- serialize(input, NULL)
    if (algo == "sha256") {
      digest::digest(serialized, algo = "sha256")
    } else {
      digest::digest(serialized, algo = "md5")
    }
  }
}

#' Compute combined hash of multiple inputs (file + params)
#'
#' Useful for SIRIUS: hash(mgf_file) + hash(parameters) as the cache key.
#'
#' @param ... Character vectors or file paths to hash together
#' @param algo Hash algorithm
#' @return Combined hash string
compute_combined_hash <- function(..., algo = "md5") {
  inputs <- list(...)
  hashes <- sapply(inputs, function(x) compute_hash(x, algo))
  combined <- paste(hashes, collapse = "_")
  if (algo == "sha256") {
    digest::digest(combined, algo = "sha256")
  } else {
    digest::digest(combined, algo = "md5")
  }
}

#' Check if a cached result exists and is valid
#'
#' @param cache_key Unique cache key (hash string)
#' @param cache_dir Cache directory path
#' @param max_age_days Maximum age of cache entry in days (NULL = no expiry)
#' @param cfg Optional config list
#' @return Path to cached file if valid, NULL otherwise
cache_lookup <- function(cache_key, cache_dir, max_age_days = NULL, cfg = NULL) {
  if (!isTRUE(cfg$performance$enable_caching %||% TRUE)) return(NULL)
  cache_file <- file.path(cache_dir, paste0(cache_key, ".qs"))

  if (!file.exists(cache_file)) {
    return(NULL)
  }

  # Check age
  if (!is.null(max_age_days)) {
    age_days <- as.numeric(difftime(Sys.time(), file.mtime(cache_file), units = "days"))
    if (age_days > max_age_days) {
      cat(sprintf("   Cache expired: %s (%.1f days > %d days)\n",
                  cache_key, age_days, max_age_days))
      unlink(cache_file)
      return(NULL)
    }
  }

  cat(sprintf("   Cache HIT: %s\n", cache_key))
  cache_file
}

#' Store a result in the cache
#'
#' @param cache_key Unique cache key
#' @param object R object to cache
#' @param cache_dir Cache directory
#' @param cfg Optional config list (for serialization settings)
#' @return Path to cached file
cache_store <- function(cache_key, object, cache_dir, cfg = NULL) {
  if (!isTRUE(cfg$performance$enable_caching %||% TRUE)) return(NULL)
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  cache_file <- file.path(cache_dir, paste0(cache_key, ".qs"))
  qs::qsave(object, cache_file, preset = "fast")
  cat(sprintf("   Cache STORE: %s\n", cache_key))
  cache_file
}

#' Evict expired cache entries
#'
#' @param cache_dir Cache directory
#' @param max_age_days Maximum age in days
#' @param pattern File pattern to match (default: "*.qs")
cache_evict_expired <- function(cache_dir, max_age_days = 30, pattern = "\\.qs$") {
  if (!dir.exists(cache_dir)) return(invisible(0))
  cache_files <- list.files(cache_dir, pattern = pattern, full.names = TRUE)
  if (length(cache_files) == 0) return(invisible(0))

  now <- Sys.time()
  n_evicted <- 0
  for (f in cache_files) {
    age_days <- as.numeric(difftime(now, file.mtime(f), units = "days"))
    if (age_days > max_age_days) {
      unlink(f)
      n_evicted <- n_evicted + 1
    }
  }
  if (n_evicted > 0) {
    cat(sprintf("   Cache eviction: removed %d expired entries from %s\n",
                n_evicted, cache_dir))
  }
  invisible(n_evicted)
}


# ══════════════════════════════════════════════════════════════════════════════
# 4. SQLite Cache for Network/API Calls
# ══════════════════════════════════════════════════════════════════════════════

#' Initialize SQLite cache database
#'
#' Creates tables for caching KEGG/MetaboAnalystR API responses.
#' Both the request hash and the response are stored as text.
#'
#' @param db_path Path to SQLite database file
#' @return Database connection object (or NULL if RSQLite unavailable)
init_sqlite_cache <- function(db_path) {
  if (!requireNamespace("RSQLite", quietly = TRUE)) {
    cat("   RSQLite package not available, SQLite cache disabled\n")
    cat("   Install with: install.packages('RSQLite')\n")
    return(NULL)
  }
  if (!requireNamespace("DBI", quietly = TRUE)) {
    cat("   DBI package not available, SQLite cache disabled\n")
    return(NULL)
  }

  dir.create(dirname(db_path), showWarnings = FALSE, recursive = TRUE)
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)

  # Create tables if they don't exist
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS api_cache (
      cache_key TEXT PRIMARY KEY,
      response TEXT,
      created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
      expires_at TIMESTAMP
    )
  ")
  DBI::dbExecute(con, "
    CREATE INDEX IF NOT EXISTS idx_api_cache_expires
    ON api_cache(expires_at)
  ")

  # Clean expired entries
  DBI::dbExecute(con, "DELETE FROM api_cache WHERE expires_at < datetime('now')")

  cat(sprintf("   SQLite cache initialized: %s\n", db_path))
  con
}

#' Look up a cached API response
#'
#' @param con SQLite database connection
#' @param cache_key Request hash key
#' @return Cached response string, or NULL if not found
sqlite_cache_lookup <- function(con, cache_key) {
  if (is.null(con)) return(NULL)
  result <- DBI::dbGetQuery(con,
    "SELECT response FROM api_cache WHERE cache_key = ? AND
     (expires_at IS NULL OR expires_at > datetime('now'))",
    params = list(cache_key))
  if (nrow(result) == 0) return(NULL)
  cat(sprintf("   SQLite cache HIT: %s\n", substr(cache_key, 1, 16)))
  result$response[1]
}

#' Store an API response in the SQLite cache
#'
#' @param con SQLite database connection
#' @param cache_key Request hash key
#' @param response Response string to cache
#' @param ttl_seconds Time-to-live in seconds (default: 86400 = 1 day)
sqlite_cache_store <- function(con, cache_key, response, ttl_seconds = 86400) {
  if (is.null(con)) return(invisible(NULL))
  DBI::dbExecute(con,
    "INSERT OR REPLACE INTO api_cache (cache_key, response, created_at, expires_at)
     VALUES (?, ?, datetime('now'), datetime('now', ? || ' seconds'))",
    params = list(cache_key, response, as.character(ttl_seconds)))
  cat(sprintf("   SQLite cache STORE: %s\n", substr(cache_key, 1, 16)))
  invisible(TRUE)
}

#' Close SQLite cache connection
#'
#' @param con SQLite database connection
close_sqlite_cache <- function(con) {
  if (!is.null(con) && DBI::dbIsValid(con)) {
    DBI::dbDisconnect(con)
    cat("   SQLite cache connection closed\n")
  }
  invisible(TRUE)
}


# ══════════════════════════════════════════════════════════════════════════════
# 5. Memory Management Utilities
# ══════════════════════════════════════════════════════════════════════════════

#' Perform garbage collection with informative logging
#'
#' @param label Optional label describing what memory was just freed
#' @param force Force gc() even if aggressive_gc is off
perform_gc <- function(label = NULL, force = FALSE, cfg = NULL) {
  gc_freq <- if (!is.null(cfg)) cfg$performance$gc_frequency %||% "after_each_step" else "after_each_step"
  if (gc_freq == "never" && !force) return(invisible(FALSE))

  before <- gc(reset = TRUE)
  mem_used <- sum(before[, 1] * before[, 2]) / 1024^2  # MB

  if (!is.null(label)) {
    cat(sprintf("   [GC] %s: %.1f MB collected\n", label, mem_used))
  }
  invisible(mem_used)
}

#' Clear large objects from memory
#'
#' Removes and garbage-collects specified objects from the calling environment.
#'
#' @param ... Object names (as character strings) to remove
#' @param env Environment to clean (default: parent frame)
clear_objects <- function(..., env = parent.frame()) {
  obj_names <- as.character(unlist(list(...)))
  for (nm in obj_names) {
    if (exists(nm, envir = env)) {
      obj_size <- tryCatch(
        utils::object.size(get(nm, envir = env)) / 1024^2,
        error = function(e) NA
      )
      rm(list = nm, envir = env)
      if (!is.na(obj_size)) {
        cat(sprintf("   Freed: '%s' (%.1f MB)\n", nm, obj_size))
      }
    }
  }
  gc()
  invisible(TRUE)
}

#' Estimate memory footprint of MS2 spectra data
#'
#' @param object mass_dataset object
#' @param cfg Config list (for threshold)
#' @return Total size in MB, invisibly
check_ms2_memory <- function(object, cfg = NULL) {
  ms2 <- object@ms2_data
  if (is.null(ms2) || length(ms2) == 0) {
    cat("   No MS2 spectra loaded\n")
    return(invisible(0))
  }

  total_size <- 0
  n_spectra <- 0
  for (md in ms2) {
    if (isS4(md) && inherits(md, "ms2_data")) {
      for (spec in md@ms2_spectra) {
        if (!is.null(spec)) {
          total_size <- total_size + utils::object.size(spec)
          n_spectra <- n_spectra + 1
        }
      }
    }
  }
  total_size_mb <- total_size / 1024^2
  threshold <- if (!is.null(cfg)) cfg$performance$max_ms2_spectra_mb %||% 500 else 500

  cat(sprintf("   MS2 spectra: %d spectra, %.1f MB\n", n_spectra, total_size_mb))
  if (total_size_mb > threshold) {
    cat(sprintf("   WARNING: MS2 spectra exceed %.0f MB threshold.\n", threshold))
    cat("   Consider reducing the number of MS2 files or increasing\n")
    cat("   max_ms2_spectra_mb in the config.\n")
  }
  invisible(total_size_mb)
}


# ══════════════════════════════════════════════════════════════════════════════
# 6. Matrix-Optimized Numerical Operations
# ══════════════════════════════════════════════════════════════════════════════

#' Compute RSD (Relative Standard Deviation) using pure matrix operations
#'
#' 10-50× faster than apply(data, 1, calc_rsd) for large matrices.
#'
#' @param mat Numeric matrix (features x samples)
#' @param min_values Minimum non-NA values required (default: 3)
#' @return Numeric vector of RSD values (NA for features with insufficient data)
matrix_rsd <- function(mat, min_values = 3) {
  # Ensure matrix
  if (!is.matrix(mat)) mat <- as.matrix(mat)

  # Compute mean and sd using colMeans/colSds for speed
  # Transpose: we want per-row stats, but colMeans is faster
  # Use matrixStats if available
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    means <- matrixStats::rowMeans2(mat, na.rm = TRUE)
    sds <- matrixStats::rowSds(mat, na.rm = TRUE)
    n_valid <- matrixStats::rowCounts(mat, value = NA, na.rm = TRUE)  # count valid
    n_valid <- ncol(mat) - n_valid
  } else {
    # Fallback: still much faster than apply with a custom function
    means <- rowMeans(mat, na.rm = TRUE)
    # Fast row SD via sweep
    sds <- sqrt(rowSums((mat - means)^2, na.rm = TRUE) / (ncol(mat) - 1))
    n_valid <- rowSums(!is.na(mat))
  }

  rsd <- sds / means * 100
  rsd[n_valid < min_values | means == 0 | is.na(means)] <- NA
  rsd
}

#' Fast blank filter using matrix operations
#'
#' @param expression_data Feature matrix (features x samples)
#' @param blank_cols Indices/names of blank sample columns
#' @param bio_cols Indices/names of biological sample columns
#' @param fold_change Minimum fold-change threshold
#' @return Logical vector: TRUE = retain feature
fast_blank_filter <- function(expression_data, blank_cols, bio_cols, fold_change = 3) {
  expr_mat <- as.matrix(expression_data)

  blank_idx <- which(colnames(expr_mat) %in% blank_cols)
  bio_idx <- which(colnames(expr_mat) %in% bio_cols)

  if (length(blank_idx) < 2) {
    return(rep(TRUE, nrow(expr_mat)))
  }

  # Median across blanks (row-wise), mean across biological samples
  blank_med <- matrixStats::rowMedians(expr_mat[, blank_idx, drop = FALSE], na.rm = TRUE)
  bio_mean <- rowMeans(expr_mat[, bio_idx, drop = FALSE], na.rm = TRUE)

  pass <- bio_mean >= fold_change * blank_med
  pass[is.na(pass) | blank_med == 0] <- TRUE
  pass
}

#' Fast missing value filter using matrix operations
#'
#' @param expression_data Feature matrix (features x samples)
#' @param class_vec Named vector of class assignments (names = sample IDs)
#' @param max_missing Maximum allowed fraction of missing values per group
#' @return Logical vector: TRUE = retain feature
fast_missing_filter <- function(expression_data, class_vec, max_missing = 0.8) {
  expr_mat <- as.matrix(expression_data)
  groups <- unique(class_vec)
  n_features <- nrow(expr_mat)

  # For each group, compute per-feature missing fraction
  keep <- rep(TRUE, n_features)
  for (grp in groups) {
    grp_samps <- names(class_vec)[class_vec == grp]
    grp_idx <- which(colnames(expr_mat) %in% grp_samps)
    if (length(grp_idx) == 0) next

    grp_data <- expr_mat[, grp_idx, drop = FALSE]
    missing_frac <- rowSums(is.na(grp_data)) / ncol(grp_data)
    keep <- keep & (missing_frac <= max_missing)
  }
  keep
}


# ══════════════════════════════════════════════════════════════════════════════
# 7. Refactored: Parallel Differential Analysis
# ══════════════════════════════════════════════════════════════════════════════

#' Parallel differential analysis (limma) across multiple comparisons
#'
#' Uses future_lapply to run independent limma models in parallel.
#' If parallel fails, falls back to sequential loop.
#'
#' @param object2 mass_dataset object
#' @param sample_info Sample information data.frame
#' @param comparisons List of comparison objects (from generate_comparisons)
#' @param cfg Validated config list
#' @return List of differential results (same structure as original run_differential)
run_differential_parallel <- function(object2, sample_info, comparisons, cfg) {
  cat("\n===== Step 4: Differential analysis (limma, parallel) =====\n")

  out_dir <- file.path(cfg$project$output_dir, "04_Differential")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  diff_cfg <- cfg$differential
  alpha <- diff_cfg$p_value_cutoff
  fc_cut <- diff_cfg$fc_threshold
  logfc_cutoff <- log2(fc_cut)

  sig_metric <- diff_cfg$significance_metric %||% "p_value"
  p_col <- if (identical(sig_metric, "adj_p_value")) "adj.P.Val" else "P.Value"
  p_lab <- if (identical(sig_metric, "adj_p_value")) "adj. p-value (FDR)" else "p-value"
  cat(sprintf("   Significance: |log2FC| >= %.3f & %s < %.3g\n", logfc_cutoff, p_lab, alpha))

  Nomalize_data <- object2@expression_data
  log_data <- log2(Nomalize_data)
  core_group_vec <- sample_info$class[match(colnames(Nomalize_data), sample_info$sample_id)]
  all_samples <- colnames(Nomalize_data)
  viz_labels <- cfg$visualization$group_labels

  # Pre-filter comparisons that have sufficient samples
  valid_comparisons <- list()
  for (cmp in comparisons) {
    ctrl_samples  <- all_samples[core_group_vec == cmp$ctrl]
    treat_samples <- all_samples[core_group_vec == cmp$treat]
    if (length(ctrl_samples) >= 2 && length(treat_samples) >= 2) {
      valid_comparisons[[length(valid_comparisons) + 1]] <- cmp
    } else {
      cat(sprintf("   ⚠️  Skipping %s: insufficient samples (ctrl=%d, treat=%d)\n",
                  cmp$label %||% cmp$name, length(ctrl_samples), length(treat_samples)))
    }
  }

  cat(sprintf("   Processing %d comparisons in parallel...\n", length(valid_comparisons)))

  # Run one comparison (to be used by future_lapply)
  run_one_comparison <- function(cmp) {
    cmp_name  <- cmp$name
    ctrl_grp  <- cmp$ctrl
    treat_grp <- cmp$treat
    label     <- cmp$label

    ctrl_samples  <- all_samples[core_group_vec == ctrl_grp]
    treat_samples <- all_samples[core_group_vec == treat_grp]

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

    # Volcano plot data
    Diff$.plot_p <- Diff[[p_col]]
    Significant <- ifelse(
      Diff$.plot_p < alpha & abs(Diff$logFC) > logfc_cutoff,
      ifelse(Diff$logFC > logfc_cutoff, "Up", "Down"), "Not"
    )

    list(
      name          = cmp_name,
      label         = label,
      all           = Diff,
      sig           = Diff_sig,
      ctrl_samples  = ctrl_samples,
      treat_samples = treat_samples,
      n_up          = sum(Diff_sig$logFC > 0, na.rm = TRUE),
      n_down        = sum(Diff_sig$logFC < 0, na.rm = TRUE),
      comparison    = cmp,
      volcano_data  = data.frame(logFC = Diff$logFC, plot_p = Diff$.plot_p,
                                 Significant = Significant, stringsAsFactors = FALSE)
    )
  }

  # Execute in parallel with sequential fallback
  results <- tryCatch({
    future_lapply(valid_comparisons, run_one_comparison,
                  future.scheduling = 1.0,
                  future.seed = TRUE)
  }, error = function(e) {
    cat(sprintf("   WARNING: Parallel differential analysis failed: %s\n", conditionMessage(e)))
    cat("   Falling back to sequential\n")
    lapply(valid_comparisons, run_one_comparison)
  })

  # Build named list and write outputs
  diff_results <- list()
  for (res in results) {
    if (is.null(res)) next
    cmp_name <- res$name
    diff_results[[cmp_name]] <- list(
      all           = res$all,
      sig           = res$sig,
      ctrl_samples  = res$ctrl_samples,
      treat_samples = res$treat_samples,
      n_up          = res$n_up,
      n_down        = res$n_down,
      comparison    = res$comparison
    )

    cat(sprintf("   %s: %d significant (up: %d, down: %d)\n",
                res$label, nrow(res$sig), res$n_up, res$n_down))

    # Write CSV
    fast_write(res$sig,
               file.path(out_dir, paste0(cmp_name, "_differential_peaks.csv")),
               row.names = TRUE)

    # Volcano plot
    vd <- res$volcano_data
    p <- ggplot(vd, aes(x = logFC, y = -log10(plot_p))) +
      geom_point(alpha = 0.4, size = 3.5, aes(color = Significant)) +
      ylab(sprintf("-log10(%s)", p_lab)) +
      scale_color_manual(values = c("blue4", "grey", "red3")) +
      geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff),
                 lty = 4, col = "black", lwd = 0.8) +
      geom_hline(yintercept = -log10(alpha), lty = 4, col = "black", lwd = 0.8) +
      labs(title = res$label) +
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
# 8. Refactored: Parallel MetFrag Analysis
# ══════════════════════════════════════════════════════════════════════════════

#' Parallel MetFrag in silico fragmentation validation
#'
#' The biggest bottleneck in the original pipeline. Refactored to:
#'   1. Pre-compute the MS2 spectra lookup table (avoids repeated scans)
#'   2. Use future_map for parallel feature processing
#'   3. Implement file-hash caching per feature
#'   4. Fall back to sequential if parallel fails
#'
#' @param app3 Current annotation data.frame
#' @param object mass_dataset object with MS2 spectra
#' @param cfg Validated config list
#' @return Updated annotation data.frame with MetFrag columns
run_metfrag_analysis_parallel <- function(app3, object, cfg) {
  cat("\n===== Step 6c: MetFrag in silico fragmentation (parallel + cached) =====\n")

  if (!isTRUE(cfg$metfrag$enabled)) {
    cat("   MetFrag disabled in config, skipping\n")
    return(app3)
  }
  if (!requireNamespace("metfRag", quietly = TRUE)) {
    cat("   metfRag package not installed. Skipping.\n")
    return(app3)
  }
  suppressPackageStartupMessages(library(metfRag))

  if (is.null(app3) || nrow(app3) == 0) {
    cat("   No identifications to validate, skipping\n")
    return(app3)
  }

  spectra_data <- object@ms2_data
  if (is.null(spectra_data) || length(spectra_data) == 0) {
    cat("   No MS2 spectra in object, cannot run MetFrag\n")
    return(app3)
  }

  cat(sprintf("   Validating %d identified features with MetFrag...\n", nrow(app3)))

  ann_cfg <- cfg$annotation
  metfrag_cfg <- cfg$metfrag
  cache_dir <- file.path(cfg$project$working_dir,
                         cfg$performance$cache_dir %||% "cache", "metfrag")
  max_age <- cfg$performance$cache_max_age_days %||% 30
  enable_cache <- isTRUE(cfg$performance$enable_caching %||% TRUE)

  # Pre-build MS2 lookup table: variable_id -> list(ms2_spec, precursor_mz, precursor_rt)
  # This avoids O(n*m) repeated scans of the spectra_data list
  cat("   Building MS2 spectra lookup table...\n")
  ms2_lookup <- new.env(parent = emptyenv())
  for (md in spectra_data) {
    if (!isS4(md) || !inherits(md, "ms2_data")) next
    for (j in seq_along(md@variable_id)) {
      vid <- md@variable_id[j]
      spec <- md@ms2_spectra[[j]]
      if (is.null(spec) || is.null(vid) || is.na(vid)) next
      ms2_lookup[[vid]] <- list(
        spec = spec,
        mz = if (length(md@ms2_mz) >= j) md@ms2_mz[j] else NA,
        rt = if (length(md@ms2_rt) >= j) md@ms2_rt[j] else NA
      )
    }
  }
  cat(sprintf("   Lookup table: %d features with MS2 spectra\n", length(ls(ms2_lookup))))

  # Evict expired cache entries
  if (enable_cache) cache_evict_expired(cache_dir, max_age)

  # Identify features that have both annotation and MS2
  features_to_process <- which(
    app3$variable_id %in% ls(ms2_lookup) &
    !is.na(app3$mz) &
    !is.na(app3$Compound.name) & app3$Compound.name != ""
  )
  cat(sprintf("   Features to process: %d\n", length(features_to_process)))

  if (length(features_to_process) == 0) {
    cat("   No features with both annotation and MS2 spectra, skipping\n")
    return(app3)
  }

  # Initialize output columns
  app3$metfrag_score <- NA_real_
  app3$metfrag_num_explained_peaks <- NA_integer_
  app3$metfrag_fragmenter_score <- NA_real_

  # Process one feature (for parallel execution)
  process_one_feature <- function(i) {
    var_id <- app3$variable_id[i]
    lookup <- ms2_lookup[[var_id]]
    if (is.null(lookup)) return(NULL)

    precursor_mz <- app3$mz[i]
    if (is.na(precursor_mz)) return(NULL)

    compound_name <- app3$Compound.name[i]
    formula <- if ("Formula" %in% colnames(app3)) {
      f <- as.character(app3$Formula[i])
      if (is.na(f)) "" else f
    } else ""

    ms2_spec <- lookup$spec
    if (is.null(ms2_spec)) return(NULL)
    ms2_df <- as.data.frame(ms2_spec)
    if (nrow(ms2_df) == 0) return(NULL)

    # Check cache first
    if (enable_cache) {
      param_hash <- compute_combined_hash(
        var_id, precursor_mz, compound_name, formula,
        as.matrix(ms2_df),
        metfrag_cfg$ppm_tol %||% 15,
        metfrag_cfg$mzabs %||% 0.005,
        algo = cfg$performance$hash_algorithm %||% "md5"
      )
      cached <- cache_lookup(param_hash, cache_dir, max_age, cfg)
      if (!is.null(cached)) {
        return(qs::qread(cached))
      }
    }

    # Build MetFrag settings
    s <- create.settings.sample()
    s[["MetFragDatabaseType"]] <- "PubChem"
    s[["DatabaseSearchRelativeMassDeviation"]] <- metfrag_cfg$ppm_tol %||% ann_cfg$ms1_match_ppm
    s[["FragmentPeakMatchAbsoluteMassDeviation"]] <- metfrag_cfg$mzabs %||% 0.005
    s[["FragmentPeakMatchRelativeMassDeviation"]] <- metfrag_cfg$ppm_tol %||% 15
    s[["NeutralPrecursorMass"]] <- precursor_mz
    s[["SampleName"]] <- compound_name
    if (formula != "") s[["NeutralPrecursorMolecularFormula"]] <- formula
    s[["PeakList"]] <- cbind(mz = ms2_df$mz, intensity = ms2_df$intensity)

    # Number of threads for MetFrag
    s[["NumberThreads"]] <- cfg$performance$metfrag_threads %||% 2

    result <- tryCatch(run.metfrag(s), error = function(e) NULL)

    if (!is.null(result) && nrow(result) > 0) {
      best <- result[1, , drop = FALSE]
      out <- list(
        score = if ("Score" %in% colnames(best)) as.numeric(best$Score[1]) else NA,
        n_explained = if ("NumberExplainedPeaks" %in% colnames(best)) as.integer(best$NumberExplainedPeaks[1]) else NA,
        fragmenter_score = if ("FragmenterScore" %in% colnames(best)) as.numeric(best$FragmenterScore[1]) else NA
      )
      # Cache the result
      if (enable_cache) {
        cache_store(param_hash, out, cache_dir, cfg)
      }
      return(out)
    }
    NULL
  }

  # Run in parallel with sequential fallback
  n_workers <- nbrOfWorkers()
  if (n_workers > 1 && length(features_to_process) >= 10) {
    cat(sprintf("   Processing %d features in parallel (%d workers)...\n",
                length(features_to_process), n_workers))

    results <- tryCatch({
      future_lapply(features_to_process, process_one_feature,
                    future.scheduling = 1.0, future.seed = TRUE)
    }, error = function(e) {
      cat(sprintf("   WARNING: Parallel MetFrag failed: %s\n", conditionMessage(e)))
      cat("   Falling back to sequential\n")
      lapply(features_to_process, process_one_feature)
    })
  } else {
    cat(sprintf("   Processing %d features sequentially...\n", length(features_to_process)))
    results <- lapply(features_to_process, process_one_feature)
  }

  # Merge results back
  n_scored <- 0
  for (k in seq_along(features_to_process)) {
    i <- features_to_process[k]
    res <- results[[k]]
    if (is.null(res)) next
    app3$metfrag_score[i] <- res$score
    app3$metfrag_num_explained_peaks[i] <- res$n_explained
    app3$metfrag_fragmenter_score[i] <- res$fragmenter_score
    if (!is.na(res$fragmenter_score) && res$fragmenter_score > 0) n_scored <- n_scored + 1
  }

  cat(sprintf("   MetFrag complete: %d processed, %d scored\n",
              length(features_to_process), n_scored))

  # Update confidence levels
  if (n_scored > 0) {
    app3 <- app3 %>%
      mutate(
        confidence_level = ifelse(
          !is.na(metfrag_fragmenter_score) & metfrag_fragmenter_score > 0,
          "Level 1 — Confirmed", confidence_level
        )
      )
    cat(sprintf("   Upgraded %d features to Level 1 based on MetFrag validation\n", n_scored))
  }

  app3
}


# ══════════════════════════════════════════════════════════════════════════════
# 9. Refactored: Parallel Boxplot Generation
# ══════════════════════════════════════════════════════════════════════════════

#' Parallel boxplot generation using future_map
#'
#' @param object2 mass_dataset object
#' @param sample_info Sample information data.frame
#' @param app3 Annotation table
#' @param diff_results Output of run_differential()
#' @param cfg Validated config list
#' @return Invisible NULL
run_boxplot_parallel <- function(object2, sample_info, app3, diff_results, cfg) {
  cat("\n===== Step 9: Boxplot (parallel) =====\n")

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

  # Collect top-N metabolites
  metabolite_set <- character(0)
  for (dr in diff_results) {
    if (is.null(dr) || is.null(dr$sig) || nrow(dr$sig) == 0) next
    top_sig <- dr$sig %>% arrange(P.Value) %>% head(5)
    metabolite_set <- c(metabolite_set, rownames(top_sig))
  }
  metabolite_set <- unique(metabolite_set)
  if (length(metabolite_set) == 0) {
    cat("   No significant metabolites to plot\n")
    return(invisible(NULL))
  }

  # Deduplicate by compound name
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
  plot_metabolites <- head(plot_metabolites, viz_cfg$boxplot_top_n %||% 30)
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

  # Draw one boxplot (for parallel)
  draw_one_boxplot <- function(vid) {
    meta_row <- app3_expr[app3_expr$variable_id == vid, ]
    if (nrow(meta_row) == 0) return(NULL)

    compound_name <- ifelse(
      !is.na(meta_row$Compound.name[1]) && meta_row$Compound.name[1] != "",
      meta_row$Compound.name[1], vid)

    box_values <- as.numeric(meta_row[1, all_samples])
    box_df <- data.frame(value = box_values, group = core_group_vec, stringsAsFactors = FALSE)
    box_df <- box_df[!is.na(box_df$value), ]
    if (nrow(box_df) < 4 || length(unique(box_df$group)) < 2) return(NULL)

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

    safe_name <- gsub("[^A-Za-z0-9]", "_", compound_name)
    save_plot(file.path(out_dir, paste0("boxplot_", safe_name)),
              plot = p, width = 8, height = 7, dpi = viz_cfg$dpi)
    TRUE
  }

  # Execute in parallel with sequential fallback
  n_workers <- nbrOfWorkers()
  if (n_workers > 1 && length(plot_metabolites) >= 6) {
    cat(sprintf("   Generating %d boxplots in parallel (%d workers)...\n",
                length(plot_metabolites), n_workers))
    tryCatch({
      future_map(plot_metabolites, draw_one_boxplot, .progress = FALSE,
                 .options = furrr_options(seed = TRUE))
    }, error = function(e) {
      cat(sprintf("   WARNING: Parallel boxplot failed: %s\n", conditionMessage(e)))
      cat("   Falling back to sequential\n")
      for (vid in plot_metabolites) draw_one_boxplot(vid)
    })
  } else {
    for (vid in plot_metabolites) draw_one_boxplot(vid)
  }

  cat("   Boxplots complete\n")
  invisible(NULL)
}


# ══════════════════════════════════════════════════════════════════════════════
# 10. Refactored: Parallel Annotated Volcano Plots
# ══════════════════════════════════════════════════════════════════════════════

#' Parallel annotated volcano plot generation
#'
#' @param diff_results Output of run_differential()
#' @param app3 Annotation table
#' @param cfg Validated config list
#' @return Invisible NULL
run_annotated_volcano_parallel <- function(diff_results, app3, cfg) {
  cat("\n===== Step 6.5: Annotated volcano plots (parallel) =====\n")

  if (is.null(app3) || nrow(app3) == 0) {
    cat("   No annotations available, skipping\n")
    return(invisible(NULL))
  }

  out_dir <- file.path(cfg$project$output_dir, "04_Differential")
  diff_cfg <- cfg$differential
  alpha <- diff_cfg$p_value_cutoff
  logfc_cutoff <- log2(diff_cfg$fc_threshold)
  sig_metric <- diff_cfg$significance_metric %||% "p_value"
  p_col <- if (identical(sig_metric, "adj_p_value")) "adj.P.Val" else "P.Value"
  p_lab <- if (identical(sig_metric, "adj_p_value")) "adj. p-value (FDR)" else "p-value"

  annot_lookup <- app3 %>%
    dplyr::select(variable_id, Compound.name) %>%
    filter(!is.na(Compound.name) & Compound.name != "")

  # Draw one volcano (for parallel)
  draw_one_volcano <- function(cmp_name) {
    dr <- diff_results[[cmp_name]]
    if (is.null(dr)) return(NULL)

    diff_all <- dr$all
    label <- dr$comparison$label

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
      geom_hline(yintercept = -log10(alpha), lty = 4, col = "black", lwd = 0.8) +
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
    TRUE
  }

  cmp_names <- names(diff_results)
  n_workers <- nbrOfWorkers()
  if (n_workers > 1 && length(cmp_names) >= 3) {
    cat(sprintf("   Generating %d volcano plots in parallel (%d workers)...\n",
                length(cmp_names), n_workers))
    tryCatch({
      future_map(cmp_names, draw_one_volcano, .progress = FALSE,
                 .options = furrr_options(seed = TRUE))
    }, error = function(e) {
      cat(sprintf("   WARNING: Parallel volcano failed: %s\n", conditionMessage(e)))
      cat("   Falling back to sequential\n")
      for (nm in cmp_names) draw_one_volcano(nm)
    })
  } else {
    for (nm in cmp_names) draw_one_volcano(nm)
  }

  invisible(NULL)
}


# ══════════════════════════════════════════════════════════════════════════════
# 11. Refactored: SIRIUS with Hash-Based Caching
# ══════════════════════════════════════════════════════════════════════════════

#' SIRIUS analysis with file-hash caching
#'
#' Before running SIRIUS, computes an MD5 hash of the input MGF file and
#' all parameters. If the output directory exists and the hash matches,
#' skips the computation entirely (memoization). If SIRIUS is re-run with
#' different parameters, the hash changes and a new computation is triggered.
#'
#' @param object mass_dataset object with MS2 spectra loaded
#' @param app3 Current annotation data.frame
#' @param cfg Validated config list
#' @return Updated annotation data.frame with SIRIUS columns
run_sirius_analysis_cached <- function(object, app3, cfg) {
  cat("\n===== Step 6b: SIRIUS (cached + multi-threaded) =====\n")

  if (!isTRUE(cfg$sirius$enabled)) {
    cat("   SIRIUS disabled in config, skipping\n")
    return(app3)
  }

  sirius_path <- cfg$sirius$path
  if (is.null(sirius_path) || !file.exists(sirius_path)) {
    cat(sprintf("   SIRIUS not found at: %s, skipping\n", sirius_path))
    return(app3)
  }
  cat(sprintf("   SIRIUS found at: %s\n", sirius_path))

  # ── Export MS2 spectra to MGF ──────────────────────────────────────────
  mgf_dir <- file.path(cfg$project$output_dir, "sirius_input")
  dir.create(mgf_dir, showWarnings = FALSE, recursive = TRUE)
  mgf_file <- file.path(mgf_dir, "ms2_spectra.mgf")

  spectra_data <- object@ms2_data
  if (is.null(spectra_data) || length(spectra_data) == 0) {
    cat("   No MS2 spectra in object, cannot run SIRIUS\n")
    return(app3)
  }

  # Write MGF using fast concatenation (avoid repeated cat() calls)
  cat("   Exporting MS2 spectra to MGF...\n")
  mgf_lines <- character(0)
  for (i in seq_along(spectra_data)) {
    md <- spectra_data[[i]]
    if (is.null(md) || length(md@ms2_spectra) == 0) next
    for (j in seq_along(md@ms2_spectra)) {
      spec <- md@ms2_spectra[[j]]
      if (is.null(spec) || nrow(as.data.frame(spec)) == 0) next
      var_id <- md@variable_id[j]
      pmz <- md@ms2_mz[j]
      prt <- md@ms2_rt[j]
      if (is.null(var_id) || is.na(var_id) || is.null(pmz) || is.na(pmz)) next

      ms2_df <- as.data.frame(spec)
      charge <- ifelse(cfg$feature_extraction$polarity == "positive", "1+", "1-")

      block <- c(
        "BEGIN IONS",
        sprintf("TITLE=%s_scan%d", var_id, j),
        sprintf("PEPMASS=%.6f", pmz),
        sprintf("RTINSECONDS=%.2f", if (is.null(prt) || is.na(prt)) 0 else prt),
        sprintf("CHARGE=%s", charge),
        sprintf("%.6f %.4f", ms2_df$mz, ms2_df$intensity),
        "END IONS"
      )
      mgf_lines <- c(mgf_lines, block)
    }
  }
  writeLines(mgf_lines, mgf_file)
  cat(sprintf("   Exported %d spectra to: %s\n",
              sum(grepl("BEGIN IONS", mgf_lines)), mgf_file))
  rm(mgf_lines); gc()

  # ── Compute cache hash ──────────────────────────────────────────────────
  cache_dir <- file.path(cfg$project$working_dir,
                         cfg$performance$cache_dir %||% "cache", "sirius")
  enable_cache <- isTRUE(cfg$performance$enable_caching %||% TRUE)
  max_age <- cfg$performance$cache_max_age_days %||% 30
  sirius_threads <- cfg$performance$sirius_threads %||% 4

  if (enable_cache) {
    cache_evict_expired(cache_dir, max_age)

    # Hash the MGF + all parameters that affect SIRIUS output
    cache_key <- compute_combined_hash(
      mgf_file,
      list(
        ppm_max = cfg$sirius$ppm_max %||% 15,
        polarity = cfg$feature_extraction$polarity,
        sirius_path = sirius_path,
        sirius_version = tryCatch(system2(sirius_path, "--version", stdout = TRUE),
                                   error = function(e) "unknown"),
        threads = sirius_threads
      ),
      algo = cfg$performance$hash_algorithm %||% "md5"
    )

    # Check if output already exists and hash matches
    project_space <- cfg$sirius$project_space
    sirius_output <- ifelse(grepl("\\.sirius$", project_space),
                            project_space,
                            file.path(project_space, "project.sirius"))
    hash_file <- file.path(dirname(sirius_output), ".sirius_cache_hash")

    if (file.exists(sirius_output) && file.exists(hash_file)) {
      stored_hash <- readLines(hash_file, warn = FALSE)[1]
      if (identical(stored_hash, cache_key)) {
        cat(sprintf("   Cache HIT: SIRIUS output unchanged (hash: %s)\n",
                    substr(cache_key, 1, 16)))
        cat("   Skipping SIRIUS computation, loading cached results...\n")
        return(merge_sirius_results(app3, sirius_output, cfg))
      } else {
        cat("   Cache KEY MISMATCH: parameters changed, re-running SIRIUS\n")
      }
    }
  }

  # ── Run SIRIUS with multi-threading ────────────────────────────────────
  project_space <- cfg$sirius$project_space
  sirius_output <- ifelse(grepl("\\.sirius$", project_space),
                          project_space,
                          file.path(project_space, "project.sirius"))
  dir.create(dirname(sirius_output), showWarnings = FALSE, recursive = TRUE)

  ppm_max <- cfg$sirius$ppm_max %||% 15
  threads_arg <- sprintf("--threads=%d", sirius_threads)

  # Formula prediction
  cat("   Launching SIRIUS formula prediction...\n")
  rc_formula <- system2(sirius_path,
    args = c("--log", "WARNING", threads_arg,
             "-i", mgf_file, "-o", sirius_output,
             "formula", "-p", "qtof", "--ppm-max", ppm_max),
    wait = TRUE)
  if (rc_formula != 0) {
    cat(sprintf("   WARNING: SIRIUS formula returned exit code %d\n", rc_formula))
  }

  # Structure database search
  cat("   Launching SIRIUS structure database search...\n")
  rc_structure <- system2(sirius_path,
    args = c("--log", "WARNING", threads_arg,
             "-i", mgf_file, "-o", sirius_output,
             "structure", "-d", "HMDB,CHEBI"),
    wait = TRUE)
  if (rc_structure != 0) {
    cat(sprintf("   WARNING: SIRIUS structure returned exit code %d\n", rc_structure))
  }

  # Fingerprint prediction
  cat("   Launching SIRIUS fingerprint prediction...\n")
  rc_fp <- system2(sirius_path,
    args = c("--log", "WARNING", threads_arg,
             "-i", mgf_file, "-o", sirius_output, "fingerprint"),
    wait = TRUE)
  if (rc_fp != 0) {
    cat(sprintf("   WARNING: SIRIUS fingerprint returned exit code %d\n", rc_fp))
  }

  cat(sprintf("   SIRIUS analysis complete. Results saved to: %s\n", sirius_output))

  # Store cache hash
  if (enable_cache) {
    writeLines(cache_key, hash_file)
    cat(sprintf("   Cache hash stored: %s\n", hash_file))
  }

  # Merge results
  merge_sirius_results(app3, sirius_output, cfg)
}

#' Merge SIRIUS results from summary files into the annotation table
#'
#' @param app3 Annotation data.frame
#' @param sirius_output Path to SIRIUS project output
#' @param cfg Config list
#' @return Updated annotation data.frame
merge_sirius_results <- function(app3, sirius_output, cfg) {
  sirius_result_dir <- dirname(sirius_output)
  summary_files <- list.files(sirius_result_dir,
    pattern = "summary.*\\.tsv$", recursive = TRUE, full.names = TRUE)

  if (length(summary_files) == 0) {
    cat("   No SIRIUS summary files found. Results may need manual review.\n")
    return(app3)
  }

  cat(sprintf("   Found %d SIRIUS summary file(s)\n", length(summary_files)))
  for (sf in summary_files) {
    sirius_results <- tryCatch(fast_read(sf), error = function(e) NULL)
    if (is.null(sirius_results) || nrow(sirius_results) == 0) next

    cat(sprintf("   SIRIUS results: %d rows from %s\n", nrow(sirius_results), basename(sf)))

    sirius_results$variable_id <- gsub("_.*$", "", sirius_results$title)
    sirius_cols <- c("variable_id")
    if ("molecularFormula" %in% colnames(sirius_results)) {
      sirius_cols <- c(sirius_cols, "molecularFormula")
    }
    if ("score" %in% colnames(sirius_results)) {
      sirius_cols <- c(sirius_cols, "score")
    }
    for (extra_col in c("InChI", "InChIkey", "smiles", "xlogp", "name")) {
      if (extra_col %in% colnames(sirius_results)) {
        sirius_cols <- c(sirius_cols, extra_col)
      }
    }

    sirius_merge <- sirius_results %>%
      dplyr::select(dplyr::any_of(sirius_cols)) %>%
      dplyr::rename(
        sirius_formula = dplyr::any_of("molecularFormula"),
        sirius_score   = dplyr::any_of("score"),
        sirius_InChI   = dplyr::any_of("InChI"),
        sirius_InChIkey = dplyr::any_of("InChIkey"),
        sirius_smiles  = dplyr::any_of("smiles"),
        sirius_name    = dplyr::any_of("name")
      )

    app3 <- app3 %>% left_join(sirius_merge, by = "variable_id")
    n_sirius <- sum(!is.na(app3$sirius_formula))
    cat(sprintf("   SIRIUS results merged: %d features with formula/structure\n", n_sirius))

    if (n_sirius > 0) {
      app3 <- app3 %>%
        mutate(
          confidence_level = ifelse(
            !is.na(sirius_formula) &
              confidence_level %in% c("Level 3 — Putatively annotated", "Level 4 — Unknown"),
            "Level 2 — Putatively identified",
            confidence_level
          )
        )
      cat("   Upgraded confidence levels for SIRIUS-validated features\n")
    }
  }

  app3
}


# ══════════════════════════════════════════════════════════════════════════════
# 12. Refactored: Fast Preprocessing
# ══════════════════════════════════════════════════════════════════════════════

#' Fast preprocessing pipeline with matrix-optimized operations
#'
#' Replaces the original run_preprocessing with:
#'   - fread for peak table loading
#'   - matrix-based RSD calculation (matrix_rsd)
#'   - matrix-based blank filter (fast_blank_filter)
#'   - Explicit memory cleanup after MS2 loading
#'
#' @param cfg Validated config list
#' @param meta_list Output of load_metadata()
#' @param sample_info_all Output of build_sample_info_all()
#' @return List with same structure as run_preprocessing
run_preprocessing_fast <- function(cfg, meta_list, sample_info_all) {
  cat("\n===== Step 2: Data preprocessing (optimized) =====\n")

  # ── Load peak table with fread ─────────────────────────────────────────
  peak_dir <- cfg$data$peak_table_dir
  peak_file <- file.path(peak_dir, "Peak_table_for_cleaning.csv")
  if (!file.exists(peak_file)) {
    peak_file <- file.path(peak_dir, "peak_table_for_cleaning.csv")
  }
  if (!file.exists(peak_file)) {
    stop(sprintf("Peak table not found: %s", peak_file))
  }

  raw_data <- fast_read(peak_file, row.names = 1)
  cat(sprintf("=> Peak table: %d features x %d samples\n", nrow(raw_data), ncol(raw_data)))

  expression_data <- as.matrix(raw_data[, -1:-2, drop = FALSE])
  variable_info <- raw_data[, 1:2, drop = FALSE]
  variable_info$variable_id <- rownames(raw_data)
  variable_info <- variable_info[, c(3, 1, 2)]
  rownames(variable_info) <- NULL

  # ── Match group membership ─────────────────────────────────────────────
  col_map <- meta_list$col_map
  fn_col <- col_map$filename

  match_group <- function(sid, meta, fn_col, grp_col) {
    idx <- which(meta[[fn_col]] == paste0(sid, ".mzML"))
    if (length(idx) > 0) return(meta[[grp_col]][idx[1]])
    idx <- which(gsub("\\.mzML$", "", meta[[fn_col]]) == sid)
    if (length(idx) > 0) return(meta[[grp_col]][idx[1]])
    gsub("[0-9-].*$", "", sid)
  }

  sample_ids_in_data <- colnames(expression_data)
  classes <- vapply(sample_ids_in_data, function(sid) {
    match_group(sid, meta_list$all, fn_col, col_map$group_id)
  }, character(1))

  # Tag MS2 samples
  ms_col <- col_map$ms_level
  if (!is.null(ms_col) && ms_col %in% colnames(meta_list$all)) {
    ms2_fns <- meta_list$all[[fn_col]][meta_list$all[[ms_col]] == "MSMS"]
    ms2_ids <- gsub("\\.mzML$", "", ms2_fns, ignore.case = TRUE)
    classes[colnames(expression_data) %in% ms2_ids] <- paste0(classes[colnames(expression_data) %in% ms2_ids], "-MSMS")
  }

  sample_info_all <- data.frame(
    sample_id = sample_ids_in_data,
    class = classes,
    stringsAsFactors = FALSE
  )
  cat("=> Sample types:\n")
  print(table(sample_info_all$class))

  f <- cfg$filtering
  blank_name <- cfg$sample_roles$blank
  qc_name <- cfg$sample_roles$qc_ms1

  # ── 2a: Blank filter (matrix-optimized) ────────────────────────────────
  cat("\n----- Step 2a: Blank background filter (fast) -----\n")
  blank_samples <- sample_info_all$sample_id[
    grepl(paste0("^", blank_name), sample_info_all$class) &
      !grepl("-MSMS$", sample_info_all$class)]
  bio_samples <- sample_info_all$sample_id[
    !(sample_info_all$class %in% c(blank_name, paste0(blank_name, "-MSMS"),
                                   qc_name, paste0(qc_name, "-MSMS")))]

  cat(sprintf("   Blank samples: %d, Biological samples: %d\n",
              length(blank_samples), length(bio_samples)))

  if (length(blank_samples) >= 2) {
    blank_pass <- fast_blank_filter(expression_data, blank_samples, bio_samples,
                                     f$blank_fold_change)
    cat(sprintf("   Retained: %d/%d (removed %d background features)\n",
                sum(blank_pass), length(blank_pass), sum(!blank_pass)))
    expression_data <- expression_data[blank_pass, , drop = FALSE]
    variable_info <- variable_info[blank_pass, , drop = FALSE]
  } else {
    cat("   Insufficient blank samples, skipping blank filter\n")
  }

  # ── 2b: Missing value filter (fast) ────────────────────────────────────
  cat("\n----- Step 2b: Missing value filter (fast) -----\n")
  bio_class <- sample_info_all$class[match(bio_samples, sample_info_all$sample_id)]
  names(bio_class) <- bio_samples

  mv_pass <- fast_missing_filter(expression_data, bio_class, f$max_missing_per_group)
  cat(sprintf("   Retained: %d/%d (removed %d high-missing features)\n",
              sum(mv_pass), length(mv_pass), sum(!mv_pass)))
  expression_data <- expression_data[mv_pass, , drop = FALSE]
  variable_info <- variable_info[mv_pass, , drop = FALSE]

  # ── 2c: QC-RSD filter (matrix-optimized) ───────────────────────────────
  cat("\n----- Step 2c: QC-RSD filter (fast) -----\n")
  qc_samples <- sample_info_all$sample_id[sample_info_all$class == qc_name]
  cat(sprintf("   QC samples (MS1 only): %d\n", length(qc_samples)))

  if (length(qc_samples) >= 3) {
    qc_data <- expression_data[, qc_samples, drop = FALSE]
    qc_rsd <- matrix_rsd(qc_data, min_values = 3)
    rsd_pass <- is.na(qc_rsd) | qc_rsd <= f$rsd_threshold
    cat(sprintf("   Retained: %d/%d (removed %d RSD > %.0f%% features)\n",
                sum(rsd_pass), length(rsd_pass), sum(!rsd_pass), f$rsd_threshold))
    expression_data <- expression_data[rsd_pass, , drop = FALSE]
    variable_info <- variable_info[rsd_pass, , drop = FALSE]
  } else {
    cat("   Insufficient QC samples, skipping RSD filter\n")
  }

  # ── 2d: Imputation + SVR correction ────────────────────────────────────
  # (identical to original — tidymass normalize_data is the bottleneck here)
  cat("\n----- Step 2d: Imputation + normalization -----\n")
  norm_cfg <- cfg$normalization

  sample_info_svr <- sample_info_all
  sample_info_svr$injection.order <- seq_len(nrow(sample_info_svr))
  sample_info_svr$class[!sample_info_svr$class %in%
    c(blank_name, paste0(blank_name, "-MSMS"),
      qc_name, paste0(qc_name, "-MSMS"))] <- "Subject"
  sample_info_svr$class[sample_info_svr$class == qc_name] <- "QC"

  object <- create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info_svr,
    variable_info = variable_info
  )

  mv_filled <- impute_mv(object = object, method = "minimum")
  cat("   Missing value imputation complete (minimum method)\n")

  isOK <- apply(mv_filled@expression_data, 1, function(row) {
    !any(row < f$min_intensity, na.rm = TRUE)
  })
  cat(sprintf("   Low-signal filter (< %g): retained %d/%d features\n",
              f$min_intensity, sum(isOK), length(isOK)))
  mv_filled@expression_data <- mv_filled@expression_data[isOK, , drop = FALSE]
  mv_filled@variable_info <- mv_filled@variable_info[isOK, , drop = FALSE]

  # Clear intermediate objects
  rm(raw_data, object); gc()

  qc_count <- sum(sample_info_all$class == qc_name)
  if (qc_count >= norm_cfg$qc_min_samples) {
    cat(sprintf("   Running QC drift correction (%d QC samples)...\n", qc_count))
    data_ok <- apply(mv_filled@expression_data, 1, function(r) all(is.finite(r)))
    if (sum(data_ok) < nrow(mv_filled@expression_data)) {
      cat(sprintf("   Removing %d non-finite features\n", sum(!data_ok)))
      mv_filled@expression_data <- mv_filled@expression_data[data_ok, , drop = FALSE]
      mv_filled@variable_info <- mv_filled@variable_info[data_ok, , drop = FALSE]
    }

    qc_in_obj <- mv_filled@sample_info$sample_id[mv_filled@sample_info$class == "QC"]
    qc_data2 <- mv_filled@expression_data[, qc_in_obj, drop = FALSE]
    qc_var <- Rfast::colVars(qc_data2)  # fast column variance

    # Use matrixStats::colVars if Rfast not available
    if (any(is.na(qc_var)) || length(qc_var) == 0) {
      qc_var <- apply(qc_data2, 1, function(x) var(as.numeric(x)))
    }
    zero_var <- qc_var == 0 | is.na(qc_var)
    if (sum(zero_var) > 0) {
      cat(sprintf("   Removing %d zero-variance QC features\n", sum(zero_var)))
      mv_filled@expression_data <- mv_filled@expression_data[!zero_var, , drop = FALSE]
      mv_filled@variable_info <- mv_filled@variable_info[!zero_var, , drop = FALSE]
    }

    if (nrow(mv_filled@expression_data) > 10) {
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
    } else {
      cat("   Insufficient features, using median normalization\n")
      object2 <- normalize_data(mv_filled, method = norm_cfg$fallback_method)
    }
    cat("   Normalization complete\n")
  } else {
    cat(sprintf("   Insufficient QC samples, using median normalization\n"))
    object2 <- normalize_data(mv_filled, method = norm_cfg$fallback_method)
  }

  # Clean up imputation objects
  rm(mv_filled, isOK, data_ok, qc_data, qc_data2, qc_var, zero_var); gc()

  # Restore class labels
  orig_class <- sample_info_all$class
  names(orig_class) <- sample_info_all$sample_id
  object2@sample_info$class <- orig_class[object2@sample_info$sample_id]

  # ── 2e: Extract core groups ────────────────────────────────────────────
  cat("\n----- Step 2e: Extract core groups -----\n")
  core_classes <- setdiff(unique(sample_info_all$class),
    c(blank_name, paste0(blank_name, "-MSMS"),
      qc_name, paste0(qc_name, "-MSMS")))
  core_samples <- sample_info_all$sample_id[sample_info_all$class %in% core_classes]
  cat(sprintf("   Core samples: %d\n", length(core_samples)))

  sample_info <- object2@sample_info
  sample_info <- sample_info[sample_info$sample_id %in% core_samples, , drop = FALSE]
  rownames(sample_info) <- NULL

  expression_data <- object2@expression_data[, core_samples, drop = FALSE]
  variable_info <- object2@variable_info

  cat(sprintf("   Final: %d features x %d samples\n", nrow(expression_data), ncol(expression_data)))
  print(table(sample_info$class))

  object2 <- create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

  # ── 2f: Load MS2 spectra ──────────────────────────────────────────────
  cat("\n----- Step 2f: Load MS2 spectra -----\n")
  ms2_loaded <- FALSE
  meta_ms2 <- meta_list$ms2

  if (nrow(meta_ms2) >= 1) {
    blank_name <- cfg$sample_roles$blank
    grp_col_ms2 <- meta_list$col_map$group_id
    is_blank_ms2 <- meta_ms2[[grp_col_ms2]] == blank_name
    if (any(is_blank_ms2)) {
      cat(sprintf("   Excluding %d MS2 blank run(s)\n", sum(is_blank_ms2)))
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

      # Check MS2 memory footprint
      if (ms2_loaded && !is.null(cfg$performance)) {
        check_ms2_memory(object2, cfg)
      }
    }
  } else {
    cat("   No MS2 samples in metadata, skipping\n")
  }

  # Final GC
  rm(orig_class, sample_info_all, meta_ms2, ms2_paths, ms2_files, ms2_exists); gc()

  list(
    object2 = object2,
    sample_info = sample_info,
    variable_info = variable_info,
    expression_data = expression_data,
    core_classes = core_classes,
    ms2_loaded = ms2_loaded
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# 13. Override helper: Wrap fast_* functions into the main pipeline
# ══════════════════════════════════════════════════════════════════════════════

#' Check if the performance module is active and return the config
#'
#' @param cfg Full config list
#' @return TRUE if performance optimization is enabled
is_performance_enabled <- function(cfg) {
  perf <- cfg$performance
  if (is.null(perf)) return(FALSE)
  io_engine <- perf$io_engine %||% "base"
  io_engine != "base"
}

#' Replace all I/O functions in the calling environment with fast versions
#'
#' This should be called at the top of main() to replace read.csv/write.csv
#' globally with fread/fwrite across the entire pipeline.
#'
#' @param cfg Config list
install_fast_io <- function(cfg) {
  if (!is_performance_enabled(cfg)) {
    cat("   High-performance I/O disabled (io_engine = 'base')\n")
    return(invisible(FALSE))
  }

  if (!requireNamespace("data.table", quietly = TRUE)) {
    cat("   data.table not available, cannot install fast I/O\n")
    return(invisible(FALSE))
  }

  cat("   Fast I/O available: fast_read() / fast_write() ready\n")
  cat("   Use fread/fwrite for all CSV operations\n")
  invisible(TRUE)
}