###############################################################################
# cck8_utils.R — shared helpers for the CCK8 x metabolomics analysis
#
# Conventions follow code/this_project (base R + jsonlite + ggplot2).
# features_of() / align_features() are adapted from cross_batch_integration.R
# so that group IDs (G#####_z..._rt...) match the stored integration exactly.
###############################################################################

suppressPackageStartupMessages({
  library(jsonlite)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a

fmt_time <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
logmsg <- function(...) cat(sprintf("[%s] ", fmt_time()), ..., "\n", sep = "")

readj <- function(f) {
  if (is.null(f) || !file.exists(f)) return(NULL)
  tryCatch(fromJSON(f, simplifyVector = FALSE), error = function(e) NULL)
}

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
read_cck8_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  cfg$paths$project_root <- normalizePath(cfg$paths$project_root, mustWork = TRUE)
  for (pk in c("cck8_xlsx", "browse_table", "file_lists_dir", "results_dir", "integration_dir")) {
    v <- cfg$paths[[pk]]
    if (!is.null(v) && !startsWith(v, "/")) cfg$paths[[pk]] <- file.path(cfg$paths$project_root, v)
  }
  if (startsWith(cfg$paths$output_dir_rel, "/")) {
    cfg$paths$output_dir <- normalizePath(cfg$paths$output_dir_rel, mustWork = FALSE)
  } else {
    cfg$paths$output_dir <- file.path(cfg$paths$project_root, cfg$paths$output_dir_rel)
  }
  for (sub in c("data", "figures", "reports", "web_export")) {
    dir.create(file.path(cfg$paths$output_dir, sub), showWarnings = FALSE, recursive = TRUE)
  }
  cfg
}

# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------
spearman_pair <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(unique(x)) < 2 || length(unique(y)) < 2 || length(x) < 5) {
    return(c(rho = NA_real_, p = NA_real_, n = length(x)))
  }
  ct <- suppressWarnings(stats::cor.test(x, y, method = "spearman"))
  c(rho = as.numeric(ct$estimate), p = as.numeric(ct$p.value), n = length(x))
}

bootstrap_spearman_ci <- function(x, y, n_boot = 1000, seed = 42) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 5 || length(unique(x)) < 2 || length(unique(y)) < 2) return(c(NA_real_, NA_real_))
  set.seed(seed)
  idx <- matrix(sample.int(n, n * n_boot, replace = TRUE), ncol = n_boot)
  rhos <- apply(idx, 2, function(ii) {
    if (length(unique(x[ii])) < 2 || length(unique(y[ii])) < 2) return(NA_real_)
    stats::cor(x[ii], y[ii], method = "spearman")
  })
  rhos <- rhos[is.finite(rhos)]
  if (length(rhos) == 0) return(c(NA_real_, NA_real_))
  as.numeric(stats::quantile(rhos, c(0.025, 0.975), names = FALSE))
}

# ---------------------------------------------------------------------------
# Feature alignment (adapted from cross_batch_integration.R — keep in sync)
# ---------------------------------------------------------------------------
features_of <- function(features, batch_id) {
  if (is.null(features) || length(features) == 0) return(NULL)
  id <- vapply(features, function(f) if (is.null(f$variable_id)) "" else as.character(f$variable_id), character(1))
  mz <- vapply(features, function(f) if (is.null(f$mz)) NA_real_ else as.numeric(as.character(f$mz)), numeric(1))
  rt <- vapply(features, function(f) if (is.null(f$rt)) NA_real_ else as.numeric(as.character(f$rt)), numeric(1))
  name <- vapply(features, function(f) if (is.null(f$compound_name)) "" else as.character(f$compound_name), character(1))
  hmdb <- vapply(features, function(f) if (is.null(f$hmdb_id)) NA_character_ else as.character(f$hmdb_id), character(1))
  kegg <- vapply(features, function(f) if (is.null(f$kegg_id)) NA_character_ else as.character(f$kegg_id), character(1))
  formula <- vapply(features, function(f) if (is.null(f$formula)) "" else as.character(f$formula), character(1))

  df <- data.frame(
    batch_id = as.integer(batch_id),
    variable_id = id, mz = mz, rt = rt,
    compound_name = name, hmdb_id = hmdb, kegg_id = kegg, formula = formula,
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$mz) & !is.na(df$rt), , drop = FALSE]
  df
}

align_features <- function(feature_df, mz_ppm = 12, rt_tol = 20) {
  n <- nrow(feature_df)
  if (n == 0) return(feature_df)

  ord <- order(feature_df$mz)
  df <- feature_df[ord, ]
  rownames(df) <- NULL

  mz_vec <- df$mz; rt_vec <- df$rt; b_vec <- df$batch_id

  parent <- seq_len(n)
  find_root <- function(i) {
    curr <- i
    while (parent[curr] != curr) { parent[curr] <<- parent[parent[curr]]; curr <- parent[curr] }
    curr
  }
  union_nodes <- function(i, j) {
    ri <- find_root(i); rj <- find_root(j)
    if (ri != rj) parent[rj] <<- ri
  }

  for (i in seq_len(n - 1)) {
    mzi <- mz_vec[i]; rti <- rt_vec[i]; bi <- b_vec[i]
    max_mz <- mzi * (1 + mz_ppm * 1e-6)
    j <- i + 1L
    while (j <= n && mz_vec[j] <= max_mz) {
      if (b_vec[j] != bi && abs(rt_vec[j] - rti) <= rt_tol) union_nodes(i, j)
      j <- j + 1L
    }
  }

  roots <- vapply(seq_len(n), find_root, integer(1))
  df$root_id <- roots

  root_stats <- aggregate(cbind(mz, rt) ~ root_id, data = df, FUN = median)
  root_mzs <- setNames(root_stats$mz, root_stats$root_id)
  root_rts <- setNames(root_stats$rt, root_stats$root_id)

  df$group_id <- sprintf("G%05d_z%.3f_rt%d", df$root_id,
                         root_mzs[as.character(df$root_id)],
                         round(root_rts[as.character(df$root_id)]))
  df$canonical_mz <- root_mzs[as.character(df$root_id)]
  df$canonical_rt <- round(root_rts[as.character(df$root_id)])

  orig_order <- match(seq_len(n), ord)
  df[orig_order, ]
}

# Build the group -> annotation table for all valid batches (same inputs and
# parameters as cross_batch_integration.R so group IDs are identical).
build_group_annotation <- function(cfg, verbose = TRUE) {
  feats_list <- list()
  for (bid in cfg$integration$batches) {
    web <- file.path(cfg$paths$results_dir, paste0("Batch", bid), "detail", "web_export")
    m <- readj(file.path(web, "metabolites.json"))
    if (is.null(m)) { logmsg("  WARNING: no metabolites.json for batch ", bid); next }
    f <- features_of(m$metabolites, bid)
    if (!is.null(f) && nrow(f) > 0) feats_list[[as.character(bid)]] <- f
    if (verbose) logmsg("  batch ", bid, ": ", if (is.null(f)) 0L else nrow(f), " features")
  }
  combined <- do.call(rbind, feats_list)
  aligned <- align_features(combined, mz_ppm = cfg$integration$mz_ppm, rt_tol = cfg$integration$rt_tol)

  # one row per group: pick the best-annotated member (non-empty name preferred)
  has_name <- ifelse(aligned$compound_name != "", 0L, 1L)
  ord <- order(aligned$root_id, has_name, aligned$variable_id)
  grp <- aligned[ord, ]
  grp <- grp[!duplicated(grp$group_id), , drop = FALSE]
  grp <- grp[, c("group_id", "canonical_mz", "canonical_rt", "compound_name",
                 "formula", "hmdb_id", "kegg_id"), drop = FALSE]
  rownames(grp) <- NULL
  list(groups = grp, aligned = aligned)
}

# Parse integrated_logFC_matrix.csv columns "B{batch}:{drug}_{High|Low}_vs_CT{n}"
parse_comparison_columns <- function(colnames_vec, valid_batches) {
  keep <- vapply(colnames_vec, function(cn) {
    m <- regmatches(cn, regexec("^B([0-9]+):([0-9]+)_(High|Low)_vs_CT[0-9]+$", cn))[[1]]
    length(m) == 4 && as.integer(m[2]) %in% valid_batches
  }, logical(1))
  cols <- colnames_vec[keep]
  parsed <- do.call(rbind, lapply(cols, function(cn) {
    m <- regmatches(cn, regexec("^B([0-9]+):([0-9]+)_(High|Low)_vs_CT[0-9]+$", cn))[[1]]
    data.frame(batch_id = as.integer(m[2]), drug_id = as.integer(m[3]),
               concentration = m[4], colname = cn, stringsAsFactors = FALSE)
  }))
  parsed
}

# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------
save_plot <- function(p, base_path) {
  png_file <- paste0(base_path, ".png")
  pdf_file <- paste0(base_path, ".pdf")
  grDevices::png(png_file, width = 1400, height = 900, res = 110)
  print(p)
  grDevices::dev.off()
  grDevices::pdf(pdf_file, width = 7, height = 4.6)
  print(p)
  grDevices::dev.off()
  invisible(c(png_file, pdf_file))
}

md_table <- function(df, colnames_override = NULL) {
  if (!is.null(colnames_override)) colnames(df) <- colnames_override
  hdr <- paste0("| ", paste(colnames(df), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  rows <- apply(df, 1, function(r) paste0("| ", paste(vapply(r, as.character, character(1)), collapse = " | "), " |"))
  c(hdr, sep, rows)
}
