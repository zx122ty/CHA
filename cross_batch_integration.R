#!/usr/bin/env Rscript
###############################################################################
# cross_batch_integration.R — Cross-Batch Integration for herbMetabo
#
# Robust cross-batch integration for untargeted metabolomics:
#   1. Feature Alignment across batches using m/z (ppm) and RT tolerance window
#   2. Batch-effect correction with ComBat (empirical Bayes) and limma (removeBatchEffect)
#   3. Dual-stage Dimensionality Reduction (Before: Raw vs After: Corrected)
#   4. Merged differential results across all comparisons and pathways
#   5. Web-ready exports (integrated_pca.json, matrices, cross_batch_rows.json)
#
# Dependencies: base R, jsonlite, sva, limma.
###############################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("package 'jsonlite' is required")
  }
  library(jsonlite)
})

# -----------------------------------------------------------------------------
# 0. Minimal CLI parsing (base R only; "--key value" or "--key=value")
# -----------------------------------------------------------------------------
parse_CLI <- function(argv = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  i <- 1L
  n <- length(argv)
  while (i <= n) {
    a <- argv[i]
    if (grepl("^--[^=]+=.*$", a)) {
      key <- sub("^--([^=]+)=.*$", "\\1", a)
      out[[gsub("-", "_", key)]] <- sub("^--[^=]+=", "", a)
      i <- i + 1L
    } else if (startsWith(a, "--")) {
      key <- gsub("-", "_", sub("^--", "", a))
      i <- i + 1L
      val <- if (i <= n && !startsWith(argv[i], "--")) {
        v <- argv[i]; i <- i + 1L; v
      } else TRUE
      out[[key]] <- val
    } else {
      i <- i + 1L
    }
  }
  out
}

opts <- parse_CLI()
defaults <- list(
  input_dir = NULL,
  output_dir = NULL,
  batches = NULL,
  mz_ppm = 12,
  rt_tol = 20,
  min_overlap = 200,
  no_log2 = FALSE
)
for (nm in names(defaults)) if (is.null(opts[[nm]])) opts[[nm]] <- defaults[[nm]]

if (is.null(opts$input_dir)) stop("--input-dir is required")
if (is.null(opts$output_dir)) stop("--output-dir is required")
opts$mz_ppm <- as.numeric(opts$mz_ppm)
opts$rt_tol <- as.numeric(opts$rt_tol)
opts$min_overlap <- as.integer(max(1, as.integer(opts$min_overlap)))
opts$no_log2 <- isTRUE(opts$no_log2)
if (!is.null(opts$batches)) {
  opts$batches <- trimws(strsplit(opts$batches, ",")[[1]])
}
dir.create(opts$output_dir, showWarnings = FALSE, recursive = TRUE)

fmt <- function(t) format(Sys.time(), "%Y-%m-%d %H:%M:%S")
logmsg <- function(...) cat(sprintf("[%s] ", fmt(Sys.time())), ..., "\n")

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a

readj <- function(f) {
  if (!file.exists(f)) return(NULL)
  tryCatch(fromJSON(f, simplifyVector = FALSE), error = function(e) NULL)
}

# -----------------------------------------------------------------------------
# Stage 1 helpers: locate + read per-batch inputs
# -----------------------------------------------------------------------------

find_expr_file <- function(root) {
  if (is.null(root) || !dir.exists(root)) return(NULL)
  cand <- list.files(root, recursive = TRUE, full.names = TRUE, pattern = "\\.csv$")
  cand <- cand[grepl("normaliz", cand, ignore.case = TRUE)]
  if (length(cand) == 0) {
    cand <- list.files(root, recursive = TRUE, full.names = TRUE, pattern = "\\.csv$")
    cand <- cand[grepl("cleaned_expression|corrected|norm", cand, ignore.case = TRUE)]
  }
  if (length(cand) == 0) return(NULL)
  norm_idx <- grep("03_normalized|03_normal|03_PCA", cand)
  if (length(norm_idx) > 0) cand <- cand[norm_idx]
  sort(cand)[1]
}

read_expr <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  m <- read.csv(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  m <- as.matrix(m)
  storage.mode(m) <- "double"
  m
}

maybe_log2 <- function(m, force = FALSE) {
  if (is.null(m)) return(m)
  # Check if values appear to be on raw scale (> 1000)
  num_vals <- m[!is.na(m) & is.finite(m)]
  if (force || (length(num_vals) > 0 && max(num_vals) >= 1000)) {
    m[m <= 0] <- NA
    m <- log2(m)
  }
  m
}

features_of <- function(features, batch_id) {
  if (is.null(features) || length(features) == 0) return(NULL)
  id <- vapply(features, function(f) if (is.null(f$variable_id)) "" else as.character(f$variable_id), character(1))
  mz <- as.numeric(vapply(features, function(f) if (is.null(f$mz)) NA else f$mz, numeric(1)))
  rt <- as.numeric(vapply(features, function(f) if (is.null(f$rt)) NA else f$rt, numeric(1)))
  name <- vapply(features, function(f) if (is.null(f$compound_name)) "" else as.character(f$compound_name), character(1))
  hmdb <- vapply(features, function(f) if (is.null(f$hmdb_id)) NA_character_ else as.character(f$hmdb_id), character(1))
  kegg <- vapply(features, function(f) if (is.null(f$kegg_id)) NA_character_ else as.character(f$kegg_id), character(1))
  formula <- vapply(features, function(f) if (is.null(f$formula)) "" else as.character(f$formula), character(1))
  db_src <- vapply(features, function(f) if (is.null(f$database_source)) "" else as.character(f$database_source), character(1))
  conf <- vapply(features, function(f) if (is.null(f$confidence_level)) "" else as.character(f$confidence_level), character(1))

  df <- data.frame(
    batch_id = as.integer(batch_id),
    variable_id = id,
    mz = mz,
    rt = rt,
    compound_name = name,
    hmdb_id = hmdb,
    kegg_id = kegg,
    formula = formula,
    database_source = db_src,
    confidence_level = conf,
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$mz) & !is.na(df$rt), , drop = FALSE]
  df
}

parse_volcano <- function(j) {
  pts <- j$points
  if (is.null(pts) || length(pts) == 0) return(NULL)
  id  <- vapply(pts, function(p) if (is.null(p$name)) "" else as.character(p$name), character(1))
  lfc <- as.numeric(vapply(pts, function(p) if (is.null(p$logFC)) NA_character_ else as.character(p$logFC), character(1)))
  pv  <- as.numeric(vapply(pts, function(p) if (is.null(p$p_value)) NA_character_ else as.character(p$p_value), character(1)))
  fd  <- as.numeric(vapply(pts, function(p) if (is.null(p$fdr)) NA_character_ else as.character(p$fdr), character(1)))
  grp <- vapply(pts, function(p) if (is.null(p$group)) "Not" else as.character(p$group), character(1))
  hm  <- vapply(pts, function(p) if (is.null(p$hmdb_id)) NA_character_ else as.character(p$hmdb_id), character(1))
  data.frame(
    variable_id = id,
    logFC = lfc,
    p_value = pv,
    fdr = fd,
    dir = grp,
    hmdb = hm,
    stringsAsFactors = FALSE
  )
}

load_batch_diffs <- function(web_dir) {
  ff <- list.files(web_dir, pattern = "^volcano_.*\\.json$", full.names = TRUE)
  if (length(ff) == 0) return(list())
  out <- list()
  for (f in ff) {
    cmp <- sub("^volcano_", "", sub("\\.json$", "", basename(f)))
    j <- readj(f)
    df <- parse_volcano(j)
    if (!is.null(df) && nrow(df) > 0) out[[cmp]] <- df
  }
  out
}

load_batch_pathways <- function(web) {
  files <- list.files(web, pattern = "^pathway_enrichment_.*\\.json$", full.names = TRUE)
  if (length(files) == 0) return(list())
  out <- list()
  for (f in files) {
    cmp <- sub("^pathway_enrichment_|\\.json$", "", basename(f))
    j <- readj(f)
    pws <- j$pathways
    if (is.null(pws) || length(pws) == 0) next
    rows <- do.call(rbind, lapply(pws, function(p) {
      data.frame(
        pathway = if (is.null(p$pathway_name)) "" else as.character(p$pathway_name),
        p_value = as.numeric(if (is.null(p$p_value)) NA else p$p_value),
        mapped  = as.numeric(if (is.null(p$mapped_count)) 0 else p$mapped_count),
        background = as.numeric(if (is.null(p$background_count)) 0 else p$background_count),
        impact  = as.numeric(if (is.null(p$impact_score)) NA else p$impact_score),
        stringsAsFactors = FALSE
      )
    }))
    out[[cmp]] <- rows
  }
  out
}

# -----------------------------------------------------------------------------
# Stage 2: Cross-batch feature alignment (Union-Find with sliding window on m/z)
# -----------------------------------------------------------------------------
align_features <- function(feature_df, mz_ppm = 12, rt_tol = 20) {
  n <- nrow(feature_df)
  if (n == 0) return(feature_df)

  ord <- order(feature_df$mz)
  df <- feature_df[ord, ]
  rownames(df) <- NULL

  mz_vec <- df$mz
  rt_vec <- df$rt
  b_vec  <- df$batch_id

  parent <- seq_len(n)
  find_root <- function(i) {
    curr <- i
    while (parent[curr] != curr) {
      parent[curr] <<- parent[parent[curr]]
      curr <- parent[curr]
    }
    curr
  }
  union_nodes <- function(i, j) {
    ri <- find_root(i)
    rj <- find_root(j)
    if (ri != rj) parent[rj] <<- ri
  }

  for (i in seq_len(n - 1)) {
    mzi <- mz_vec[i]
    rti <- rt_vec[i]
    bi  <- b_vec[i]
    max_mz <- mzi * (1 + mz_ppm * 1e-6)

    j <- i + 1L
    while (j <= n && mz_vec[j] <= max_mz) {
      if (b_vec[j] != bi && abs(rt_vec[j] - rti) <= rt_tol) {
        union_nodes(i, j)
      }
      j <- j + 1L
    }
  }

  roots <- vapply(seq_len(n), find_root, integer(1))
  df$root_id <- roots

  # Group summary: canonical mz and rt (median of members)
  root_stats <- aggregate(cbind(mz, rt) ~ root_id, data = df, FUN = median)
  root_mzs <- setNames(root_stats$mz, root_stats$root_id)
  root_rts <- setNames(root_stats$rt, root_stats$root_id)

  # Format canonical group ID e.g. G0001_z447.123
  df$group_id <- sprintf("G%05d_z%.3f_rt%d", df$root_id, root_mzs[as.character(df$root_id)], round(root_rts[as.character(df$root_id)]))
  df$canonical_mz <- root_mzs[as.character(df$root_id)]
  df$canonical_rt <- root_rts[as.character(df$root_id)]

  # Map back to original order
  orig_order <- match(seq_len(n), ord)
  df[orig_order, ]
}

# -----------------------------------------------------------------------------
# Stage 3: Statistical Batch Correction (ComBat, limma, median)
# -----------------------------------------------------------------------------
run_combat <- function(M, batch_vec) {
  if (!requireNamespace("sva", quietly = TRUE)) return(NULL)
  logmsg("  Running ComBat (Empirical Bayes)...")
  res <- tryCatch({
    sva::ComBat(dat = M, batch = as.factor(batch_vec), mod = NULL,
                par.prior = TRUE, prior.plots = FALSE)
  }, error = function(e) {
    logmsg("    parametric ComBat failed, trying mean.only version...")
    tryCatch({
      sva::ComBat(dat = M, batch = as.factor(batch_vec), mod = NULL,
                  par.prior = FALSE, mean.only = TRUE, prior.plots = FALSE)
    }, error = function(e2) {
      logmsg("    ComBat failed completely:", e2$message)
      NULL
    })
  })
  res
}

run_limma <- function(M, batch_vec) {
  if (!requireNamespace("limma", quietly = TRUE)) return(NULL)
  logmsg("  Running limma (removeBatchEffect)...")
  res <- tryCatch({
    limma::removeBatchEffect(M, batch = as.factor(batch_vec))
  }, error = function(e) {
    logmsg("    limma failed:", e$message)
    NULL
  })
  res
}

run_median_center <- function(M, batch_vec) {
  logmsg("  Running per-batch median centering...")
  M2 <- M
  for (bb in unique(batch_vec)) {
    cols <- which(batch_vec == bb)
    b_meds <- apply(M2[, cols, drop = FALSE], 1, median, na.rm = TRUE)
    M2[, cols] <- M2[, cols] - b_meds
  }
  M2
}

# Calculate batch effect metric: F-statistic & R2 of PC1 on batch
batch_mixing_stat <- function(pca_scores, batch_vec) {
  if (is.null(pca_scores) || length(unique(batch_vec)) < 2) return(list(F = 0, R2 = 0))
  fit <- lm(pca_scores[, 1] ~ as.factor(batch_vec))
  s <- summary(fit)
  f_val <- if (!is.null(s$fstatistic)) round(s$fstatistic[1], 2) else 0
  r2_val <- round(s$r.squared, 4)
  list(F = f_val, R2 = r2_val)
}

# PCA computation helper returning formatted list for JSON
compute_pca_details <- function(M, batch_vec, sample_ids, stage_name, method_name) {
  if (is.null(M) || nrow(M) < 3 || ncol(M) < 3) return(NULL)
  pca <- prcomp(t(M), center = TRUE, scale. = FALSE)
  nc <- min(6, ncol(pca$x))
  pc <- pca$x[, 1:nc, drop = FALSE]
  var_exp <- round(100 * summary(pca)$importance[2, 1:nc], 2)
  names(var_exp) <- paste0("PC", 1:nc)

  scores <- lapply(seq_len(nrow(pc)), function(i) list(
    sample_id = sample_ids[i],
    batch_id = as.integer(batch_vec[i]),
    PC1 = round(as.numeric(pc[i, 1]), 4),
    PC2 = round(as.numeric(pc[i, 2]), 4),
    PC3 = if (nc >= 3) round(as.numeric(pc[i, 3]), 4) else 0,
    PC4 = if (nc >= 4) round(as.numeric(pc[i, 4]), 4) else 0,
    PC5 = if (nc >= 5) round(as.numeric(pc[i, 5]), 4) else 0,
    PC6 = if (nc >= 6) round(as.numeric(pc[i, 6]), 4) else 0
  ))

  stat <- batch_mixing_stat(pc, batch_vec)

  list(
    stage = stage_name,
    method = method_name,
    variance_explained = as.list(var_exp),
    batch_F_stat = stat$F,
    batch_R2 = stat$R2,
    scores = scores
  )
}

# -----------------------------------------------------------------------------
# Stage 4: Per-metabolite meta-analysis stats
# -----------------------------------------------------------------------------
fisher_p <- function(ps) {
  ps <- ps[!is.na(ps)]
  ps <- ps[ps > 0 & ps < 1]
  if (length(ps) == 0) return(NA_real_)
  stat <- -2 * sum(log(ps))
  pchisq(stat, df = 2 * length(ps), lower.tail = FALSE)
}

stouffer_p <- function(ps) {
  ps <- ps[!is.na(ps)]
  ps <- ps[ps > 0 & ps < 1]
  if (length(ps) == 0) return(NA_real_)
  z <- sum(qnorm(1 - ps)) / sqrt(length(ps))
  pnorm(z, lower.tail = FALSE)
}

integrate_metabolite <- function(gid, rows) {
  n <- nrow(rows)
  if (n == 0) return(NULL)
  n_sig <- sum(rows$dir %in% c("Up", "Down"), na.rm = TRUE)
  n_up  <- sum(rows$dir == "Up", na.rm = TRUE)
  n_dn  <- sum(rows$dir == "Down", na.rm = TRUE)
  med_fc <- stats::median(rows$logFC, na.rm = TRUE)
  min_p  <- min(rows$p_value, na.rm = TRUE)
  fisher <- if (sum(!is.na(rows$p_value)) >= 2) fisher_p(rows$p_value) else min_p
  stoufp <- if (sum(!is.na(rows$p_value)) >= 2) stouffer_p(rows$p_value) else NA

  vid_tab <- table(rows$variable_id)
  rep_vid <- names(vid_tab)[which.max(vid_tab)]
  dir_consistency <- if (n_sig > 0) max(n_up, n_dn) / n_sig else NA_real_

  list(
    group_id = as.character(gid),
    variable_id = rep_vid,
    n_treatments = n,
    n_batches = length(unique(rows$batch_id)),
    n_batches_sig = length(unique(rows$batch_id[rows$dir %in% c("Up", "Down")])),
    n_significant = n_sig,
    n_up = n_up,
    n_down = n_dn,
    median_logFC = round(med_fc, 4),
    min_p_value = round(min_p, 6),
    fisher_p = if (!is.na(fisher)) format(fisher, scientific = TRUE) else NA,
    stouffer_p = if (!is.na(stoufp)) format(stoufp, scientific = TRUE) else NA,
    direction_consistency = round(dir_consistency, 3)
  )
}

summarize_feature_groups <- function(diff_long) {
  if (is.null(diff_long) || nrow(diff_long) == 0) return(list())
  ks <- split(seq_len(nrow(diff_long)), diff_long$group_id)
  res <- lapply(names(ks), function(k) {
    r <- diff_long[ks[[k]], , drop = FALSE]
    integrate_metabolite(k, r)
  })
  res
}

# -----------------------------------------------------------------------------
# Main Driver
# -----------------------------------------------------------------------------
main <- function() {
  logmsg("cross_batch_integration.R — start")
  manifest <- readj(file.path(opts$input_dir, "integration_manifest.json"))
  if (is.null(manifest)) stop("no integration_manifest.json in --input-dir")
  inc <- manifest$included_batches
  if (length(inc) == 0) stop("manifest contains no batches")

  if (!is.null(opts$batches)) {
    inc <- inc[vapply(inc, function(b) as.character(b$batch_id) %in% opts$batches, logical(1))]
  }
  logmsg(sprintf("  manifest included batches to check: %d", length(inc)))

  batch_meta <- list()
  excluded_batches_record <- manifest$excluded_batches %||% list()

  for (b in inc) {
    bid <- as.character(b$batch_id)
    web <- b$web_dir
    logmsg(sprintf("  checking batch %s : %s", bid, web))

    man <- readj(file.path(web, "manifest.json"))
    if (is.null(man)) next

    # QC quantification gate
    qq <- man$qc_quality
    if (!is.null(qq) && !isTRUE(qq$quantification_valid)) {
      reason <- qq$exclude_reason %||% "quantification_valid = false"
      logmsg(sprintf("  batch %s EXCLUDED from quantitative integration: %s", bid, reason))
      if (!is.null(qq$qc_rsd_median)) {
        logmsg(sprintf("    QC-RSD median %.1f%%, QC-QC corr %s, %s features",
                       as.numeric(qq$qc_rsd_median),
                       format(qq$qc_qc_correlation %||% NA_real_),
                       format(qq$n_features %||% NA_integer_)))
      }
      excluded_batches_record[[bid]] <- list(
        batch_id = as.integer(bid),
        name = b$name %||% paste0("Batch", bid),
        reason = reason,
        qc_rsd_median = qq$qc_rsd_median,
        qc_qc_correlation = qq$qc_qc_correlation,
        n_features = qq$n_features
      )
      next
    }

    out_root <- if (!is.null(b$output_dir)) b$output_dir else dirname(web)
    cmp_meta <- man$treatments

    raw_feats <- features_of(readj(file.path(web, "metabolites.json"))$metabolites, bid)
    diffs <- load_batch_diffs(web)
    for (cmp in names(diffs)) diffs[[cmp]]$batch_id <- as.integer(bid)

    expr_path <- find_expr_file(out_root)
    expr <- NULL
    if (!is.null(expr_path)) {
      m_raw <- read_expr(expr_path)
      m_raw <- maybe_log2(m_raw, force = !opts$no_log2)
      expr <- m_raw
    } else {
      logmsg(sprintf("  no expression matrix for batch %s", bid))
    }

    batch_meta[[bid]] <- list(
      batch_id = bid, cmp = cmp_meta, web = web,
      diffs = diffs, expr = expr, features = raw_feats
    )
  }

  valid_bids <- names(batch_meta)
  logmsg(sprintf("  batches passing QC gate for integration: %s (count: %d)",
                 paste(valid_bids, collapse = ", "), length(valid_bids)))

  if (length(valid_bids) < 2) {
    stop("Fewer than 2 batches passed QC gate for quantitative integration.")
  }

  # ---------------- Step 2. Feature Alignment ----------------
  logmsg("Aligning features across valid batches...")
  all_feats_list <- lapply(batch_meta, function(b) b$features)
  combined_feats <- do.call(rbind, all_feats_list)
  aligned_feats <- align_features(combined_feats, mz_ppm = opts$mz_ppm, rt_tol = opts$rt_tol)

  # Attach group_id mapping back to each batch
  for (bid in valid_bids) {
    b_af <- aligned_feats[aligned_feats$batch_id == as.integer(bid), ]
    v2g <- setNames(b_af$group_id, b_af$variable_id)
    batch_meta[[bid]]$v2g <- v2g
  }

  # ---------------- Step 3. Assemble Unified Matrix ----------------
  logmsg("Assembling unified cross-batch matrix...")
  sample_names <- unlist(lapply(valid_bids, function(bid) {
    sprintf("B%s_%s", bid, colnames(batch_meta[[bid]]$expr))
  }))
  sample_batches <- sub("^B([0-9]+)_.*$", "\\1", sample_names)

  # Count presence of feature groups across batches
  gp_batch_tab <- table(aligned_feats$group_id, aligned_feats$batch_id)
  gp_n_batches <- rowSums(gp_batch_tab > 0)

  # For cross-batch matrix & batch correction, use feature groups present in >= 2 batches
  shared_groups <- names(gp_n_batches)[gp_n_batches >= 2]
  logmsg(sprintf("  total unified groups: %d | shared in >=2 batches: %d",
                 length(gp_n_batches), length(shared_groups)))

  M_raw <- matrix(NA_real_, nrow = length(shared_groups), ncol = length(sample_names),
                  dimnames = list(shared_groups, sample_names))

  for (bid in valid_bids) {
    b_em <- batch_meta[[bid]]$expr
    b_v2g <- batch_meta[[bid]]$v2g
    b_cols <- sprintf("B%s_%s", bid, colnames(b_em))

    # map rows of b_em to group_ids
    matched_vids <- intersect(rownames(b_em), names(b_v2g))
    if (length(matched_vids) == 0) next

    sub_em <- b_em[matched_vids, , drop = FALSE]
    sub_gids <- unname(b_v2g[matched_vids])

    # aggregate duplicate features mapping to same group
    for (gid in intersect(unique(sub_gids), shared_groups)) {
      rows_idx <- which(sub_gids == gid)
      if (length(rows_idx) == 1) {
        M_raw[gid, b_cols] <- sub_em[rows_idx, ]
      } else {
        M_raw[gid, b_cols] <- colMeans(sub_em[rows_idx, , drop = FALSE], na.rm = TRUE)
      }
    }
  }

  # Imputation for batch correction algorithms:
  # Replace NA with minimum observed per feature - 1.0 (log2 scale) + small Gaussian noise
  set.seed(42)
  M_imputed <- M_raw
  for (r in seq_len(nrow(M_imputed))) {
    rv <- M_imputed[r, ]
    na_idx <- which(is.na(rv))
    if (length(na_idx) > 0) {
      min_val <- min(rv, na.rm = TRUE)
      fill_val <- min_val - 1.0
      M_imputed[r, na_idx] <- rnorm(length(na_idx), mean = fill_val, sd = 0.08)
    }
  }

  # ---------------- Step 4. Batch Correction & PCA ----------------
  logmsg("Running batch effect corrections and PCA...")

  # 1. Before Correction (Raw)
  pca_before <- compute_pca_details(M_imputed, sample_batches, sample_names, "before", "raw")
  logmsg(sprintf("  [Before] Batch PC1 F-stat: %.1f, R2: %.3f",
                 pca_before$batch_F_stat, pca_before$batch_R2))

  # 2. ComBat
  M_combat <- run_combat(M_imputed, sample_batches)
  pca_combat <- if (!is.null(M_combat)) {
    compute_pca_details(M_combat, sample_batches, sample_names, "after", "combat")
  } else NULL
  if (!is.null(pca_combat)) {
    logmsg(sprintf("  [After/ComBat] Batch PC1 F-stat: %.1f, R2: %.3f",
                   pca_combat$batch_F_stat, pca_combat$batch_R2))
  }

  # 3. limma
  M_limma <- run_limma(M_imputed, sample_batches)
  pca_limma <- if (!is.null(M_limma)) {
    compute_pca_details(M_limma, sample_batches, sample_names, "after", "limma")
  } else NULL
  if (!is.null(pca_limma)) {
    logmsg(sprintf("  [After/limma] Batch PC1 F-stat: %.1f, R2: %.3f",
                   pca_limma$batch_F_stat, pca_limma$batch_R2))
  }

  # 4. Median Centering
  M_median <- run_median_center(M_imputed, sample_batches)
  pca_median <- compute_pca_details(M_median, sample_batches, sample_names, "after", "median_centering")
  logmsg(sprintf("  [After/Median] Batch PC1 F-stat: %.1f, R2: %.3f",
                 pca_median$batch_F_stat, pca_median$batch_R2))

  # Primary corrected matrix for export (prefer ComBat, then limma, then median)
  primary_method <- if (!is.null(M_combat)) "combat" else if (!is.null(M_limma)) "limma" else "median_centering"
  M_corrected_primary <- if (primary_method == "combat") M_combat else if (primary_method == "limma") M_limma else M_median
  pca_primary <- if (primary_method == "combat") pca_combat else if (primary_method == "limma") pca_limma else pca_median

  # ---------------- Step 5. Merge Differential Results ----------------
  logmsg("Merging differential analysis comparisons...")
  merged_diffs <- lapply(batch_meta, function(b) {
    if (length(b$diffs) == 0) return(NULL)
    v2g <- b$v2g
    res <- lapply(names(b$diffs), function(cmp) {
      d <- b$diffs[[cmp]]
      gid <- unname(v2g[d$variable_id])
      d$group_id <- ifelse(is.na(gid), d$variable_id, gid)
      d$cmp <- cmp
      d$batch_id <- as.character(b$batch_id)
      d
    })
    do.call(rbind, res)
  })
  diffs_rows <- do.call(rbind, merged_diffs)
  integrated_stats <- summarize_feature_groups(diffs_rows)
  logmsg(sprintf("  %d unified feature groups statistically integrated", length(integrated_stats)))

  # ---------------- Step 6. Merge Pathway Enrichments ----------------
  logmsg("Merging pathway enrichment results...")
  pw_rows <- do.call(rbind, lapply(batch_meta, function(b) {
    if (is.null(b$web)) return(NULL)
    pws <- load_batch_pathways(b$web)
    if (length(pws) == 0) return(NULL)
    do.call(rbind, lapply(names(pws), function(cmp) {
      p <- pws[[cmp]]
      p$cmp <- cmp
      p$batch_id <- as.character(b$batch_id)
      p
    }))
  }))
  integrated_paths <- list()
  if (!is.null(pw_rows) && nrow(pw_rows) > 0) {
    gpw <- split(seq_len(nrow(pw_rows)), pw_rows$pathway)
    integrated_paths <- lapply(names(gpw), function(pn) {
      r <- pw_rows[gpw[[pn]], , drop = FALSE]
      ps <- suppressWarnings(as.numeric(r$p_value))
      list(
        pathway = pn,
        n_comparisons = nrow(r),
        n_batches = length(unique(r$batch_id)),
        min_p = round(min(ps, na.rm = TRUE), 6),
        fisher_p = ifelse(sum(!is.na(ps)) >= 2,
                          format(fisher_p(ps), scientific = TRUE), NA),
        total_mapped = sum(r$mapped, na.rm = TRUE)
      )
    })
  }

  # ---------------- Step 7. Write Web Artifacts ----------------
  wdir <- file.path(opts$output_dir, "web_export")
  dir.create(wdir, showWarnings = FALSE, recursive = TRUE)

  # Multi-stage PCA JSON payload
  pca_json_data <- list(
    pca = TRUE,
    primary_method = primary_method,
    available_methods = c("combat", "limma", "median_centering"),
    stages = list(
      before = pca_before,
      after = list(
        combat = pca_combat,
        limma = pca_limma,
        median_centering = pca_median
      )
    ),
    # Backwards-compatible top-level keys
    scores = pca_primary$scores,
    variance_explained = pca_primary$variance_explained,
    n_samples = ncol(M_raw),
    n_features = nrow(M_raw)
  )
  write_json(pca_json_data, file.path(wdir, "integrated_pca.json"), pretty = TRUE, auto_unbox = TRUE)

  # Expression matrices
  write.csv(M_raw, file.path(wdir, "raw_expression_matrix.csv"))
  write.csv(M_corrected_primary, file.path(wdir, "integrated_expression_matrix.csv"))

  # DB-import rows
  dbrows <- lapply(integrated_stats, function(e) {
    list(
      variable_id = e$variable_id,
      key = e$group_id,
      integration_method = "cross_batch",
      overall_logFC = e$median_logFC,
      overall_p_value = e$fisher_p,
      n_batches_sig = e$n_batches_sig,
      direction_consistency = e$direction_consistency
    )
  })
  write_json(list(rows = dbrows), file.path(wdir, "cross_batch_rows.json"), pretty = TRUE, auto_unbox = TRUE)

  # Merged logFC matrix
  keys <- vapply(integrated_stats, function(e) e$group_id, character(1))
  cmps <- unique(diffs_rows[, c("cmp", "batch_id")])
  cmps <- cmps[order(cmps$batch_id, cmps$cmp), , drop = FALSE]
  colnam <- sprintf("B%s:%s", cmps$batch_id, cmps$cmp)
  mat <- matrix(NA_real_, nrow = length(keys), ncol = length(colnam), dimnames = list(keys, colnam))
  mi <- match(diffs_rows$group_id, keys)
  mj <- match(sprintf("B%s:%s", diffs_rows$batch_id, diffs_rows$cmp), colnam)
  keep <- !is.na(mi) & !is.na(mj)
  mat[cbind(mi[keep], mj[keep])] <- diffs_rows$logFC[keep]
  write.csv(mat, file.path(wdir, "integrated_logFC_matrix.csv"))

  # Updated integration manifest
  int_manifest_out <- list(
    integration = list(
      name = "Cross-Batch Integration",
      date = fmt(Sys.time()),
      n_batches = length(valid_bids),
      batches = as.list(valid_bids),
      primary_method = primary_method,
      methods_available = list(
        combat = !is.null(M_combat),
        limma = !is.null(M_limma),
        median_centering = TRUE,
        meta_analysis = TRUE
      ),
      qc_metrics = list(
        raw_F = pca_before$batch_F_stat,
        raw_R2 = pca_before$batch_R2,
        corrected_F = pca_primary$batch_F_stat,
        corrected_R2 = pca_primary$batch_R2
      )
    ),
    included_batches = lapply(batch_meta, function(b) list(
      batch_id = as.integer(b$batch_id),
      name = paste0("Batch", b$batch_id),
      n_treatments = length(b$cmp),
      web_dir = b$web
    )),
    excluded_batches = excluded_batches_record
  )
  write_json(int_manifest_out, file.path(wdir, "integration_manifest.json"), pretty = TRUE, auto_unbox = TRUE)

  # Report summary text
  report_lines <- c(
    "Cross-batch integration report",
    sprintf("  generated: %s", fmt(Sys.time())),
    sprintf("  included batches: %s", paste(valid_bids, collapse = ", ")),
    sprintf("  excluded batches: %s", paste(names(excluded_batches_record), collapse = ", ")),
    sprintf("  primary correction method: %s", primary_method),
    sprintf("  batch mixing metric (F-stat on PC1): before=%.1f -> after=%.1f",
            pca_before$batch_F_stat, pca_primary$batch_F_stat),
    sprintf("  comparisons: %d", nrow(cmps)),
    sprintf("  unified feature groups: %d", length(integrated_stats)),
    sprintf("  shared pathways: %d", length(integrated_paths))
  )
  writeLines(report_lines, file.path(opts$output_dir, "integration_report.txt"))

  logmsg("cross_batch_integration.R — finished successfully!")
  invisible(list(pca = pca_json_data, diffs = diffs_rows, pathways = integrated_paths))
}

if (sys.nframe() == 0) {
  main()
}
