#!/usr/bin/env Rscript
###############################################################################
# conc_utils.R — shared utilities for the concentration-response module.
# Conventions follow cck8_analysis/cck8_utils.R (config reading, group
# annotation via features_of()/align_features() adapted from
# cross_batch_integration.R — keep those in sync if integration params change).
#
# Design notes (two-concentration extensibility, plan §2.8 / A5):
#   - All dose handling goes through cfg$doses (a list), never hard-coded
#     "Low"/"High". Adding a third concentration = extend the list in
#     concentration_config.yaml + provide the new columns upstream.
#   - pairwise_dose_gradient() computes one-difference gradients for any pair;
#     fit_dose_slope() switches to OLS-on-log-dose when >=3 levels exist.
###############################################################################

suppressPackageStartupMessages({ library(jsonlite); library(dplyr) })

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
this_dir_conc <- function() {
  fa <- commandArgs(trailingOnly = FALSE)
  farg <- sub("^--file=", "", fa[grep("^--file=", fa)][1])
  dirname(normalizePath(farg, mustWork = FALSE))
}

read_conc_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  cfg$paths$project_root <- normalizePath(cfg$paths$project_root, mustWork = TRUE)
  for (pk in c("browse_table", "results_dir", "cck8_dir", "integration_dir")) {
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
  # dose levels: ordered list; validate
  doses <- cfg$doses
  stopifnot(length(doses) >= 2)
  for (d in doses) {
    stopifnot(!is.null(d$name), !is.null(d$mg_ml), length(d$cck8_cols) >= 1, !is.null(d$ms_dose))
  }
  cfg$doses_logconc <- vapply(doses, function(d) log2(as.numeric(d$mg_ml)), numeric(1))
  # KEGG cache paths (relative to this file's dir)
  for (ck in c("kegg_link_cache", "kegg_name_cache")) {
    if (!is.null(cfg$pathway[[ck]]) && !startsWith(cfg$pathway[[ck]], "/")) {
      cfg$pathway[[ck]] <- file.path(this_dir_conc(), cfg$pathway[[ck]])
    }
  }
  cfg
}

# CLI parsing shared by the runner and every step script: --config <path> and
# --steps <csv> are recognized anywhere in argv (not only as the first flag).
conc_cli_args <- function(steps_default = "0,1,2,3,4,5,6") {
  args <- commandArgs(trailingOnly = TRUE)
  pick <- function(flag) {
    i <- which(args == flag)
    if (length(i) >= 1 && length(args) > i[1]) args[i[1] + 1] else NULL
  }
  steps_arg <- pick("--steps") %||% steps_default
  list(config = pick("--config"),
       steps = as.integer(strsplit(steps_arg, ",")[[1]]),
       args = args)
}

# ---------------------------------------------------------------------------
# Statistics (same conventions as cck8_utils.R)
# ---------------------------------------------------------------------------
spearman_pair <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(unique(x)) < 2 || length(unique(y)) < 2 || length(x) < 5) {
    return(c(rho = NA_real_, p = NA_real_))
  }
  st <- tryCatch(stats::cor.test(x, y, method = "spearman", exact = FALSE),
                 error = function(e) NULL)
  if (is.null(st)) return(c(rho = NA_real_, p = NA_real_))
  c(rho = as.numeric(st$estimate), p = as.numeric(st$p.value))
}

mad_scale <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  1.4826 * stats::mad(x, na.rm = TRUE)
}

# ---------------------------------------------------------------------------
# Dose-response primitives (extensible to >=3 concentrations)
# ---------------------------------------------------------------------------
# Per-herb phenotype slope from replicate values at each dose level.
# x_by_dose: named numeric vector of per-dose means; sd_by_dose / n_by_dose:
# per-dose sample sd and replicate count (for the SE of the difference).
# 2 levels: simple difference with SE = sqrt(sd1^2/n1 + sd2^2/n2) (plan §2.1);
# z = slope/SE is the two-point test statistic (t-like, df ~ n1+n2-2).
# >=3 levels: OLS slope on log2(dose), SE from the fit (A5 extensibility).
fit_dose_slope <- function(x_by_dose, sd_by_dose = NULL, n_by_dose = NULL, doses_cfg) {
  dnames <- vapply(doses_cfg, function(d) d$name, character(1))
  means <- as.numeric(x_by_dose[match(dnames, names(x_by_dose))])
  if (any(!is.finite(means))) return(c(slope = NA_real_, z = NA_real_, se = NA_real_))
  logc <- vapply(doses_cfg, function(d) log2(as.numeric(d$mg_ml)), numeric(1))
  if (length(logc) == 2) {
    slope <- means[2] - means[1]
    sds <- as.numeric(sd_by_dose[match(dnames, names(sd_by_dose))])
    ns  <- as.integer(n_by_dose[match(dnames, names(n_by_dose))])
    se <- NA_real_
    if (all(is.finite(sds)) && all(!is.na(ns)) && all(ns >= 2)) {
      se <- sqrt(sum(sds^2 / ns))
    }
    z <- if (is.finite(se) && se > 0) slope / se else NA_real_
  } else {
    fit <- tryCatch(stats::lm(means ~ logc), error = function(e) NULL)
    if (is.null(fit)) return(c(slope = NA_real_, z = NA_real_, se = NA_real_))
    co <- stats::coef(fit); se <- summary(fit)$coefficients[2, 2]
    slope <- co[2]; z <- if (se > 0) co[2] / se else NA_real_
  }
  c(slope = slope, z = z, se = se)
}

# Per-herb metabolic gradient for one dose pair (plan §2.2).
pairwise_dose_gradient <- function(mat_high, mat_low) {
  # both: numeric vectors over feature groups; NA where either side missing
  g <- mat_high - mat_low
  g[!is.finite(g)] <- NA_real_
  g
}

# Global robust z-standardization of a gradient vector (plan §2.2).
robust_z <- function(x) {
  med <- stats::median(x, na.rm = TRUE)
  sc <- mad_scale(x)
  if (!is.finite(sc) || sc == 0) return(rep(NA_real_, length(x)))
  z <- (x - med) / sc
  z[!is.finite(z)] <- NA_real_
  z
}

# Bootstrap CI for a pooled Spearman correlation (plan §2.8: n=1000, seed=42,
# resample herbs with replacement). x: groups x herbs matrix; y: per-herb vector.
bootstrap_spearman_ci <- function(x, y, B = 1000, seed = 42) {
  set.seed(seed)
  nh <- ncol(x)
  rhos <- rep(NA_real_, B)   # NA (not 0) for draws where the test fails
  for (b in seq_len(B)) {
    idx <- sample.int(nh, replace = TRUE)
    zres <- matrix(y[idx], nrow = nrow(x), ncol = nh, byrow = TRUE)
    ok <- is.finite(x[, idx, drop = FALSE]) & is.finite(zres)
    xs <- x[ok]; ys <- zres[ok]
    if (length(xs) < 10 || length(unique(xs)) < 2 || length(unique(ys)) < 2) next
    st <- tryCatch(stats::cor.test(xs, ys, method = "spearman", exact = FALSE),
                   error = function(e) NULL)
    if (!is.null(st)) rhos[b] <- as.numeric(st$estimate)
  }
  rhos <- rhos[is.finite(rhos)]
  if (length(rhos) < 10) return(c(rho = NA_real_, lo = NA_real_, hi = NA_real_))
  q <- stats::quantile(rhos, c(0.025, 0.975), names = FALSE)
  c(rho = median(rhos), lo = q[1], hi = q[2])
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
  root_mzs <- setNames(root_stats$mz, as.character(root_stats$root_id))
  root_rts <- setNames(root_stats$rt, as.character(root_stats$root_id))

  df$group_id <- sprintf("G%05d_z%.3f_rt%d", df$root_id,
                         root_mzs[as.character(df$root_id)],
                         round(root_rts[as.character(df$root_id)]))
  df$canonical_mz <- root_mzs[as.character(df$root_id)]
  df$canonical_rt <- round(root_rts[as.character(df$root_id)])

  orig_order <- match(seq_len(n), ord)
  df[orig_order, ]
}

# ---------------------------------------------------------------------------
# Group annotation: one row per integrated group with inherited member info.
# Handles both multi-member groups (G#####_z..._rt...) and single-member
# fallback keys (M..._NEG etc.) via the cross_batch_rows.json mapping.
# Returns data.frame(group_id, canonical_mz, canonical_rt, compound_name,
#                    formula, hmdb_id, kegg_id, n_members)
# ---------------------------------------------------------------------------
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

  # best-annotated member per group (non-empty name preferred)
  has_name <- ifelse(aligned$compound_name != "", 0L, 1L)
  ord <- order(aligned$root_id, has_name, aligned$variable_id)
  grp <- aligned[ord, ]
  grp <- grp[!duplicated(grp$group_id), , drop = FALSE]

  # true member count per group (all features aligned into it across batches);
  # bare-key rows below are single-member fallbacks.
  nm_tab <- table(aligned$group_id)
  grp$n_members <- as.integer(nm_tab[as.character(grp$group_id)])

  # single-member groups: their key in the matrix is the variable_id itself.
  # Map them via cross_batch_rows.json (key -> variable_id), then look up the
  # member feature row (first match; a variable_id can appear in several batches).
  rows_j <- readj(file.path(cfg$paths$integration_dir, "cross_batch_rows.json"))
  if (!is.null(rows_j) && !is.null(rows_j$rows)) {
    rws <- do.call(rbind, lapply(rows_j$rows, function(x) {
      data.frame(key = x$key, variable_id = x$variable_id, stringsAsFactors = FALSE)
    }))
    single <- rws[!grepl("^G", rws$key), , drop = FALSE]
    if (nrow(single) > 0) {
      feat_idx <- split(seq_len(nrow(combined)), combined$variable_id)
      single_rows <- lapply(seq_len(nrow(single)), function(i) {
        idx <- feat_idx[[single$variable_id[i]]]
        if (is.null(idx)) return(NULL)
        h <- combined[idx[1], ]
        data.frame(group_id = single$key[i], canonical_mz = h$mz, canonical_rt = round(h$rt),
                   compound_name = h$compound_name, formula = h$formula, hmdb_id = h$hmdb_id,
                   kegg_id = h$kegg_id, n_members = 1L, stringsAsFactors = FALSE)
      })
      single_rows <- do.call(rbind, single_rows[!vapply(single_rows, is.null, logical(1))])
      grp <- rbind(grp[, c("group_id", "canonical_mz", "canonical_rt", "compound_name",
                           "formula", "hmdb_id", "kegg_id", "n_members")], single_rows)
    }
  }

  rownames(grp) <- NULL
  grp <- grp[, c("group_id", "canonical_mz", "canonical_rt", "compound_name",
                 "formula", "hmdb_id", "kegg_id", "n_members"), drop = FALSE]
  # dedupe (a single-member key could collide in theory)
  grp <- grp[!duplicated(grp$group_id), , drop = FALSE]
  grp
}

# variable_id -> integrated-matrix row key, per batch.
# Uses the GLOBAL alignment (same call as build_group_annotation / the
# integration module): every aligned feature maps to its group G-key. Local
# per-batch alignment must NOT be used here — its G-keys are numbered within
# the batch and do not match the matrix. A few matrix rows are keyed by a bare
# variable_id (features absent from metabolites.json); for those, fall back to
# the variable_id itself.
build_feature_key_map <- function(cfg) {
  feats_list <- list()
  for (bid in cfg$integration$batches) {
    web <- file.path(cfg$paths$results_dir, paste0("Batch", bid), "detail", "web_export")
    m <- readj(file.path(web, "metabolites.json"))
    if (is.null(m)) next
    f <- features_of(m$metabolites, bid)
    if (!is.null(f) && nrow(f) > 0) feats_list[[as.character(bid)]] <- f
  }
  combined <- do.call(rbind, feats_list)
  aligned <- align_features(combined, mz_ppm = cfg$integration$mz_ppm, rt_tol = cfg$integration$rt_tol)
  split(setNames(aligned$group_id, aligned$variable_id), as.character(aligned$batch_id))
}

# ---------------------------------------------------------------------------
# KEGG compound -> pathway map from cached bulk table (plan §2.6)
# ---------------------------------------------------------------------------
load_kegg_link_table <- function(cfg, verbose = TRUE) {
  f <- cfg$pathway$kegg_link_cache
  if (!file.exists(f)) {
    logmsg("KEGG link cache missing: ", f)
    logmsg("Fetching https://rest.kegg.jp/link/pathway/cpd (one-time, may take minutes)...")
    dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
    ok <- tryCatch({
      con <- url("https://rest.kegg.jp/link/pathway/cpd", open = "r", blocking = TRUE)
      on.exit(close(con))
      n <- readLines(con, warn = FALSE)
      cat(n, file = f, sep = "\n")
      length(n) > 1000
    }, error = function(e) { logmsg("  download failed:", e$message); FALSE })
    if (!ok) stop("Could not fetch KEGG compound-pathway table; re-run with network access.")
  }
  tab <- read.delim(f, header = FALSE, col.names = c("cpd", "pathway"),
                    stringsAsFactors = FALSE, check.names = FALSE)
  if (verbose) logmsg("KEGG link table: ", nrow(tab), " links,", length(unique(tab$pathway)), " pathways")
  tab
}

# ---------------------------------------------------------------------------
# Integrated logFC matrix loading + comparison parsing
# ---------------------------------------------------------------------------
load_logfc_matrix <- function(cfg, doses_cfg = NULL) {
  mat_path <- file.path(cfg$paths$integration_dir, "integrated_logFC_matrix.csv")
  logmsg("Reading ", mat_path)
  M <- read.csv(mat_path, row.names = 1, check.names = FALSE)
  M <- as.matrix(M); storage.mode(M) <- "double"

  # columns: B{batch}:{drug}_{Dose}_vs_CT{n}
  doses_cfg <- doses_cfg %||% cfg$doses
  dose_pat <- paste(vapply(doses_cfg, function(d) d$ms_dose, character(1)), collapse = "|")
  pat <- sprintf("^B([0-9]+):([0-9]+)_(%s)_vs_CT[0-9]+$", dose_pat)
  keep <- vapply(colnames(M), function(cn) {
    m <- regmatches(cn, regexec(pat, cn))[[1]]
    length(m) == 4 && as.integer(m[2]) %in% cfg$integration$batches
  }, logical(1))
  cols <- colnames(M)[keep]
  cmp_info <- do.call(rbind, lapply(cols, function(cn) {
    m <- regmatches(cn, regexec(pat, cn))[[1]]
    data.frame(batch_id = as.integer(m[2]), drug_id = as.integer(m[3]),
               dose = m[4], colname = cn, stringsAsFactors = FALSE)
  }))
  rownames(cmp_info) <- NULL
  list(matrix = M[, keep, drop = FALSE], cmp_info = cmp_info)
}

# ---------------------------------------------------------------------------
# Plotting / reporting helpers (same as cck8 module)
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
  # readable precision for reports (CSV outputs keep full precision);
  # integers and non-finite values pass through unchanged
  for (cn in names(df)) {
    if (is.numeric(df[[cn]])) {
      v <- df[[cn]]
      df[[cn]] <- ifelse(is.finite(v) & v != round(v), signif(v, 4), v)
    }
  }
  hdr <- paste0("| ", paste(colnames(df), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  rows <- apply(df, 1, function(r) paste0("| ", paste(vapply(r, as.character, character(1)), collapse = " | "), " |"))
  c(hdr, sep, rows)
}

# Per-category 2x2 Fisher enrichment of selected classes vs the rest.
# tab: categories x classes contingency table (from table()). For each
# (class, category): exact two-sided Fisher test of
#   [in class & category] vs [outside class & category], given the class
#   margin (proper per-category test — an omnibus r x 2 table yields ONE
#   p-value for the whole table, not one per category).
# odds_ratio = conditional MLE (fisher.test estimate on the 2x2);
# BH-FDR across ALL rows of the returned table.
fisher_category_enrichment <- function(tab, classes, class_label = "class",
                                       category_label = "category", count_label = "n_in_class") {
  rows <- list()
  for (cl in classes) {
    if (!cl %in% colnames(tab)) next
    a <- as.integer(tab[, cl]); names(a) <- rownames(tab)
    rest <- rowSums(tab) - a
    keep <- (a + rest) > 0 & a > 0
    if (!any(keep)) next
    a <- a[keep]; rest <- rest[keep]
    n_cl <- sum(a); n_out <- sum(rest)
    for (k in seq_along(a)) {
      m <- matrix(c(a[k], n_cl - a[k], rest[k], n_out - rest[k]), nrow = 2,
                  dimnames = list(category = c(names(a)[k], "other"),
                                  set = c("in_class", "out_class")))
      ft <- tryCatch(stats::fisher.test(m), error = function(e) NULL)
      if (is.null(ft)) next
      rows[[length(rows) + 1]] <- data.frame(
        class = cl, category = names(a)[k],
        n_in = unname(a[k]), n_total = unname(a[k] + rest[k]),
        p_fisher = as.numeric(ft$p.value),
        odds_ratio = unname(as.numeric(ft$estimate)),
        stringsAsFactors = FALSE)
    }
  }
  enr <- do.call(rbind, rows)
  if (is.null(enr)) return(NULL)
  enr$fdr <- as.numeric(stats::p.adjust(enr$p_fisher, method = "BH"))
  enr <- enr[order(enr$p_fisher), c("class", "category", "n_in", "n_total",
                                    "p_fisher", "fdr", "odds_ratio"), drop = FALSE]
  colnames(enr) <- c(class_label, category_label, count_label, "n_total", "p_fisher", "fdr", "odds_ratio")
  rownames(enr) <- NULL
  enr
}

write_report <- function(lines, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  writeLines(lines, path)
  logmsg("Report: ", path)
}
