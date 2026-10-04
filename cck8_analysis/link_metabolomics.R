#!/usr/bin/env Rscript
###############################################################################
# link_metabolomics.R — Step 2: herb-level association between CCK8 viability
# and metabolomic perturbation magnitude (per comparison, from stored volcano
# + expression JSONs). Writes data/metabolomics_metrics.csv,
# data/herb_linkage.csv, reports/linkage_report.md, figures/linkage_*.
###############################################################################
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

this_dir <- {
  fa <- commandArgs(trailingOnly = FALSE)
  farg <- sub("^--file=", "", fa[grep("^--file=", fa)][1])
  dirname(normalizePath(farg, mustWork = FALSE))
}
source(file.path(this_dir, "cck8_utils.R"))

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a
fmt_time <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
logmsg <- function(...) cat(sprintf("[%s] ", fmt_time()), ..., "\n", sep = "")

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1 && args[[1]] == "--config") args[[2]] else "cck8_config.yaml"
cfg <- read_cck8_config(cfg_path)
out_dir <- cfg$paths$output_dir

norm <- read.csv(file.path(out_dir, "data", "toxicity_class.csv"), stringsAsFactors = FALSE)
valid_batches <- cfg$integration$batches

# --- drug -> batch map from per-batch file lists -----------------------------
drug_batch <- do.call(rbind, lapply(valid_batches, function(bid) {
  f <- file.path(cfg$paths$file_lists_dir,
                 sprintf("mzml_file_list_Batch%d_convert.csv", bid))
  df <- read.csv(f, stringsAsFactors = FALSE)
  m <- df[df$drug_id %in% as.character(1:999), c("drug_id"), drop = FALSE]
  m$batch_id <- bid
  m$drug_id <- as.integer(m$drug_id)
  m[!duplicated(m), ]
}))

# --- per-comparison metrics ---------------------------------------------------
parse_volcano_metrics <- function(j, drug_id) {
  pts <- j$points
  id  <- vapply(pts, function(p) if (is.null(p$name)) "" else as.character(p$name), character(1))
  lfc <- vapply(pts, function(p) if (is.null(p$logFC)) NA_real_ else as.numeric(as.character(p$logFC)), numeric(1))
  pv  <- vapply(pts, function(p) if (is.null(p$p_value)) NA_real_ else as.numeric(as.character(p$p_value)), numeric(1))
  fd  <- vapply(pts, function(p) if (is.null(p$fdr)) NA_real_ else as.numeric(as.character(p$fdr)), numeric(1))
  sig <- !is.na(fd) & fd < cfg$linkage$fdr_cutoff
  data.frame(
    batch_id = drug_batch$batch_id[drug_batch$drug_id == drug_id][1],
    drug_id = drug_id,
    comparison = j$comparison,
    concentration = sub(".*_(High|Low)_vs_CT[0-9]+$", "\\1", j$comparison),
    n_total = length(id),
    n_significant = sum(sig),
    n_up = sum(sig & lfc > 0),
    n_down = sum(sig & lfc < 0),
    median_abs_logFC_sig = if (sum(sig) > 0) median(abs(lfc[sig])) else NA_real_,
    max_neglog10p = -log10(min(pv, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}

# global perturbation distance: mean |z(treat) - z(ctrl)| across features
perturbation_distance <- function(j, drug_id, concentration) {
  if (is.null(j) || is.null(j$data) || is.null(j$samples)) return(NA_real_)
  samples <- as.character(unlist(j$samples))
  dose_sfx <- if (concentration == "High") "H" else "L"
  treat <- samples[grepl(paste0("^", drug_id, "_[0-9]+_", dose_sfx, "$"), samples)]
  ctrl  <- samples[grepl("^CT[0-9]+_", samples)]
  if (length(treat) == 0 || length(ctrl) == 0) return(NA_real_)
  expected_n <- length(treat) + length(ctrl)
  vals_list <- lapply(j$data, function(row) {
    v <- unlist(row$values)[c(treat, ctrl)]
    as.numeric(v)
  })
  if (any(vapply(vals_list, length, integer(1)) != expected_n)) return(NA_real_)
  mat <- do.call(rbind, vals_list)   # features x samples
  mat <- mat[stats::complete.cases(mat), , drop = FALSE]
  if (nrow(mat) < 10 || ncol(mat) != expected_n) return(NA_real_)
  z <- scale(mat)
  treat_idx <- seq_len(length(treat))
  ctrl_idx <- (length(treat) + 1L):ncol(z)
  d <- abs(rowMeans(z[, treat_idx, drop = FALSE]) - rowMeans(z[, ctrl_idx, drop = FALSE]))
  mean(d, na.rm = TRUE)
}

logmsg("Scanning volcano + expression JSONs for batches ", paste(valid_batches, collapse = ","))
metric_rows <- list()
for (bid in valid_batches) {
  web <- file.path(cfg$paths$results_dir, paste0("Batch", bid), "detail", "web_export")
  vf <- list.files(web, pattern = "^volcano_[0-9]+_(High|Low)_vs_CT[0-9]+\\.json$", full.names = TRUE)
  for (f in vf) {
    cmp <- sub("^volcano_", "", sub("\\.json$", "", basename(f)))
    drug_id <- as.integer(sub("^([0-9]+)_(High|Low)_vs_CT[0-9]+$", "\\1", cmp))
    j <- readj(f)
    if (is.null(j)) next
    row <- parse_volcano_metrics(j, drug_id)
    ej <- readj(file.path(web, sprintf("expression_%s.json", cmp)))
    row$perturbation_distance <- perturbation_distance(ej, drug_id, row$concentration)
    metric_rows[[length(metric_rows) + 1]] <- row
  }
}
metrics <- do.call(rbind, metric_rows)
write.csv(metrics, file.path(out_dir, "data", "metabolomics_metrics.csv"), row.names = FALSE)
logmsg("Comparison-level metrics: ", nrow(metrics), " rows")

# --- herb x dose linkage table -------------------------------------------------
linkage <- norm %>%
  select(match_id, drug_name_zh, batch_id, toxicity_class, low_mean, high_mean,
         dose_delta, no_quant_ms, qc_flag) %>%
  filter(batch_id %in% valid_batches)

link_low <- dplyr::left_join(
  linkage %>% dplyr::mutate(drug_id = match_id, concentration = "Low", viability = low_mean),
  metrics, by = c("drug_id", "concentration")
)
link_high <- dplyr::left_join(
  linkage %>% dplyr::mutate(drug_id = match_id, concentration = "High", viability = high_mean),
  metrics, by = c("drug_id", "concentration")
)
herb_linkage <- rbind(link_low, link_high)
write.csv(herb_linkage, file.path(out_dir, "data", "herb_linkage.csv"), row.names = FALSE)

# --- Spearman correlations ------------------------------------------------------
n_boot <- cfg$linkage$bootstrap_n; seed <- cfg$linkage$seed
cor_vars <- c("n_significant", "median_abs_logFC_sig", "max_neglog10p", "perturbation_distance")
res_rows <- lapply(cor_vars, function(v) {
  r_low  <- spearman_pair(herb_linkage$viability[herb_linkage$concentration == "Low"],
                          herb_linkage[[v]][herb_linkage$concentration == "Low"])
  ci_low <- bootstrap_spearman_ci(herb_linkage$viability[herb_linkage$concentration == "Low"],
                                  herb_linkage[[v]][herb_linkage$concentration == "Low"], n_boot, seed)
  r_high <- spearman_pair(herb_linkage$viability[herb_linkage$concentration == "High"],
                          herb_linkage[[v]][herb_linkage$concentration == "High"])
  ci_high <- bootstrap_spearman_ci(herb_linkage$viability[herb_linkage$concentration == "High"],
                                   herb_linkage[[v]][herb_linkage$concentration == "High"], n_boot, seed)
  data.frame(
    metric = v,
    low_rho = r_low["rho"], low_p = r_low["p"],
    low_ci_lo = ci_low[1], low_ci_hi = ci_low[2],
    high_rho = r_high["rho"], high_p = r_high["p"],
    high_ci_lo = ci_high[1], high_ci_hi = ci_high[2],
    n = r_low["n"]
  )
})
cor_res <- do.call(rbind, res_rows)
write.csv(cor_res, file.path(out_dir, "data", "linkage_correlations.csv"), row.names = FALSE)

# --- report ---------------------------------------------------------------------
fmt_rho <- function(x) if (is.na(x)) "NA" else sprintf("%.3f", x)
fmt_p <- function(x) if (is.na(x)) "NA" else format.pval(x, digits = 2)
lines <- c(
  "# CCK8 x Metabolomics Perturbation Linkage", "",
  paste0("_Generated: ", fmt_time(), "_"), "",
  sprintf("n = %d herbs with quantitative MS data (batch 4 excluded). Each herb contributes one Low and one High dose point.", nrow(linkage)), "",
  "Metrics per comparison (from stored volcano / expression JSONs):", "",
  "- `n_significant`: FDR < 0.05 features vs batch control (CT*)",
  "- `median_abs_logFC_sig`: median |log2FC| among significant features",
  "- `max_neglog10p`: strongest single-feature p-value",
  "- `perturbation_distance`: mean |z(treated) - z(control)| across all features (global profile shift)", "",
  "Note: OPLS-DA results are unavailable in the stored exports (`available: false` for all comparisons), so R2Y/Q2 are not used.", "",
  "## Spearman correlation: viability vs perturbation magnitude", ""
)
lines <- c(lines, md_table(cor_res[, c("metric","low_rho","low_p","low_ci_lo","low_ci_hi","high_rho","high_p","high_ci_lo","high_ci_hi","n")],
                           c("metric","rho (Low)","p (Low)","CI lo","CI hi","rho (High)","p (High)","CI lo","CI hi","n")), "",
  "Interpretation: a **negative** rho means herbs that perturb the metabolic profile more strongly tend to have lower viability (more toxic).", "")

dir.create(file.path(out_dir, "reports"), showWarnings = FALSE, recursive = TRUE)
writeLines(lines, file.path(out_dir, "reports", "linkage_report.md"))

# --- figures ---------------------------------------------------------------------
make_scatter <- function(metric, label) {
  d <- herb_linkage %>%
    filter(!is.na(.data[[metric]])) %>%
    mutate(cls = factor(toxicity_class, levels = c("toxic_low","toxic_high","neutral","proliferative")))
  p <- ggplot(d, aes(x = .data[[metric]], y = viability, color = cls)) +
    geom_point(alpha = 0.8, size = 1.7) +
    geom_smooth(method = "loess", se = TRUE, color = "grey35", fill = "grey70") +
    scale_color_brewer(palette = "Set2") +
    labs(x = label, y = "Mean CCK8 viability (control-normalized)",
         title = sprintf("Viability vs %s (n=%d herb-dose points)", label, nrow(d))) +
    theme_minimal(base_size = 13)
  save_plot(p, file.path(out_dir, "figures", sprintf("linkage_%s", metric)))
}
make_scatter("n_significant", "N significant features (FDR<0.05)")
make_scatter("perturbation_distance", "Global perturbation distance")

# dose-response within herb: delta viability vs delta n_sig
delta_tab <- dplyr::inner_join(
  link_high %>% select(drug_id, high_mean, n_significant = n_significant),
  link_low %>% select(drug_id, low_mean, n_significant_low = n_significant),
  by = "drug_id") %>%
  mutate(d_viability = high_mean - low_mean, d_nsig = n_significant - n_significant_low)
r_delta <- spearman_pair(delta_tab$d_viability, delta_tab$d_nsig)
logmsg("Dose-delta Spearman (dViability vs dNsig): rho=", fmt_rho(r_delta["rho"]), " p=", fmt_p(r_delta["p"]))
lines2 <- c("", "## Dose contrast", "",
  sprintf("Within-herb change from Low to High: Spearman(dViability, dN_sig) = **%s** (p = %s, n=%d).",
          fmt_rho(r_delta["rho"]), fmt_p(r_delta["p"]), r_delta["n"]))
cat(lines2, sep = "\n")
con <- file(file.path(out_dir, "reports", "linkage_report.md"), open = "a")
writeLines(lines2, con); close(con)
logmsg("Done: linkage report + figures written")
