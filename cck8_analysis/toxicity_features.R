#!/usr/bin/env Rscript
###############################################################################
# toxicity_features.R — Step 3 (core): feature-group level association between
# integrated cross-batch logFC and CCK8 viability.
#   A) Spearman per group vs high_mean / dose_delta, BH-FDR
#   B) two-group contrast: toxic vs neutral (Wilcoxon on group logFC)
#   C) overlap of A/B candidate sets with per-herb significant features
# Writes data/toxicity_features.csv, reports/toxicity_features_report.md,
# figures/toxicity_volcano_{high,delta}.
###############################################################################
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

this_dir <- {
  fa <- commandArgs(trailingOnly = FALSE)
  farg <- sub("^--file=", "", fa[grep("^--file=", fa)][1])
  dirname(normalizePath(farg, mustWork = FALSE))
}
source(file.path(this_dir, "cck8_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1 && args[[1]] == "--config") args[[2]] else "cck8_config.yaml"
cfg <- read_cck8_config(cfg_path)
out_dir <- cfg$paths$output_dir

# --- load integrated logFC matrix ---------------------------------------------
mat_path <- file.path(cfg$paths$integration_dir, "integrated_logFC_matrix.csv")
logmsg("Reading ", mat_path)
M <- read.csv(mat_path, row.names = 1, check.names = FALSE)
M <- as.matrix(M); storage.mode(M) <- "double"

cmp_info <- parse_comparison_columns(colnames(M), cfg$integration$batches)
stopifnot(nrow(cmp_info) == ncol(M))
logmsg("Comparisons parsed: ", nrow(cmp_info))

# --- herb viability ------------------------------------------------------------
norm <- read.csv(file.path(out_dir, "data", "toxicity_class.csv"), stringsAsFactors = FALSE)
herbs <- norm %>% filter(batch_id %in% cfg$integration$batches)
cmp_info$match_id <- cmp_info$drug_id   # drug_id == id_match (verified for all 5 batches)
cmp_info$high_mean <- herbs$high_mean[match(cmp_info$match_id, herbs$match_id)]
cmp_info$dose_delta <- herbs$dose_delta[match(cmp_info$match_id, herbs$match_id)]
cmp_info$toxicity_class <- herbs$toxicity_class[match(cmp_info$match_id, herbs$match_id)]

# herb-level vectors: one value per comparison column
high_by_col <- cmp_info$high_mean
delta_by_col <- cmp_info$dose_delta

# --- A) Spearman per group ------------------------------------------------------
logmsg("Per-group Spearman vs viability...")
n_groups <- nrow(M)
rho_high <- numeric(n_groups); p_high <- numeric(n_groups)
rho_delta <- numeric(n_groups); p_delta <- numeric(n_groups)
for (g in seq_len(n_groups)) {
  x <- M[g, ]
  rh <- spearman_pair(x, high_by_col)
  rho_high[g] <- rh["rho"]; p_high[g] <- rh["p"]
  rd <- spearman_pair(x, delta_by_col)
  rho_delta[g] <- rd["rho"]; p_delta[g] <- rd["p"]
}

fdr_high <- as.numeric(p.adjust(p_high, method = "BH"))
fdr_delta <- as.numeric(p.adjust(p_delta, method = "BH"))

cut_rho <- cfg$linkage$min_abs_rho; cut_fdr <- cfg$linkage$fdr_cutoff
res <- data.frame(
  group_id = rownames(M),
  rho_high = rho_high, p_high = p_high, fdr_high = fdr_high,
  rho_delta = rho_delta, p_delta = p_delta, fdr_delta = fdr_delta,
  sig_high = !is.na(fdr_high) & fdr_high < cut_fdr & abs(rho_high) >= cut_rho,
  sig_delta = !is.na(fdr_delta) & fdr_delta < cut_fdr & abs(rho_delta) >= cut_rho,
  stringsAsFactors = FALSE
)

# --- B) two-group contrast: toxic vs neutral ------------------------------------
toxic_cols <- cmp_info$colname[cmp_info$toxicity_class %in% c("toxic_high", "toxic_low")]
neutral_cols <- cmp_info$colname[cmp_info$toxicity_class == "neutral"]
logmsg("Two-group contrast: ", length(toxic_cols), " toxic vs ", length(neutral_cols), " neutral comparisons")

wilcox_rows <- lapply(seq_len(n_groups), function(g) {
  xt <- M[g, toxic_cols]; xn <- M[g, neutral_cols]
  xt <- xt[is.finite(xt)]; xn <- xn[is.finite(xn)]
  if (length(xt) < 3 || length(xn) < 3) return(c(p = NA_real_, n_t = length(xt), n_n = length(xn),
                                                 med_t = NA_real_, med_n = NA_real_))
  wt <- suppressWarnings(stats::wilcox.test(xt, xn))
  c(p = as.numeric(wt$p.value), n_t = length(xt), n_n = length(xn),
    med_t = median(xt), med_n = median(xn))
})
wres <- as.data.frame(do.call(rbind, wilcox_rows), stringsAsFactors = FALSE)
colnames(wres) <- c("p_wilcox", "n_toxic_obs", "n_neutral_obs", "median_logFC_toxic", "median_logFC_neutral")
res$p_wilcox <- wres$p_wilcox
res$median_logFC_toxic <- wres$median_logFC_toxic
res$median_logFC_neutral <- wres$median_logFC_neutral
res$fdr_wilcox <- as.numeric(p.adjust(res$p_wilcox, method = "BH"))
res$sig_two_group <- !is.na(res$fdr_wilcox) & res$fdr_wilcox < cut_fdr

# --- annotation ------------------------------------------------------------------
logmsg("Rebuilding group annotation (must match stored integration IDs)...")
grp_build <- build_group_annotation(cfg)
grp_anno <- grp_build$groups
aligned_all <- grp_build$aligned
matched_ids <- intersect(res$group_id, grp_anno$group_id)
logmsg("Group IDs matched to stored matrix: ", length(matched_ids), " / ", nrow(res))
res <- merge(res, grp_anno, by = "group_id", all.x = TRUE)

out_csv <- file.path(out_dir, "data", "toxicity_features.csv")
write.csv(res, out_csv, row.names = FALSE)
logmsg("Wrote ", out_csv)

# --- C) overlap with per-herb significant features --------------------------------
# For toxic herbs: their batch's volcano significant variable_ids -> group ids (deduplicated)
flagged_gids <- res$group_id[res$sig_high | res$sig_two_group]
toxic_sig_gids <- character(0)
sig_pts_total <- 0L; sig_pts_unmapped <- 0L
for (bid in cfg$integration$batches) {
  web <- file.path(cfg$paths$results_dir, paste0("Batch", bid), "detail", "web_export")
  vf <- list.files(web, pattern = "^volcano_[0-9]+_(High|Low)_vs_CT[0-9]+\\.json$", full.names = TRUE)
  for (f in vf) {
    cmp <- sub("^volcano_", "", sub("\\.json$", "", basename(f)))
    drug_id <- as.integer(sub("^([0-9]+)_(High|Low)_vs_CT[0-9]+$", "\\1", cmp))
    if (!(drug_id %in% herbs$match_id)) next
    cls <- herbs$toxicity_class[herbs$match_id == drug_id]
    if (!(cls %in% c("toxic_high", "toxic_low"))) next
    j <- readj(f)
    pts <- j$points
    ids <- vapply(pts, function(p) if (is.null(p$name)) "" else as.character(p$name), character(1))
    fd  <- vapply(pts, function(p) if (is.null(p$fdr)) NA_real_ else as.numeric(as.character(p$fdr)), numeric(1))
    sig_ids <- ids[!is.na(fd) & fd < cut_fdr]
    b_aligned <- aligned_all[aligned_all$batch_id == as.integer(bid), ]
    if (nrow(b_aligned) == 0) next
    v2g <- setNames(b_aligned$group_id, b_aligned$variable_id)
    sig_pts_total <- sig_pts_total + length(sig_ids)
    sig_pts_unmapped <- sig_pts_unmapped + sum(is.na(v2g[sig_ids]))
    toxic_sig_gids <- c(toxic_sig_gids, unique(v2g[sig_ids]))
  }
}
toxic_sig_gids <- unique(toxic_sig_gids)
overlap_hits <- sum(toxic_sig_gids %in% flagged_gids)
# chance baseline: prevalence of flagged groups among all annotated groups
baseline_rate <- length(flagged_gids) / max(length(unique(res$group_id)), 1)


# --- report ------------------------------------------------------------------------
top_high <- res %>% filter(sig_high) %>% arrange(fdr_high) %>% head(30)
top_delta <- res %>% filter(sig_delta) %>% arrange(fdr_delta) %>% head(30)
top_2g <- res %>% filter(sig_two_group) %>% arrange(fdr_wilcox) %>% head(30)

lines <- c(
  "# Toxicity-Associated Metabolic Features (Integrated Groups)", "",
  paste0("_Generated: ", fmt_time(), "_"), "",
  sprintf("Universe: **%d** integrated feature groups x %d herb-dose comparisons (batches %s; batch 4 excluded).",
          n_groups, nrow(cmp_info), paste(cfg$integration$batches, collapse = "/")), "",
  "## A) Per-group Spearman association with CCK8 viability", "",
  sprintf("Cutoff: FDR < %.2f AND |rho| >= %.2f.", cut_fdr, cut_rho), "",
  sprintf("- vs **high-dose mean viability**: %d groups significant (top 30 below)", sum(res$sig_high)),
  ""
)
if (nrow(top_high) > 0) lines <- c(lines, md_table(top_high[, c("group_id","compound_name","formula","rho_high","fdr_high")]))
lines <- c(lines, "", sprintf("- vs **dose delta** (high - low viability): %d groups significant", sum(res$sig_delta)), "")
if (nrow(top_delta) > 0) lines <- c(lines, md_table(top_delta[, c("group_id","compound_name","formula","rho_delta","fdr_delta")]))

lines <- c(lines, "", "## B) Two-group contrast: toxic vs neutral", "",
  sprintf("Wilcoxon test per group on logFC across %d toxic vs %d neutral comparisons; BH-FDR < %.2f.",
          length(toxic_cols), length(neutral_cols), cut_fdr), "",
  sprintf("- **%d groups** significant (top 30 below)", sum(res$sig_two_group)), "")
if (nrow(top_2g) > 0) lines <- c(lines, md_table(top_2g[, c("group_id","compound_name","formula","median_logFC_toxic","median_logFC_neutral","fdr_wilcox")]))

lines <- c(lines, "", "## C) Robustness: overlap with per-herb significant features", "",
  sprintf("Across all toxic herbs' comparisons (deduplicated), %d unique significant feature groups; **%d (%.1f%%)** are also flagged by A or B.",
          length(toxic_sig_gids), overlap_hits, 100 * overlap_hits / max(length(toxic_sig_gids), 1)),
  sprintf("Chance baseline: %.1f%% of all annotated groups are flagged by A or B; an overlap well above this rate indicates the two independent screens agree.",
          100 * baseline_rate), "",
  "## Caveats", "",
  sprintf("- %d of %d matrix rows could be annotated (the remainder are single-batch features without mz/RT, i.e. not part of any aligned group).", length(matched_ids), n_groups),
  sprintf("- n = %d herb-dose comparisons (%d Low + %d High); per-group correlations are screening-level evidence.",
          nrow(cmp_info), sum(cmp_info$concentration == "Low"), sum(cmp_info$concentration == "High")),
  "- Groups present in only one batch contribute a single observation; cross-batch shared groups are the more reliable ones.",
  sprintf("- Step C overlap: %d of %d toxic-herb FDR-significant features (%.1f%%) have no variable_id in the batch metabolites.json and thus no group mapping; they are excluded from the deduplicated count (slight undercount).",
          sig_pts_unmapped, sig_pts_total, 100 * sig_pts_unmapped / max(sig_pts_total, 1L)),
  "- Direction: negative rho vs high_mean / positive rho vs dose_delta indicates 'more toxic herbs show lower/higher logFC'.")

dir.create(file.path(out_dir, "reports"), showWarnings = FALSE, recursive = TRUE)
writeLines(lines, file.path(out_dir, "reports", "toxicity_features_report.md"))

# --- volcano figures ------------------------------------------------------------------
vol_plot <- function(rho, fdr, title, base) {
  d <- res %>% mutate(sig = sig_high | sig_delta | sig_two_group)
  p <- ggplot(d, aes(x = rho, y = -log10(pmax(fdr, 1e-300)), color = sig)) +
    geom_point(alpha = 0.5, size = 0.9) +
    geom_vline(xintercept = c(-cut_rho, cut_rho), linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("TRUE" = "#d62728", "FALSE" = "grey70")) +
    labs(x = sprintf("Spearman rho (group logFC vs %s)", title), y = "-log10(FDR)",
         title = sprintf("Toxicity-associated integrated feature groups (n=%d)", nrow(d))) +
    theme_minimal(base_size = 13) + guides(color = "none")
  save_plot(p, file.path(out_dir, "figures", base))
}
vol_plot(rho_high, fdr_high, "high-dose viability", "toxicity_volcano_high")
vol_plot(rho_delta, fdr_delta, "dose delta (high-low)", "toxicity_volcano_delta")
logmsg("Done: toxicity features report + volcanos")
