#!/usr/bin/env Rscript
###############################################################################
# feature_viability_linkage.R — Step 3 (Q2): concentration-sensitive markers.
# Per integrated feature group g: Spearman( grad(g,d), z_viab(d) ) across the
# herbs for which g is dose-measured; BH-FDR over all groups; direction
# consistency; robustness = overlap with per-herb significant features
# (volcano FDR < cutoff, both doses).
# Writes data/feature_dose_linkage.csv, figures/feature_dose_linkage_*.
###############################################################################
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
this_dir <- {
  fa <- commandArgs(trailingOnly = FALSE)
  farg <- sub("^--file=", "", fa[grep("^--file=", fa)][1])
  dirname(normalizePath(farg, mustWork = FALSE))
}
source(file.path(this_dir, "conc_utils.R"))

main <- function() {
  cfg_path <- conc_cli_args()$config %||% file.path(this_dir, "concentration_config.yaml")
  cfg <- read_conc_config(cfg_path)
  out_dir <- cfg$paths$output_dir

  norm <- read.csv(file.path(out_dir, "data", "conc_inputs.csv"), stringsAsFactors = FALSE)
  gr <- readRDS(file.path(out_dir, "data", "gradient_z.rds"))
  G <- gr$G; Z <- gr$Z
  herb_ids <- colnames(G)
  z_viab <- norm$z_viab[match(herb_ids, norm$match_id)]

  # Primary test: POOLED Spearman per group — pool all (herb, gradient) pairs of a
  # group across herbs. Each group contributes its own n_herbs_used (~15-80), so the
  # pooled test has far more power than a per-group mean-based test at this sparsity.
  # BH-FDR over all groups. (Per-group gradients are only measurable where both doses
  # of that herb measured that group's features — median 29 herbs.)
  n_groups <- nrow(G)
  logmsg("Pooled Spearman per group vs standardized viability slope (", n_groups, " groups)...")

  rho <- pval <- frac_same_sign <- n_used <- rep(NA_real_, n_groups)  # NA (not 0) for untested groups
  for (g in seq_len(n_groups)) {
    x <- G[g, ]; z <- z_viab
    ok <- is.finite(x) & is.finite(z)
    if (sum(ok) < 8) next
    sp <- spearman_pair(x[ok], z[ok])
    rho[g] <- sp["rho"]; pval[g] <- sp["p"]
    n_used[g] <- sum(ok)
    frac_same_sign[g] <- mean(sign(x[ok]) == sign(z[ok]))
  }

  fdr <- as.numeric(p.adjust(pval, method = "BH"))
  cut_rho <- cfg$metabolic$min_abs_rho
  cut_fdr <- cfg$metabolic$fdr_cutoff
  # direction_consistency is ALIGNED to the association direction (for negative
  # rho it is 1 - frac_same_sign): high = the group's gradients follow the
  # viability slope as the association predicts. frac_same_sign (raw, unsigned)
  # is kept alongside for transparency.
  dir_consistency <- ifelse(is.na(rho), NA_real_,
                     ifelse(rho < 0, 1 - frac_same_sign, frac_same_sign))
  res <- data.frame(
    group_id = rownames(G),
    rho_vs_zviab = round(rho, 5), p_value = pval, fdr = fdr,
    n_herbs_used = n_used,
    direction_consistency = round(dir_consistency, 4),
    frac_same_sign = round(frac_same_sign, 4),
    direction_consistent = !is.na(dir_consistency) & dir_consistency >= cfg$metabolic$direction_consistency_min,
    sig_dose_sensitive = !is.na(fdr) & fdr < cut_fdr & abs(rho) >= cut_rho,
    stringsAsFactors = FALSE
  )

  # Global pooled association + bootstrap CI (plan §2.8: n=1000, seed from config)
  zmat <- matrix(z_viab, nrow = nrow(G), ncol = ncol(G), byrow = TRUE)
  pool_mask <- is.finite(G) & is.finite(zmat)
  global_sp <- spearman_pair(G[pool_mask], zmat[pool_mask])
  boot_ci <- bootstrap_spearman_ci(G, z_viab, B = cfg$metabolic$bootstrap_n, seed = cfg$metabolic$seed)
  logmsg(sprintf("Global pooled Spearman rho=%.3f (p=%.3g); bootstrap 95%% CI [%.3f, %.3f] over %d pairs",
                 global_sp["rho"], global_sp["p"], boot_ci["lo"], boot_ci["hi"], sum(pool_mask)))

  # --- annotation (inherited from member features) ----------------------------------
  logmsg("Building group annotation...")
  grp_anno <- build_group_annotation(cfg, verbose = TRUE)
  res <- merge(res, grp_anno, by = "group_id", all.x = TRUE)

  out_csv <- file.path(out_dir, "data", "feature_dose_linkage.csv")
  write.csv(res, out_csv, row.names = FALSE)
  logmsg("Wrote ", out_csv, "; sig_dose_sensitive: ", sum(res$sig_dose_sensitive))

  # --- robustness: overlap with per-herb significant features --------------------------
  # Volcano points are mapped to integrated-matrix row keys via the GLOBAL
  # alignment (build_feature_key_map); rows keyed by a bare variable_id map to
  # themselves, mirroring cross_batch_integration.R step 5.
  logmsg("Overlap check vs per-herb significant features...")
  flagged <- res$group_id[res$sig_dose_sensitive]
  key_map <- build_feature_key_map(cfg)
  vol_fdr_cut <- cfg$metabolic$volcano_fdr_cutoff
  dose_alt <- paste(vapply(cfg$doses, function(d) d$ms_dose, character(1)), collapse = "|")
  herb_sig_gids <- character(0); sig_pts_total <- 0L; unmapped <- 0L
  for (bid in cfg$integration$batches) {
    v2g <- key_map[[as.character(bid)]]
    if (is.null(v2g)) next
    web <- file.path(cfg$paths$results_dir, paste0("Batch", bid), "detail", "web_export")
    vf <- list.files(web, pattern = sprintf("^volcano_[0-9]+_(%s)_vs_CT[0-9]+\\.json$", dose_alt), full.names = TRUE)
    for (f in vf) {
      cmp <- sub("^volcano_", "", sub("\\.json$", "", basename(f)))
      drug_id <- as.integer(sub(sprintf("^([0-9]+)_(%s)_vs_CT[0-9]+$", dose_alt), "\\1", cmp))
      if (!(drug_id %in% herb_ids)) next
      j <- readj(f)
      if (is.null(j) || is.null(j$points)) next
      ids <- vapply(j$points, function(p) if (is.null(p$name)) "" else as.character(p$name), character(1))
      # fdr can be the literal string "NA" (upstream writer) — coercion then
      # yields NA intentionally; suppress the per-point coercion warning
      fd  <- vapply(j$points, function(p) if (is.null(p$fdr)) NA_real_ else suppressWarnings(as.numeric(as.character(p$fdr))), numeric(1))
      sig_ids <- ids[!is.na(fd) & fd < vol_fdr_cut]
      gids <- v2g[sig_ids]
      unmapped <- unmapped + sum(is.na(gids))
      gids[is.na(gids)] <- sig_ids[is.na(gids)]   # bare-variable_id matrix rows
      sig_pts_total <- sig_pts_total + length(sig_ids)
      herb_sig_gids <- c(herb_sig_gids, unique(gids))
    }
  }
  herb_sig_gids <- unique(herb_sig_gids)
  overlap_hits <- sum(flagged %in% herb_sig_gids)
  baseline_rate <- length(herb_sig_gids) / max(nrow(res), 1)

  # Enrichment robustness: do the TOP |rho| groups (rank-based, no p-value
  # threshold) over-represent per-herb significant features? Hypergeometric test.
  topk <- res[order(-abs(res$rho_vs_zviab))[seq_len(min(500, sum(!is.na(res$rho_vs_zviab))))], "group_id"]
  n_pop <- nrow(res)
  n_flagged_pop <- length(intersect(herb_sig_gids, res$group_id))
  n_topk_hit <- length(intersect(topk, herb_sig_gids))
  # one-sided upper-tail: P(X >= observed hits) under random draw of top-k size
  if (n_flagged_pop == 0) {
    hyper <- NA_real_   # nothing to test: no per-herb-significant groups at all
  } else {
    hyper <- tryCatch(stats::phyper(n_topk_hit - 1, n_flagged_pop, n_pop - n_flagged_pop,
                                    length(topk), lower.tail = FALSE),
                      error = function(e) NA_real_)
  }

  # --- figures --------------------------------------------------------------------------
  df <- res %>% filter(!is.na(rho_vs_zviab))
  p1 <- ggplot(df, aes(x = rho_vs_zviab, y = -log10(p_value))) +
    geom_point(aes(color = sig_dose_sensitive), alpha = 0.5, size = 1.2) +
    geom_vline(xintercept = c(-cut_rho, cut_rho), linetype = "dashed", color = "grey40") +
    scale_color_manual(values = c("TRUE" = "#d62728", "FALSE" = "grey70")) +
    labs(x = sprintf("Spearman rho: dose gradient vs viability slope"), y = expression(-log[10](p)),
         title = "Q2 - Feature-group dose gradients vs concentration-dependent viability") +
    theme_minimal() + theme(legend.position = "none")
  save_plot(p1, file.path(out_dir, "figures", "feature_dose_linkage_scatter"))

  # --- report -------------------------------------------------------------------------------
  top <- res %>% filter(sig_dose_sensitive) %>% arrange(fdr) %>% head(50)
  lines <- c(
    "# Q2 - Concentration-Sensitive Metabolic Markers", "",
    sprintf("_Generated: %s_", fmt_time()), "",
    "Per group: POOLED Spearman correlation of all herb dose-gradient values (High logFC - Low logFC,",
    "paired NA masking) with the standardized viability dose slope z_viab; BH-FDR over all integrated groups.",
    "(Pooled because each group's gradient is only measurable in ~15-80 herbs due to per-batch feature coverage.)",
    sprintf("Global pooled Spearman (all herb-group pairs): rho = %.3f, p = %s; bootstrap 95%% CI [%.3f, %.3f] (B=%d, seed=%d).",
            global_sp["rho"], format(global_sp["p"], scientific = TRUE, digits = 3),
            boot_ci["lo"], boot_ci["hi"], cfg$metabolic$bootstrap_n, cfg$metabolic$seed),
    "Note: the pooled p treats all herb-group pairs as independent (pairs share herbs) and is anti-conservative; the herb-resampled bootstrap CI is the primary inference.",
    sprintf("Cutoffs: FDR < %.2f and |rho| >= %.2f -> **%d** concentration-sensitive marker groups.", cut_fdr, cut_rho, sum(res$sig_dose_sensitive)), "",
    "## Top 50 (by FDR)", ""
  )
  if (nrow(top) > 0) {
    lines <- c(lines, md_table(top[, c("group_id", "compound_name", "kegg_id", "rho_vs_zviab", "fdr", "n_herbs_used", "direction_consistency")]))
  } else lines <- c(lines, "_none_")
  lines <- c(lines, "", "## Robustness: overlap with per-herb significant features", "",
    sprintf("- Per-herb significant points (volcano FDR < %.2f, all dose levels, %d herbs): %d (bare-variable_id rows: %d)",
            vol_fdr_cut, length(herb_ids), sig_pts_total, unmapped),
    sprintf("- Unique groups flagged by any herb: %d (chance baseline for random sets of this size: %.1f%%)",
            length(herb_sig_gids), 100 * baseline_rate),
    sprintf("- Overlap with Q2 candidate set (%d groups): **%d** (%.1f%%)",
            length(flagged), overlap_hits, 100 * overlap_hits / max(length(flagged), 1)),
    "",
    "## Enrichment robustness (rank-based)", "",
    sprintf("- Top %d groups by |rho| contain **%d** per-herb-significant groups", length(topk), n_topk_hit),
    sprintf("  (expected under random draw: %.1f; hypergeometric p = %s)",
            length(topk) * n_flagged_pop / n_pop,
            if (is.na(hyper)) "NA" else format(hyper, scientific = TRUE, digits = 3)), "",
    "Interpretation: an overlap well above the baseline rate indicates that dose-gradient concordance is not an artifact of generic per-herb significance.", "")
  write_report(lines, file.path(out_dir, "reports", "feature_dose_linkage_report.md"))
  invisible(NULL)
}

if (!interactive()) main()
