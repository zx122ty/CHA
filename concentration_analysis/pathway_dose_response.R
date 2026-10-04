#!/usr/bin/env Rscript
###############################################################################
# pathway_dose_response.R — Step 5 (Q4): dose-dependent pathway rewiring.
# group -> KEGG compound (all member features, mz/RT window) -> pathway
# (cached KEGGREST bulk table). Per pathway: pooled Spearman of the per-herb
# pathway gradient vs z_viab + toxic_high-vs-neutral Wilcoxon on the gradient.
# Exploratory: only ~13% of groups carry a KEGG compound ID (plan §2.6).
# Writes data/pathway_dose_response.csv, figures/pathway_dose_response_*.
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
  Z <- gr$Z                       # groups x dose-paired herbs (robust z of gradient)
  herb_ids <- colnames(Z)
  z_viab <- norm$z_viab[match(herb_ids, norm$match_id)]
  tox_class <- norm$toxicity_class[match(herb_ids, norm$match_id)]

  # --- KEGG caches -------------------------------------------------------------------
  link_tab <- load_kegg_link_table(cfg)
  name_f <- cfg$pathway$kegg_name_cache
  if (!file.exists(name_f)) {
    logmsg("Fetching KEGG pathway name list (one-time)...")
    dir.create(dirname(name_f), showWarnings = FALSE, recursive = TRUE)
    ok <- tryCatch({
      con <- url("https://rest.kegg.jp/list/pathway", open = "r", blocking = TRUE)
      on.exit(close(con))
      n <- readLines(con, warn = FALSE)
      cat(n, file = name_f, sep = "\n")
      length(n) > 100
    }, error = function(e) { logmsg("  download failed:", e$message); FALSE })
    if (!ok) stop("Could not fetch KEGG pathway name list; re-run with network access.")
  }
  names_tab <- read.delim(name_f, header = FALSE, col.names = c("pathway_id", "pathway_name"),
                          stringsAsFactors = FALSE, check.names = FALSE, sep = "\t")
  # cache may carry bare "mapXXXXX" or prefixed "path:mapXXXXX" ids depending on
  # how it was fetched — normalize so the match below always works
  names_tab$pathway_id <- sub("^path:", "", names_tab$pathway_id)

  # --- group -> pathway map (ALL member features, not just best-annotated) --------------
  logmsg("Rebuilding aligned members for KEGG mapping...")
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

  # one compound can map to several pathways -> long table (compound x pathway)
  kegg_cpd <- paste0("cpd:", aligned$kegg_id)
  cpd_idx <- match(kegg_cpd, link_tab$cpd)
  ok_cpd <- !is.na(aligned$kegg_id) & !is.na(cpd_idx)
  gp_long <- data.frame(group_id = aligned$group_id[ok_cpd],
                        pathway = link_tab$pathway[cpd_idx[ok_cpd]], stringsAsFactors = FALSE)
  gp_long <- unique(gp_long)
  logmsg("Member features mapped to a pathway: ", sum(!is.na(aligned$kegg_id)), " / ", nrow(aligned),
         "; group-pathway pairs: ", nrow(gp_long))

  # --- per-pathway gradient matrix (pathways x herbs) ------------------------------------
  pw_ids <- sort(unique(gp_long$pathway))
  # link table ids are "path:mapXXXXX" while the name cache is bare "mapXXXXX"
  pw_bare <- sub("^path:", "", pw_ids)
  pw_names <- names_tab$pathway_name[match(pw_bare, names_tab$pathway_id)]
  names(pw_names) <- pw_ids

  pw_groups <- split(gp_long$group_id, gp_long$pathway)
  min_g <- cfg$pathway$min_groups_per_pathway
  PW <- matrix(NA_real_, nrow = length(pw_ids), ncol = ncol(Z),
               dimnames = list(pw_ids, colnames(Z)))
  for (pi in seq_along(pw_ids)) {
    gs <- rownames(Z)[rownames(Z) %in% pw_groups[[pi]]]
    if (length(gs) < min_g) next
    PW[pi, ] <- colMeans(Z[gs, , drop = FALSE], na.rm = TRUE)
  }

  # --- tests -----------------------------------------------------------------------------
  res_rows <- list()
  for (pi in seq_along(pw_ids)) {
    x <- PW[pi, ]; z <- z_viab
    ok <- is.finite(x) & is.finite(z)
    if (sum(ok) < 8) next
    sp <- spearman_pair(x[ok], z[ok])
    # toxic_high vs neutral on the gradient
    xt <- x[tox_class == "toxic_high"]; xn <- x[tox_class == "neutral"]
    xt <- xt[is.finite(xt)]; xn <- xn[is.finite(xn)]
    wt <- if (length(xt) >= 3 && length(xn) >= 3) suppressWarnings(stats::wilcox.test(xt, xn)) else NULL
    res_rows[[pi]] <- data.frame(
      pathway_id = pw_ids[pi],
      pathway_name = pw_names[pw_ids[pi]],
      n_groups_mapped = length(pw_groups[[pi]]),
      n_herbs_used = sum(ok),
      rho_vs_zviab = sp["rho"], p_rho = sp["p"],
      med_gradient_toxic_high = if (!is.null(wt)) median(xt) else NA_real_,
      med_gradient_neutral = if (!is.null(wt)) median(xn) else NA_real_,
      p_wilcox = if (!is.null(wt)) as.numeric(wt$p.value) else NA_real_,
      n_toxic_high = length(xt), n_neutral = length(xn),
      stringsAsFactors = FALSE
    )
  }
  res <- do.call(rbind, res_rows)
  rownames(res) <- NULL
  if (length(res_rows) == 0) stop("No pathway had >=8 finite herb-gradient pairs; check KEGG mapping.")
  res$fdr_rho <- as.numeric(p.adjust(res$p_rho, "BH"))
  res$fdr_wilcox <- as.numeric(p.adjust(res$p_wilcox, "BH"))
  cut_fdr <- cfg$pathway$fdr_cutoff; cut_rho <- cfg$pathway$min_abs_rho
  res$sig_dose_dep <- !is.na(res$fdr_rho) & res$fdr_rho < cut_fdr & abs(res$rho_vs_zviab) >= cut_rho
  res$sig_toxic_contrast <- !is.na(res$fdr_wilcox) & res$fdr_wilcox < cut_fdr

  out_csv <- file.path(out_dir, "data", "pathway_dose_response.csv")
  write.csv(res, out_csv, row.names = FALSE)
  logmsg("Wrote ", out_csv, "; sig_dose_dep: ", sum(res$sig_dose_dep), "; sig_toxic_contrast: ", sum(res$sig_toxic_contrast))

  # --- figure ---------------------------------------------------------------------------
  top <- res %>% filter(!is.na(rho_vs_zviab)) %>% arrange(-abs(rho_vs_zviab)) %>% head(15)
  p1 <- ggplot(top, aes(x = reorder(pathway_name, rho_vs_zviab), y = rho_vs_zviab, fill = rho_vs_zviab < 0)) +
    geom_col(color = "white") + coord_flip() +
    scale_fill_manual(name = NULL, values = c("TRUE" = "#d62728", "FALSE" = "#3182bd")) +
    labs(x = NULL, y = "Spearman rho (pathway gradient vs viability slope)",
         title = "Q4 - Top pathways by dose-response concordance") + theme_minimal()
  save_plot(p1, file.path(out_dir, "figures", "pathway_dose_response_top"))

  # --- report ------------------------------------------------------------------------------
  lines <- c(
    "# Q4 - Dose-Dependent Pathway Rewiring (KEGG)", "",
    sprintf("_Generated: %s_  (EXPLORATORY: KEGG compound coverage ~13%% of groups)", fmt_time()), "",
    "Per pathway: mean per-herb gradient z over its mapped groups; pooled Spearman vs z_viab;",
    "plus toxic_high-vs-neutral Wilcoxon on the same gradient. BH-FDR across all pathways.", "",
    sprintf("- Pathways tested: **%d** (min %d mapped groups); sig_dose_dep (FDR<%.2f, |rho|>=%.2f): **%d**; sig_toxic_contrast: **%d**",
            nrow(res), min_g, cut_fdr, cut_rho, sum(res$sig_dose_dep), sum(res$sig_toxic_contrast)), "",
    "## Top 30 by |rho|", ""
  )
  top30 <- res %>% filter(!is.na(rho_vs_zviab)) %>% arrange(-abs(rho_vs_zviab)) %>% head(30)
  lines <- c(lines, md_table(top30[, c("pathway_id", "pathway_name", "n_groups_mapped", "rho_vs_zviab", "fdr_rho",
                                       "med_gradient_toxic_high", "med_gradient_neutral", "fdr_wilcox")]))
  if (sum(res$sig_dose_dep) > 0) {
    sig <- res %>% filter(sig_dose_dep) %>% arrange(fdr_rho)
    lines <- c(lines, "", "## Significant dose-dependent pathways (FDR)", "")
    lines <- c(lines, md_table(sig[, c("pathway_id", "pathway_name", "n_groups_mapped", "rho_vs_zviab", "fdr_rho")]))
  }
  write_report(lines, file.path(out_dir, "reports", "pathway_dose_response_report.md"))
  invisible(NULL)
}

if (!interactive()) main()
