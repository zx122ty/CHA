#!/usr/bin/env Rscript
###############################################################################
# metabolic_dose_gradient.R — Step 2: per-herb metabolic dose gradients.
#   grad(g,d) = logFC_{d,High} - logFC_{d,Low}  (paired NA masking, plan §2.2/A6)
#   global robust z-standardization per group across herbs
#   per-herb magnitude metab_gradient_mag(d) = mean_g |z_g|
#   classification: metab_sensitive (>=P75) / metab_insensitive (<=P25) / intermediate
# Writes data/metabolic_gradient.csv, data/gradient_z.rds (for steps 3 & 5).
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
  lf <- readRDS(file.path(out_dir, "data", "logfc_parsed.rds"))
  M <- lf$matrix; cmp <- lf$cmp_info
  doses_cfg <- cfg$doses

  # dose-paired herbs (A4/A6)
  paired <- norm[norm$ms_dose_paired, ]
  logmsg("Dose-paired herbs: ", nrow(paired))

  d1 <- doses_cfg[[1]]$ms_dose   # Low
  d2 <- doses_cfg[[2]]$ms_dose   # High

  # column index per (drug, dose)
  cmp_idx <- split(seq_len(nrow(cmp)), paste(cmp$drug_id, cmp$dose))
  herb_ids <- paired$match_id
  G <- matrix(NA_real_, nrow = nrow(M), ncol = length(herb_ids),
              dimnames = list(rownames(M), as.character(herb_ids)))

  for (i in seq_along(herb_ids)) {
    h <- herb_ids[i]
    ch <- cmp_idx[[paste(h, d2)]]   # High column
    cl <- cmp_idx[[paste(h, d1)]]   # Low column
    if (is.null(ch) || is.null(cl)) next
    G[, i] <- pairwise_dose_gradient(M[, ch], M[, cl])
  }

  # coverage per group (how many herbs have a finite gradient)
  cov_g <- rowSums(is.finite(G))
  logmsg("Group coverage: median ", as.integer(median(cov_g)), " / ", length(herb_ids), " herbs")

  # min herbs with a finite gradient for z-standardization (also used in the
  # report below — single threshold so code and text cannot drift apart)
  thr <- max(5L, as.integer(round(0.1 * ncol(G))))
  # global robust z-standardization per group (plan §2.2)
  Z <- t(apply(G, 1, function(x) {
    if (sum(is.finite(x)) < thr) return(rep(NA_real_, length(x)))
    robust_z(x)
  }))

  colnames(Z) <- herb_ids   # t() drops column names; restore for downstream steps

  # per-herb magnitude: mean |z| across groups (robust to outlier features)
  mag <- colMeans(abs(Z), na.rm = TRUE)
  names(mag) <- herb_ids
  n_groups_used <- colSums(is.finite(Z))

  res <- data.frame(
    match_id = as.integer(names(mag)),
    drug_name_zh = paired$drug_name_zh,
    batch_id = paired$batch_id,
    metab_gradient_mag = round(as.numeric(mag), 5),
    n_groups_with_gradient = n_groups_used,
    stringsAsFactors = FALSE
  )
  # groups excluded from z-standardization: below coverage threshold or constant
  # gradient (MAD=0 — mostly exact-zero High/Low logFC pairs; cannot be scaled)
  const_rows <- vapply(seq_len(nrow(G)), function(i) {
    v <- G[i, ]; v <- v[is.finite(v)]
    length(v) >= thr && stats::mad(v) == 0
  }, logical(1))
  n_low_cov <- sum(cov_g < thr)

  q_hi <- quantile(res$metab_gradient_mag, cfg$metabolic$quantile_high, na.rm = TRUE)
  q_lo <- quantile(res$metab_gradient_mag, cfg$metabolic$quantile_low, na.rm = TRUE)
  res$metab_class <- ifelse(res$metab_gradient_mag >= q_hi, "metab_sensitive",
                     ifelse(res$metab_gradient_mag <= q_lo, "metab_insensitive", "intermediate"))

  out_csv <- file.path(out_dir, "data", "metabolic_gradient.csv")
  write.csv(res, out_csv, row.names = FALSE)
  saveRDS(list(Z = Z, G = G, cov_g = cov_g), file.path(out_dir, "data", "gradient_z.rds"))
  logmsg("Wrote ", out_csv, " and gradient_z.rds")

  # --- figure: magnitude distribution with class bands --------------------------------
  df <- res %>% mutate(metab_class = factor(metab_class, levels = c("metab_insensitive", "intermediate", "metab_sensitive")))
  p <- ggplot(df, aes(x = metab_gradient_mag, fill = metab_class)) +
    geom_histogram(bins = 40, color = "white") +
    geom_vline(xintercept = c(q_lo, q_hi), linetype = "dashed", color = "grey30") +
    scale_fill_manual(values = c("metab_insensitive" = "#9ecae1", "intermediate" = "#deebf7", "metab_sensitive" = "#3182bd")) +
    labs(x = "metabolic dose-gradient magnitude (mean |z| across groups)", y = "n herbs",
         title = sprintf("Q2/Q3 - Metabolic dose-response magnitude (%d dose-paired herbs)", nrow(res))) +
    theme_minimal() + theme(legend.position = "none")
  save_plot(p, file.path(out_dir, "figures", "metabolic_gradient_dist"))

  # --- report ----------------------------------------------------------------------------
  lines <- c(
    "# Metabolic Dose Gradients (per herb x feature group)", "",
    sprintf("_Generated: %s_  (A6: paired NA masking; ComBat-corrected integrated logFC)", fmt_time()), "",
    sprintf("- Dose-paired herbs: **%d** / %d (batch 4 excluded)", nrow(res), nrow(norm)),
    sprintf("- Groups with gradient in >=50%% of herbs: %d / %d", sum(cov_g >= 0.5 * length(herb_ids)), nrow(M)),
    sprintf("- Excluded from z-standardization: %d constant-gradient groups (MAD=0), %d below coverage threshold (<%d herbs)", sum(const_rows), n_low_cov, thr),
    sprintf("- magnitude quantiles: P25=%.3f, median=%.3f, P75=%.3f", q_lo, median(res$metab_gradient_mag), q_hi), "",
    "## Class counts (metab)", ""
  )
  lines <- c(lines, md_table(as.data.frame(table(res$metab_class))))
  top_sens <- res %>% filter(metab_class == "metab_sensitive") %>% arrange(-metab_gradient_mag) %>% head(20)
  top_ins  <- res %>% filter(metab_class == "metab_insensitive") %>% arrange(metab_gradient_mag) %>% head(20)
  lines <- c(lines, "", "## Most metabolically dose-sensitive herbs (top 20)", "")
  lines <- c(lines, md_table(top_sens[, c("match_id", "drug_name_zh", "batch_id", "metab_gradient_mag")]))
  lines <- c(lines, "", "## Least metabolically dose-responsive herbs (bottom 20)", "")
  lines <- c(lines, md_table(top_ins[, c("match_id", "drug_name_zh", "batch_id", "metab_gradient_mag")]))
  write_report(lines, file.path(out_dir, "reports", "metabolic_gradient_report.md"))
  invisible(NULL)
}

if (!interactive()) main()
