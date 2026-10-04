#!/usr/bin/env Rscript
###############################################################################
# phenotype_dose_response.R — Step 1 (Q1): classify herbs by concentration-
# response shape on the CCK8 phenotype side (all 161 herbs, A4).
# Classes (plan §2.3): linear_sensitive / threshold_high / threshold_low /
# proliferative (A3: active) / insensitive. Descriptive (A5); QC flags carried.
# Writes data/phenotype_response.csv, figures/response_shape_*.
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
  p <- cfg$phenotype

  # --- classification (plan §2.3 rules, refined on observed distribution) ---------
  # Order matters: proliferative first (A3), then inhibitory subclasses by shape.
  # Effect columns are derived from the config dose names (created in step 0), so
  # renaming/adding dose levels needs no code change here.
  dose_names <- vapply(cfg$doses, function(d) d$name, character(1))
  eff_hi_col <- paste0("effect_", tolower(dose_names[length(dose_names)]))   # top dose
  eff_lo_col <- paste0("effect_", tolower(dose_names[1]))                    # reference dose
  z <- norm$z_viab; hi_eff <- norm[[eff_hi_col]]; lo_eff <- norm[[eff_lo_col]]
  cls <- rep("insensitive", nrow(norm))
  is_prolif <- (1 - norm[[eff_hi_col]]) > p$proliferative_threshold   # == high_mean > threshold
  inhib <- (z <= -p$z_cut) & !is_prolif
  eff_ratio <- ifelse(lo_eff > 0, hi_eff / lo_eff, NA_real_)
  is_linear   <- inhib & (lo_eff >= p$effect_cut) &
                 (eff_ratio >= p$linear_ratio_min) & (eff_ratio <= p$linear_ratio_max)
  is_thr_high <- inhib & !is_linear & (lo_eff < p$effect_cut) & (hi_eff >= p$effect_cut)
  is_saturat  <- inhib & !is_linear & !is_thr_high & (lo_eff >= p$effect_cut) &
                 (hi_eff < lo_eff * p$saturating_factor)
  is_weak     <- inhib & !is_linear & !is_thr_high & !is_saturat   # significant slope, small effects
  cls[is_linear]   <- "linear_sensitive"
  cls[is_thr_high] <- "threshold_high"
  cls[is_saturat]  <- "saturating_low"
  cls[is_weak]     <- "weak_inhibitory"
  cls[is_prolif]   <- "proliferative"

  # activity flag (A3): toxic directions OR proliferative
  norm$pheno_class <- cls
  norm$active_pheno <- is_prolif | inhib
  norm$inhibitory <- inhib

  out_csv <- file.path(out_dir, "data", "phenotype_response.csv")
  write.csv(norm, out_csv, row.names = FALSE)
  logmsg("Wrote ", out_csv)

  # --- figures ----------------------------------------------------------------------
  df <- norm %>%
    mutate(pheno_class = factor(pheno_class, levels = c("linear_sensitive", "threshold_high",
                                                        "saturating_low", "weak_inhibitory",
                                                        "proliferative", "insensitive")))
  p1 <- ggplot(df, aes(x = z_viab, y = dose_ratio, color = pheno_class)) +
    geom_point(alpha = 0.85, size = 2) +
    geom_vline(xintercept = -p$z_cut, linetype = "dashed", color = "grey40") +
    labs(x = "standardized viability dose slope  z (negative = inhibition)",
         y = "dose_ratio (share of effect at High dose)",
         title = sprintf("Q1 - Phenotype concentration-response shape (%d herbs)", nrow(df))) +
    theme_minimal() + theme(legend.position = "bottom")
  save_plot(p1, file.path(out_dir, "figures", "response_shape_scatter"))

  tab <- df %>% count(pheno_class) %>% as.data.frame()
  p2 <- ggplot(tab, aes(x = reorder(pheno_class, n), y = n, fill = pheno_class)) +
    geom_col(color = "white") + coord_flip() +
    labs(x = NULL, y = "n herbs", title = "Q1 - response-shape class counts") + theme_minimal()
  save_plot(p2, file.path(out_dir, "figures", "response_shape_counts"))

  # --- TCM category enrichment per class (exact 2x2 Fisher per category vs rest) ------
  cat_tab <- table(df$category_major_zh, df$pheno_class)
  enr <- fisher_category_enrichment(
    cat_tab, classes = setdiff(as.character(unique(df$pheno_class)), "insensitive"),
    class_label = "pheno_class", category_label = "category_major_zh", count_label = "n_in_class")
  if (!is.null(enr)) write.csv(enr, file.path(out_dir, "data", "phenotype_category_enrichment.csv"), row.names = FALSE)

  # --- report -------------------------------------------------------------------------
  lines <- c(
    "# Q1 — Phenotype Concentration-Response Shape", "",
    sprintf("_Generated: %s_  (descriptive per A5; thresholds in `concentration_config.yaml`)", fmt_time()), "",
    "Rules: z = standardized viability dose slope (SE from replicate variance); effect = 1 − mean viability;",
    sprintf("z_cut=%.2f, effect_cut=%.2f, linear window [%.2f, %.2f], proliferative > %.2f.",
            p$z_cut, p$effect_cut, p$linear_ratio_min, p$linear_ratio_max, p$proliferative_threshold), "",
    "## Class counts", ""
  )
  lines <- c(lines, md_table(tab))
  lines <- c(lines, "", "## Notable herbs", "")
  for (cl in c("linear_sensitive", "threshold_high", "saturating_low", "weak_inhibitory", "proliferative")) {
    sub_df <- norm[norm$pheno_class == cl, c("match_id", "drug_name_zh", "low_mean", "high_mean", "z_viab", eff_lo_col, eff_hi_col, "ms_linked")]
    lines <- c(lines, sprintf("### %s (n=%d)", cl, nrow(sub_df)), "")
    lines <- c(lines, if (nrow(sub_df) > 0) md_table(head(sub_df[order(-abs(sub_df$z_viab)), ], 25)) else "_none_")
    lines <- c(lines, "")
  }
  if (!is.null(enr)) {
    top_enr <- enr %>% filter(fdr < 0.10) %>% arrange(fdr) %>% head(20)
    lines <- c(lines, "## TCM major-category enrichment (Fisher vs rest, FDR<0.10)", "")
    lines <- c(lines, if (nrow(top_enr) > 0) md_table(top_enr) else "_none at FDR < 0.10_")
  }
  write_report(lines, file.path(out_dir, "reports", "phenotype_response_report.md"))
  invisible(NULL)
}

if (!interactive()) main()
