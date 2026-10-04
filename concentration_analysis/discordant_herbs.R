#!/usr/bin/env Rscript
###############################################################################
# discordant_herbs.R — Step 4 (Q3): phenotype-vs-metabolism discordance.
# Two-axis plot: x = z_viab (phenotype dose slope), y = metab_gradient_mag.
# Quadrants (plan §2.5):
#   active_both            |z|>=z_cut & mag >= P75
#   active_phenotype_only  z<=-z_cut & mag <= P25   (toxic, metabolism flat; QC suspect)
#   active_metabolism_only |z|<z_cut & mag >= P75   (metabolically active candidates)
#   inactive               the rest
# TCM major-category Fisher enrichment per quadrant.
# Writes data/discordant_herbs.csv, figures/discordant_quadrants_*.
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

  pheno <- read.csv(file.path(out_dir, "data", "phenotype_response.csv"), stringsAsFactors = FALSE)
  metab <- read.csv(file.path(out_dir, "data", "metabolic_gradient.csv"), stringsAsFactors = FALSE)
  p <- cfg$phenotype

  # effect columns are derived from the config dose names (created in step 0)
  dose_names <- vapply(cfg$doses, function(d) d$name, character(1))
  eff_hi_col <- paste0("effect_", tolower(dose_names[length(dose_names)]))
  eff_lo_col <- paste0("effect_", tolower(dose_names[1]))
  df <- merge(pheno[, c("match_id", "category_major_zh", "z_viab",
                         "pheno_class", "low_mean", "high_mean", eff_lo_col, eff_hi_col)],
              metab, by = "match_id")
  logmsg("Dose-paired herbs for two-axis analysis: ", nrow(df))

  q_hi <- quantile(metab$metab_gradient_mag, cfg$metabolic$quantile_high, na.rm = TRUE)
  q_lo <- quantile(metab$metab_gradient_mag, cfg$metabolic$quantile_low, na.rm = TRUE)
  strong_pheno <- is.finite(df$z_viab) & (abs(df$z_viab) >= p$z_cut)
  inhib <- is.finite(df$z_viab) & df$z_viab <= -p$z_cut
  flat_meto <- df$metab_gradient_mag <= q_lo
  df$quadrant <- ifelse(strong_pheno & df$metab_gradient_mag >= q_hi, "active_both",
                ifelse(inhib & flat_meto, "active_phenotype_only",
                ifelse(!strong_pheno & df$metab_gradient_mag >= q_hi, "active_metabolism_only", "inactive")))

  out_csv <- file.path(out_dir, "data", "discordant_herbs.csv")
  write.csv(df, out_csv, row.names = FALSE)
  logmsg("Wrote ", out_csv)

  # --- figure -----------------------------------------------------------------------
  df$quadrant <- factor(df$quadrant, levels = c("active_both", "active_phenotype_only",
                                                 "active_metabolism_only", "inactive"))
  p1 <- ggplot(df, aes(x = z_viab, y = metab_gradient_mag, color = quadrant)) +
    geom_point(alpha = 0.85, size = 2) +
    geom_vline(xintercept = c(-p$z_cut, p$z_cut), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = q_hi, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = q_lo, linetype = "dotted", color = "grey55") +
    scale_color_manual(values = c("active_both" = "#d62728", "active_phenotype_only" = "#ff7f0e",
                                  "active_metabolism_only" = "#2ca02c", "inactive" = "grey65")) +
    labs(x = "phenotype dose slope z (negative = inhibition)", y = "metabolic gradient magnitude (mean |z|)",
         title = sprintf("Q3 - Phenotype vs metabolic concentration response (%d herbs)", nrow(df))) +
    theme_minimal() + theme(legend.position = "bottom")
  save_plot(p1, file.path(out_dir, "figures", "discordant_quadrants"))

  # --- TCM enrichment per quadrant (exact 2x2 Fisher per category vs rest) --------------
  cat_tab <- table(df$category_major_zh, df$quadrant)
  enr <- fisher_category_enrichment(
    cat_tab, classes = c("active_both", "active_phenotype_only", "active_metabolism_only"),
    class_label = "quadrant", category_label = "category_major_zh", count_label = "n_in_quad")
  if (!is.null(enr)) write.csv(enr, file.path(out_dir, "data", "discordant_category_enrichment.csv"), row.names = FALSE)

  # --- report ---------------------------------------------------------------------------
  lines <- c(
    "# Q3 - Phenotype-Metabolism Discordant Herbs", "",
    sprintf("_Generated: %s_  (%d dose-paired herbs; thresholds in config)", fmt_time(), nrow(df)), "",
    "Quadrants: x = standardized viability dose slope z; y = metabolic gradient magnitude.",
    sprintf("Rules (plan §2.5): active_both |z|>=%.1f & mag>=P75; active_phenotype_only z<=-%.1f & mag<=P25;", p$z_cut, p$z_cut),
    sprintf("(active_metabolism_only |z|<%.1f & mag>=P75; inactive = rest).", p$z_cut),
    sprintf("Cutoffs: |z| >= %.2f (strong phenotype); mag >= P75 = %.4f (metabolically active).", p$z_cut, q_hi), "",
    "## Quadrant counts", ""
  )
  lines <- c(lines, md_table(as.data.frame(table(df$quadrant))))
  for (cl in c("active_metabolism_only", "active_phenotype_only")) {
    sub_df <- df[df$quadrant == cl, c("match_id", "drug_name_zh", "category_major_zh", "z_viab", "metab_gradient_mag", "pheno_class")]
    lines <- c(lines, "", sprintf("## %s (n=%d) — full list", cl, nrow(sub_df)), "")
    ord <- order(-abs(sub_df$z_viab), -sub_df$metab_gradient_mag)
    lines <- c(lines, if (nrow(sub_df) > 0) md_table(sub_df[ord, ]) else "_none_")
  }
  if (!is.null(enr)) {
    top_enr <- enr %>% filter(fdr < 0.10) %>% arrange(fdr) %>% head(20)
    lines <- c(lines, "", "## TCM major-category enrichment (Fisher, FDR<0.10)", "")
    lines <- c(lines, if (nrow(top_enr) > 0) md_table(top_enr) else "_none at FDR < 0.10_")
  }
  write_report(lines, file.path(out_dir, "reports", "discordant_herbs_report.md"))
  invisible(NULL)
}

if (!interactive()) main()
