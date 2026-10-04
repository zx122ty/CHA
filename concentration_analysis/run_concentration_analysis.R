#!/usr/bin/env Rscript
###############################################################################
# run_concentration_analysis.R — entry point for the concentration-response
# module. Usage:
#   Rscript run_concentration_analysis.R --config concentration_config.yaml [--steps 0,1,2,3,4,5,6]
# Steps (plan §3):
#   0 load_inputs              CCK8 + metadata + logFC parsing + QC
#   1 phenotype_dose_response  Q1 response-shape classification (161 herbs)
#   2 metabolic_dose_gradient  per-herb x group gradients + magnitude (133 herbs)
#   3 feature_viability_linkage Q2 concentration-sensitive markers
#   4 discordant_herbs         Q3 phenotype-metabolism quadrants
#   5 pathway_dose_response    Q4 KEGG dose-dependent pathways (exploratory)
#   6 export_conc_web          Q5 DRI + web contract
###############################################################################
suppressPackageStartupMessages({ library(dplyr) })
this_dir <- {
  fa <- commandArgs(trailingOnly = FALSE)
  farg <- sub("^--file=", "", fa[grep("^--file=", fa)][1])
  dirname(normalizePath(farg, mustWork = FALSE))
}
source(file.path(this_dir, "conc_utils.R"))

main <- function() {
  ca <- conc_cli_args()
  cfg_path <- ca$config %||% file.path(this_dir, "concentration_config.yaml")
  steps <- ca$steps

  cfg <- read_conc_config(cfg_path)
  out_dir <- cfg$paths$output_dir
  t_start <- Sys.time()
  logmsg("Concentration-response analysis; steps: ", paste(steps, collapse = ","))

  run_step <- function(n, name) {
    if (!(n %in% steps)) return(invisible(NULL))
    logmsg(sprintf("=== Step %d: %s ===", n, name))
    t0 <- Sys.time()
    source(file.path(this_dir, name))   # each script runs main() when non-interactive
    logmsg(sprintf("Step %d done in %.1f s", n, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }

  run_step(0, "load_inputs.R")
  run_step(1, "phenotype_dose_response.R")
  run_step(2, "metabolic_dose_gradient.R")
  run_step(3, "feature_viability_linkage.R")
  run_step(4, "discordant_herbs.R")
  run_step(5, "pathway_dose_response.R")
  run_step(6, "export_conc_web.R")

  # --- summary report (always written when all steps ran) --------------------------
  if (setequal(steps, c(0,1,2,3,4,5,6))) {
    dri <- read.csv(file.path(out_dir, "data", "dri.csv"), stringsAsFactors = FALSE)
    feat <- read.csv(file.path(out_dir, "data", "feature_dose_linkage.csv"), stringsAsFactors = FALSE)
    pw <- read.csv(file.path(out_dir, "data", "pathway_dose_response.csv"), stringsAsFactors = FALSE)
    disc <- read.csv(file.path(out_dir, "data", "discordant_herbs.csv"), stringsAsFactors = FALSE)

    fin_rho <- is.finite(feat$rho_vs_zviab)
    top_rho_line <- if (any(fin_rho)) {
      k <- which(fin_rho)[which.max(abs(feat$rho_vs_zviab[fin_rho]))]
      sprintf("- Top |rho|: %.2f (group %s); global pooled rho with bootstrap CI is reported in `feature_dose_linkage_report.md`.",
              abs(feat$rho_vs_zviab[k]), feat$group_id[k])
    } else "- Top |rho|: no group had a finite rho."

    lines <- c(
      "# Concentration-Response Analysis — Summary", "",
      sprintf("_Generated: %s_  (runtime %.1f s; config: `concentration_config.yaml`)", fmt_time(),
              as.numeric(difftime(Sys.time(), t_start, units = "secs"))), "",
      "## Design in one line", "",
      sprintf("%d-level dose design (%s; values assumed, A1): per-herb quantities are dose contrasts with replicate-variance SE; formal inference is cross-herb (bootstrap/BH-FDR). See plan §2.",
              length(cfg$doses),
              paste(vapply(cfg$doses, function(d) sprintf("%s=%.1f mg/mL", d$name, as.numeric(d$mg_ml)), character(1)),
                    collapse = ", ")), "",
      sprintf("## Q1 — Phenotype response shape (%d herbs)", nrow(dri)), ""
    )
    lines <- c(lines, md_table(count(dri, pheno_class) %>% as.data.frame()))
    lines <- c(lines, "",
      sprintf("- %d herbs show a significant inhibitory dose slope (z <= -%.1f); of those, %d are threshold-type (effect concentrated at High), only %d linear.",
              sum(dri$z_viab <= -cfg$phenotype$z_cut & is.finite(dri$z_viab)), cfg$phenotype$z_cut,
              sum(dri$pheno_class == "threshold_high"), sum(dri$pheno_class == "linear_sensitive")), "",
      "## Q2 — Concentration-sensitive markers", "",
      sprintf("- **%d** of %d integrated groups pass FDR<%.2f & |rho|>=%.2f (pooled per-group Spearman vs viability slope).",
              sum(feat$sig_dose_sensitive), nrow(feat), cfg$metabolic$fdr_cutoff, cfg$metabolic$min_abs_rho),
      top_rho_line, "",
      sprintf("## Q3 — Discordant herbs (%d dose-paired)", nrow(disc)), ""
    )
    lines <- c(lines, md_table(count(disc, quadrant) %>% as.data.frame()))
    lines <- c(lines, "",
      sprintf("- **active_metabolism_only** (metabolically active, phenotype-flat) = pharmacologically active candidates: %d herbs.",
              sum(disc$quadrant == "active_metabolism_only")), "",
      "## Q4 — Pathway level (exploratory)", "",
      sprintf("- %d KEGG pathways tested; sig_dose_dep: %d; sig_toxic_contrast: %d (KEGG coverage ~13%% of groups).",
              nrow(pw), sum(pw$sig_dose_dep), sum(pw$sig_toxic_contrast)), "",
      "## Q5 — DRI (website)", "",
      sprintf("- DRI = %.1f·|z_viab capped at %.0f| + %.1f·(metab magnitude / P95); exported to `web_export/concentration.json` (%d herbs).",
              cfg$dri$weight_pheno, cfg$dri$z_cap, cfg$dri$weight_meto, nrow(dri)), "",
      "## Caveats (must accompany any citation)", "",
      "- Dose values 0.2/1 mg/mL are **assumed** (A1) — confirm with the lab; standardized z is unaffected, absolute slopes are not.",
      "- Single CCK8 plate assumed (A2); n=2 replicates per dose → per-herb classes are descriptive (A5).",
      sprintf("- Batch 4 herbs (%d) have phenotype-only entries (ms_linked=false); their DRI is NULL (composite index needs both sides) but quadrant labels use the phenotype side only.", sum(!dri$ms_linked)),
      "- Q2/Q4 negative FDR results are informative: at two concentrations, global concentration-dependent signatures are weak; the rank-based and quadrant views carry the signal.", "")
    write_report(lines, file.path(out_dir, "reports", "concentration_summary.md"))
  }
  logmsg("All done.")
}

main()
