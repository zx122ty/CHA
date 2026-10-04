#!/usr/bin/env Rscript
###############################################################################
# load_inputs.R — Step 0: assemble the concentration-response input tables.
#   - CCK8 replicate values + toxicity class (from cck8_analysis module)
#   - Browse_Table metadata (categories, taste/meridian)
#   - per-herb phenotype dose slope z (fit_dose_slope on replicate values)
#   - integrated logFC matrix parsed into per-dose comparison tables
# Writes data/conc_inputs.csv (one row per herb), reports/qc_report.md.
###############################################################################
suppressPackageStartupMessages({ library(dplyr) })
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
  doses_cfg <- cfg$doses

  # --- CCK8 normalized values ------------------------------------------------------
  cck8_csv <- file.path(cfg$paths$cck8_dir, cfg$cck8$normalized_csv)
  norm <- read.csv(cck8_csv, stringsAsFactors = FALSE, check.names = FALSE)
  logmsg("CCK8 herbs: ", nrow(norm))

  tox <- read.csv(file.path(cfg$paths$cck8_dir, cfg$cck8$toxicity_class_csv), stringsAsFactors = FALSE)
  norm$toxicity_class <- tox$toxicity_class[match(norm$match_id, tox$match_id)]

  # --- metadata sanity check -----------------------------------------------------------
  # cck8_normalized.csv already carries Browse_Table fields (pinyin, drug_name_en,
  # batch_id, category_*, taste_meridian_zh) from the cck8 module's join.
  bt <- readr::read_csv(cfg$paths$browse_table, show_col_types = FALSE)
  n_missing <- sum(!norm$match_id %in% bt$match_id)
  if (n_missing > 0) logmsg("WARNING: ", n_missing, " herbs missing from Browse_Table")
  norm$batch_id <- as.integer(norm$batch_id)

  # --- per-herb phenotype dose slope ---------------------------------------------------
  dose_names <- vapply(doses_cfg, function(d) d$name, character(1))
  means_vec <- sds_vec <- ns_vec <- vector("list", length(dose_names))
  for (i in seq_along(doses_cfg)) {
    cols <- doses_cfg[[i]]$cck8_cols
    vals <- as.matrix(norm[, cols])
    means_vec[[i]] <- rowMeans(vals, na.rm = TRUE)
    sds_vec[[i]]   <- apply(vals, 1, stats::sd)          # NA when n<2
    ns_vec[[i]]    <- rowSums(!is.na(vals))
  }
  names(means_vec) <- names(sds_vec) <- names(ns_vec) <- dose_names

  herb_i <- function(lst, i) setNames(vapply(lst, function(v) v[i], numeric(1)), dose_names)
  slope_res <- t(vapply(seq_len(nrow(norm)), function(i) {
    fit_dose_slope(herb_i(means_vec, i), herb_i(sds_vec, i), herb_i(ns_vec, i), doses_cfg)
  }, numeric(3)))
  colnames(slope_res) <- c("dose_slope_viab", "z_viab", "se_z")
  norm <- cbind(norm, slope_res)

  # effect magnitudes vs control (per dose) + dose_ratio (plan §2.1)
  for (i in seq_along(doses_cfg)) {
    eff <- 1 - unlist(means_vec[[i]])
    norm[[paste0("effect_", tolower(dose_names[i]))]] <- eff
  }
  if (length(doses_cfg) == 2) {
    hi_eff <- norm[[paste0("effect_", tolower(dose_names[2]))]]   # top dose
    lo_eff <- norm[[paste0("effect_", tolower(dose_names[1]))]]   # reference dose
    denom <- hi_eff + lo_eff
    norm$dose_ratio <- ifelse(is.finite(denom) & denom > 0, hi_eff / denom, NA_real_)
  } else {
    effs <- vapply(seq_along(doses_cfg), function(i) unlist(means_vec[[i]]), numeric(nrow(norm)))
    effs <- sweep(effs, 2, 1, "-")
    norm$dose_ratio <- ifelse(rowSums(effs) > 0, effs[, ncol(effs)] / rowSums(effs), NA_real_)
  }

  # QC flags (A5: low-confidence slopes kept but flagged).
  # Same convention as the cck8 module (load_cck8.R): pairwise absolute
  # difference between replicates > rep_discordance => flag.
  rep_disc_ok <- vapply(seq_len(nrow(norm)), function(i) {
    all(vapply(doses_cfg, function(d) {
      v <- as.numeric(as.matrix(norm[i, d$cck8_cols]))
      ok <- is.finite(v)
      if (sum(ok) < 2) return(TRUE)
      v <- v[ok]
      max(abs(outer(v, v, "-")), na.rm = TRUE) <= cfg$cck8$rep_discordance
    }, logical(1)))
  }, logical(1))
  norm$low_confidence_slope <- is.na(norm$z_viab) | !rep_disc_ok

  # MS-linked flag (A4): herb's batch is in the integrated set
  norm$ms_linked <- norm$batch_id %in% cfg$integration$batches

  # --- integrated logFC matrix -----------------------------------------------------------
  lf <- load_logfc_matrix(cfg, doses_cfg)
  logmsg("logFC matrix: ", nrow(lf$matrix), " groups x ", ncol(lf$matrix), " comparisons")
  saveRDS(lf, file.path(out_dir, "data", "logfc_parsed.rds"))

  # per-herb dose-pair availability (A6: paired NA masking)
  cmp <- lf$cmp_info
  pair_avail <- unique(cmp[, c("drug_id", "dose")])
  # cmp$dose holds ms_dose suffixes (not display names) — compare like with like,
  # same as metabolic_dose_gradient.R
  d1 <- doses_cfg[[1]]$ms_dose; d2 <- doses_cfg[[2]]$ms_dose
  has1 <- pair_avail$drug_id[pair_avail$dose == d1]
  has2 <- pair_avail$drug_id[pair_avail$dose == d2]
  norm[[paste0("ms_has_", tolower(d1))]] <- norm$match_id %in% has1
  norm[[paste0("ms_has_", tolower(d2))]] <- norm$match_id %in% has2
  norm$ms_dose_paired <- norm$ms_linked & norm[[paste0("ms_has_", tolower(d1))]] & norm[[paste0("ms_has_", tolower(d2))]]

  out_csv <- file.path(out_dir, "data", "conc_inputs.csv")
  write.csv(norm, out_csv, row.names = FALSE)
  logmsg("Wrote ", out_csv, " (", nrow(norm), " herbs; ms_linked: ", sum(norm$ms_linked),
         "; dose-paired: ", sum(norm$ms_dose_paired), ")")

  # --- QC report -----------------------------------------------------------------------
  lines <- c(
    "# Concentration-Response Input QC", "",
    sprintf("_Generated: %s_", fmt_time()), "",
    sprintf("- Herbs: **%d** (CCK8); MS-linked (valid batches %s): **%d**; batch 4 only: **%d**",
            nrow(norm), paste(cfg$integration$batches, collapse = "/"), sum(norm$ms_linked), sum(!norm$ms_linked)),
    sprintf("- Dose levels (assumed, A1): %s", paste(vapply(doses_cfg, function(d) sprintf("%s=%.1f mg/mL", d$name, d$mg_ml), character(1)), collapse = ", ")),
    "", "## Phenotype dose-slope distribution", "",
    sprintf("- finite z_viab: %d / %d; low_confidence_slope flags: %d",
            sum(is.finite(norm$z_viab)), nrow(norm), sum(norm$low_confidence_slope)),
    sprintf("- median z_viab = %.2f; P10/P90 = %.2f / %.2f",
            median(norm$z_viab, na.rm = TRUE),
            quantile(norm$z_viab, c(0.10, 0.90), na.rm = TRUE)[1],
            quantile(norm$z_viab, c(0.10, 0.90), na.rm = TRUE)[2]),
    "", "## Herbs with low-confidence slopes (replicate discordance or NA)", ""
  )
  lc <- norm[norm$low_confidence_slope, c("match_id", "drug_name_zh", "batch_id", "z_viab")]
  lines <- c(lines, if (nrow(lc) > 0) md_table(head(lc, 40)) else "_none_")
  write_report(lines, file.path(out_dir, "reports", "qc_report.md"))
  invisible(NULL)
}

if (!interactive()) main()
