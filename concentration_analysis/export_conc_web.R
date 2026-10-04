#!/usr/bin/env Rscript
###############################################################################
# export_conc_web.R — Step 6 (Q5): Dose-Response Index + web contract.
# DRI combines the phenotype dose slope (|z| capped) and metabolic gradient
# magnitude (P95-scaled); quadrant labels for the website facet.
# Writes data/dri.csv, web_export/concentration.json (flat, stable keys —
# same contract style as cck8.json; for a future `ingest_concentration`).
###############################################################################
suppressPackageStartupMessages({ library(jsonlite) })
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
  d <- cfg$dri

  df <- merge(pheno[, c("match_id", "drug_name_zh", "drug_name_en", "pinyin", "batch_id",
                         "category_major_zh", "category_minor_zh", "low_mean", "high_mean",
                         "z_viab", "pheno_class", "active_pheno", "toxicity_class",
                         "ms_linked", "ms_dose_paired", "low_confidence_slope")],
              metab[, c("match_id", "metab_gradient_mag", "metab_class")],
              by = "match_id", all.x = TRUE)

  # DRI components (plan §2.7)
  df$dri_pheno <- sign(df$z_viab) * pmin(abs(df$z_viab), d$z_cap) / d$z_cap   # NA-safe: NA propagates
  meto_scale <- if (d$meto_p95) quantile(df$metab_gradient_mag, 0.95, na.rm = TRUE) else max(df$metab_gradient_mag, na.rm = TRUE)
  if (!is.finite(meto_scale) || meto_scale <= 0) meto_scale <- NA_real_   # degenerate scale -> dri_meto NA
  df$dri_meto <- df$metab_gradient_mag / meto_scale
  df$dri <- d$weight_pheno * abs(df$dri_pheno) + d$weight_meto * df$dri_meto

  # quadrant labels — same rules as discordant_herbs.R (plan §2.5). Herbs without
  # metabolic data (batch 4) count as "metabolism flat" (nothing measured), so a
  # strongly inhibitory one lands on active_phenotype_only; their DRI stays NULL
  # because the composite index requires both sides.
  q_hi <- quantile(metab$metab_gradient_mag, cfg$metabolic$quantile_high, na.rm = TRUE)
  q_lo <- quantile(metab$metab_gradient_mag, cfg$metabolic$quantile_low, na.rm = TRUE)
  z_cut <- cfg$phenotype$z_cut
  strong_pheno <- is.finite(df$z_viab) & abs(df$z_viab) >= z_cut
  inhib <- is.finite(df$z_viab) & df$z_viab <= -z_cut
  hi_meto <- !is.na(df$metab_gradient_mag) & df$metab_gradient_mag >= q_hi
  flat_meto <- is.na(df$metab_gradient_mag) | df$metab_gradient_mag <= q_lo
  df$dri_quadrant <- ifelse(strong_pheno & hi_meto, "active_both",
                    ifelse(inhib & flat_meto, "active_phenotype_only",
                    ifelse(!strong_pheno & hi_meto, "active_metabolism_only", "inactive")))

  out_csv <- file.path(out_dir, "data", "dri.csv")
  write.csv(df, out_csv, row.names = FALSE)
  logmsg("Wrote ", out_csv)

  # --- web contract ------------------------------------------------------------------
  # row-object lists (same flat style as cck8.json's top_toxicity_features;
  # a future ingest_concentration migration can parse these directly)
  df_rows <- function(x) {
    lapply(seq_len(nrow(x)), function(i) {
      r <- as.list(x[i, , drop = FALSE])
      for (cn in names(r)) {
        if (is.numeric(r[[cn]])) {
          v <- as.numeric(r[[cn]])
          if (!is.finite(v)) {
            r[[cn]] <- NULL
          } else {
            # whole numbers (counts, ids) stay integers; fractional values keep
            # 3 significant digits so tiny fdr values don't collapse to 0
            r[[cn]] <- if (v == round(v)) as.integer(v) else signif(v, 3)
          }
        } else if (is.logical(r[[cn]])) {
          r[[cn]] <- isTRUE(r[[cn]])
        } else if (is.na(r[[cn]])) {
          r[[cn]] <- NULL
        }
      }
      r
    })
  }

  top_feat <- tryCatch({
    f <- read.csv(file.path(out_dir, "data", "feature_dose_linkage.csv"), stringsAsFactors = FALSE)
    df_rows(f[order(-abs(f$rho_vs_zviab))[seq_len(min(50, nrow(f)))],
      c("group_id", "compound_name", "kegg_id", "rho_vs_zviab", "fdr", "n_herbs_used", "direction_consistency")])
  }, error = function(e) NULL)

  top_pw <- tryCatch({
    pw <- read.csv(file.path(out_dir, "data", "pathway_dose_response.csv"), stringsAsFactors = FALSE)
    df_rows(pw[order(-abs(pw$rho_vs_zviab))[seq_len(min(20, nrow(pw)))],
        c("pathway_id", "pathway_name", "n_groups_mapped", "rho_vs_zviab", "fdr_rho")])
  }, error = function(e) NULL)

  herb_json <- lapply(seq_len(nrow(df)), function(i) {
    r <- df[i, ]
    list(
      match_id = as.integer(r$match_id),
      drug_name_zh = r$drug_name_zh,
      drug_name_en = if (is.na(r$drug_name_en)) NULL else r$drug_name_en,
      pinyin = if (is.na(r$pinyin)) NULL else r$pinyin,
      batch_id = as.integer(r$batch_id),
      category_major_zh = if (is.na(r$category_major_zh)) NULL else r$category_major_zh,
      low_mean = round(as.numeric(r$low_mean), 4),
      high_mean = round(as.numeric(r$high_mean), 4),
      z_viab = if (is.finite(r$z_viab)) round(as.numeric(r$z_viab), 3) else NULL,
      pheno_class = r$pheno_class,
      toxicity_class = if (is.na(r$toxicity_class)) NULL else r$toxicity_class,
      metab_gradient_mag = if (!is.na(r$metab_gradient_mag)) round(as.numeric(r$metab_gradient_mag), 5) else NULL,
      metab_class = if (is.na(r$metab_class)) NULL else r$metab_class,
      dri_pheno = if (is.finite(r$dri_pheno)) round(as.numeric(r$dri_pheno), 4) else NULL,
      dri_meto = if (is.finite(r$dri_meto)) round(as.numeric(r$dri_meto), 4) else NULL,
      dri = if (is.finite(r$dri)) round(as.numeric(r$dri), 4) else NULL,
      dri_quadrant = r$dri_quadrant,
      ms_linked = isTRUE(r$ms_linked),
      ms_dose_paired = isTRUE(r$ms_dose_paired),
      low_confidence_slope = isTRUE(r$low_confidence_slope)
    )
  })

  payload <- list(
    meta = list(
      module = "concentration_response",
      generated = fmt_time(),
      n_herbs = nrow(df),
      n_dose_paired = sum(!is.na(df$metab_gradient_mag)),
      doses_mg_ml = vapply(cfg$doses, function(x) x$mg_ml, numeric(1)),
      dose_names = vapply(cfg$doses, function(x) x$name, character(1)),
      assumptions = list(dose_confirmed = FALSE, single_plate = TRUE, proliferative_active = TRUE),
      config_ref = "code/this_project/concentration_analysis/concentration_config.yaml"
    ),
    herbs = herb_json,
    top_dose_sensitive_features = if (is.null(top_feat)) NULL else top_feat,
    top_pathways = if (is.null(top_pw)) NULL else top_pw
  )

  out_json <- file.path(out_dir, "web_export", "concentration.json")
  # null="null": jsonlite's default renders NULL list elements as {} — the web
  # contract requires JSON null (ingest code treats {} as a present value).
  write_json(payload, out_json, pretty = TRUE, auto_unbox = TRUE, na = "null", null = "null")
  logmsg("Wrote ", out_json)
  invisible(NULL)
}

if (!interactive()) main()
