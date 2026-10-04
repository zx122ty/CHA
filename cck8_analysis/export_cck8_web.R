#!/usr/bin/env Rscript
###############################################################################
# export_cck8_web.R — Step 6: write web_export/cck8.json following the existing
# annotation-contract style (stable schema for the Django ingest + Vue pages).
###############################################################################

suppressPackageStartupMessages({ library(dplyr) })

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

norm <- read.csv(file.path(out_dir, "data", "toxicity_class.csv"), stringsAsFactors = FALSE)

# replicate columns as stored in the CSV (read.csv applies make.names: "-" -> ".")
low_cols_csv  <- gsub("-", ".", cfg$cck8$low_cols, fixed = TRUE)
high_cols_csv <- gsub("-", ".", cfg$cck8$high_cols, fixed = TRUE)
stopifnot(all(c(low_cols_csv, high_cols_csv) %in% names(norm)))

# optional: per-herb linkage metrics from Step 2 (for the website scatter)
linkage_metrics <- NULL
lk_path <- file.path(out_dir, "data", "herb_linkage.csv")
if (file.exists(lk_path)) {
  lk <- read.csv(lk_path, stringsAsFactors = FALSE)
  keep <- intersect(c("match_id", "concentration", "viability", "n_significant",
                      "median_abs_logFC_sig", "perturbation_distance"), names(lk))
  linkage_metrics <- lk[, keep, drop = FALSE]
}

# optional: attach top toxicity-associated features per herb (from Step 3)
feat_path <- file.path(out_dir, "data", "toxicity_features.csv")
top_feats <- NULL
if (file.exists(feat_path)) {
  tf <- read.csv(feat_path, stringsAsFactors = FALSE)
  top_feats <- tf %>% filter(sig_high | sig_two_group) %>% arrange(fdr_high) %>% head(50)
}

lk_metric <- function(mid, conc, field) {
  if (is.null(linkage_metrics)) return(NULL)
  sub_df <- linkage_metrics[linkage_metrics$match_id == mid &
                             linkage_metrics$concentration == conc, , drop = FALSE]
  if (nrow(sub_df) == 0 || is.na(sub_df[[field]][1])) return(NULL)
  as.numeric(sub_df[[field]][1])
}

herb_rows <- lapply(seq_len(nrow(norm)), function(i) {
  r <- norm[i, ]
  list(
    match_id = as.integer(r$match_id),
    drug_name_zh = r$drug_name_zh,
    batch_id = as.integer(r$batch_id),
    low_reps = round(c(as.numeric(r[[low_cols_csv[1]]]), as.numeric(r[[low_cols_csv[2]]])), 4),
    high_reps = round(c(as.numeric(r[[high_cols_csv[1]]]), as.numeric(r[[high_cols_csv[2]]])), 4),
    low_mean = round(as.numeric(r$low_mean), 4),
    high_mean = round(as.numeric(r$high_mean), 4),
    dose_delta = round(as.numeric(r$dose_delta), 4),
    toxicity_class = r$toxicity_class,
    qc_flag = isTRUE(r$qc_flag),
    no_quant_ms = isTRUE(r$no_quant_ms),
    n_sig_low = lk_metric(r$match_id, "Low", "n_significant"),
    n_sig_high = lk_metric(r$match_id, "High", "n_significant"),
    perturbation_distance_low = lk_metric(r$match_id, "Low", "perturbation_distance"),
    perturbation_distance_high = lk_metric(r$match_id, "High", "perturbation_distance")
  )
})

payload <- list(
  assay = "CCK8",
  generated_at = fmt_time(),
  source_file = cfg$paths$cck8_xlsx,
  assumptions = list(
    values_normalized_to_control = TRUE,
    assumed_low_conc_mg_ml = cfg$cck8$assumed_low_conc_mg_ml,
    assumed_high_conc_mg_ml = cfg$cck8$assumed_high_conc_mg_ml,
    cell_system = "rat primary hepatocytes (same as metabolomics experiment)",
    toxic_threshold = cfg$classification$toxic_threshold,
    proliferative_threshold = cfg$classification$proliferative_threshold
  ),
  n_herbs = nrow(norm),
  class_counts = list(
    toxic_high = sum(norm$toxicity_class == "toxic_high"),
    toxic_low = sum(norm$toxicity_class == "toxic_low"),
    proliferative = sum(norm$toxicity_class == "proliferative"),
    neutral = sum(norm$toxicity_class == "neutral")
  ),
  herbs = herb_rows,
  top_toxicity_features = if (is.null(top_feats)) NULL else lapply(seq_len(nrow(top_feats)), function(i) {
    r <- top_feats[i, ]
    list(group_id = r$group_id, compound_name = r$compound_name,
         rho_high = round(r$rho_high, 3), fdr_high = r$fdr_high)
  })
)

out_json <- file.path(out_dir, "web_export", "cck8.json")
# null="null": jsonlite's default renders NULL list elements as {} — the web
# contract requires JSON null (ingest code treats {} as a present value).
write_json(payload, out_json, pretty = TRUE, auto_unbox = TRUE, na = "null", null = "null")
logmsg("Wrote ", out_json, " (", nrow(norm), " herbs)")
