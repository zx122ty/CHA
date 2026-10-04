#!/usr/bin/env Rscript
###############################################################################
# load_cck8.R — Step 0a: read CCK8 xlsx, join drug metadata, derive per-herb
# summary columns and QC flags. Writes data/cck8_normalized.csv.
###############################################################################
suppressPackageStartupMessages({ library(readxl); library(dplyr); library(readr) })

this_dir <- {
  fa <- commandArgs(trailingOnly = FALSE)
  farg <- sub("^--file=", "", fa[grep("^--file=", fa)][1])
  dirname(normalizePath(farg, mustWork = FALSE))
}
source(file.path(this_dir, "cck8_utils.R"))

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a
fmt_time <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
logmsg <- function(...) cat(sprintf("[%s] ", fmt_time()), ..., "\n", sep = "")

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1 && args[[1]] == "--config") args[[2]] else "cck8_config.yaml"
cfg <- read_cck8_config(cfg_path)
out_dir <- cfg$paths$output_dir
dir.create(file.path(out_dir, "data"), showWarnings = FALSE, recursive = TRUE)

# --- read CCK8 ---------------------------------------------------------------
cck8_path <- cfg$paths$cck8_xlsx
logmsg("Reading CCK8: ", cck8_path)
cck8 <- as.data.frame(read_excel(cck8_path))
nm <- names(cck8)
nm[nm == "低浓度.1"] <- "低浓度-1"; nm[nm == "低浓度.2"] <- "低浓度-2"
nm[nm == "高浓度.1"] <- "高浓度-1"; nm[nm == "高浓度.2"] <- "高浓度-2"
names(cck8) <- nm
stopifnot(all(c("match_id", "drug_name_zh") %in% names(cck8)))

low_cols <- cfg$cck8$low_cols; high_cols <- cfg$cck8$high_cols
stopifnot(all(c(low_cols, high_cols) %in% names(cck8)))

cck8 <- cck8 %>%
  mutate(across(all_of(c(low_cols, high_cols)), as.numeric)) %>%
  mutate(
    low_mean = rowMeans(across(all_of(low_cols))),
    high_mean = rowMeans(across(all_of(high_cols))),
    rep_disc_low = abs(.data[[low_cols[1]]] - .data[[low_cols[2]]]),
    rep_disc_high = abs(.data[[high_cols[1]]] - .data[[high_cols[2]]]),
    dose_delta = high_mean - low_mean   # negative => more toxic at high dose
  ) %>%
  mutate(match_id = as.integer(match_id))

# --- join metadata -----------------------------------------------------------
bt <- readr::read_csv(cfg$paths$browse_table,
                      show_col_types = FALSE)
keep_cols <- intersect(c("match_id", "pinyin", "drug_name_en",
                         "batch_id", "category_major_zh", "category_minor_zh",
                         "taste_meridian_zh"), names(bt))
meta <- bt %>% select(all_of(keep_cols))
meta$match_id <- as.integer(meta$match_id)

norm <- cck8 %>%
  left_join(meta, by = "match_id") %>%
  mutate(
    no_quant_ms = batch_id == 4L,   # batch 4: QC gate excluded it from quantification
    flag_rep_disc = rep_disc_low > cfg$classification$rep_discordance |
                    rep_disc_high > cfg$classification$rep_discordance,
    flag_extreme_low = high_mean < cfg$classification$extreme_low | low_mean < cfg$classification$extreme_low,
    flag_extreme_high = high_mean > cfg$classification$extreme_high | low_mean > cfg$classification$extreme_high,
    qc_flag = flag_rep_disc | flag_extreme_low | flag_extreme_high
  )

n_bad <- sum(is.na(norm$batch_id))
if (n_bad > 0) logmsg("WARNING: ", n_bad, " herbs did not join metadata")

out_csv <- file.path(out_dir, "data", "cck8_normalized.csv")
write.csv(norm, out_csv, row.names = FALSE)
logmsg("Wrote ", out_csv, " (", nrow(norm), " herbs)")
