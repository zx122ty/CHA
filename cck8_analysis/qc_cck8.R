#!/usr/bin/env Rscript
###############################################################################
# qc_cck8.R — Step 0b: QC report for the CCK8 dataset (no plate layout is
# available, so QC relies on replicate discordance + extreme values).
# Writes reports/qc_report.md.
###############################################################################
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

norm <- read.csv(file.path(out_dir, "data", "cck8_normalized.csv"), stringsAsFactors = FALSE)
thr <- cfg$classification$rep_discordance

flagged <- norm[norm$qc_flag, ]
by_type <- data.frame(
  flag = c(sprintf("replicate discordance > %.2f (either dose)", thr),
           sprintf("extreme low < %.2f", cfg$classification$extreme_low),
           sprintf("extreme high > %.2f", cfg$classification$extreme_high)),
  n = c(sum(norm$flag_rep_disc), sum(norm$flag_extreme_low), sum(norm$flag_extreme_high))
)

lines <- c(
  "# CCK8 QC Report",
  "",
  paste0("_Generated: ", fmt_time(), " from `", cfg$paths$cck8_xlsx, "`_"),
  "",
  "## Overview",
  "",
  sprintf("- Herbs: **%d** (2 doses x 2 replicates each; values assumed control-normalized)", nrow(norm)),
  sprintf("- Low-dose mean viability: median %.3f (range %.3f-%.3f)",
          median(norm$low_mean), min(norm$low_mean), max(norm$low_mean)),
  sprintf("- High-dose mean viability: median %.3f (range %.3f-%.3f)",
          median(norm$high_mean), min(norm$high_mean), max(norm$high_mean)),
  "",
  "## Flags (marked, not removed)",
  ""
)
lines <- c(lines, md_table(by_type))
lines <- c(lines, "")

if (nrow(flagged) > 0) {
  lines <- c(lines, "### Flagged herbs", "")
  fl <- flagged[, c("match_id", "drug_name_zh", "batch_id", "low_mean", "high_mean",
                    "rep_disc_low", "rep_disc_high")]
  lines <- c(lines, md_table(fl))
} else {
  lines <- c(lines, "No herbs flagged.")
}

lines <- c(lines, "", "## Caveats", "",
  "- No control-well raw values (OD450) are available in the v3 file, so plate-level QC (edge effect, plate drift) cannot be checked. Assumption: values are already normalized to control (=1).",
  "- Only 2 replicates per dose; replicate discordance is the only within-herb consistency check.",
  paste0("- ", sum(norm$no_quant_ms), " herbs belong to batch 4, which is excluded from quantitative metabolomics (QC gate). Their CCK8 values are retained but cannot be linked to MS data."),
  "- Assumed dosing: low = 0.2 mg/mL, high = 1.0 mg/mL (same as the metabolomics experiment, protocol_v3) — unconfirmed.")

dir.create(file.path(out_dir, "reports"), showWarnings = FALSE, recursive = TRUE)
out_md <- file.path(out_dir, "reports", "qc_report.md")
writeLines(lines, out_md)
logmsg("Wrote ", out_md)
