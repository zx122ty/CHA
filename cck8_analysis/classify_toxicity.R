#!/usr/bin/env Rscript
###############################################################################
# classify_toxicity.R — Step 1: threshold-based toxicity classification with
# sensitivity analysis. Writes data/toxicity_class.csv, reports/toxicity_summary.md,
# figures/toxicity_classification.{png,pdf}.
###############################################################################
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

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
toxic_thr <- cfg$classification$toxic_threshold
proliferative_thr <- cfg$classification$proliferative_threshold

classify_at <- function(df, thr) {
  cls <- ifelse(df$high_mean < thr & df$low_mean >= thr, "toxic_high",
        ifelse(df$low_mean < thr, "toxic_low",
        ifelse(df$high_mean > proliferative_thr, "proliferative", "neutral")))
  cls
}

norm$toxicity_class <- classify_at(norm, toxic_thr)
n_toxic <- sum(norm$toxicity_class %in% c("toxic_high", "toxic_low"))

# sensitivity table: class membership per threshold
sens_rows <- lapply(cfg$classification$sensitivity_thresholds, function(thr) {
  cls <- classify_at(norm, thr)
  data.frame(
    threshold = thr,
    n_toxic_high = sum(cls == "toxic_high"),
    n_toxic_low = sum(cls == "toxic_low"),
    n_toxic_total = sum(cls %in% c("toxic_high", "toxic_low")),
    n_proliferative = sum(cls == "proliferative"),
    n_neutral = sum(cls == "neutral")
  )
})
sens <- do.call(rbind, sens_rows)

# which toxic herbs are in batch 4 (no quantitative MS data)
toxic_df <- norm[norm$toxicity_class %in% c("toxic_high", "toxic_low"), ]
n_toxic_b4 <- sum(toxic_df$no_quant_ms)

out_csv <- file.path(out_dir, "data", "toxicity_class.csv")
write.csv(norm, out_csv, row.names = FALSE)

# --- report ------------------------------------------------------------------
lines <- c(
  "# Toxicity Classification Summary", "",
  paste0("_Generated: ", fmt_time(), "_"), "",
  sprintf("Rules (default threshold **< %.2f** viability; proliferative **> %.2f**):", toxic_thr, proliferative_thr), "",
  "- `toxic_high`: high-dose mean < thr AND low-dose mean >= thr",
  "- `toxic_low`: low-dose mean < thr (toxic at the lower dose)",
  "- `proliferative`: high-dose mean > prolif. thr (kept separate; excluded from toxic/neutral contrasts)",
  "- `neutral`: everything else", "",
  "## Class counts (default threshold)", ""
)
cls_levels <- c("toxic_high", "toxic_low", "proliferative", "neutral")
cls_tab <- data.frame(
  class = cls_levels,
  n = vapply(cls_levels, function(cl) sum(norm$toxicity_class == cl), integer(1)),
  of_which_batch4_no_ms = vapply(cls_levels,
                                 function(cl) sum(norm$no_quant_ms & norm$toxicity_class == cl), integer(1))
)
lines <- c(lines, md_table(cls_tab), "",
  "## Sensitivity to threshold", "")
lines <- c(lines, md_table(sens), "",
  sprintf("**%d** of the %d toxic herbs sit in batch 4 (no quantitative MS data) — CCK8-only for them.", n_toxic_b4, n_toxic), "",
  "## Toxic herbs (default threshold)", "")
lines <- c(lines, md_table(toxic_df[, c("match_id","drug_name_zh","batch_id","toxicity_class","low_mean","high_mean","dose_delta","qc_flag")]))

dir.create(file.path(out_dir, "reports"), showWarnings = FALSE, recursive = TRUE)
writeLines(lines, file.path(out_dir, "reports", "toxicity_summary.md"))
logmsg("Wrote toxicity summary; toxic n=", n_toxic)

# --- figure ------------------------------------------------------------------
plot_df <- norm %>%
  mutate(class = factor(toxicity_class, levels = c("toxic_low","toxic_high","neutral","proliferative")))
p <- ggplot(plot_df, aes(x = low_mean, y = high_mean, color = class)) +
  geom_point(alpha = 0.85, size = 1.6) +
  geom_abline(intercept = toxic_thr, slope = 0, linetype = "dashed", color = "grey40") +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey60") +
  scale_color_brewer(palette = "Set2") +
  labs(x = sprintf("Low-dose mean viability (assumed %.1f mg/mL)", cfg$cck8$assumed_low_conc_mg_ml),
       y = sprintf("High-dose mean viability (assumed %.1f mg/mL)", cfg$cck8$assumed_high_conc_mg_ml),
       title = "CCK8 toxicity classification (n=161 herbs, rat primary hepatocytes)") +
  theme_minimal(base_size = 13) +
  guides(color = guide_legend(title = "Class"))
save_plot(p, file.path(out_dir, "figures", "toxicity_classification"))
logmsg("Wrote figure toxicity_classification")
