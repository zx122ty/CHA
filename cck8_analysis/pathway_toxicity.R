#!/usr/bin/env Rscript
###############################################################################
# pathway_toxicity.R — Step 4: compare KEGG/disease pathway enrichment between
# toxic and neutral herbs (stored pathway_enrichment_*.json; sparse: only some
# comparisons have hits). Writes data/pathway_toxicity.csv,
# reports/pathway_toxicity_report.md, figures/pathway_toxicity.{png,pdf}.
###############################################################################
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

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
herbs <- norm %>% filter(batch_id %in% cfg$integration$batches)

# --- collect pathway hits per comparison ---------------------------------------
rows <- list()
for (bid in cfg$integration$batches) {
  web <- file.path(cfg$paths$results_dir, paste0("Batch", bid), "detail", "web_export")
  pf <- list.files(web, pattern = "^pathway_enrichment_[0-9]+_(High|Low)_vs_CT[0-9]+\\.json$", full.names = TRUE)
  for (f in pf) {
    cmp <- sub("^pathway_enrichment_", "", sub("\\.json$", "", basename(f)))
    drug_id <- as.integer(sub("^([0-9]+)_(High|Low)_vs_CT[0-9]+$", "\\1", cmp))
    if (!(drug_id %in% herbs$match_id)) next
    j <- readj(f)
    pw <- j$pathways
    if (is.null(pw) || length(pw) == 0) next
    df <- data.frame(
      match_id = drug_id,
      drug_name_zh = herbs$drug_name_zh[herbs$match_id == drug_id],
      batch_id = bid,
      comparison = cmp,
      concentration = sub(".*_(High|Low)_vs_CT[0-9]+$", "\\1", cmp),
      pathway_name = vapply(pw, function(x) if (is.null(x$pathway_name)) "" else as.character(x$pathway_name), character(1)),
      p_value = vapply(pw, function(x) if (is.null(x$p_value)) NA_real_ else as.numeric(as.character(x$p_value)), numeric(1)),
      mapped_count = as.integer(vapply(pw, function(x) if (is.null(x$mapped_count)) NA_integer_ else as.integer(x$mapped_count), integer(1))),
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1]] <- df
  }
}
pw_all <- do.call(rbind, rows)
if (is.null(pw_all) || nrow(pw_all) == 0) {
  stop("No pathway hits found in any valid-batch comparison")
}
pw_all$toxicity_class <- herbs$toxicity_class[match(pw_all$match_id, herbs$match_id)]
write.csv(pw_all, file.path(out_dir, "data", "pathway_hits_raw.csv"), row.names = FALSE)

# --- per-herb pathway sets --------------------------------------------------------
herb_pw <- pw_all %>%
  group_by(match_id, drug_name_zh, toxicity_class) %>%
  summarise(n_pathways = n_distinct(pathway_name), min_p = min(p_value, na.rm = TRUE), .groups = "drop")

# --- Fisher: pathway over-represented in toxic vs neutral herbs ---------------------
toxic_herbs <- unique(herb_pw$match_id[herb_pw$toxicity_class %in% c("toxic_high", "toxic_low")])
neutral_herbs <- unique(herb_pw$match_id[herb_pw$toxicity_class == "neutral"])

pathway_universe <- sort(unique(pw_all$pathway_name))
res_rows <- lapply(pathway_universe, function(pw) {
  # count herbs (not comparisons)
  t_n <- length(unique(pw_all$match_id[pw_all$pathway_name == pw & pw_all$match_id %in% toxic_herbs]))
  n_n <- length(unique(pw_all$match_id[pw_all$pathway_name == pw & pw_all$match_id %in% neutral_herbs]))
  tab <- matrix(c(t_n, length(toxic_herbs) - t_n, n_n, length(neutral_herbs) - n_n), nrow = 2,
                dimnames = list(c("toxic", "neutral"), c("has_pw", "no_pw")))
  ft <- suppressWarnings(stats::fisher.test(tab, simulate.p.value = TRUE, B = 10000))
  data.frame(pathway_name = pw, n_toxic_herbs = t_n, n_neutral_herbs = n_n,
             fisher_p = as.numeric(ft$p.value), stringsAsFactors = FALSE)
})
pw_res <- do.call(rbind, res_rows)
pw_res$fdr <- as.numeric(p.adjust(pw_res$fisher_p, method = "BH"))
pw_res$sig_enriched_toxic <- pw_res$fdr < cfg$linkage$fdr_cutoff & pw_res$n_toxic_herbs > 0
write.csv(pw_res, file.path(out_dir, "data", "pathway_toxicity.csv"), row.names = FALSE)

# --- report ---------------------------------------------------------------------------
sig_pw <- pw_res %>% filter(sig_enriched_toxic) %>% arrange(fisher_p)
lines <- c(
  "# Pathway Enrichment: Toxic vs Neutral Herbs", "",
  paste0("_Generated: ", fmt_time(), "_"), "",
  sprintf("Pathway hits available for **%d of %d** comparisons (enrichment is sparse in the stored exports).",
          length(unique(pw_all$comparison)), sum(vapply(cfg$integration$batches, function(bid) {
            web <- file.path(cfg$paths$results_dir, paste0("Batch", bid), "detail", "web_export")
            length(list.files(web, pattern = "^pathway_enrichment_[0-9]+_(High|Low)_vs_CT[0-9]+\\.json$"))
          }, integer(1)))),
  sprintf("Fisher exact test per pathway: herb with >=1 hit on that pathway — toxic (n=%d) vs neutral (n=%d); BH-FDR < %.2f.",
          length(toxic_herbs), length(neutral_herbs), cfg$linkage$fdr_cutoff), "",
  "## Pathways enriched in toxic herbs", ""
)
if (nrow(sig_pw) > 0) {
  lines <- c(lines, md_table(sig_pw[, c("pathway_name","n_toxic_herbs","n_neutral_herbs","fisher_p","fdr")]))
} else {
  lines <- c(lines, "_No pathway reached FDR < 0.05._", "",
    "Given the sparsity of stored pathway hits (disease-pathway mapping), treat this step as exploratory; the feature-level results (Step 3) are the primary readout.")
}
lines <- c(lines, "", "## All pathways observed", "")
lines <- c(lines, md_table(pw_res %>% arrange(fisher_p) %>% head(40) %>% select(pathway_name, n_toxic_herbs, n_neutral_herbs, fisher_p, fdr)))

dir.create(file.path(out_dir, "reports"), showWarnings = FALSE, recursive = TRUE)
writeLines(lines, file.path(out_dir, "reports", "pathway_toxicity_report.md"))

# --- figure -----------------------------------------------------------------------------
d <- pw_res %>% filter(n_toxic_herbs + n_neutral_herbs >= 3) %>% head(25) %>%
  mutate(pathway_name = factor(pathway_name, levels = rev(pathway_name)))
p <- ggplot(d, aes(x = fisher_p, y = pathway_name)) +
  geom_point(aes(color = n_toxic_herbs > n_neutral_herbs), size = 2) +
  scale_x_log10() +
  scale_color_manual(values = c("TRUE" = "#d62728", "FALSE" = "#1f77b4"),
                     labels = c("TRUE" = "more in toxic", "FALSE" = "not enriched")) +
  geom_vline(xintercept = cfg$linkage$fdr_cutoff, linetype = "dashed", color = "grey50") +
  labs(x = "Fisher p (log10)", y = NULL,
       title = "Pathway enrichment: toxic vs neutral herbs") +
  theme_minimal(base_size = 12) + theme(axis.text.x = element_text(size = 9))
save_plot(p, file.path(out_dir, "figures", "pathway_toxicity"))
logmsg("Done: pathway toxicity report; ", nrow(sig_pw), " pathways FDR-significant")
