#!/usr/bin/env Rscript
# Rebuild web_export/metabolites.json for a batch from the already-computed
# identification report, so a change to build_metabolites_json() (or the report
# that writes compound_source) can be applied without re-running the pipeline.
#
# Reads the Excel "All_Identifications" sheet written by
# generate_identification_report() — the same final candidate table (with
# SIRIUS columns and feature-level MSI levels) as all_annotated_metabolites.csv.
# features.json needs the mass_dataset and is only written by the pipeline.
#
# Usage: Rscript refresh_web_export.R <batch_output_dir>
#   e.g. Rscript refresh_web_export.R /path/results/Batch5/detail

suppressPackageStartupMessages({
  library(dplyr)
  library(jsonlite)
  library(openxlsx)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript refresh_web_export.R <batch_detail_dir>")
detail_dir <- normalizePath(args[1], mustWork = TRUE)

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
if (is.null(script_dir) || script_dir == ".") script_dir <- getwd()
source(file.path(script_dir, "export_for_web.R"))

xlsx <- file.path(detail_dir, "06_Annotation", "Compound_identification_report.xlsx")
if (!file.exists(xlsx)) stop(sprintf("Identification report not found: %s", xlsx))
app3 <- read.xlsx(xlsx, sheet = "All_Identifications", check.names = FALSE)
cat(sprintf("Loaded %d annotation rows x %d columns\n", nrow(app3), ncol(app3)))

# object2 is only used for expression data, which build_metabolites_json() does
# not touch — a NULL is fine here.
metabolites <- build_metabolites_json(app3, NULL)

out <- file.path(detail_dir, "web_export", "metabolites.json")
write_json_export(metabolites, out)
cat(sprintf("Wrote %s (%d metabolites)\n", out, length(metabolites$metabolites)))
