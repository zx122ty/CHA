#!/usr/bin/env Rscript
###############################################################################
# run_cck8_analysis.R — orchestrator for the CCK8 x metabolomics analysis.
# Usage: Rscript run_cck8_analysis.R --config cck8_config.yaml [--steps 0,1,2,3,4,6]
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- "cck8_config.yaml"
if ("--config" %in% args) cfg_path <- args[[which(args == "--config") + 1L]]
steps <- c("0", "1", "2", "3", "4", "6")
if ("--steps" %in% args) steps <- strsplit(args[[which(args == "--steps") + 1L]], ",")[[1]]

this_dir <- {
  fa <- commandArgs(trailingOnly = FALSE)
  farg <- sub("^--file=", "", fa[grep("^--file=", fa)][1])
  dirname(normalizePath(farg, mustWork = FALSE))
}
source(file.path(this_dir, "cck8_utils.R"))

cfg <- read_cck8_config(cfg_path)
logmsg("=== CCK8 x metabolomics analysis ===")
logmsg("Config: ", cfg_path)
logmsg("Output: ", cfg$paths$output_dir)

run_step <- function(name, script) {
  logmsg(sprintf("--- Step %s: %s ---", name, basename(script)))
  t0 <- Sys.time()
  ok <- system2("Rscript", c(shQuote(file.path(this_dir, script)), "--config", shQuote(cfg_path)), wait = TRUE) == 0
  if (!ok) stop(sprintf("Step %s (%s) failed", name, basename(script)))
  logmsg(sprintf("Step %s done in %.1f s", name, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

if ("0" %in% steps) {
  run_step("0a load+normalize", "load_cck8.R")
  run_step("0b qc report", "qc_cck8.R")
}
if ("1" %in% steps) run_step("1 toxicity classification", "classify_toxicity.R")
if ("2" %in% steps) run_step("2 herb-level linkage", "link_metabolomics.R")
if ("3" %in% steps) run_step("3 feature-level association", "toxicity_features.R")
if ("4" %in% steps) run_step("4 pathway comparison", "pathway_toxicity.R")
if ("6" %in% steps) run_step("6 web export", "export_cck8_web.R")

logmsg("=== All requested steps complete ===")
