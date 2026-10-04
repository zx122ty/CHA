#!/usr/bin/env Rscript
###############################################################################
# make_config.R — Metadata-driven YAML Config Generator
#
# Reads a metadata CSV, auto-detects its column structure, identifies groups,
# and generates a fully-populated YAML configuration file for the
# metabolomics_pipeline.R enhanced pipeline (with SIRIUS, MetFrag, CAMERA support).
#
# Usage:
#   Rscript make_config.R \
#     --metadata  /path/to/metadata.csv \
#     --raw-dir   /path/to/mzML \
#     --output    /path/to/config.yaml \
#     [--project-name MyProject] \
#     [--polarity  positive|negative] \
#     [--sirius-path /path/to/sirius] \
#     [--sirius-project /path/to/project_space] \
#     [--ms1-db-dir /path/to/ms1/databases] \
#     [--ms2-db-dir /path/to/ms2/databases]
###############################################################################

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
})

# ══════════════════════════════════════════════════════════════════════════════
# 0. CLI definition
# ══════════════════════════════════════════════════════════════════════════════

option_list <- list(
  make_option(c("--metadata"), type = "character", default = NULL,
              help = "Path to metadata CSV file", metavar = "FILE"),
  make_option(c("--raw-dir"), type = "character", default = NULL,
              help = "Path to raw mzML data directory", metavar = "DIR"),
  make_option(c("--output"), type = "character", default = "config.yaml",
              help = "Output YAML config file path [default: %default]",
              metavar = "FILE"),
  make_option(c("--project-name"), type = "character", default = NULL,
              help = "Project name (defaults to basename of --metadata parent dir)"),
  make_option(c("--polarity"), type = "character", default = NULL,
              help = "Override polarity: 'positive' or 'negative' (auto-detected if omitted)"),
  make_option(c("--working-dir"), type = "character", default = NULL,
              help = "Working directory (defaults to current directory)"),
  make_option(c("--finalOutput"), type = "character", default = NULL,
              help = "Final output directory (output_dir = finalOutput/detail, peak_table_dir = finalOutput/peak_table). If omitted, outputs are placed alongside the metadata file.", metavar = "DIR"),
  make_option(c("--sirius-path"), type = "character", default = NULL,
              help = "Path to SIRIUS CLI executable (e.g. /home/zhengxiao/soft/sirius/bin/sirius). If omitted, SIRIUS features will be disabled.",
              metavar = "FILE"),
  make_option(c("--sirius-project"), type = "character", default = NULL,
              help = "Path to SIRIUS project space directory. Defaults to <output_dir>/sirius_project",
              metavar = "DIR"),
  make_option(c("--ms1-db-dir"), type = "character", default = NULL,
              help = "Path to MS1 accurate-mass database directory. Overrides the default path.",
              metavar = "DIR"),
  make_option(c("--ms2-db-dir"), type = "character", default = NULL,
              help = "Path to MS2 spectral database directory. Overrides the default path.",
              metavar = "DIR"),
  make_option(c("--batch-id"), type = "integer", default = NULL,
              help = "Batch ID (1-6) for web export. If omitted, batch_id will be NA.",
              metavar = "INT"),
  make_option(c("--split-concentration"), dest = "split_concentration",
              action = "store_true", default = FALSE,
              help = "Split each drug into separate High/Low treatment groups (e.g. 10_High_vs_CT1) using the concentration column"),
  make_option(c("--oplsda-permutations"), dest = "oplsda_permutations",
              type = "integer", default = NULL,
              help = "OPLS-DA permutation count (default 1000). Use 100 for quick test runs.",
              metavar = "INT"),
  make_option(c("--cores"), type = "integer", default = NULL,
              help = "Parallel workers for OPLS-DA permutations etc. (0 or omit = auto-detect all cores). Also scales SIRIUS --threads.",
              metavar = "INT"),
  make_option(c("--peakwidth"), type = "character", default = NULL,
              help = "Chromatographic peak width range in seconds as 'min,max' (default '1,8'; observed peaks here are 0.7-3 s; old default was 5,30). Only affects NEW feature extractions — an existing Peak_table_for_cleaning.csv is reused as-is.",
              metavar = "RANGE")
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  usage = "Rscript make_config.R --metadata <FILE> --raw-dir <DIR> [options]"
))

# Validate required args
if (is.null(opt$metadata)) stop("--metadata is required. Use --help for usage.")
if (is.null(opt$`raw-dir`))  stop("--raw-dir is required. Use --help for usage.")

# ── Default paths ──────────────────────────────────────────────────────────
# SIRIUS path: use user-provided or try common default
sirius_path <- opt$`sirius-path`
if (is.null(sirius_path) || sirius_path == "") {
  # Check common install locations
  for (candidate in c("/home/data/shareData/soft/sirius/bin/sirius",
                       "/usr/local/bin/sirius",
                       "/home/zhengxiao/soft/sirius/bin/sirius")) {
    if (file.exists(candidate)) {
      sirius_path <- candidate
      break
    }
  }
  if (is.null(sirius_path) || sirius_path == "") {
    sirius_path <- "/path/to/sirius/bin/sirius"  # placeholder
  }
}

# Internal-standard reference CSV: resolve relative to this script's own
# location (<project root>/code/this_project/make_config.R → ../metadata/),
# because the pipeline setwd()s to working_dir before reading the path.
is_csv_path <- ""
{
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  candidates <- c(
    if (length(script_arg) > 0) {
      file.path(
        normalizePath(file.path(dirname(sub("^--file=", "", script_arg[1])),
                                "..", "metadata"), mustWork = FALSE),
        "internal_standards.csv")
    },
    file.path("code", "metadata", "internal_standards.csv")
  )
  for (candidate in candidates) {
    if (file.exists(candidate)) {
      is_csv_path <- normalizePath(candidate)
      break
    }
  }
  if (is_csv_path == "") is_csv_path <- candidates[1]
}

# ══════════════════════════════════════════════════════════════════════════════
# 1. Read and inspect metadata
# ══════════════════════════════════════════════════════════════════════════════

cat("===== make_config: Metadata-driven config generator =====\n")

metadata_file <- opt$metadata
if (!file.exists(metadata_file)) {
  stop(sprintf("Metadata file not found: %s", metadata_file))
}
meta <- read.csv(metadata_file, stringsAsFactors = FALSE)
cat(sprintf("=> Loaded metadata: %d rows × %d columns\n", nrow(meta), ncol(meta)))
cat(sprintf("   Columns: %s\n", paste(colnames(meta), collapse = ", ")))

# ══════════════════════════════════════════════════════════════════════════════
# 2. Auto-detect column roles
# ══════════════════════════════════════════════════════════════════════════════

cols <- tolower(colnames(meta))
col_map <- list(
  filename  = NULL,
  ms_level  = NULL,
  group_id  = NULL,
  replicate = NULL,
  mode      = NULL
)

# Filename column: contains ".mzML" or named "filename"/"file_name"/"file"
for (nm in c("file_name", "filename", "file", "sample", "sample_name")) {
  idx <- which(cols == tolower(nm))
  if (length(idx) == 1) {
    col_map$filename <- colnames(meta)[idx]
    cat(sprintf("   Detected filename column: '%s'\n", col_map$filename))
    break
  }
}
if (is.null(col_map$filename)) {
  # Heuristic: first column containing ".mzML" in values
  for (cn in colnames(meta)) {
    if (any(grepl("\\.mzML$", meta[[cn]], ignore.case = TRUE))) {
      col_map$filename <- cn
      cat(sprintf("   Detected filename column (heuristic): '%s'\n", cn))
      break
    }
  }
}

# MS level column
for (nm in c("ms_level", "mslevel", "ms", "level")) {
  idx <- which(cols == tolower(nm))
  if (length(idx) == 1) {
    col_map$ms_level <- colnames(meta)[idx]
    cat(sprintf("   Detected MS level column: '%s'\n", col_map$ms_level))
    vals <- unique(meta[[col_map$ms_level]])
    cat(sprintf("     Values: %s\n", paste(vals, collapse = ", ")))
    break
  }
}
if (is.null(col_map$ms_level)) {
  cat("   ⚠️  No MS level column detected; assuming all files are MS1\n")
}

# Group/drug ID column
for (nm in c("drug_id", "group", "group_id", "treatment", "class", "condition", "mzmine_sample_type")) {
  idx <- which(cols == tolower(nm))
  if (length(idx) == 1) {
    col_map$group_id <- colnames(meta)[idx]
    cat(sprintf("   Detected group_id column: '%s'\n", col_map$group_id))
    vals <- unique(meta[[col_map$group_id]])
    cat(sprintf("     Unique groups (%d): %s\n", length(unique(vals)),
                paste(head(unique(vals), 20), collapse = ", "),
                if (length(unique(vals)) > 20) "..." else ""))
    break
  }
}

# Replicate column
for (nm in c("replicate_id", "replicate", "rep", "repid")) {
  idx <- which(cols == tolower(nm))
  if (length(idx) == 1) {
    col_map$replicate <- colnames(meta)[idx]
    cat(sprintf("   Detected replicate column: '%s'\n", col_map$replicate))
    break
  }
}

# Mode/polarity column
for (nm in c("mode", "polarity", "ion_mode", "ionmode")) {
  idx <- which(cols == tolower(nm))
  if (length(idx) == 1) {
    col_map$mode <- colnames(meta)[idx]
    cat(sprintf("   Detected mode column: '%s'\n", col_map$mode))
    vals <- unique(meta[[col_map$mode]])
    cat(sprintf("     Values: %s\n", paste(vals, collapse = ", ")))
    break
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. Determine polarity
# ══════════════════════════════════════════════════════════════════════════════

if (!is.null(opt$polarity)) {
  polarity <- tolower(opt$polarity)
  cat(sprintf("   Polarity (from --polarity): %s\n", polarity))
} else if (!is.null(col_map$mode)) {
  mode_vals <- tolower(unique(meta[[col_map$mode]]))
  if (any(grepl("positive", mode_vals))) {
    polarity <- "positive"
  } else if (any(grepl("negative", mode_vals))) {
    polarity <- "negative"
  } else {
    polarity <- "positive"  # default
  }
  cat(sprintf("   Polarity (detected from mode): %s\n", polarity))
} else {
  polarity <- "positive"
  cat(sprintf("   Polarity (default): %s\n", polarity))
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. Identify special sample roles
# ══════════════════════════════════════════════════════════════════════════════

group_col <- col_map$group_id
all_groups <- unique(meta[[group_col]])
all_groups <- all_groups[!is.na(all_groups) & all_groups != ""]

# Special role detection heuristics
# Control groups in this study are named CT1..CT6 (one per batch), so match
# bare "CT"/"control"/... AND the CT+digits form. Without the digits form the
# control is misclassified as a treatment group and reference_group falls back
# to "CT", which matches no sample — every comparison is then skipped.
detect_role <- function(groups) {
  roles <- list(control = NULL, qc_ms1 = NULL, blank = NULL, treatment = c())
  for (g in groups) {
    gl <- tolower(as.character(g))
    if (gl %in% c("ct", "c", "control", "ctrl", "nc", "normal") || grepl("^ct[0-9]*$", gl)) {
      roles$control <- g
    } else if (gl %in% c("qc", "q")) {
      roles$qc_ms1 <- g
    } else if (gl %in% c("blank", "bk", "blk", "bl")) {
      roles$blank <- g
    } else if (grepl("^blank", gl) && !grepl("msms", gl)) {
      roles$blank <- g
    } else if (gl %in% c("msms", "qc-msms", "qc_msms")) {
      # MSMS pattern — handled automatically, not a separate role
    } else {
      roles$treatment <- c(roles$treatment, g)
    }
  }
  roles
}

roles <- detect_role(all_groups)

# ── Optional: split treatment groups by concentration (High/Low) ──────────
# Sets metadata_columns.group_suffix so build_sample_info_all() appends the
# concentration to every non-special sample's class (CT1/QC/Blank unchanged).
group_suffix_col <- NULL
if (isTRUE(opt$split_concentration)) {
  for (nm in c("drug_concentration", "concentration", "dose", "dose_level")) {
    idx <- which(cols == tolower(nm))
    if (length(idx) == 1) {
      group_suffix_col <- colnames(meta)[idx]
      break
    }
  }
  if (is.null(group_suffix_col)) {
    stop("--split-concentration: no concentration column found (looked for drug_concentration/concentration/dose/dose_level)")
  }
  special <- c(roles$control, roles$qc_ms1, roles$blank)
  suffix_vals <- as.character(meta[[group_suffix_col]])
  trt_rows <- !(meta[[group_col]] %in% special) & !is.na(suffix_vals) & suffix_vals != ""
  roles$treatment <- sort(unique(paste0(meta[[group_col]][trt_rows], "_",
                                        suffix_vals[trt_rows])))
  cat(sprintf("   Split by '%s': %d treatment groups (e.g. %s)\n",
              group_suffix_col, length(roles$treatment),
              paste(head(roles$treatment, 4), collapse = ", ")))
}

cat(sprintf("   Control:    %s\n", if (is.null(roles$control)) "⚠️ NOT FOUND" else roles$control))
cat(sprintf("   QC (MS1):   %s\n", if (is.null(roles$qc_ms1)) "⚠️ NOT FOUND" else roles$qc_ms1))
cat(sprintf("   Blank:      %s\n", if (is.null(roles$blank)) "⚠️ NOT FOUND" else roles$blank))
cat(sprintf("   Treatment groups (%d): %s\n", length(roles$treatment),
            paste(head(roles$treatment, 10), collapse = ", "),
            if (length(roles$treatment) > 10) "..." else ""))

# Also detect MSMS samples
if (!is.null(col_map$ms_level)) {
  ms_level_col <- col_map$ms_level
  ms1_rows <- meta[[ms_level_col]] == "MS"
  ms2_rows <- meta[[ms_level_col]] == "MSMS"
  cat(sprintf("   MS1 samples: %d, MS2 samples: %d\n", sum(ms1_rows), sum(ms2_rows)))
}

# ══════════════════════════════════════════════════════════════════════════════
# 5. Generate group colors
# ══════════════════════════════════════════════════════════════════════════════

# Built-in palette for common roles + auto-assignment for treatments
generate_group_colors <- function(roles) {
  colors <- list()
  if (!is.null(roles$control)) colors[[roles$control]] <- "#58CDD9"
  if (!is.null(roles$qc_ms1))  colors[[roles$qc_ms1]]  <- "#BEBEBE"
  if (!is.null(roles$blank))   colors[[roles$blank]]   <- "#3C5488FF"

  # For treatments, cycle through ggsci-like colors
  trt_palette <- c("#D20A13", "#088247", "#FFD121", "#7E6148FF",
                   "#5BC0EB", "#F39B7FFF", "#BC3C29FF", "#0072B5FF",
                   "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF",
                   "#FFDC91FF", "#EE4C97FF", "#00A087FF", "#8491B4FF",
                   "#CC3333", "#8C510A", "#01665E", "#5C88DA",
                   "#A6CEE3", "#B2DF8A", "#FB9A99", "#FF7F00",
                   "#CAB2D6", "#33A02C", "#E31A1C", "#1F78B4",
                   "#B15928", "#FDBF6F", "#6A3D9A")
  for (i in seq_along(roles$treatment)) {
    g <- roles$treatment[[i]]
    if (is.null(colors[[g]])) {
      colors[[g]] <- trt_palette[(i - 1) %% length(trt_palette) + 1]
    }
  }
  colors
}

group_colors <- generate_group_colors(roles)

# ══════════════════════════════════════════════════════════════════════════════
# 6. Generate group labels
# ══════════════════════════════════════════════════════════════════════════════

generate_group_labels <- function(roles) {
  labels <- list()
  if (!is.null(roles$control)) labels[[roles$control]] <- "Control"
  if (!is.null(roles$qc_ms1))  labels[[roles$qc_ms1]]  <- "QC (MS1)"
  if (!is.null(roles$blank))   labels[[roles$blank]]   <- "Blank"
  for (g in roles$treatment) {
    labels[[g]] <- g  # generic: use group ID as label
  }
  labels
}

group_labels <- generate_group_labels(roles)

# ══════════════════════════════════════════════════════════════════════════════
# 7. Detect working directory and project name
# ══════════════════════════════════════════════════════════════════════════════

wd <- opt$`working-dir`
if (is.null(wd)) wd <- getwd()
cat(sprintf("=> Working directory: %s\n", wd))

project_name <- opt$`project-name`
if (is.null(project_name)) {
  # Derive from metadata parent dir name
  project_name <- basename(dirname(normalizePath(metadata_file)))
  if (project_name == "." || project_name == "") {
    project_name <- "Untitled_Metabolomics_Project"
  }
}
cat(sprintf("=> Project name: %s\n", project_name))

# ══════════════════════════════════════════════════════════════════════════════
# 8. Build and write YAML config
# ══════════════════════════════════════════════════════════════════════════════

# Determine output dir relative paths
raw_dir_rel <- opt$`raw-dir`
metadata_rel <- metadata_file
# If paths are absolute and under wd, make them relative
# (simple approach: just store as-is; pipeline resolves relative to working_dir)

# Determine output directories based on --finalOutput or metadata location
if (!is.null(opt$finalOutput)) {
  final_out <- opt$finalOutput
  output_dir_val  <- file.path(final_out, "detail")
  peak_table_val  <- file.path(final_out, "peak_table")
} else {
  output_dir_val  <- file.path(dirname(metadata_rel), paste0("Result_", project_name))
  peak_table_val  <- file.path(dirname(metadata_rel), "Result")
}

# Peak width range (seconds) for centWave/massprocesser. The old default
# c(5,30) was wider than every observed peak (0.7-3 s), risking merged
# isomers and missed narrow peaks. Override with --peakwidth "min,max".
peakwidth_val <- list(1, 8)
if (!is.null(opt$peakwidth)) {
  pw_parts <- suppressWarnings(as.numeric(strsplit(opt$peakwidth, ",")[[1]]))
  if (length(pw_parts) == 2 && all(is.finite(pw_parts)) &&
      pw_parts[1] > 0 && pw_parts[1] < pw_parts[2]) {
    peakwidth_val <- as.list(pw_parts)
    cat(sprintf("=> peakwidth override: %g-%g s\n", pw_parts[1], pw_parts[2]))
  } else {
    cat("WARNING: --peakwidth expects 'min,max' seconds; keeping default 1,8\n")
  }
}

config <- list(
  project = list(
    name        = project_name,
    batch_id    = if (!is.null(opt$`batch-id`)) opt$`batch-id` else NA_integer_,
    working_dir = wd,
    output_dir  = output_dir_val
  ),
  data = list(
    metadata_file  = metadata_rel,
    raw_mzml_dir   = raw_dir_rel,
    peak_table_dir = peak_table_val
  ),
  metadata_columns = list(
    filename    = if (!is.null(col_map$filename))  col_map$filename  else "file_name",
    ms_level    = if (!is.null(col_map$ms_level))  col_map$ms_level  else "MS_level",
    group_id    = if (!is.null(col_map$group_id))  col_map$group_id  else "drug_id",
    replicate   = if (!is.null(col_map$replicate)) col_map$replicate else "replicate_id",
    mode        = if (!is.null(col_map$mode))      col_map$mode      else "mode",
    # When set, the pipeline appends this column's value to each treatment
    # sample's class (10 → 10_High / 10_Low); control/QC/blank are unchanged.
    group_suffix = group_suffix_col
  ),
  sample_roles = list(
    control = if (!is.null(roles$control)) roles$control else "CT",
    qc_ms1  = if (!is.null(roles$qc_ms1))  roles$qc_ms1  else "QC",
    blank   = if (!is.null(roles$blank))   roles$blank   else "Blank"
  ),
  feature_extraction = list(
    polarity             = polarity,
    ppm                  = 15,
    peakwidth            = if (!is.null(opt$peakwidth)) as.numeric(strsplit(opt$peakwidth, ",")[[1]]) else c(10, 60),
    snthresh             = 5,
    noise                = 500,
    threads              = 10,
    output_tic           = FALSE,
    output_bpc           = FALSE,
    output_rt_correction_plot = FALSE,
    min_fraction         = 0.5,
    fill_peaks           = FALSE
  ),
  filtering = list(
    blank_fold_change   = 3.0,
    max_missing_per_group = 0.8,
    min_detected_per_group = 2,     # absolute floor: NAs beyond this are not rescued by imputation
    rsd_threshold       = 30.0,
    min_intensity       = 1000.0
  ),
  normalization = list(
    qc_min_samples  = 3,
    svr_threads     = 4,
    fallback_method = "median"
  ),
  performance = list(
    cores              = if (!is.null(opt$cores)) opt$cores else 0,  # 0 = auto-detect all cores
    parallel_backend   = "multisession",   # multisession | multicore | sequential
    parallel_setup_timeout = 30,
    io_engine          = "data.table",    # data.table | base
    serialize_engine   = "qs",             # qs | fst | rds
    enable_caching     = TRUE,
    cache_dir          = "cache",
    # --threads passed to the SIRIUS CLI. SIRIUS saturates whatever it is given,
    # and the batches run sequentially, so use most of the machine (cap 16 to
    # leave headroom for the R process and the OS).
    sirius_threads     = if (!is.null(opt$cores) && opt$cores > 0)
                           as.integer(opt$cores)
                         else max(1L, min(16L, parallel::detectCores() - 2L)),
    aggressive_gc      = TRUE
  ),
  differential = list(
    fc_threshold               = 1.2,
    p_value_cutoff             = 0.05,
    significance_metric        = "adj_p_value",  # "p_value" (raw) or "adj_p_value" (FDR)
    p_adjust_method            = "fdr",
    reference_group            = if (!is.null(roles$control)) roles$control else "CT",
    auto_generate_comparisons  = TRUE,
    extra_comparisons          = list()
  ),
  oplsda = list(
    enabled          = TRUE,
    n_permutations   = if (!is.null(opt$oplsda_permutations)) opt$oplsda_permutations else 1000,
    vip_threshold    = 1.0,
    comparison_subset = list()  # empty = all
  ),
  overlap_analysis = list(
    enabled                    = TRUE,
    venn_max_sets              = 5,
    upset_max_sets             = 20,
    min_comparisons_per_metabolite = 2
  ),
  kegg = list(
    # tidymass enrich_hmdb runs against the HMDB/SMPDB collection; the
    # "kegg" folder name is legacy. Filter on BH-adjusted q-values.
    q_value_cutoff             = 0.05
  ),
  metaboanalyst_advanced = list(
    # Requires the MetaboAnalystR package (installed here: 4.2.0).
    # MSEA/topology download reference data from the MetaboAnalyst server on
    # first use; every submodule is wrapped in tryCatch so failures degrade
    # to skips, not pipeline aborts. Set enabled=false to skip Step 10 entirely.
    enabled = TRUE,
    msea = list(
      enabled = TRUE, method = "globaltest", stat_type = "tstat",
      pathway_db = "kegg", n_top_sets = 20
    ),
    pathway_topology = list(
      enabled = TRUE, metric = "rbc", pathway_db = "kegg", combine_queries = TRUE
    ),
    biomarker = list(
      enabled = TRUE, rf_ntree = 500, rf_mtry = NULL, rf_importance = "accuracy",
      n_top_features = 20, roc_threshold = 0.7, roc_ci_method = "delong",
      roc_n_bootstrap = 2000
    ),
    chem_class_enrichment = list(
      enabled = TRUE, taxon_level = "super_class", method = "fisher",
      p_cutoff = 0.05, p_adjust_method = "fdr"
    )
  ),
  annotation = list(
    column           = "rp",
    ms1_db_dir      = if (!is.null(opt$`ms1-db-dir`)) opt$`ms1-db-dir` else "/home/data/shareData/16tdisk/Project/tidymass2/database/MS1_Database",
    ms2_db_dir      = if (!is.null(opt$`ms2-db-dir`)) opt$`ms2-db-dir` else "/home/data/shareData/16tdisk/Project/tidymass2/database/MS2_Compound_database",
    inhouse_db_path = "/home/data/shareData/16tdisk/Project/tidymass2/database/MS1_Database/inhouse_Metabolite.database",
    # The in-house DB's mz column does not match its own Formula column
    # (3/295 rows consistent, off by up to +69 Da), so as priority #1 it would
    # overwrite correct HMDB/KEGG hits with wrong names. Keep it disabled until
    # the mz/mz.neg/mz.pos columns are rebuilt from Formula.
    inhouse_enabled = FALSE,
    ms1_match_ppm   = 10,
    ms2_match_ppm   = 30,
    ms2_match_tol   = 0.5,
    rt_tol_inhouse  = 60,
    rt_tol_external = 180,
    # Scoring weights for MS1+MS2 joint matching (only used for MS2 databases)
    ms1_match_weight = 0.25,
    rt_match_weight  = 0.25,
    ms2_match_weight = 0.5,
    total_score_tol  = 0.35,
    candidate_num    = 3,
    db_priority     = list("MS2_HMDB", "MS2_GNPS", "MS2_GNPS_IH",
                           "MS2_GNPS_NIH", "MS2_MassBank", "MS2_MoNA",
                           "HMDB", "KEGG", "ChEBI", "FooDB", "PubChem")
  ),
  # ── CAMERA: Adduct & isotope grouping ──────────────────────────────────
  camera = list(
    enabled     = FALSE,    # Requires Bioconductor packages: CAMERA, xcms
    polarity    = polarity,
    perfwhm     = 0.6,
    ppm         = 15,
    cor_eic_th  = 0.75,
    graphMethod = "lpc",
    pval        = 0.05,
    calcCiS     = TRUE,
    calcCaS     = FALSE
  ),
  # ── SIRIUS: Molecular formula & structure prediction ──────────────────
  sirius = list(
    # SIRIUS 6.3 CLI requires a licensed login (`sirius login`). Without it the
    # binary prints "Login ERROR", writes no project, and exits 0 — so the
    # pipeline reports success while contributing zero annotations. A valid
    # academic token is cached in ~/.sirius-6.3 on this machine (verified
    # 2026-09-23), so this is enabled by default. Re-run `sirius login` if the
    # token expires.
    enabled       = TRUE,
    path          = sirius_path,
    project_space = if (!is.null(opt$`sirius-project`)) opt$`sirius-project` else file.path(output_dir_val, "sirius_project"),
    ppm_max       = 15,
    # CSI:FingerID structure databases. Natural-product collections
    # (PLANTCYC, KNAPSACK, GNPS) matter for TCM herbs; HMDB/CHEBI are
    # human-centric and cover few plant secondary metabolites.
    #
    # NOTE: the plant collections are NOT bundled with SIRIUS 6.3 — they must be
    # downloaded once into the workspace (`sirius custom-db-downloader`), or
    # `structures -d PLANTCYC` finds zero compounds and, with `--recompute`, can
    # stall for hours. See sirius_available_databases() in the pipeline: the
    # requested list is filtered against what is actually installed.
    structure_db  = "HMDB,CHEBI",
    # SIRIUS skips already-computed results, so changing structure_db alone has
    # no effect. Set true for one run to force the search to be redone.
    recompute_structure = FALSE,
    # Hard wall-clock cap per SIRIUS invocation. 6.3.7 sometimes finishes writing
    # the project and then never exits; without this the batch hangs for hours.
    timeout_sec = 3600
  ),
  # ── MetFrag: In silico fragmentation validation ───────────────────────
  metfrag = list(
    enabled  = FALSE,     # Requires metfRag R package + compatible Java (≤17)
    ppm_tol  = 15,
    mzabs    = 0.005
  ),
  # ── Internal standards: isotope-labeled amino acid mix (2 uM, extraction-solvent spike) ──
  # CSV lives next to the batch metadata under <project root>/code/metadata/;
  # resolve it here (at config-gen time, relative to make_config.R's own dir)
  # because the pipeline setwd()s to working_dir before reading it.
  internal_standards = list(
    enabled     = TRUE,
    file        = is_csv_path,
    mz_tol_ppm  = 15,
    rt_tol_sec  = 20
  ),
  visualization = list(
    dpi              = 300,
    group_colors     = group_colors,
    group_labels     = group_labels,
    heatmap_top_n    = 20,
    boxplot_top_n    = 30
  )
)

# Handle yaml package differences: as.yaml needs column.major = FALSE for clean output.
# The logical handler must return a "verbatim"-classed string so booleans are emitted
# UNQUOTED (true / false). A plain character return gets quoted ('true'), which read_yaml
# then parses back as a STRING — and isTRUE("true") is FALSE, silently disabling every
# isTRUE()-gated step (OPLS-DA, auto_generate_comparisons, overlap_analysis).
yaml_str <- as.yaml(config, column.major = FALSE, indent = 2,
                    handlers = list(
                      logical = function(x) {
                        v <- ifelse(x, "true", "false")
                        class(v) <- "verbatim"
                        v
                      }
                    ))

# Write output
output_path <- opt$output
dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
writeLines(c(
  "# =============================================================================",
  sprintf("# Metabolomics Pipeline Configuration — %s", project_name),
  sprintf("# Auto-generated by make_config.R on %s", Sys.Date()),
  "# =============================================================================",
  "# Edit this file before running the pipeline.",
  "# Lines with null values need your manual input.",
  "# =============================================================================",
  ""
), output_path)
write(yaml_str, file = output_path, append = TRUE)

cat(sprintf("\n===== Config written to: %s =====\n", normalizePath(output_path)))
cat(sprintf("  Project: %s\n", project_name))
cat(sprintf("  Groups:  %d total (%d treatment + control/qc/blank)\n",
            length(all_groups), length(roles$treatment)))
cat(sprintf("  Polarity: %s\n", polarity))
cat("\nFeatures enabled:\n")
cat(sprintf("  SIRIUS:     %s (%s)\n",
            if (isTRUE(config$sirius$enabled)) "YES" else "NO", config$sirius$path))
cat(sprintf("  MetFrag:    %s\n", if (isTRUE(config$metfrag$enabled)) "YES" else "NO (requires metfRag + Java <=17)"))
cat(sprintf("  CAMERA:     %s (requires CAMERA/xcms Bioconductor)\n", if (isTRUE(config$camera$enabled)) "YES" else "NO"))
cat(sprintf("  OPLS-DA:    %s\n", if (isTRUE(config$oplsda$enabled)) "YES" else "NO"))
cat(sprintf("  MSI Levels: YES (Level 1-4 confidence scoring)\n"))
cat("\nNext steps:\n")
cat(sprintf("  1. Review and edit: %s\n", output_path))
cat(sprintf("  2. Run pipeline: Rscript metabolomics_pipeline.R --config %s\n", output_path))
