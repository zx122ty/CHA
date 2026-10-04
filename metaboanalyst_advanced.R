# ══════════════════════════════════════════════════════════════════════════════
# 9.5. MetaboAnalystR Advanced Analysis (Step 10)
# ══════════════════════════════════════════════════════════════════════════════
#
# Integrates four advanced downstream analytical methods from MetaboAnalystR:
#   1. MSEA (Metabolite Set Enrichment Analysis) — QSEA/GlobalTest using
#      all-feature continuous statistics, not just significant features
#   2. Pathway Topology Analysis — KEGG graph topology (Relative Betweenness
#      Centrality or Impact Centrality) combined with enrichment p-values
#   3. Biomarker Analysis — Random Forest feature selection + ROC curves
#   4. Chemical Class Enrichment — HMDB Super Class / ClassyFire taxonomy
#      enrichment
#
# All wrapped in tryCatch for graceful degradation, using the same logging
# style and output conventions as the parent pipeline.
# ══════════════════════════════════════════════════════════════════════════════

#' Convert mass_dataset to MetaboAnalystR mSet object
#'
#' CRITICAL STATE MANAGEMENT:
#' The input data (object2) has ALREADY been through:
#'   1. Blank background filtering         (Step 2a)
#'   2. Missing value filter                (Step 2b)
#'   3. QC-RSD filter                       (Step 2c)
#'   4. Minimum-value imputation            (Step 2d)
#'   5. SVR (or LOESS/median) normalization (Step 2d)
#'   6. log2 transformation                 (Step 4, in run_differential)
#'
#' Therefore we MUST initialize the mSet object with the pre-processed matrix
#' and tell MetaboAnalystR to skip ALL of its internal normalization and
#' filtering steps. This is achieved by:
#'   - Using InitDataObjects("pktable", "stat") — pktable mode bypasses raw
#'     spectral processing
#'   - Pre-populating mSet$dataSet$norm with the log2-transformed matrix
#'     and setting mSet$dataSet$preproc = TRUE so Read.TextData skips
#'     normalization
#'   - Setting mSet$dataSet$filt to bypass feature filtering
#'   - After Read.TextData, we OVERWRITE mSet$dataSet$norm with our matrix
#'     and set mSet$dataSet$proc and mSet$dataSet$preproc flags to TRUE
#'
#' @param object2 mass_dataset object (filtered, imputed, SVR-normalized)
#' @param sample_info data.frame with sample_id and class columns
#' @param diff_results Output of run_differential(), for all-feature statistics
#' @param cfg Validated config list
#' @return List with mSet (MetaboAnalystR object), success (logical), and
#'         the mapping table used for feature-to-ID conversion
convert_mass_dataset_to_mSet <- function(object2, sample_info, diff_results, cfg) {
  cat("\n----- Converting mass_dataset to MetaboAnalystR mSetObj -----\n")

  # ── Check MetaboAnalystR availability ──────────────────────────────────
  if (!requireNamespace("MetaboAnalystR", quietly = TRUE)) {
    cat("   MetaboAnalystR package not installed. Skipping.\n")
    cat("   Install with: install.packages('MetaboAnalystR',\n")
    cat("     repos = c('https://cran.r-project.org', 'https://bioconductor.org'))\n")
    return(list(mSet = NULL, success = FALSE, id_map = NULL))
  }
  suppressPackageStartupMessages(library(MetaboAnalystR))

  # ── Extract the log2-transformed expression matrix ─────────────────────
  # The expression data is already SVR-normalized in object2.
  # We log2-transform it (consistent with run_differential) and use it
  # as the "normalized" matrix that MetaboAnalystR will operate on.
  expr_mat <- object2@expression_data
  # Ensure no zeros or negatives before log2 (imputation with "minimum" in
  # Step 2d guarantees all values > 0, but double-check for safety)
  if (any(expr_mat <= 0, na.rm = TRUE)) {
    cat("   Warning: found <=0 values in expression data, adding pseudo-count\n")
    min_val <- min(expr_mat[expr_mat > 0], na.rm = TRUE)
    expr_mat[expr_mat <= 0] <- min_val * 0.5
  }
  log2_mat <- log2(as.matrix(expr_mat))
  # Transpose: rows = samples, cols = features (MetaboAnalystR convention)
  log2_mat_t <- t(log2_mat)
  n_samples <- nrow(log2_mat_t)
  n_features <- ncol(log2_mat_t)
  cat(sprintf("   Expression matrix: %d samples x %d features\n", n_samples, n_features))

  # ── Build class labels ─────────────────────────────────────────────────
  sample_order <- colnames(object2@expression_data)
  class_vec <- sample_info$class[match(sample_order, sample_info$sample_id)]
  # Replace any NA class with "Unknown"
  class_vec[is.na(class_vec)] <- "Unknown"
  # Convert to MetaboAnalystR format: clean, factorized
  class_vec <- make.names(class_vec)              # remove special chars
  class_vec <- factor(class_vec, levels = unique(class_vec))
  cat(sprintf("   Sample classes: %s\n", paste(levels(class_vec), collapse = ", ")))

  # ── Build feature ID mapping table ─────────────────────────────────────
  # Map variable_id -> HMDB.ID, KEGG.ID, Compound.name from annotation
  # This is used later by MSEA and topology analysis
  cat("   Building feature-ID mapping from annotation table...\n")
  id_map <- data.frame(
    variable_id = rownames(log2_mat),
    stringsAsFactors = FALSE
  )
  # Try to load app3 annotation if available (passed from parent pipeline)
  # We'll populate enrichment IDs later in the calling function

  # ── Initialize mSet object ─────────────────────────────────────────────
  # Use "pktable" type to bypass raw spectral processing
  mSet <- tryCatch({
    MetaboAnalystR::InitDataObjects("pktable", "stat", FALSE, default.dpi = 100)
  }, error = function(e) {
    cat(sprintf("   ERROR: InitDataObjects failed: %s\n", conditionMessage(e)))
    return(NULL)
  })
  if (is.null(mSet)) return(list(mSet = NULL, success = FALSE, id_map = id_map))

  # ── Prepare the data matrix for MetaboAnalystR ─────────────────────────
  # MetaboAnalystR's Read.TextData expects a CSV-like format:
  #   - Column 1: sample labels
  #   - Column 2: class labels
  #   - Remaining columns: feature intensities
  # We'll write a temporary CSV and use Read.TextData, then override
  # the internal normalization flags.
  require(utils)

  # Build the data frame: rows = samples, first col = label, second = class
  df_for_ma <- data.frame(
    label = sample_order,
    class = as.character(class_vec),
    log2_mat_t,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Write to temp file
  tmp_csv <- tempfile(pattern = "metaboanalyst_input_", fileext = ".csv")
  write.csv(df_for_ma, file = tmp_csv, row.names = FALSE, quote = FALSE)
  on.exit(unlink(tmp_csv), add = TRUE)

  # ── Read data into mSet ────────────────────────────────────────────────
  cat("   Reading data into MetaboAnalystR...\n")
  mSet <- tryCatch({
    MetaboAnalystR::Read.TextData(mSet, tmp_csv, "colu", "disc")
  }, error = function(e) {
    cat(sprintf("   ERROR: Read.TextData failed: %s\n", conditionMessage(e)))
    return(NULL)
  })
  if (is.null(mSet)) return(list(mSet = NULL, success = FALSE, id_map = id_map))

  # ── CRITICAL: Override normalization to prevent double-normalization ──
  # MetaboAnalystR's Read.TextData with "colu" format performs some
  # internal processing. We OVERWRITE the normalized data matrix with
  # our pre-processed log2 matrix and set flags to skip all normalization.
  # This is the key trick to feed pre-processed data into MetaboAnalystR.
  cat("   Overriding mSet normalization state (data already pre-processed)...\n")

  mSet <- tryCatch({
    # Sanity check: the data dimensions should match
    if (nrow(mSet$dataSet$orig) != n_samples) {
      cat(sprintf("   WARNING: sample count mismatch (orig=%d, expected=%d)\n",
                  nrow(mSet$dataSet$orig), n_samples))
    }

    # Overwrite the normalized data with our log2 matrix
    # mSet$dataSet$norm is the slot MetaboAnalystR uses for all downstream
    mSet$dataSet$norm <- log2_mat_t
    colnames(mSet$dataSet$norm) <- rownames(log2_mat)
    rownames(mSet$dataSet$norm) <- sample_order

    # Set the preprocessing flags to indicate data is already normalized
    mSet$dataSet$preproc <- TRUE
    mSet$dataSet$proc <- TRUE
    mSet$dataSet$filt <- TRUE
    mSet$dataSet$row.sel <- TRUE
    mSet$dataSet$col.sel <- TRUE

    # Copy the class labels
    mSet$dataSet$cls <- class_vec
    mSet$dataSet$orig.cls <- class_vec

    # Set the number of samples and features
    mSet$dataSet$mu.num <- n_features
    mSet$dataSet$filt.num <- n_features
    mSet$dataSet$proc.num <- n_features

    cat("   mSetObj normalization state overridden successfully\n")
    mSet
  }, error = function(e) {
    cat(sprintf("   ERROR overriding normalization state: %s\n", conditionMessage(e)))
    return(NULL)
  })

  if (is.null(mSet)) return(list(mSet = NULL, success = FALSE, id_map = id_map))

  cat("   mSetObj conversion complete\n")
  list(mSet = mSet, success = TRUE, id_map = id_map)
}


#' Run MSEA (Metabolite Set Enrichment Analysis) via MetaboAnalystR
#'
#' Implements Quantitative Enrichment Analysis (QSEA) using GlobalTest or
#' globalANCOVA, which uses continuous statistics (limma t-statistics or
#' log2FC) for ALL detected features, not just significant ones. This is
#' more powerful than the simple ORA in Step 7.
#'
#' @param mSet A MetaboAnalystR mSet object (from convert_mass_dataset_to_mSet)
#' @param id_map Data.frame mapping variable_id -> HMDB/KEGG IDs
#' @param app3 Annotation table from run_annotation()
#' @param diff_results Output of run_differential() for t-statistics
#' @param cfg Config list (metaboanalyst_advanced$msea sub-block)
#' @param out_dir Output directory for this step
#' @return Updated mSet (or NULL on failure)
run_msea_analysis <- function(mSet, id_map, app3, diff_results, cfg, out_dir) {
  cat("\n----- MSEA: Metabolite Set Enrichment Analysis -----\n")

  if (!requireNamespace("MetaboAnalystR", quietly = TRUE)) {
    cat("   MetaboAnalystR not available, skipping MSEA\n")
    return(NULL)
  }

  msea_cfg <- cfg$metaboanalyst_advanced$msea
  method <- msea_cfg$method %||% "globaltest"  # globaltest, ora
  stat_type <- msea_cfg$stat_type %||% "tstat" # tstat, log2fc
  pathway_db <- msea_cfg$pathway_db %||% "kegg"
  n_top <- msea_cfg$n_top_sets %||% 20
  p_cutoff <- cfg$differential$p_value_cutoff %||% 0.05

  # ── Build the compound-to-pathway ID mapping ───────────────────────────
  # We need a mapping from the features in the mSet to either HMDB or KEGG IDs
  cat("   Building compound-to-pathway ID mapping...\n")

  feature_ids <- rownames(mSet$dataSet$norm)
  if (is.null(app3) || nrow(app3) == 0) {
    cat("   No annotation table available, cannot map features to pathway IDs\n")
    cat("   MSEA requires HMDB or KEGG IDs from annotation\n")
    return(NULL)
  }

  # Build the mapping: variable_id -> HMDB.ID (preferred) or KEGG.ID
  id_mapping <- app3 %>%
    dplyr::select(variable_id, HMDB.ID, KEGG.ID, Compound.name) %>%
    filter(!is.na(HMDB.ID) & HMDB.ID != "" | !is.na(KEGG.ID) & KEGG.ID != "")

  if (nrow(id_mapping) == 0) {
    cat("   No features have HMDB or KEGG IDs, cannot run MSEA\n")
    cat("   MSEA requires at least some annotated features with pathway IDs\n")
    return(NULL)
  }

  cat(sprintf("   %d features with pathway IDs (HMDB/KEGG)\n", nrow(id_mapping)))

  # ── Prepare the compound list file for MetaboAnalystR ──────────────────
  # MetaboAnalystR's MSEA expects a compound list (one per line) with
  # a statistic (t-stat or log2FC) for each compound.
  # We need to get the t-statistics or log2FC from the differential results.
  # For multi-comparison studies, we take the first comparison or average.

  cat(sprintf("   Computing all-feature statistics (%s)...\n", stat_type))

  all_stats <- NULL
  if (length(diff_results) > 0) {
    # Use the first available comparison's statistics
    first_cmp <- diff_results[[1]]
    if (!is.null(first_cmp$all)) {
      diff_all <- first_cmp$all
      diff_all$variable_id <- rownames(diff_all)

      if (stat_type == "tstat") {
        # Use the moderated t-statistic: t = B / (s * sqrt(v))
        # In limma's topTable, this is the 't' column
        t_col <- if ("t" %in% colnames(diff_all)) "t" else "B"
        diff_all$stat <- diff_all[[t_col]]
      } else {
        # Use log2 fold change
        diff_all$stat <- diff_all$logFC
      }

      # Merge with ID mapping
      all_stats <- id_mapping %>%
        left_join(diff_all %>% dplyr::select(variable_id, stat),
                  by = "variable_id") %>%
        filter(!is.na(stat) & is.finite(stat))
    }
  }

  if (is.null(all_stats) || nrow(all_stats) == 0) {
    cat("   Could not compute differential statistics, falling back to ORA\n")
    method <- "ora"
  }

  # ── Run MSEA for each comparison ───────────────────────────────────────
  for (cmp_name in names(diff_results)) {
    dr <- diff_results[[cmp_name]]
    if (is.null(dr)) next

    label <- dr$comparison$label %||% cmp_name
    cat(sprintf("   Processing: %s\n", label))

    # Get the t-statistics for this specific comparison
    diff_all <- dr$all
    diff_all$variable_id <- rownames(diff_all)

    if (stat_type == "tstat") {
      t_col <- if ("t" %in% colnames(diff_all)) "t" else "B"
      diff_all$stat <- diff_all[[t_col]]
    } else {
      diff_all$stat <- diff_all$logFC
    }

    cmp_stats <- id_mapping %>%
      left_join(diff_all %>% dplyr::select(variable_id, stat),
                by = "variable_id") %>%
      filter(!is.na(stat) & is.finite(stat))

    if (nrow(cmp_stats) == 0) {
      cat("     No features with statistics for this comparison, skipping\n")
      next
    }

    # ── Write the compound list file for MetaboAnalystR ──────────────────
    # For GlobalTest/QSEA: we need a two-column file:
    #   CompoundID \t statistic
    # For ORA: we need a one-column list of compound IDs (significant only)
    cmp_dir <- file.path(out_dir, cmp_name)
    dir.create(cmp_dir, showWarnings = FALSE, recursive = TRUE)

    if (method == "globaltest" || method == "qsea") {
      # Use all features with continuous statistics
      # Preference: HMDB.ID > KEGG.ID
      cmp_stats <- cmp_stats %>%
        mutate(
          pathway_id = ifelse(!is.na(HMDB.ID) & HMDB.ID != "",
                              HMDB.ID, KEGG.ID)
        ) %>%
        filter(!is.na(pathway_id) & pathway_id != "")

      # For GlobalTest, we need unique compound-by-statistic
      # Take the median stat for compounds appearing multiple times
      compound_stats <- cmp_stats %>%
        group_by(pathway_id) %>%
        summarise(stat = median(stat, na.rm = TRUE), .groups = "drop") %>%
        filter(!is.na(stat))

      if (nrow(compound_stats) < 5) {
        cat("     Too few mapped compounds (<5), skipping\n")
        next
      }

      # Write as two-column format
      cmpd_file <- file.path(cmp_dir, "compound_statistics.txt")
      write.table(compound_stats, file = cmpd_file,
                  sep = "\t", row.names = FALSE, col.names = FALSE,
                  quote = FALSE)
      cat(sprintf("     %d compounds with statistics written to: %s\n",
                  nrow(compound_stats), basename(cmp_dir)))

      # ── Run GlobalTest via MetaboAnalystR ─────────────────────────────
      mSet <- tryCatch({
        mSet <- MetaboAnalystR::InitDataObjects("conc", "mset.qsea", FALSE, default.dpi = 100)
        mSet <- MetaboAnalystR::Setup.MapData(mSet, cmpd_file)
        mSet <- MetaboAnalystR::CrossReferencing(mSet, pathway_db)
        mSet <- MetaboAnalystR::CreateMappingResultTable(mSet)

        # Run QSEA/GlobalTest
        mSet <- MetaboAnalystR::SetMsetType(mSet, pathway_db)
        mSet <- MetaboAnalystR::SetMetabolomeFilter(mSet, FALSE)
        mSet <- MetaboAnalystR::CalculateGlobalTest(mSet)
        mSet
      }, error = function(e) {
        cat(sprintf("     ERROR in GlobalTest for %s: %s\n", label, conditionMessage(e)))
        return(NULL)
      })

    } else {
      # ── ORA mode (fallback) ───────────────────────────────────────────
      # Use only significant features
      sig_features <- cmp_stats %>%
        filter(abs(stat) >= log2(cfg$differential$fc_threshold %||% 1.2)) %>%
        mutate(
          pathway_id = ifelse(!is.na(HMDB.ID) & HMDB.ID != "",
                              HMDB.ID, KEGG.ID)
        ) %>%
        filter(!is.na(pathway_id) & pathway_id != "")

      # Take unique compounds
      sig_compounds <- unique(sig_features$pathway_id)

      if (length(sig_compounds) < 3) {
        cat("     Too few significant compounds (<3), skipping\n")
        next
      }

      cmpd_file <- file.path(cmp_dir, "significant_compounds.txt")
      writeLines(sig_compounds, cmpd_file)
      cat(sprintf("     %d significant compounds written\n", length(sig_compounds)))

      mSet <- tryCatch({
        mSet <- MetaboAnalystR::InitDataObjects("conc", "mset.ora", FALSE, default.dpi = 100)
        mSet <- MetaboAnalystR::Setup.MapData(mSet, cmpd_file)
        mSet <- MetaboAnalystR::CrossReferencing(mSet, pathway_db)
        mSet <- MetaboAnalystR::CreateMappingResultTable(mSet)
        mSet <- MetaboAnalystR::SetMsetType(mSet, pathway_db)
        mSet <- MetaboAnalystR::SetMetabolomeFilter(mSet, FALSE)
        mSet <- MetaboAnalystR::CalculateOra(mSet)
        mSet
      }, error = function(e) {
        cat(sprintf("     ERROR in ORA for %s: %s\n", label, conditionMessage(e)))
        return(NULL)
      })
    }

    if (is.null(mSet)) next

    # ── Extract and save results ─────────────────────────────────────────
    msea_result <- tryCatch({
      MetaboAnalystR::GetMsetResult(mSet)
    }, error = function(e) NULL)

    if (!is.null(msea_result) && nrow(msea_result) > 0) {
      # Add hit ratio column
      msea_result$hit_ratio <- msea_result$Hit / msea_result$Total
      msea_result <- msea_result[order(msea_result$pval), , drop = FALSE]

      write.csv(msea_result,
                file.path(cmp_dir, "msea_results.csv"),
                row.names = FALSE)
      cat(sprintf("     MSEA results: %d pathway sets, %d significant (p<%.3g)\n",
                  nrow(msea_result),
                  sum(msea_result$pval < p_cutoff, na.rm = TRUE),
                  p_cutoff))

      # ── Plot top N sets ───────────────────────────────────────────────
      plot_data <- msea_result %>%
        head(n_top) %>%
        mutate(pathway = factor(Pathway, levels = rev(unique(Pathway))))

      p <- ggplot(plot_data, aes(x = hit_ratio, y = pathway)) +
        geom_point(aes(size = -log10(pval), color = -log10(pval))) +
        scale_color_gradientn(colours = c("#46bac2", "#f7ca64", "#e34a33"),
                              name = "-log10(p)") +
        scale_size_continuous(name = "-log10(p)") +
        labs(x = "Hit Ratio", y = "",
             title = sprintf("MSEA — %s (%s)", label,
                             ifelse(method == "globaltest", "GlobalTest", "ORA"))) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              plot.title = element_text(hjust = 0.5, size = 12),
              axis.text.y = element_text(size = rel(1.1)))

      save_plot(file.path(cmp_dir, "msea_bubble"),
                plot = p, width = 10, height = max(5, nrow(plot_data) * 0.4))
    } else {
      cat("     No MSEA results returned\n")
    }
  }

  cat("   MSEA analysis complete\n")
  mSet
}


#' Run Pathway Topology Analysis via MetaboAnalystR
#'
#' Combines enrichment p-values with KEGG pathway topology measures
#' (Relative Betweenness Centrality or Impact Centrality) to identify
#' the most biologically relevant pathways.
#'
#' @param mSet A MetaboAnalystR mSet object (from MSEA step, or fresh)
#' @param app3 Annotation table
#' @param diff_results Output of run_differential()
#' @param cfg Config list
#' @param out_dir Output directory for this step
#' @return NULL (results saved to files)
run_pathway_topology <- function(mSet, app3, diff_results, cfg, out_dir) {
  cat("\n----- Pathway Topology Analysis -----\n")

  if (!requireNamespace("MetaboAnalystR", quietly = TRUE)) {
    cat("   MetaboAnalystR not available, skipping\n")
    return(NULL)
  }

  pt_cfg <- cfg$metaboanalyst_advanced$pathway_topology
  metric <- pt_cfg$metric %||% "rbc"
  pathway_db <- pt_cfg$pathway_db %||% "kegg"
  p_cutoff <- cfg$differential$p_value_cutoff %||% 0.05

  if (is.null(app3) || nrow(app3) == 0) {
    cat("   No annotations available, skipping\n")
    return(NULL)
  }

  # ── Build compound ID list from annotation ─────────────────────────────
  # For KEGG pathway topology, we need KEGG compound IDs
  cat("   Building KEGG compound ID list from annotations...\n")

  for (cmp_name in names(diff_results)) {
    dr <- diff_results[[cmp_name]]
    if (is.null(dr)) next

    label <- dr$comparison$label %||% cmp_name
    cat(sprintf("   Processing: %s\n", label))

    # Get significant features and merge with annotations
    diff_sig <- dr$sig
    if (is.null(diff_sig) || nrow(diff_sig) == 0) {
      cat("     No significant features, skipping\n")
      next
    }
    diff_sig$variable_id <- rownames(diff_sig)

    annotated <- app3 %>%
      inner_join(diff_sig %>% dplyr::select(variable_id),
                 by = "variable_id") %>%
      distinct(Compound.name, .keep_all = TRUE)

    # Extract KEGG IDs (prefer KEGG over HMDB for topology)
    kegg_ids <- annotated$KEGG.ID
    kegg_ids <- kegg_ids[!is.na(kegg_ids) & kegg_ids != ""]
    # Clean KEGG IDs: ensure they start with "C" (KEGG compound prefix)
    kegg_ids <- unique(gsub("^.*?(C\\d+)", "\\1", kegg_ids))
    kegg_ids <- kegg_ids[grepl("^C\\d+", kegg_ids)]

    if (length(kegg_ids) < 3) {
      cat(sprintf("     Only %d valid KEGG IDs, falling back to HMDB IDs\n",
                  length(kegg_ids)))
      # Fall back to HMDB IDs
      hmdb_ids <- annotated$HMDB.ID
      hmdb_ids <- unique(hmdb_ids[!is.na(hmdb_ids) & hmdb_ids != ""])
      if (length(hmdb_ids) < 3) {
        cat("     Too few pathway IDs, skipping\n")
        next
      }
      compound_ids <- hmdb_ids
      id_type <- "hmdb"
    } else {
      compound_ids <- kegg_ids
      id_type <- "kegg"
    }

    # ── Write compound list and run topology analysis ────────────────────
    cmp_dir <- file.path(out_dir, cmp_name)
    dir.create(cmp_dir, showWarnings = FALSE, recursive = TRUE)
    cmpd_file <- file.path(cmp_dir, "compound_ids.txt")
    writeLines(compound_ids, cmpd_file)

    cat(sprintf("     %d compounds for topology analysis (%s)\n",
                length(compound_ids), id_type))

    # Initialize MetaboAnalystR for pathway analysis
    mSet_pt <- tryCatch({
      pt <- MetaboAnalystR::InitDataObjects("conc", "pathqea", FALSE, default.dpi = 100)
      pt <- MetaboAnalystR::Setup.MapData(pt, cmpd_file)

      # Cross-reference with the specified pathway database
      pt <- MetaboAnalystR::CrossReferencing(pt, pathway_db)  # "kegg" or "smpdb"
      pt <- MetaboAnalystR::CreateMappingResultTable(pt)

      # Set the pathway analysis parameters
      pt <- MetaboAnalystR::SetMsetType(pt, pathway_db)
      pt <- MetaboAnalystR::SetMetabolomeFilter(pt, FALSE)

      # Run the pathway analysis with topology measure
      pt <- MetaboAnalystR::CalculatePathwayScore(
        pt, method = metric   # "rbc" or "impact"
      )
      pt
    }, error = function(e) {
      cat(sprintf("     ERROR: %s\n", conditionMessage(e)))
      return(NULL)
    })

    if (is.null(mSet_pt)) next

    # ── Extract and save results ─────────────────────────────────────────
    path_result <- tryCatch({
      MetaboAnalystR::GetPathwayResult(mSet_pt)
    }, error = function(e) NULL)

    if (!is.null(path_result) && nrow(path_result) > 0) {
      # Add a combined score: -log10(p) * topology
      if ("pval" %in% colnames(path_result) &&
          metric %in% colnames(path_result)) {
        path_result$combined_score <- -log10(path_result$pval + 1e-10) *
          path_result[[metric]]
      }

      path_result <- path_result[order(path_result$pval), , drop = FALSE]
      write.csv(path_result,
                file.path(cmp_dir, "pathway_topology_results.csv"),
                row.names = FALSE)
      cat(sprintf("     %d pathways with topology scores\n", nrow(path_result)))

      # ── Plot: enrichment vs topology scatter ───────────────────────────
      if (nrow(path_result) >= 3) {
        plot_data <- path_result %>%
          mutate(
            significant = pval < p_cutoff,
            label = ifelse(rank(pval) <= 10 | rank(rev(sort(.data[[metric]]))) <= 5,
                           Pathway, "")
          )

        p <- ggplot(plot_data, aes(x = .data[[metric]],
                                    y = -log10(pval + 1e-10))) +
          geom_point(aes(size = Hit, color = significant), alpha = 0.7) +
          geom_text_repel(aes(label = label), size = 3,
                          max.overlaps = 15, box.padding = 0.4) +
          scale_color_manual(values = c("TRUE" = "#D20A13", "FALSE" = "#5BC0EB")) +
          labs(x = sprintf("Pathway Topology (%s)", metric),
               y = "-log10(p-value)",
               title = sprintf("Pathway Topology — %s", label)) +
          theme_bw() +
          theme(panel.grid = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 12))

        save_plot(file.path(cmp_dir, "pathway_topology_scatter"),
                  plot = p, width = 10, height = 8)

        # ── Topology bar plot (top 15) ──────────────────────────────────
        top15 <- plot_data %>%
          arrange(desc(combined_score)) %>%
          head(15) %>%
          mutate(Pathway = factor(Pathway, levels = rev(unique(Pathway))))

        p2 <- ggplot(top15, aes(x = combined_score, y = Pathway)) +
          geom_bar(stat = "identity", aes(fill = -log10(pval + 1e-10))) +
          scale_fill_gradientn(colours = c("#46bac2", "#f7ca64", "#e34a33"),
                               name = "-log10(p)") +
          labs(x = "Combined Score (-log10(p) × topology)", y = "",
               title = sprintf("Top Pathways — %s", label)) +
          theme_bw() +
          theme(panel.grid = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 12))

        save_plot(file.path(cmp_dir, "pathway_topology_bar"),
                  plot = p2, width = 10, height = max(5, nrow(top15) * 0.4))
      }
    } else {
      cat("     No pathway topology results returned\n")
    }
  }

  cat("   Pathway topology analysis complete\n")
  invisible(NULL)
}


#' Run Biomarker Analysis (Random Forest + ROC curves) via MetaboAnalystR
#'
#' Implements:
#'   1. Random Forest classification for feature importance ranking
#'   2. ROC curve analysis for top N biomarker candidates
#'   3. AUC calculation with confidence intervals
#'
#' For each comparison, the top-ranked features from RF are evaluated
#' individually for diagnostic potential via ROC curves.
#'
#' @param mSet MetaboAnalystR mSet object (from conversion)
#' @param object2 mass_dataset object (for expression data)
#' @param sample_info Sample info data.frame
#' @param diff_results Output of run_differential()
#' @param cfg Config list
#' @param out_dir Output directory
#' @return NULL (results saved to files)
run_biomarker_analysis <- function(mSet, object2, sample_info, diff_results,
                                    cfg, out_dir) {
  cat("\n----- Biomarker Analysis (Random Forest + ROC) -----\n")

  if (!requireNamespace("MetaboAnalystR", quietly = TRUE)) {
    cat("   MetaboAnalystR not available, cannot run biomarker analysis\n")
    return(NULL)
  }
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    cat("   randomForest not available, cannot run biomarker analysis\n")
    return(NULL)
  }
  if (!requireNamespace("pROC", quietly = TRUE)) {
    cat("   pROC package not available, ROC curve analysis will be skipped\n")
    cat("   Install with: install.packages('pROC')\n")
    # Allow RF to proceed, ROC will be skipped
  }

  bio_cfg <- cfg$metaboanalyst_advanced$biomarker
  rf_ntree <- bio_cfg$rf_ntree %||% 500
  rf_mtry <- bio_cfg$rf_mtry  # NULL = default
  rf_importance <- bio_cfg$rf_importance %||% "accuracy"
  n_top <- bio_cfg$n_top_features %||% 20
  roc_threshold <- bio_cfg$roc_threshold %||% 0.7
  roc_ci_method <- bio_cfg$roc_ci_method %||% "delong"
  roc_n_bootstrap <- bio_cfg$roc_n_bootstrap %||% 2000

  # ── Extract expression data ────────────────────────────────────────────
  expr_data <- t(object2@expression_data)  # samples x features
  feature_names <- colnames(expr_data)

  # ── Process each comparison ────────────────────────────────────────────
  for (cmp_name in names(diff_results)) {
    dr <- diff_results[[cmp_name]]
    if (is.null(dr) || is.null(dr$sig) || nrow(dr$sig) == 0) {
      cat(sprintf("   Skipping %s: no significant features\n", cmp_name))
      next
    }

    label <- dr$comparison$label %||% cmp_name
    cat(sprintf("   Processing: %s\n", label))

    # Get samples for this comparison
    ctrl_samples <- dr$ctrl_samples
    treat_samples <- dr$treat_samples
    all_samples <- c(ctrl_samples, treat_samples)

    if (length(ctrl_samples) < 3 || length(treat_samples) < 3) {
      cat(sprintf("     Insufficient samples (ctrl=%d, treat=%d), need >=3 each\n",
                  length(ctrl_samples), length(treat_samples)))
      next
    }

    # Subset expression data to these samples
    sample_idx <- intersect(all_samples, rownames(expr_data))
    if (length(sample_idx) < 4) {
      cat("     Insufficient overlapping samples, skipping\n")
      next
    }

    # Focus on significant features only (dimension reduction)
    sig_features <- intersect(rownames(dr$sig), feature_names)
    if (length(sig_features) < 3) {
      cat("     Too few significant features, skipping\n")
      next
    }

    rf_data <- expr_data[sample_idx, sig_features, drop = FALSE]
    response <- factor(ifelse(sample_idx %in% treat_samples, "Treatment", "Control"),
                       levels = c("Control", "Treatment"))

    # Convert to data.frame for randomForest
    rf_df <- as.data.frame(rf_data)
    rf_df$class <- response

    cmp_dir <- file.path(out_dir, cmp_name)
    dir.create(cmp_dir, showWarnings = FALSE, recursive = TRUE)

    # ── Random Forest feature selection ──────────────────────────────────
    cat(sprintf("     Running Random Forest (%d trees, %d samples, %d features)...\n",
                rf_ntree, nrow(rf_df), length(sig_features)))

    rf_result <- tryCatch({
      mtry_val <- if (is.null(rf_mtry)) max(1, floor(sqrt(length(sig_features)))) else rf_mtry
      rf <- randomForest::randomForest(
        class ~ .,
        data = rf_df,
        ntree = rf_ntree,
        mtry = mtry_val,
        importance = TRUE,
        proximity = FALSE,
        na.action = na.omit
      )
      rf
    }, error = function(e) {
      cat(sprintf("     ERROR: Random Forest failed: %s\n", conditionMessage(e)))
      return(NULL)
    })

    if (is.null(rf_result)) {
      cat(sprintf("     Random Forest failed for %s, skipping\n", label))
      next
    }

    cat(sprintf("     RF OOB error rate: %.2f%%\n",
                rf_result$err.rate[nrow(rf_result$err.rate), 1] * 100))

    # ── Extract variable importance ──────────────────────────────────────
    importance_mat <- randomForest::importance(rf_result)
    # Determine which importance column to use
    if (rf_importance == "accuracy" && "MeanDecreaseAccuracy" %in% colnames(importance_mat)) {
      imp_col <- "MeanDecreaseAccuracy"
    } else if ("MeanDecreaseGini" %in% colnames(importance_mat)) {
      imp_col <- "MeanDecreaseGini"
    } else {
      imp_col <- colnames(importance_mat)[1]
    }

    importance_df <- data.frame(
      variable_id = rownames(importance_mat),
      importance = importance_mat[, imp_col],
      stringsAsFactors = FALSE
    )
    importance_df <- importance_df[order(importance_df$importance, decreasing = TRUE), ]

    write.csv(importance_df,
              file.path(cmp_dir, "rf_importance.csv"),
              row.names = FALSE)
    cat(sprintf("     Top features by %s:\n", imp_col))
    for (i in 1:min(5, nrow(importance_df))) {
      cat(sprintf("       %d. %s (importance: %.3f)\n",
                  i, importance_df$variable_id[i], importance_df$importance[i]))
    }

    # ── Plot variable importance ─────────────────────────────────────────
    top_n_imp <- importance_df %>% head(min(n_top * 2, nrow(importance_df)))
    top_n_imp$variable_id <- factor(top_n_imp$variable_id,
                                    levels = rev(top_n_imp$variable_id))

    p <- ggplot(top_n_imp, aes(x = importance, y = variable_id)) +
      geom_bar(stat = "identity", fill = "#088247", alpha = 0.8) +
      labs(x = sprintf("Importance (%s)", imp_col), y = "",
           title = sprintf("RF Variable Importance — %s", label)) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 12),
            axis.text.y = element_text(size = rel(0.8)))

    save_plot(file.path(cmp_dir, "rf_importance"),
              plot = p, width = 10, height = max(4, nrow(top_n_imp) * 0.3))

    # ── ROC Curve Analysis for top N biomarker candidates ────────────────
    if (!requireNamespace("pROC", quietly = TRUE)) {
      cat("     pROC not available, skipping ROC curve analysis\n")
      # Save the RF model for later use
      saveRDS(rf_result, file.path(cmp_dir, "random_forest_model.rds"))
      next
    }

    cat(sprintf("     Computing ROC curves for top %d features...\n",
                min(n_top, nrow(importance_df))))

    top_features <- importance_df$variable_id[1:min(n_top, nrow(importance_df))]
    roc_results <- list()

    for (feat in top_features) {
      if (!feat %in% colnames(rf_data)) next
      values <- rf_data[, feat]
      if (length(unique(values)) < 2) next

      roc_obj <- tryCatch({
        pROC::roc(response, values, quiet = TRUE, ci = TRUE,
                  ci.method = roc_ci_method)
      }, error = function(e) {
        # Fallback without CI
        tryCatch({
          pROC::roc(response, values, quiet = TRUE)
        }, error = function(e2) NULL)
      })

      if (!is.null(roc_obj)) {
        auc_val <- as.numeric(pROC::auc(roc_obj))
        ci_val <- tryCatch(as.numeric(pROC::ci.auc(roc_obj)), error = function(e) NULL)

        roc_results[[feat]] <- list(
          roc = roc_obj,
          auc = auc_val,
          ci = ci_val
        )
      }
    }

    # Filter to features with AUC >= threshold
    roc_results <- roc_results[sapply(roc_results, function(r) r$auc >= roc_threshold)]

    if (length(roc_results) == 0) {
      cat(sprintf("     No features met AUC threshold >= %.2f\n", roc_threshold))
      # Save RF model anyway
      saveRDS(rf_result, file.path(cmp_dir, "random_forest_model.rds"))
      next
    }

    cat(sprintf("     %d features with AUC >= %.2f\n", length(roc_results), roc_threshold))

    # ── Save ROC summary table ──────────────────────────────────────────
    roc_summary <- data.frame(
      variable_id = names(roc_results),
      auc = sapply(roc_results, function(r) r$auc),
      ci_lower = sapply(roc_results, function(r) if (!is.null(r$ci)) r$ci[1] else NA),
      ci_upper = sapply(roc_results, function(r) if (!is.null(r$ci)) r$ci[3] else NA),
      stringsAsFactors = FALSE
    )
    roc_summary <- roc_summary[order(roc_summary$auc, decreasing = TRUE), ]

    write.csv(roc_summary,
              file.path(cmp_dir, "roc_summary.csv"),
              row.names = FALSE)

    # ── Plot individual ROC curves (top 9) ──────────────────────────────
    top_n_roc <- head(names(roc_results), min(9, length(roc_results)))

    # Individual ROC plots
    for (feat in top_n_roc) {
      rr <- roc_results[[feat]]
      auc_label <- sprintf("AUC = %.3f", rr$auc)
      if (!is.null(rr$ci) && length(rr$ci) == 3) {
        auc_label <- sprintf("AUC = %.3f (%.3f - %.3f)",
                             rr$auc, rr$ci[1], rr$ci[3])
      }

      # Create a data.frame for ggplot
      roc_df <- data.frame(
        specificity = rev(rr$roc$specificities),
        sensitivity = rev(rr$roc$sensitivities)
      )

      p <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
        geom_line(color = "#D20A13", linewidth = 1.2) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                     color = "grey50", linewidth = 0.8) +
        annotate("text", x = 0.75, y = 0.25, label = auc_label,
                 size = 4.5, hjust = 0, fontface = "italic") +
        labs(x = "1 - Specificity", y = "Sensitivity",
             title = sprintf("ROC: %s", feat)) +
        coord_fixed() +
        theme_bw() +
        theme(panel.grid = element_blank(),
              plot.title = element_text(hjust = 0.5, size = 10))

      save_plot(file.path(cmp_dir, paste0("roc_", gsub("[^A-Za-z0-9]", "_", feat))),
                plot = p, width = 6, height = 6)
    }

    # ── Multi-ROC overlay plot (all features with AUC >= threshold) ─────
    if (length(roc_results) >= 2) {
      # Build combined ROC plot data
      roc_plot_data <- do.call(rbind, lapply(names(roc_results), function(feat) {
        rr <- roc_results[[feat]]
        data.frame(
          variable_id = feat,
          specificity = rev(rr$roc$specificities),
          sensitivity = rev(rr$roc$sensitivities),
          stringsAsFactors = FALSE
        )
      }))

      # Add AUC labels
      auc_labels <- data.frame(
        variable_id = names(roc_results),
        auc = sapply(roc_results, function(r) r$auc),
        stringsAsFactors = FALSE
      )
      auc_labels$label <- sprintf("%s (AUC=%.3f)", auc_labels$variable_id, auc_labels$auc)
      roc_plot_data <- roc_plot_data %>%
        left_join(auc_labels, by = "variable_id")

      # Use a qualitative color palette
      n_curves <- length(unique(roc_plot_data$variable_id))
      palette <- c("#D20A13", "#088247", "#0072B5FF", "#FFD121",
                   "#7E6148FF", "#BC3C29FF", "#20854EFF", "#7876B1FF",
                   "#E18727FF", "#F39B7FFF", "#5BC0EB", "#EE4C97FF")

      p <- ggplot(roc_plot_data, aes(x = 1 - specificity, y = sensitivity,
                                      color = label)) +
        geom_line(linewidth = 0.8) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                     color = "grey50", linewidth = 0.6) +
        scale_color_manual(values = rep(palette, length.out = n_curves)) +
        labs(x = "1 - Specificity", y = "Sensitivity",
             color = "Feature (AUC)",
             title = sprintf("ROC Curves — %s", label)) +
        coord_fixed() +
        theme_bw() +
        theme(panel.grid = element_blank(),
              plot.title = element_text(hjust = 0.5, size = 12),
              legend.text = element_text(size = 8))

      save_plot(file.path(cmp_dir, "roc_multi_overlay"),
                plot = p, width = 9, height = 7)
    }

    # ── Save RF model ───────────────────────────────────────────────────
    saveRDS(rf_result, file.path(cmp_dir, "random_forest_model.rds"))
    cat(sprintf("     Random Forest model saved\n"))
  }

  cat("   Biomarker analysis complete\n")
  invisible(NULL)
}


#' Run Chemical Class Enrichment Analysis
#'
#' Performs enrichment analysis based on chemical taxonomy (HMDB Super Class,
#' Class, or Sub Class) rather than biochemical pathways. Uses the annotation
#' table's HMDB chemical taxonomy annotations (or ClassyFire classifications)
#' to test for over-representation of chemical classes among significant
#' features.
#'
#' This is a self-contained implementation that does not require MetaboAnalystR
#' (it uses Fisher's exact test or globaltest-style approach). If MetaboAnalystR
#' is available, it can also use its MSEA infrastructure with custom metabolite
#' sets derived from chemical taxonomy.
#'
#' @param app3 Annotation table with HMDB taxonomy columns
#' @param diff_results Output of run_differential()
#' @param cfg Config list
#' @param out_dir Output directory
#' @return Invisible data.frame with enrichment results
run_chem_class_enrichment <- function(app3, diff_results, cfg, out_dir) {
  cat("\n----- Chemical Class Enrichment Analysis -----\n")

  cc_cfg <- cfg$metaboanalyst_advanced$chem_class_enrichment
  taxon_level <- cc_cfg$taxon_level %||% "super_class"
  method <- cc_cfg$method %||% "fisher"
  p_cutoff <- cc_cfg$p_cutoff %||% 0.05
  p_adjust <- cc_cfg$p_adjust_method %||% "fdr"

  # ── Check for taxonomy annotations ─────────────────────────────────────
  # Look for HMDB taxonomy columns in app3
  # HMDB taxonomy columns typically follow the pattern:
  #   HMDB.super.class, HMDB.class, HMDB.sub.class
  # or ClassyFire columns:
  #   ClassyFire.superclass, ClassyFire.class, ClassyFire.subclass
  # We also try the tidymass convention: super_class, main_class, sub_class

  if (is.null(app3) || nrow(app3) == 0) {
    cat("   No annotations available, skipping\n")
    return(invisible(NULL))
  }

  # Detect taxonomy columns
  taxon_cols <- c(
    "super_class", "main_class", "sub_class",
    "HMDB.super.class", "HMDB.class", "HMDB.sub.class",
    "ClassyFire.superclass", "ClassyFire.class", "ClassyFire.subclass",
    "Super.class", "Class", "Sub.class"
  )
  available_taxon <- intersect(taxon_cols, colnames(app3))

  if (length(available_taxon) == 0) {
    cat("   No chemical taxonomy columns found in annotation table.\n")
    cat("   To enable chemical class enrichment, ensure HMDB taxonomy\n")
    cat("   columns are present in the annotation database.\n")
    cat("   Available columns: ", paste(head(colnames(app3), 30), collapse = ", "),
        if (ncol(app3) > 30) "...", "\n")
    return(invisible(NULL))
  }

  cat(sprintf("   Using taxonomy level: '%s'\n", taxon_level))

  # Find the column that best matches the requested level
  taxon_col <- NULL
  taxon_patterns <- switch(taxon_level,
    super_class = c("super_class", "Super.class", "HMDB.super.class",
                    "ClassyFire.superclass", "superclass"),
    class = c("main_class", "Class", "HMDB.class", "ClassyFire.class",
              "chem_class"),
    sub_class = c("sub_class", "Sub.class", "HMDB.sub.class",
                  "ClassyFire.subclass", "subclass")
  )
  for (pat in taxon_patterns) {
    matches <- grep(pat, colnames(app3), ignore.case = TRUE, value = TRUE)
    if (length(matches) > 0) {
      taxon_col <- matches[1]
      break
    }
  }

  if (is.null(taxon_col)) {
    cat(sprintf("   No column matching '%s' found, trying first available taxonomy column\n",
                taxon_level))
    cat(sprintf("   Available taxonomy columns: %s\n",
                paste(available_taxon, collapse = ", ")))
    taxon_col <- available_taxon[1]
  }

  cat(sprintf("   Using taxonomy column: '%s'\n", taxon_col))

  # ── Build the taxonomy universe ────────────────────────────────────────
  # Get the taxonomy for all annotated features
  taxonomy_data <- app3 %>%
    dplyr::select(variable_id, taxonomy = dplyr::all_of(taxon_col)) %>%
    filter(!is.na(taxonomy) & taxonomy != "" & taxonomy != "NA")

  if (nrow(taxonomy_data) == 0) {
    cat("   No features have taxonomy annotations, skipping\n")
    return(invisible(NULL))
  }

  cat(sprintf("   %d features with taxonomy annotations\n", nrow(taxonomy_data)))

  # ── Process each comparison ────────────────────────────────────────────
  all_results <- list()

  for (cmp_name in names(diff_results)) {
    dr <- diff_results[[cmp_name]]
    if (is.null(dr) || is.null(dr$sig) || nrow(dr$sig) == 0) next

    label <- dr$comparison$label %||% cmp_name
    cat(sprintf("   Processing: %s\n", label))

    # Get significant features
    sig_vars <- rownames(dr$sig)

    # Merge with taxonomy
    sig_tax <- taxonomy_data %>%
      filter(variable_id %in% sig_vars)

    if (nrow(sig_tax) < 3) {
      cat("     Too few significant features with taxonomy, skipping\n")
      next
    }

    # ── Compute enrichment ───────────────────────────────────────────────
    # Build the contingency table for each taxon class
    all_classes <- unique(taxonomy_data$taxonomy)
    sig_classes <- unique(sig_tax$taxonomy)

    if (method == "fisher") {
      # Fisher's exact test for each class
      n_total <- nrow(taxonomy_data)
      n_sig <- nrow(sig_tax)

      enrichment_list <- lapply(all_classes, function(cls) {
        in_class_total <- sum(taxonomy_data$taxonomy == cls)
        in_class_sig <- sum(sig_tax$taxonomy == cls)

        # Create 2x2 contingency table
        #                  | In class | Not in class
        #   Significant    | a        | b
        #   Not significant| c        | d
        a <- in_class_sig
        b <- n_sig - in_class_sig
        c <- in_class_total - in_class_sig
        d <- n_total - in_class_total - b

        # Small expected count check
        if (a + b + c + d < 3 || a + c < 2) return(NULL)

        ft <- tryCatch(
          fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater"),
          error = function(e) NULL
        )

        if (is.null(ft)) return(NULL)

        data.frame(
          taxonomy = cls,
          n_total = in_class_total,
          n_sig = in_class_sig,
          expected = round(in_class_total * n_sig / n_total, 2),
          fold_enrichment = if (in_class_total > 0 && in_class_sig > 0)
                              round((in_class_sig / n_sig) / (in_class_total / n_total), 2)
                            else NA,
          p_value = ft$p.value,
          stringsAsFactors = FALSE
        )
      })

      enrichment_df <- do.call(rbind, enrichment_list)
      if (is.null(enrichment_df) || nrow(enrichment_df) == 0) {
        cat("     No enrichment results, skipping\n")
        next
      }

      # Adjust p-values
      enrichment_df$p_adjusted <- p.adjust(enrichment_df$p_value, method = p_adjust)
      enrichment_df <- enrichment_df[order(enrichment_df$p_value), , drop = FALSE]

    } else if (method == "globaltest") {
      # Use a globaltest-inspired approach: compute the sum of t-statistics
      # for features in each class vs. the overall distribution
      # This is a simplified version; full GlobalTest requires the raw data
      cat("     GlobalTest method requires raw expression data,\n")
      cat("     falling back to Fisher's exact test\n")

      # Fall back to Fisher
      n_total <- nrow(taxonomy_data)
      n_sig <- nrow(sig_tax)

      enrichment_list <- lapply(all_classes, function(cls) {
        in_class_total <- sum(taxonomy_data$taxonomy == cls)
        in_class_sig <- sum(sig_tax$taxonomy == cls)
        a <- in_class_sig; b <- n_sig - in_class_sig
        c <- in_class_total - in_class_sig; d <- n_total - in_class_total - b
        if (a + b + c + d < 3 || a + c < 2) return(NULL)
        ft <- tryCatch(fisher.test(matrix(c(a, b, c, d), nrow = 2),
                                    alternative = "greater"), error = function(e) NULL)
        if (is.null(ft)) return(NULL)
        data.frame(taxonomy = cls, n_total = in_class_total, n_sig = in_class_sig,
                   expected = round(in_class_total * n_sig / n_total, 2),
                   fold_enrichment = if (in_class_total > 0 && in_class_sig > 0)
                     round((in_class_sig / n_sig) / (in_class_total / n_total), 2) else NA,
                   p_value = ft$p.value, stringsAsFactors = FALSE)
      })
      enrichment_df <- do.call(rbind, enrichment_list)
      if (is.null(enrichment_df) || nrow(enrichment_df) == 0) next
      enrichment_df$p_adjusted <- p.adjust(enrichment_df$p_value, method = p_adjust)
      enrichment_df <- enrichment_df[order(enrichment_df$p_value), , drop = FALSE]
    }

    # ── Save results ────────────────────────────────────────────────────
    cmp_dir <- file.path(out_dir, cmp_name)
    dir.create(cmp_dir, showWarnings = FALSE, recursive = TRUE)

    write.csv(enrichment_df,
              file.path(cmp_dir, "chem_class_enrichment.csv"),
              row.names = FALSE)
    cat(sprintf("     %d chemical classes tested, %d significant (p<%.3g)\n",
                nrow(enrichment_df),
                sum(enrichment_df$p_adjusted < p_cutoff, na.rm = TRUE),
                p_cutoff))

    # ── Plot ─────────────────────────────────────────────────────────────
    sig_enrich <- enrichment_df %>%
      filter(p_adjusted < p_cutoff) %>%
      head(20)

    if (nrow(sig_enrich) > 0) {
      sig_enrich <- sig_enrich %>%
        arrange(desc(fold_enrichment)) %>%
        mutate(taxonomy = factor(taxonomy, levels = taxonomy))

      p <- ggplot(sig_enrich, aes(x = fold_enrichment, y = taxonomy)) +
        geom_bar(stat = "identity", aes(fill = -log10(p_adjusted + 1e-10)),
                 alpha = 0.85) +
        geom_text(aes(label = n_sig), hjust = -0.2, size = 3) +
        scale_fill_gradientn(colours = c("#46bac2", "#f7ca64", "#e34a33"),
                             name = "-log10(p.adj)") +
        labs(x = "Fold Enrichment", y = "",
             title = sprintf("Chemical Class Enrichment — %s (%s)",
                             label, taxon_level)) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              plot.title = element_text(hjust = 0.5, size = 12),
              axis.text.y = element_text(size = rel(0.9)))

      save_plot(file.path(cmp_dir, "chem_class_enrichment"),
                plot = p, width = 10, height = max(4, nrow(sig_enrich) * 0.35))
    }

    # Show top classes
    top5 <- head(enrichment_df, 5)
    for (i in seq_len(nrow(top5))) {
      cat(sprintf("       %s: n_sig=%d, fold=%.2f, p=%.4g\n",
                  top5$taxonomy[i], top5$n_sig[i],
                  top5$fold_enrichment[i], top5$p_adjusted[i]))
    }

    all_results[[cmp_name]] <- enrichment_df
  }

  cat("   Chemical class enrichment complete\n")
  invisible(all_results)
}


#' Run MetaboAnalystR Advanced Analyses (Step 10)
#'
#' Master orchestrator for all four MetaboAnalystR downstream modules:
#'   1. MSEA (QSEA/GlobalTest)
#'   2. Pathway Topology Analysis
#'   3. Biomarker Analysis (RF + ROC)
#'   4. Chemical Class Enrichment
#'
#' Each module is wrapped in tryCatch for graceful degradation — if a module
#' fails, the pipeline continues with the next one.
#'
#' @param object2 mass_dataset object (pre-processed)
#' @param sample_info Sample information data.frame
#' @param diff_results Output of run_differential()
#' @param app3 Annotation table from run_annotation()
#' @param cfg Validated config list
#' @return Invisible list with module results
run_metaboanalyst_advanced <- function(object2, sample_info, diff_results,
                                       app3, cfg) {
  cat("\n===== Step 10: MetaboAnalystR Advanced Analysis =====\n")

  ma_cfg <- cfg$metaboanalyst_advanced
  if (!isTRUE(ma_cfg$enabled)) {
    cat("   MetaboAnalystR advanced analysis disabled in config, skipping\n")
    return(invisible(NULL))
  }

  # ── Create output directory ────────────────────────────────────────────
  out_dir <- file.path(cfg$project$output_dir, "10_MetaboAnalystR_Advanced")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  cat(sprintf("   Output directory: %s\n", out_dir))

  # ── Step 10a: Convert mass_dataset to mSetObj ──────────────────────────
  conversion <- convert_mass_dataset_to_mSet(object2, sample_info,
                                              diff_results, cfg)
  if (!conversion$success || is.null(conversion$mSet)) {
    cat("   WARNING: Could not initialize MetaboAnalystR mSet object.\n")
    cat("   All MetaboAnalystR modules will be skipped.\n")
    cat("   Check that MetaboAnalystR is installed and the data is valid.\n")
    return(invisible(NULL))
  }
  mSet <- conversion$mSet
  id_map <- conversion$id_map

  # ── Step 10b: MSEA (Metabolite Set Enrichment Analysis) ────────────────
  msea_dir <- file.path(out_dir, "MSEA")
  dir.create(msea_dir, showWarnings = FALSE, recursive = TRUE)
  msea_result <- NULL
  if (isTRUE(ma_cfg$msea$enabled)) {
    tryCatch({
      msea_result <- run_msea_analysis(mSet, id_map, app3, diff_results,
                                        cfg, msea_dir)
    }, error = function(e) {
      cat(sprintf("   WARNING: MSEA failed: %s\n", conditionMessage(e)))
      cat("   Continuing with next module...\n")
    })
  } else {
    cat("   MSEA disabled in config, skipping\n")
  }

  # ── Step 10c: Pathway Topology Analysis ────────────────────────────────
  topo_dir <- file.path(out_dir, "Pathway_Topology")
  dir.create(topo_dir, showWarnings = FALSE, recursive = TRUE)
  topo_result <- NULL
  if (isTRUE(ma_cfg$pathway_topology$enabled)) {
    tryCatch({
      topo_result <- run_pathway_topology(mSet, app3, diff_results, cfg, topo_dir)
    }, error = function(e) {
      cat(sprintf("   WARNING: Pathway Topology failed: %s\n", conditionMessage(e)))
      cat("   Continuing with next module...\n")
    })
  } else {
    cat("   Pathway Topology disabled in config, skipping\n")
  }

  # ── Step 10d: Biomarker Analysis (RF + ROC) ────────────────────────────
  bio_dir <- file.path(out_dir, "Biomarker_Analysis")
  dir.create(bio_dir, showWarnings = FALSE, recursive = TRUE)
  bio_result <- NULL
  if (isTRUE(ma_cfg$biomarker$enabled)) {
    tryCatch({
      bio_result <- run_biomarker_analysis(mSet, object2, sample_info, diff_results,
                                           cfg, bio_dir)
    }, error = function(e) {
      cat(sprintf("   WARNING: Biomarker Analysis failed: %s\n", conditionMessage(e)))
      cat("   Continuing with next module...\n")
    })
  } else {
    cat("   Biomarker Analysis disabled in config, skipping\n")
  }

  # ── Step 10e: Chemical Class Enrichment ────────────────────────────────
  chem_dir <- file.path(out_dir, "Chemical_Class_Enrichment")
  dir.create(chem_dir, showWarnings = FALSE, recursive = TRUE)
  chem_result <- NULL
  if (isTRUE(ma_cfg$chem_class_enrichment$enabled)) {
    tryCatch({
      chem_result <- run_chem_class_enrichment(app3, diff_results, cfg, chem_dir)
    }, error = function(e) {
      cat(sprintf("   WARNING: Chemical Class Enrichment failed: %s\n", conditionMessage(e)))
      cat("   Continuing...\n")
    })
  } else {
    cat("   Chemical Class Enrichment disabled in config, skipping\n")
  }

  cat("   Step 10: MetaboAnalystR advanced analysis complete\n")

  # ── Store results globally for web export ──────────────────────────────
  msea_results_global <<- msea_result
  topology_results_global <<- topo_result
  biomarker_results_global <<- bio_result
  chem_class_results_global <<- chem_result

  invisible(list(
    mSet = mSet,
    msea_result = msea_result,
    chem_result = chem_result
  ))
}