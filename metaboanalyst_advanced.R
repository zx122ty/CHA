# ══════════════════════════════════════════════════════════════════════════════
# 9.5. Advanced downstream analysis (Step 10)
# ══════════════════════════════════════════════════════════════════════════════
#
# Three per-comparison analyses on top of Steps 4-7:
#   1. Pathway ORA on identified compounds — significant vs all MSI Level 1-2
#      compounds with a KEGG ID, over the organism's KEGG pathways (the
#      mummichog library of Step 7)
#   2. Biomarker analysis — Random Forest importance + per-feature ROC curves
#   3. Chemical class enrichment — one class per feature (the identified
#      compound's, CANOPUS, or the class all accurate-mass candidates share)
#
# MetaboAnalystR is not used here any more. Its MSEA/topology read compound
# IDs from the annotation, which for Level 3 features are accurate-mass
# candidates (about half of them chance matches); its 4.x API takes vectors
# where this module passed files; and its reference libraries now ship in qs2
# format, which 4.2 cannot read. Pathway topology needs KEGG pathway graphs
# and is not computed. Step 7 (mummichog on all features) is the main pathway
# analysis; this ORA covers only the few identified compounds.
#
# Every analysis is wrapped in tryCatch: a failure is reported and skipped.
# Results are returned per comparison in the shape export_for_web.R reads.
# ══════════════════════════════════════════════════════════════════════════════

#' One row per feature with its MSI level, display name and identified IDs
#'
#' @param app3 Candidate table from finalize_annotation()
#' @return data.frame (variable_id, level, display_name, identified_kegg_id)
advanced_feature_table <- function(app3) {
  if (is.null(app3) || nrow(app3) == 0) {
    return(data.frame(variable_id = character(0), level = character(0),
                      display_name = character(0), identified_kegg_id = character(0)))
  }
  ord <- if ("candidate_rank" %in% colnames(app3)) order(app3$variable_id, app3$candidate_rank) else seq_len(nrow(app3))
  f <- app3[ord, , drop = FALSE]
  f <- f[!duplicated(f$variable_id), , drop = FALSE]
  get <- function(col) if (col %in% colnames(f)) as.character(f[[col]]) else rep(NA_character_, nrow(f))
  level <- get("feature_level")
  level[is.na(level)] <- get("confidence_level")[is.na(level)]
  data.frame(variable_id = as.character(f$variable_id), level = level,
             display_name = get("display_name"),
             identified_kegg_id = get("identified_kegg_id"),
             stringsAsFactors = FALSE)
}

#' Significant features of one comparison (same p column and cutoff as Step 4)
significant_ids <- function(dr, cfg) {
  if (is.null(dr$all)) return(character(0))
  p_col <- if (identical(cfg$differential$significance_metric, "adj_p_value")) "adj.P.Val" else "P.Value"
  p <- dr$all[[p_col]]
  rownames(dr$all)[!is.na(p) & p < cfg$differential$p_value_cutoff]
}

#' Pathway over-representation of significant identified compounds
#'
#' Background = every identified (MSI Level 1-2) compound of the batch with a
#' KEGG ID that is in at least one pathway of the organism's library; the
#' query = those of them significant in the comparison. Fisher's exact test
#' per pathway with at least one query hit, BH over the tested pathways.
#' A comparison needs >= msea$min_significant (default 3) significant
#' identified compounds in the library, otherwise it is skipped.
#'
#' @param app3 Candidate table from finalize_annotation()
#' @param diff_results Output of run_differential()
#' @param cfg Validated config list
#' @param out_dir Output directory (MSEA/)
#' @return Named list per comparison: list(result, method), or NULL
run_identified_pathway_ora <- function(app3, diff_results, cfg, out_dir) {
  cat("\n----- Pathway ORA on identified compounds (KEGG) -----\n")
  lib_file <- mummichog_library_file(cfg)
  if (is.null(lib_file)) {
    cat("   KEGG pathway library unavailable, skipping\n")
    return(NULL)
  }
  lib <- qs::qread(lib_file)
  pw_cpds <- lib$pathways$cpds
  names(pw_cpds) <- lib$pathways$name
  in_lib <- unique(unlist(pw_cpds))

  feat <- advanced_feature_table(app3)
  feat$identified_kegg_id <- sub("[;,| ].*$", "", feat$identified_kegg_id)   # first ID if several
  ident <- feat[grepl("Level [12]", feat$level) & !is.na(feat$identified_kegg_id) &
                  feat$identified_kegg_id %in% in_lib, , drop = FALSE]
  universe <- unique(ident$identified_kegg_id)
  min_sig <- as.integer(cfg$metaboanalyst_advanced$msea$min_significant %||% 3)
  cat(sprintf("   Identified compounds in the KEGG %s library: %d (background)\n",
              cfg$kegg$organism %||% "mmu", length(universe)))
  if (length(universe) < min_sig) {
    cat(sprintf("   Fewer than %d identified compounds; pathway enrichment is Step 7 (mummichog)\n", min_sig))
    return(NULL)
  }

  results <- list()
  n_skipped <- 0
  for (cmp_name in names(diff_results)) {
    dr <- diff_results[[cmp_name]]
    if (is.null(dr)) next
    hits <- unique(ident$identified_kegg_id[ident$variable_id %in% significant_ids(dr, cfg)])
    if (length(hits) < min_sig) {
      n_skipped <- n_skipped + 1
      next
    }
    rows <- lapply(names(pw_cpds), function(pw) {
      members <- intersect(pw_cpds[[pw]], universe)
      k <- length(intersect(members, hits))
      if (k == 0) return(NULL)
      tab <- matrix(c(k, length(members) - k,
                      length(hits) - k, length(universe) - length(members) - length(hits) + k), 2)
      data.frame(set_name = pw, hit_count = k, total_count = length(members),
                 p_value = stats::fisher.test(tab, alternative = "greater")$p.value,
                 hits = paste(intersect(members, hits), collapse = ";"),
                 stringsAsFactors = FALSE)
    })
    res <- do.call(rbind, rows)
    if (is.null(res) || nrow(res) == 0) {
      n_skipped <- n_skipped + 1
      next
    }
    res$fdr <- stats::p.adjust(res$p_value, "BH")
    res <- res[order(res$p_value), , drop = FALSE]
    cmp_dir <- file.path(out_dir, cmp_name)
    dir.create(cmp_dir, showWarnings = FALSE, recursive = TRUE)
    write.csv(res, file.path(cmp_dir, "identified_pathway_ora.csv"), row.names = FALSE)
    results[[cmp_name]] <- list(result = res,
                                method = "ORA, identified compounds (Fisher, BH)")
  }
  cat(sprintf("   %d comparisons tested, %d skipped (< %d significant identified compounds)\n",
              length(results), n_skipped, min_sig))
  if (length(results) == 0) NULL else results
}

#' Biomarker analysis: Random Forest importance + ROC curves
#'
#' Per comparison: a Random Forest (control vs treatment) on the significant
#' features, then a ROC curve for each of the top n_top_features by
#' importance. Features are labelled with their display names.
#'
#' @param object2 mass_dataset (pre-processed)
#' @param app3 Candidate table (for display names)
#' @param diff_results Output of run_differential()
#' @param cfg Validated config list
#' @param out_dir Output directory (Biomarker_Analysis/)
#' @return Named list per comparison: list(feature_importance, roc_summary, roc_auc)
run_biomarker_analysis <- function(object2, app3, diff_results, cfg, out_dir) {
  cat("\n----- Biomarker Analysis (Random Forest + ROC) -----\n")
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    cat("   randomForest not available, skipping\n")
    return(NULL)
  }
  has_proc <- requireNamespace("pROC", quietly = TRUE)
  if (!has_proc) cat("   pROC not available, ROC curves skipped\n")

  bio_cfg <- cfg$metaboanalyst_advanced$biomarker
  rf_ntree <- as.integer(bio_cfg$rf_ntree %||% 500)
  rf_mtry <- bio_cfg$rf_mtry
  rf_importance <- bio_cfg$rf_importance %||% "accuracy"
  n_top <- as.integer(bio_cfg$n_top_features %||% 20)
  roc_threshold <- as.numeric(bio_cfg$roc_threshold %||% 0.7)
  roc_ci_method <- bio_cfg$roc_ci_method %||% "delong"

  names_map <- advanced_feature_table(app3)
  label_of <- function(ids) {
    lab <- names_map$display_name[match(ids, names_map$variable_id)]
    ifelse(is.na(lab), ids, lab)
  }
  expr_data <- log2(t(as.matrix(object2@expression_data)))   # samples x features

  results <- list()
  for (cmp_name in names(diff_results)) {
    dr <- diff_results[[cmp_name]]
    if (is.null(dr)) next
    label <- dr$comparison$label %||% cmp_name
    sig_features <- intersect(significant_ids(dr, cfg), colnames(expr_data))
    samples <- intersect(c(dr$ctrl_samples, dr$treat_samples), rownames(expr_data))
    if (length(sig_features) < 3 || length(intersect(dr$ctrl_samples, samples)) < 3 ||
        length(intersect(dr$treat_samples, samples)) < 3) next

    rf_data <- expr_data[samples, sig_features, drop = FALSE]
    response <- factor(ifelse(samples %in% dr$treat_samples, "Treatment", "Control"),
                       levels = c("Control", "Treatment"))
    set.seed(1)
    rf <- tryCatch(randomForest::randomForest(
      x = rf_data, y = response, ntree = rf_ntree, importance = TRUE,
      mtry = if (is.null(rf_mtry)) max(1, floor(sqrt(ncol(rf_data)))) else rf_mtry),
      error = function(e) {
        cat(sprintf("   %s: Random Forest failed: %s\n", label, conditionMessage(e)))
        NULL
      })
    if (is.null(rf)) next

    imp <- randomForest::importance(rf)
    imp_col <- if (rf_importance == "accuracy" && "MeanDecreaseAccuracy" %in% colnames(imp)) {
      "MeanDecreaseAccuracy"
    } else "MeanDecreaseGini"
    fi <- data.frame(variable_id = rownames(imp), display_name = label_of(rownames(imp)),
                     importance = imp[, imp_col], stringsAsFactors = FALSE)
    fi <- fi[order(fi$importance, decreasing = TRUE), , drop = FALSE]
    rownames(fi) <- fi$variable_id

    cmp_dir <- file.path(out_dir, cmp_name)
    dir.create(cmp_dir, showWarnings = FALSE, recursive = TRUE)
    write.csv(fi, file.path(cmp_dir, "rf_importance.csv"), row.names = FALSE)

    top_imp <- head(fi, n_top)
    top_imp$display_name <- factor(top_imp$display_name, levels = rev(unique(top_imp$display_name)))
    p <- ggplot(top_imp, aes(x = importance, y = display_name)) +
      geom_bar(stat = "identity", fill = "#088247", alpha = 0.8) +
      labs(x = sprintf("Importance (%s)", imp_col), y = "",
           title = sprintf("RF variable importance — %s (OOB error %.0f%%)", label,
                           100 * rf$err.rate[nrow(rf$err.rate), "OOB"])) +
      theme_bw() +
      theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, size = 11),
            axis.text.y = element_text(size = rel(0.8)))
    save_plot(file.path(cmp_dir, "rf_importance"), plot = p, width = 10,
              height = max(4, nrow(top_imp) * 0.3))

    roc_summary <- NULL
    if (has_proc) {
      roc_rows <- lapply(head(fi$variable_id, n_top), function(v) {
        r <- tryCatch(pROC::roc(response, rf_data[, v], quiet = TRUE, ci = TRUE,
                                ci.method = roc_ci_method), error = function(e) NULL)
        if (is.null(r)) return(NULL)
        ci <- as.numeric(r$ci)
        data.frame(variable_id = v, display_name = label_of(v), auc = as.numeric(pROC::auc(r)),
                   ci_lower = ci[1], ci_upper = ci[3], stringsAsFactors = FALSE)
      })
      roc_summary <- do.call(rbind, roc_rows)
      if (!is.null(roc_summary)) {
        roc_summary <- roc_summary[order(roc_summary$auc, decreasing = TRUE), , drop = FALSE]
        write.csv(roc_summary, file.path(cmp_dir, "roc_summary.csv"), row.names = FALSE)
      }
    }
    results[[cmp_name]] <- list(
      feature_importance = head(fi, n_top),
      roc_summary = roc_summary,
      roc_auc = if (is.null(roc_summary)) NA_real_ else max(roc_summary$auc),
      n_auc_above_threshold = if (is.null(roc_summary)) NA_integer_ else sum(roc_summary$auc >= roc_threshold)
    )
  }
  cat(sprintf("   %d comparisons analysed\n", length(results)))
  if (length(results) == 0) NULL else results
}

#' Chemical class of each feature
#'
#' Identified features (Level 1-2): the class of the identified compound. Other
#' features: the CANOPUS (ClassyFire, from MS2) class, else the class every
#' accurate-mass candidate of the feature shares (class-level annotation is MSI
#' Level 3 by definition; features whose candidates disagree get none).
#'
#' @param app3 Candidate table from finalize_annotation()
#' @param level "super_class" or "class"
#' @return data.frame(variable_id, taxonomy, source)
feature_chemical_class <- function(app3, level = "super_class") {
  hmdb_col <- if (level == "class") "Class" else "Super_class"
  canopus_col <- if (level == "class") "canopus_cf_class" else "canopus_cf_superclass"
  get <- function(col) if (col %in% colnames(app3)) as.character(app3[[col]]) else rep(NA_character_, nrow(app3))
  d <- data.frame(variable_id = as.character(app3$variable_id),
                  level = get("feature_level"), evidence = get("evidence"),
                  hmdb_class = get(hmdb_col), canopus = get(canopus_col),
                  is_rank1 = if ("candidate_rank" %in% colnames(app3)) app3$candidate_rank == 1 else TRUE,
                  stringsAsFactors = FALSE)
  d$hmdb_class[!nzchar(d$hmdb_class)] <- NA
  d$canopus[!nzchar(d$canopus)] <- NA

  out <- lapply(split(d, d$variable_id), function(x) {
    ident <- grepl("Level [12]", x$level[1])
    rank1 <- x[which(x$is_rank1)[1], , drop = FALSE]
    if (ident && !is.na(rank1$hmdb_class)) return(c(rank1$hmdb_class, "identified compound"))
    can <- stats::na.omit(x$canopus)
    if (length(can) > 0) return(c(can[1], "CANOPUS"))
    ms1 <- x$hmdb_class[x$evidence %in% "MS1_endogenous"]
    if (length(ms1) > 0 && !anyNA(ms1) && length(unique(ms1)) == 1) return(c(ms1[1], "all mass candidates"))
    c(NA_character_, NA_character_)
  })
  res <- data.frame(variable_id = names(out),
                    taxonomy = vapply(out, `[`, character(1), 1),
                    source = vapply(out, `[`, character(1), 2),
                    stringsAsFactors = FALSE)
  res[!is.na(res$taxonomy), , drop = FALSE]
}

#' Chemical class enrichment of significant features
#'
#' One class per feature (feature_chemical_class()); Fisher's exact test
#' (one-sided, over-representation) of significant vs all classified
#' features per class with >= 2 features, BH over the tested classes.
#'
#' @param app3 Candidate table from finalize_annotation()
#' @param diff_results Output of run_differential()
#' @param cfg Validated config list
#' @param out_dir Output directory (Chemical_Class_Enrichment/)
#' @return Named list per comparison: list(result, taxon_level)
run_chem_class_enrichment <- function(app3, diff_results, cfg, out_dir) {
  cat("\n----- Chemical Class Enrichment Analysis -----\n")
  cc_cfg <- cfg$metaboanalyst_advanced$chem_class_enrichment
  taxon_level <- cc_cfg$taxon_level %||% "super_class"
  p_cutoff <- as.numeric(cc_cfg$p_cutoff %||% 0.05)

  tax <- feature_chemical_class(app3, taxon_level)
  if (nrow(tax) == 0) {
    cat("   No feature has a chemical class, skipping\n")
    return(NULL)
  }
  cat(sprintf("   %d features with a %s (%s)\n", nrow(tax), taxon_level,
              paste(sprintf("%s: %d", names(table(tax$source)), table(tax$source)), collapse = ", ")))
  classes <- names(which(table(tax$taxonomy) >= 2))

  results <- list()
  for (cmp_name in names(diff_results)) {
    dr <- diff_results[[cmp_name]]
    if (is.null(dr)) next
    label <- dr$comparison$label %||% cmp_name
    sig <- tax$variable_id %in% significant_ids(dr, cfg)
    if (sum(sig) < 3) next
    res <- do.call(rbind, lapply(classes, function(cls) {
      in_cls <- tax$taxonomy == cls
      tab <- matrix(c(sum(sig & in_cls), sum(sig & !in_cls),
                      sum(!sig & in_cls), sum(!sig & !in_cls)), 2)
      data.frame(class_name = cls, hit_count = sum(sig & in_cls), total_count = sum(in_cls),
                 expected = round(sum(in_cls) * sum(sig) / nrow(tax), 2),
                 p_value = stats::fisher.test(tab, alternative = "greater")$p.value,
                 stringsAsFactors = FALSE)
    }))
    if (is.null(res) || nrow(res) == 0) next
    res$fold_enrichment <- round(res$hit_count / pmax(res$expected, 1e-9), 2)
    res$fdr <- stats::p.adjust(res$p_value, "BH")
    res <- res[order(res$p_value), , drop = FALSE]

    cmp_dir <- file.path(out_dir, cmp_name)
    dir.create(cmp_dir, showWarnings = FALSE, recursive = TRUE)
    write.csv(res, file.path(cmp_dir, "chem_class_enrichment.csv"), row.names = FALSE)
    sig_res <- res[res$fdr < p_cutoff, , drop = FALSE]
    if (nrow(sig_res) > 0) {
      sig_res$class_name <- factor(sig_res$class_name, levels = rev(sig_res$class_name))
      p <- ggplot(sig_res, aes(x = fold_enrichment, y = class_name)) +
        geom_bar(stat = "identity", aes(fill = -log10(fdr + 1e-10)), alpha = 0.85) +
        geom_text(aes(label = hit_count), hjust = -0.2, size = 3) +
        scale_fill_gradientn(colours = c("#46bac2", "#f7ca64", "#e34a33"), name = "-log10(FDR)") +
        labs(x = "Fold enrichment", y = "",
             title = sprintf("Chemical class enrichment — %s (%s)", label, taxon_level)) +
        theme_bw() +
        theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, size = 12))
      save_plot(file.path(cmp_dir, "chem_class_enrichment"), plot = p, width = 10,
                height = max(4, nrow(sig_res) * 0.35))
    }
    results[[cmp_name]] <- list(result = res, taxon_level = taxon_level)
  }
  n_sig <- sum(vapply(results, function(r) any(r$result$fdr < p_cutoff), logical(1)))
  cat(sprintf("   %d comparisons tested, %d with a class at FDR < %.2f\n",
              length(results), n_sig, p_cutoff))
  if (length(results) == 0) NULL else results
}

#' Run the advanced analyses (Step 10)
#'
#' @param object2 mass_dataset object (pre-processed)
#' @param sample_info Sample information data.frame (unused; kept for the caller)
#' @param diff_results Output of run_differential()
#' @param app3 Candidate table from finalize_annotation() (all features)
#' @param cfg Validated config list
#' @return Invisible list with module results (also stored in *_results_global
#'   for the web export)
run_metaboanalyst_advanced <- function(object2, sample_info, diff_results,
                                       app3, cfg) {
  cat("\n===== Step 10: Advanced analysis (identified-compound ORA, biomarkers, chemical classes) =====\n")
  ma_cfg <- cfg$metaboanalyst_advanced
  if (!isTRUE(ma_cfg$enabled)) {
    cat("   Advanced analysis disabled in config, skipping\n")
    return(invisible(NULL))
  }
  out_dir <- file.path(cfg$project$output_dir, "10_MetaboAnalystR_Advanced")
  # outputs of the MetaboAnalystR version (never complete) would mix with these
  unlink(list.files(out_dir, full.names = TRUE), recursive = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  run_module <- function(enabled, name, expr) {
    if (!isTRUE(enabled)) {
      cat(sprintf("   %s disabled in config, skipping\n", name))
      return(NULL)
    }
    tryCatch(expr, error = function(e) {
      cat(sprintf("   WARNING: %s failed: %s\n", name, conditionMessage(e)))
      NULL
    })
  }
  ora_result <- run_module(ma_cfg$msea$enabled, "Identified-compound pathway ORA",
                           run_identified_pathway_ora(app3, diff_results, cfg, file.path(out_dir, "MSEA")))
  if (isTRUE(ma_cfg$pathway_topology$enabled)) {
    cat("   Pathway topology needs KEGG pathway graphs; not computed (pathways: Step 7)\n")
  }
  bio_result <- run_module(ma_cfg$biomarker$enabled, "Biomarker analysis",
                           run_biomarker_analysis(object2, app3, diff_results, cfg,
                                                  file.path(out_dir, "Biomarker_Analysis")))
  chem_result <- run_module(ma_cfg$chem_class_enrichment$enabled, "Chemical class enrichment",
                            run_chem_class_enrichment(app3, diff_results, cfg,
                                                      file.path(out_dir, "Chemical_Class_Enrichment")))
  cat("   Step 10 complete\n")

  msea_results_global <<- ora_result
  topology_results_global <<- NULL
  biomarker_results_global <<- bio_result
  chem_class_results_global <<- chem_result
  invisible(list(msea_result = ora_result, biomarker_result = bio_result,
                 chem_result = chem_result))
}
