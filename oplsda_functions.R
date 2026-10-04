###############################################################################
# oplsda_functions.R — OPLS-DA 分析模块
# 依赖: ropls, ggplot2
# 职责: OPLS-DA 建模 + Score Plot + S-Plot + VIP + 置换检验 (≥1000次)
#       输出所有结果供主 pipeline 调用 + 供 Step 6.5 注释更新
###############################################################################

# ══════════════════════════════════════════════════════════════════════════════
# 主入口函数
# ══════════════════════════════════════════════════════════════════════════════

#' OPLS-DA 全流程分析
#'
#' 对每个比较对执行: 数据校验 → 建模 → Score Plot → S-Plot → VIP → 置换检验
#' 结果保存为 .csv 表 + .pdf/.png 图 + .rds 模型（供 Step 6.5 注释更新）
#'
#' @param expression_data  表达矩阵 (features × samples)
#' @param sample_info      样本信息 (需含 sample_id, class 列)
#' @param comparisons      比较对列表
#' @param comparison_labels 比较对显示名向量
#' @param output_dir       输出目录
#' @param n_perm           置换检验次数 (默认 1000)
#' @param vip_threshold    VIP 高亮阈值 (默认 1.0)
#' @param seed             随机种子 (默认 42)
#'
#' @return 不可见地返回 metrics 汇总 data.frame
run_oplsda_analysis <- function(expression_data,
                                 sample_info,
                                 comparisons,
                                 comparison_labels,
                                 output_dir,
                                 n_perm           = 1000,
                                 vip_threshold    = 1.0,
                                 seed             = 42) {

  cat("\n===== Step 3.5: OPLS-DA 分析 =====\n")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  all_samples  <- colnames(expression_data)
  group_vec    <- sample_info$class[match(all_samples, sample_info$sample_id)]
  names(group_vec) <- all_samples

  # 配色（与主 pipeline 的 mycol 一致）
  group_colors <- list(
    "CT" = "#58CDD9",
    "A"  = "#D20A13",
    "B"  = "#088247",
    "AB" = "#FFD121"
  )

  all_metrics <- list()
  all_vip     <- list()

  for (cmp in comparisons) {
    cmp_name  <- cmp$name
    ctrl_grp  <- cmp$ctrl
    treat_grp <- cmp$treat
    label     <- comparison_labels[cmp_name]
    if (is.null(label) || is.na(label) || length(label) == 0) {
      label <- cmp_name
    }

    cat(sprintf("\n  => OPLS-DA: %s (%s vs %s)\n", label, treat_grp, ctrl_grp))

    # ── Step 3.5.1: 数据校验 ────────────────────────────────────────────────
    ctrl_samples  <- all_samples[group_vec == ctrl_grp]
    treat_samples <- all_samples[group_vec == treat_grp]

    # 组数校验
    n_groups <- length(unique(c(ctrl_grp, treat_grp)))
    if (n_groups != 2) {
      cat(sprintf("  ⛔ 错误: OPLS-DA 仅支持两组对比。检测到 %s 组 (%s)。\n",
                  n_groups, paste(unique(c(ctrl_grp, treat_grp)), collapse = "/")))
      cat("  ⛔ 请提供仅含两组的数据，或改用 PLS-DA 模块。跳过该比较对。\n")
      next
    }

    # 样本量检查
    if (length(ctrl_samples) < 3 || length(treat_samples) < 3) {
      cat(sprintf("  ⚠️  警告: 样本量过少 (Ctrl=%s, Treat=%s)，",
                  length(ctrl_samples), length(treat_samples)))
      cat("交叉验证可能无法稳定评估 Q²。跳过该比较对。\n")
      next
    }

    # ── Step 3.5.2: 数据准备 ────────────────────────────────────────────────
    sample_ids <- c(ctrl_samples, treat_samples)
    expr <- as.matrix(expression_data[, sample_ids, drop = FALSE])

    # 安全过滤: 移除全 NA / 常数 / 高缺失特征
    row_sds    <- apply(expr, 1, sd, na.rm = TRUE)
    na_prop    <- rowMeans(is.na(expr))
    keep       <- !is.na(row_sds) & row_sds > 0 & na_prop <= 0.5

    if (sum(keep) < 10) {
      cat(sprintf("  ⚠️  有效特征数不足 (%s < 10)，跳过 OPLS-DA\n", sum(keep)))
      next
    }
    expr <- expr[keep, , drop = FALSE]

    # 填补剩余 NA（行均值填补）
    expr <- t(apply(expr, 1, function(r) {
      r[is.na(r)] <- mean(r, na.rm = TRUE)
      r
    }))

    # 转置为 (samples × features)
    X <- t(expr)

    # 响应变量 (factor)
    Y <- factor(c(rep(ctrl_grp, length(ctrl_samples)),
                  rep(treat_grp, length(treat_samples))),
                levels = c(ctrl_grp, treat_grp))
    names(Y) <- sample_ids

    # ropls defaults to 7-fold CV, which errors when nrow(X) <= 7 — i.e. for
    # every 3-vs-3 comparison. Cap it at the largest valid value.
    cv_i <- max(2, min(7, nrow(X) - 1))

    # 颜色（保持与主 pipeline 一致）
    ctrl_color <- group_colors[[ctrl_grp]]
    treat_color <- group_colors[[treat_grp]]
    if (is.null(ctrl_color)) ctrl_color <- "#58CDD9"
    if (is.null(treat_color)) treat_color <- "#D20A13"
    cmp_colors <- c(treat_color, ctrl_color)
    names(cmp_colors) <- c(treat_grp, ctrl_grp)

    # ── Step 3.5.3: OPLS-DA 建模 ────────────────────────────────────────────
    set.seed(seed)

    oplsda_res <- tryCatch({
      ropls::opls(
        x       = X,
        y       = Y,
        crossvalI = cv_i,
        predI   = 1,
        orthoI  = NA,
        scaleC  = "standard",
        info.txtC = "none",
        fig.pdfC  = "none"
      )
    }, error = function(e) {
      cat(sprintf("  ⛔ OPLS-DA 建模失败: %s\n", conditionMessage(e)))
      return(NULL)
    })

    if (is.null(oplsda_res)) next

    # 验证模型是否有效（orthoI=NA 在部分数据上可能返回空模型）
    model_valid <- tryCatch({
      nrow(oplsda_res@summaryDF) > 0 &&
      ncol(oplsda_res@scoreMN) >= 1 &&
      length(oplsda_res@summaryDF[["R2Y(cum)"]]) > 0
    }, error = function(e) FALSE)

    if (!model_valid) {
      cat("  ⚠️  orthoI=NA 模型无效，尝试 orthoI=1 ...\n")
      oplsda_res <- tryCatch({
        ropls::opls(x = X, y = Y, predI = 1, orthoI = 1, crossvalI = cv_i,
                     scaleC = "standard", info.txtC = "none", fig.pdfC = "none")
      }, error = function(e) NULL)

      if (!is.null(oplsda_res) && nrow(oplsda_res@summaryDF) > 0) {
        cat("  ✓ orthoI=1 模型成功\n")
      } else {
        cat("  ⚠️  orthoI=1 也失败，尝试 orthoI=0 ...\n")
        oplsda_res <- tryCatch({
          ropls::opls(x = X, y = Y, predI = 1, orthoI = 0, crossvalI = cv_i,
                       scaleC = "standard", info.txtC = "none", fig.pdfC = "none")
        }, error = function(e) NULL)

        if (is.null(oplsda_res) || nrow(oplsda_res@summaryDF) == 0) {
          cat("  ⚠️  所有 OPLS-DA 模型均无效，跳过该比较对\n")
          next
        }
      }
    }

    # 提取核心指标
    r2x <- oplsda_res@summaryDF$`R2X(cum)`
    r2y <- oplsda_res@summaryDF$`R2Y(cum)`
    q2  <- oplsda_res@summaryDF$`Q2(cum)`
    cat(sprintf("  ✓ 模型: R²X(cum)=%.4f, R²Y(cum)=%.4f, Q²(cum)=%.4f\n", r2x, r2y, q2))

    # ── Step 3.5.4: Score Plot (t1 vs to1) ──────────────────────────────────
    score_plot_obj <- tryCatch({
      t1 <- oplsda_res@scoreMN[, 1]

      ortho_score_mat <- tryCatch(oplsda_res@orthoScoreMN, error = function(e) NULL)
      has_ortho <- !is.null(ortho_score_mat) && ncol(ortho_score_mat) >= 1

      scores_df <- data.frame(
        sample_id = names(t1),
        t1        = as.numeric(t1),
        group     = Y[names(t1)],
        stringsAsFactors = FALSE
      )

      y_col  <- NULL
      y_lab  <- NULL

      if (has_ortho) {
        to1     <- ortho_score_mat[, 1]
        scores_df$to1 <- as.numeric(to1)
        y_col   <- "to1"
        y_lab   <- "Orthogonal Component (to1)"
        cat("  ✓ 正交成分数 (orthoI):", ncol(ortho_score_mat), "\n")
      } else if (ncol(oplsda_res@scoreMN) >= 2) {
        scores_df$t2 <- oplsda_res@scoreMN[, 2]
        y_col <- "t2"
        y_lab <- paste0("Predictive Component (t2)")
        cat("  ✓ 无正交成分，使用第2预测成分\n")
      } else {
        cat("  · 仅含 1 个预测成分，无法绘制 2D Score Plot\n")
      }

      if (!is.null(y_col)) {
        p_score <- ggplot2::ggplot(
          scores_df,
          ggplot2::aes(x = t1, y = .data[[y_col]], color = group)
        ) +
          ggplot2::geom_point(size = 4, alpha = 0.85) +
          ggplot2::stat_ellipse(
            ggplot2::aes(fill = group),
            level = 0.95, type = "norm",
            geom = "polygon", alpha = 0.10, show.legend = FALSE
          ) +
          ggplot2::scale_color_manual(
            values = cmp_colors,
            labels = names(cmp_colors)
          ) +
          ggplot2::scale_fill_manual(values = cmp_colors, guide = "none") +
          ggplot2::labs(
            title = paste0("OPLS-DA Score Plot: ", label),
            x = "Predictive Component (t1)",
            y = y_lab,
            color = "Group"
          ) +
          ggplot2::theme_bw() +
          ggplot2::theme(
            panel.grid      = ggplot2::element_blank(),
            aspect.ratio    = 1,
            plot.title      = ggplot2::element_text(hjust = 0.5, size = 14),
            axis.title      = ggplot2::element_text(size = 13),
            axis.text       = ggplot2::element_text(size = 11),
            legend.position = "right"
          )

        save_plot(file.path(output_dir, paste0(cmp_name, "_score_plot")),
                  plot = p_score, width = 8, height = 7)
        cat("  ✓ Score Plot 已保存\n")
      }

      # 返回 Score Plot 对象（供标注更新时重用标题等信息）
      list(scores_df = scores_df, has_ortho = has_ortho, y_col = y_col, y_lab = y_lab,
           p_score = if (exists("p_score")) p_score else NULL)

    }, error = function(e) {
      cat(sprintf("  ⚠️  Score Plot 绘制失败: %s\n", conditionMessage(e)))
      return(list(scores_df = NULL, has_ortho = FALSE, y_col = NULL, y_lab = NULL, p_score = NULL))
    })

    # ── Step 3.5.5: VIP 值计算 ─────────────────────────────────────────────
    # Primary: use ropls' own VIP (@vipVn), which is populated whenever the
    # model is fitted with crossvalI (all fit branches above pass cv_i).
    # Falls back to the manual loading/weight formula only if the slot is
    # empty, and finally to NA.
    p_pred <- tryCatch(as.numeric(oplsda_res@loadingMN[, 1]),
                       error = function(e) rep(NA_real_, ncol(X)))

    vip_values <- tryCatch({
      slot_vip <- ropls::getVipVn(oplsda_res)
      if (!is.null(slot_vip) && length(slot_vip) == ncol(X)) {
        as.numeric(slot_vip)
      } else {
        cat("  · @vipVn 为空，回退手动 VIP 公式\n")
        t_pred <- oplsda_res@scoreMN[, 1]
        n_var  <- ncol(X)

        # 尝试获取正交权重和正交得分
        w_orth <- tryCatch(oplsda_res@orthoWeightMN, error = function(e) NULL)
        t_orth <- tryCatch(oplsda_res@orthoScoreMN, error = function(e) NULL)

        has_ortho_vip <- !is.null(w_orth) && !is.null(t_orth) &&
                         is.matrix(w_orth) && is.matrix(t_orth) &&
                         ncol(w_orth) >= 1 && ncol(t_orth) >= 1 &&
                         nrow(w_orth) == n_var

        SS_pred <- sum(t_pred^2)

        if (has_ortho_vip) {
          n_ortho   <- ncol(w_orth)
          SS_ortho  <- apply(t_orth^2, 2, sum)
          w_sq_norm <- sweep(w_orth^2, 2, SS_ortho, "*")
          sqrt(n_var * (p_pred^2 * SS_pred + rowSums(w_sq_norm)) /
               (SS_pred + sum(SS_ortho)))
        } else {
          # 无正交成分时退化为 PLS-DA VIP
          sqrt(n_var * p_pred^2)
        }
      }
    }, error = function(e) {
      cat(sprintf("  ⚠️  VIP 计算失败: %s\n", conditionMessage(e)))
      rep(NA_real_, ncol(X))
    })
    names(vip_values) <- colnames(X)

    # ── Step 3.5.6: S-Plot ──────────────────────────────────────────────────
    splot_obj <- tryCatch({
      # 获取 p[1] (loading) — ropls 1.42.0 的 modelDF 无 p1 列, 改用 loadingMN
      p1_vec <- tryCatch({
        if (!is.null(oplsda_res@loadingMN) && ncol(oplsda_res@loadingMN) >= 1) {
          as.numeric(oplsda_res@loadingMN[, 1])
        } else {
          NULL
        }
      }, error = function(e) NULL)

      cor_vec <- tryCatch({
        cor(t(expr), oplsda_res@scoreMN[, 1], use = "pairwise.complete.obs")[, 1]
      }, error = function(e) rep(NA_real_, nrow(expr)))

      if (!is.null(p1_vec) && length(p1_vec) == nrow(expr)) {
        s_df <- data.frame(
          variable_id       = rownames(expr),
          covariance_p1     = p1_vec,
          correlation_pcorr1 = cor_vec,
          VIP               = vip_values,
          stringsAsFactors  = FALSE
        )
        s_df$highlight <- !is.na(s_df$VIP) & s_df$VIP > vip_threshold
        n_highlight <- sum(s_df$highlight, na.rm = TRUE)

        p_splot <- ggplot2::ggplot(s_df, ggplot2::aes(x = covariance_p1, y = correlation_pcorr1)) +
          ggplot2::geom_point(
            ggplot2::aes(color = highlight, size = highlight),
            alpha = 0.65
          ) +
          ggplot2::scale_color_manual(
            values = c("FALSE" = "#808180FF", "TRUE" = "#D20A13"),
            labels = c("FALSE" = paste0("VIP ≤ ", vip_threshold),
                       "TRUE"  = paste0("VIP > ", vip_threshold)),
            name   = "Variable Importance"
          ) +
          ggplot2::scale_size_manual(values = c("FALSE" = 1.5, "TRUE" = 3), guide = "none") +
          ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
          ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
          ggplot2::labs(
            title = paste0("S-Plot: ", label),
            x = "Covariance (p[1])",
            y = "Correlation (p(corr)[1])"
          ) +
          ggplot2::theme_bw() +
          ggplot2::theme(
            panel.grid      = ggplot2::element_blank(),
            aspect.ratio    = 1,
            plot.title      = ggplot2::element_text(hjust = 0.5, size = 14),
            axis.title      = ggplot2::element_text(size = 13),
            axis.text       = ggplot2::element_text(size = 11),
            legend.position = "right"
          )

        save_plot(file.path(output_dir, paste0(cmp_name, "_splot")),
                  plot = p_splot, width = 9, height = 7)
        cat(sprintf("  ✓ S-Plot 已保存 (高亮 %s 个 VIP > %s 的变量)\n",
                    n_highlight, vip_threshold))

        list(s_df = s_df, p1_vec = p1_vec, cor_vec = cor_vec,
             p_splot = p_splot, n_highlight = n_highlight)
      } else {
        cat("  ⚠️  p[1] 向量不可用或维度不匹配，跳过 S-Plot\n")
        list(s_df = NULL, p1_vec = NULL, cor_vec = cor_vec, p_splot = NULL, n_highlight = 0)
      }
    }, error = function(e) {
      cat(sprintf("  ⚠️  S-Plot 绘制失败: %s\n", conditionMessage(e)))
      list(s_df = NULL, p1_vec = NULL, cor_vec = NULL, p_splot = NULL, n_highlight = 0)
    })

    # ── Step 3.5.7: 置换检验 ────────────────────────────────────────────────
    cat(sprintf("  · 置换检验 (n = %s)...\n", n_perm))

    Y_numeric <- as.numeric(Y)

    # One permutation fit (self-contained; safe to run on parallel workers)
    perm_one <- function(i) {
      Y_perm <- sample(Y)
      pcor <- cor(as.numeric(Y_perm), Y_numeric)

      tmp <- tryCatch({
        ropls::opls(x = X, y = Y_perm, crossvalI = cv_i,
                     predI = 1, orthoI = NA,
                     scaleC = "standard",
                     info.txtC = "none", fig.pdfC = "none")
      }, error = function(e) NULL)

      # 检查模型是否有效 (orthoI=NA 在 permuted data 上可能返回空模型)
      is_valid <- !is.null(tmp) &&
                  nrow(tmp@summaryDF) > 0 &&
                  "R2Y(cum)" %in% names(tmp@summaryDF)

      r2y <- NA_real_
      q2  <- NA_real_
      if (is_valid) {
        r2y <- tmp@summaryDF[["R2Y(cum)"]]
        q2  <- tmp@summaryDF[["Q2(cum)"]]
      } else {
        # 回退: 尝试 orthoI = 0 (PLS-DA) 用于 permuted data
        tmp2 <- tryCatch({
          ropls::opls(x = X, y = Y_perm, crossvalI = cv_i,
                       predI = 1, orthoI = 0,
                       scaleC = "standard",
                       info.txtC = "none", fig.pdfC = "none")
        }, error = function(e) NULL)

        if (!is.null(tmp2) && nrow(tmp2@summaryDF) > 0) {
          r2y <- tmp2@summaryDF[["R2Y(cum)"]]
          q2  <- tmp2@summaryDF[["Q2(cum)"]]
        }
      }
      list(cor = pcor, r2y = r2y, q2 = q2)
    }

    # Parallelize via future.apply using the backend initialized by
    # init_parallel_backend() (config: performance.cores / parallel_backend).
    #
    # NOTE: each ropls::opls() fit here costs ~1 s and ships the full X matrix
    # (thousands of features x dozens of samples) to every worker on every
    # dispatch, so fine-grained future_lapply over the permutations spends most
    # of its wall clock on serialisation. Chunking sends one X per worker per
    # chunk instead: n_perm/n_workers dispatches instead of n_perm.
    use_parallel <- requireNamespace("future.apply", quietly = TRUE) &&
                    !inherits(future::plan(), "sequential")
    if (use_parallel) {
      n_workers <- future::nbrOfWorkers()
      cat(sprintf("  · 置换检验并行: %s workers\n", n_workers))
      perm_one_chunk <- function(idx) lapply(idx, perm_one)
      chunks <- split(seq_len(n_perm),
                      cut(seq_len(n_perm), breaks = n_workers, labels = FALSE))
      perm_res <- unlist(
        future.apply::future_lapply(chunks, perm_one_chunk, future.seed = TRUE),
        recursive = FALSE)
    } else {
      perm_res <- lapply(seq_len(n_perm), perm_one)
    }

    perm_cor <- vapply(perm_res, function(r) r$cor, numeric(1))
    perm_r2y <- vapply(perm_res, function(r) r$r2y, numeric(1))
    perm_q2  <- vapply(perm_res, function(r) r$q2,  numeric(1))
    rm(perm_res)
    invisible(gc())

    # 移除失败的置换
    valid     <- !is.na(perm_r2y) & !is.na(perm_q2)
    perm_cor  <- perm_cor[valid]
    perm_r2y  <- perm_r2y[valid]
    perm_q2   <- perm_q2[valid]
    n_success <- length(perm_r2y)

    cat(sprintf("  · 有效置换: %s / %s\n", n_success, n_perm))

    # 计算截距
    q2_intercept  <- NA_real_
    r2y_intercept <- NA_real_

    if (n_success >= 3) {
      # Q² 回归: Q²_perm ~ cor(perm_Y, orig_Y)
      lm_q2  <- stats::lm(perm_q2 ~ perm_cor)
      q2_intercept  <- unname(stats::coef(lm_q2)[1])

      # R²Y 回归
      lm_r2y <- stats::lm(perm_r2y ~ perm_cor)
      r2y_intercept <- unname(stats::coef(lm_r2y)[1])

      # 过拟合检查
      if (q2_intercept > 0.05) {
        cat(sprintf("  ⚠️  过拟合风险: Q² 截距 = %.4f > 0.05\n", q2_intercept))
      } else {
        cat(sprintf("  ✓ Q² 截距 = %.4f (可接受)\n", q2_intercept))
      }

      # ── 置换检验图 ──
      # 构建数据（确保 r2y/q2 非 NULL 以保持长度一致）
      r2y_safe <- if (is.null(r2y) || length(r2y) == 0) NA_real_ else r2y
      q2_safe  <- if (is.null(q2)  || length(q2)  == 0) NA_real_ else q2

      r2y_points <- data.frame(
        cor   = c(perm_cor, 1),
        value = c(perm_r2y, r2y_safe),
        type  = c(rep("Permuted", n_success), "Original"),
        metric = "R²Y"
      )
      q2_points <- data.frame(
        cor   = c(perm_cor, 1),
        value = c(perm_q2, q2_safe),
        type  = c(rep("Permuted", n_success), "Original"),
        metric = "Q²"
      )
      plot_df <- rbind(r2y_points, q2_points)

      # 回归线预测
      cor_seq <- seq(min(c(perm_cor, 0)), max(c(perm_cor, 1)), length.out = 60)
      r2y_pred <- predict(lm_r2y, newdata = data.frame(perm_cor = cor_seq))
      q2_pred  <- predict(lm_q2,  newdata = data.frame(perm_cor = cor_seq))

      reg_df <- rbind(
        data.frame(cor = cor_seq, value = r2y_pred, metric = "R²Y"),
        data.frame(cor = cor_seq, value = q2_pred,  metric = "Q²")
      )

      # 截距标注
      intercept_ann <- sprintf("Q² intercept = %.4f", q2_intercept)
      if (q2_intercept > 0.05) intercept_ann <- paste0(intercept_ann, " ⚠️")

      p_perm <- ggplot2::ggplot() +
        ggplot2::geom_point(
          data = plot_df[plot_df$type == "Permuted", ],
          ggplot2::aes(x = cor, y = value, color = metric),
          size = 2, alpha = 0.55
        ) +
        ggplot2::geom_point(
          data = plot_df[plot_df$type == "Original", ],
          ggplot2::aes(x = cor, y = value, fill = metric),
          size = 4, shape = 23, color = "black"
        ) +
        ggplot2::geom_line(
          data = reg_df,
          ggplot2::aes(x = cor, y = value, color = metric),
          linewidth = 0.7, linetype = "dashed"
        ) +
        ggplot2::facet_wrap(~ metric, scales = "free_y", ncol = 2) +
        ggplot2::labs(
          title    = paste0("Permutation Test: ", label),
          subtitle = sprintf("n = %s permutations | %s",
                             n_success, intercept_ann),
          x        = "Correlation with original Y",
          y        = "Value"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          panel.grid      = ggplot2::element_blank(),
          plot.title      = ggplot2::element_text(hjust = 0.5, size = 14),
          plot.subtitle   = ggplot2::element_text(hjust = 0.5, size = 10, color = "grey30"),
          strip.background = ggplot2::element_rect(fill = "grey90"),
          strip.text      = ggplot2::element_text(size = 12),
          legend.position = "none"
        )

      save_plot(file.path(output_dir, paste0(cmp_name, "_permutation")),
                plot = p_perm, width = 10, height = 5.5)
      cat("  ✓ 置换检验图已保存\n")
    } else {
      cat(sprintf("  ⚠️  有效置换不足 (%s < 3)，跳过回归与绘图\n", n_success))
    }

    # ── Step 3.5.8: 构建 Feature Importance Table ──────────────────────────
    p1_val   <- if (!is.null(splot_obj$p1_vec)) splot_obj$p1_vec else NA_real_
    cor_val  <- if (!is.null(splot_obj$cor_vec)) splot_obj$cor_vec else NA_real_

    feature_df <- data.frame(
      variable_id = rownames(expr),
      VIP         = vip_values,
      p1          = p1_val,
      pcorr1      = cor_val,
      stringsAsFactors = FALSE
    )

    # ── 保存模型 RDS（供 Step 6.5 注释更新） ─────────────────────────────
    model_save <- list(
      comparison  = cmp_name,
      label       = label,
      control     = ctrl_grp,
      treat       = treat_grp,
      n_features  = nrow(expr),

      # 指标
      r2x           = r2x,
      r2y           = r2y,
      q2            = q2,
      q2_intercept  = q2_intercept,
      r2y_intercept = r2y_intercept,
      n_perm        = n_success,

      # 特征重要性
      feature_importance = feature_df,

      # 绘图用数据
      scores       = score_plot_obj$scores_df,
      has_ortho    = score_plot_obj$has_ortho,
      splot_data   = splot_obj$s_df,
      kept_features = rownames(expr),

      # 置换原始数据
      perm_cor = if (exists("perm_cor")) perm_cor else NULL,
      perm_r2y = if (exists("perm_r2y")) perm_r2y else NULL,
      perm_q2  = if (exists("perm_q2")) perm_q2 else NULL
    )

    saveRDS(model_save,
            file.path(output_dir, paste0(cmp_name, "_oplsda_model.rds")))

    # ── 收集汇总 ──────────────────────────────────────────────────────────
    tryCatch({
      # Ensure all values are single-length (not NULL, not character(0))
      safe_label <- if (is.null(label) || length(label) == 0) {
        cmp_name
      } else if (is.na(label)) {
        cmp_name
      } else {
        as.character(label)
      }
      safe_r2x <- if (length(r2x) == 0) NA_real_ else as.numeric(r2x)
      safe_r2y <- if (length(r2y) == 0) NA_real_ else as.numeric(r2y)
      safe_q2  <- if (length(q2)  == 0) NA_real_ else as.numeric(q2)
      safe_q2_intercept  <- if (length(q2_intercept)  == 0) NA_real_ else as.numeric(q2_intercept)
      safe_r2y_intercept <- if (length(r2y_intercept) == 0) NA_real_ else as.numeric(r2y_intercept)

      all_metrics[[cmp_name]] <- data.frame(
        Comparison    = cmp_name,
        Label         = safe_label,
        Control       = ctrl_grp,
        Treat         = treat_grp,
        R2X_cum       = round(safe_r2x, 4),
        R2Y_cum       = round(safe_r2y, 4),
        Q2_cum        = round(safe_q2, 4),
        Perm_N        = n_success,
        Q2_Intercept  = round(safe_q2_intercept, 4),
        R2Y_Intercept = round(safe_r2y_intercept, 4),
        Overfit_Risk  = if (!is.na(safe_q2_intercept) && safe_q2_intercept > 0.05) "YES" else "NO",
        N_Features    = nrow(expr),
        stringsAsFactors = FALSE
      )

      feat_for_vip <- feature_df
      feat_for_vip$comparison <- cmp_name
      all_vip[[cmp_name]] <- feat_for_vip

      cat("  ✓ OPLS-DA 完成\n")
    }, error = function(e) {
      cat(sprintf("  ⚠️  Metrics 收集失败 (跳过): %s\n", conditionMessage(e)))
    })
  }

  # ── 汇总输出 ─────────────────────────────────────────────────────────────
  if (length(all_metrics) > 0) {
    metrics_summary <- dplyr::bind_rows(all_metrics)
    utils::write.csv(metrics_summary,
                     file.path(output_dir, "oplsda_metrics.csv"),
                     row.names = FALSE)
    cat(sprintf("\n  ✓ Metrics 汇总表 (%s 个比较对): %s\n",
                nrow(metrics_summary),
                file.path(output_dir, "oplsda_metrics.csv")))

    vip_summary <- dplyr::bind_rows(all_vip)
    # reorder columns: comparison first
    vip_summary <- vip_summary[, c("comparison", setdiff(colnames(vip_summary), "comparison"))]
    utils::write.csv(vip_summary,
                     file.path(output_dir, "oplsda_vip.csv"),
                     row.names = FALSE)
    cat(sprintf("  ✓ VIP 汇总表 (%s 行): %s\n",
                nrow(vip_summary),
                file.path(output_dir, "oplsda_vip.csv")))
  } else {
    cat("\n  · 无成功的 OPLS-DA 模型\n")
  }

  cat("\n===== OPLS-DA 分析完成 =====\n")
  invisible(all_metrics)
}


# ══════════════════════════════════════════════════════════════════════════════
# 注释更新函数（Step 6.5 调用）
# ══════════════════════════════════════════════════════════════════════════════

#' 用代谢物注释更新 OPLS-DA 结果
#'
#' 读取 Step 3.5 保存的 .rds 文件，合并 app3 注释（含 MSI 置信度），
#' 重新生成带化合物名和置信度等级的 S-Plot 和 Feature Importance Table。
#'
#' @param oplsda_dir    OPLSDA 结果目录（含 .rds 文件）
#' @param output_dir    注释更新版输出目录
#' @param app3          注释结果 data.frame（需含 variable_id, Compound.name, confidence_level）
#' @param vip_threshold VIP 高亮阈值
update_oplsda_with_annotation <- function(oplsda_dir,
                                          output_dir,
                                          app3,
                                          vip_threshold = 1.0) {

  cat("\n  => 更新 OPLS-DA 结果（注释标注版，含置信度）...\n")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # 查找所有模型文件
  rds_files <- list.files(oplsda_dir, pattern = "_oplsda_model\\.rds$", full.names = TRUE)

  if (length(rds_files) == 0) {
    cat("  => 未找到 OPLS-DA 模型文件，跳过\n")
    return(invisible(NULL))
  }

  if (is.null(app3) || nrow(app3) == 0) {
    cat("  => 注释数据为空，跳过 OPLS-DA 注释更新\n")
    return(invisible(NULL))
  }

  # 确保 app3 有 variable_id 列
  if (!"variable_id" %in% colnames(app3)) {
    cat("  => app3 缺少 variable_id 列，跳过\n")
    return(invisible(NULL))
  }

  # 检测是否有 MSI 置信度列
  has_confidence <- "confidence_level" %in% colnames(app3)

  all_vip_annotated <- list()

  for (rds_file in rds_files) {
    model <- tryCatch(readRDS(rds_file), error = function(e) NULL)
    if (is.null(model)) {
      cat(sprintf("  ⚠️  无法读取: %s\n", basename(rds_file)))
      next
    }

    cmp_name <- model$comparison
    label    <- model$label
    cat(sprintf("  · 正在更新: %s\n", label))

    # ── 合并注释（含置信度等级）──
    fi <- model$feature_importance
    # 动态选择要 join 的列
    join_cols <- c("variable_id", "Compound.name", "HMDB.ID", "KEGG.ID")
    if (has_confidence) {
      join_cols <- c(join_cols, "confidence_level")
    }
    fi_annotated <- fi %>%
      dplyr::left_join(
        app3 %>% dplyr::select(dplyr::any_of(join_cols)),
        by = "variable_id"
      ) %>%
      dplyr::mutate(
        label_name = dplyr::if_else(
          !is.na(Compound.name) & Compound.name != "",
          Compound.name,
          variable_id
        )
      )

    # ── 保存注释版 Feature Importance Table ──
    utils::write.csv(fi_annotated,
                     file.path(output_dir, paste0(cmp_name, "_vip_annotated.csv")),
                     row.names = FALSE)
    fi_annotated$comparison <- cmp_name
    all_vip_annotated[[cmp_name]] <- fi_annotated

    n_annotated <- sum(!is.na(fi_annotated$Compound.name) & fi_annotated$Compound.name != "")
    n_confident <- if (has_confidence) {
      sum(!is.na(fi_annotated$confidence_level) &
            grepl("Level [12]", fi_annotated$confidence_level))
    } else 0
    cat(sprintf("    Feature Importance: %s / %s 已注释",
                n_annotated, nrow(fi_annotated)))
    if (has_confidence && n_confident > 0) {
      cat(sprintf("（其中 %d 个为 Level 1-2 高置信度）", n_confident))
    }
    cat("\n")

    # ── 注释版 S-Plot ──
    if (!is.null(model$splot_data) && nrow(model$splot_data) > 0) {
      s_df <- model$splot_data
      # 合并注释和置信度
      splot_join_cols <- c("variable_id", "Compound.name")
      if (has_confidence) {
        splot_join_cols <- c(splot_join_cols, "confidence_level")
      }
      s_df <- s_df %>%
        dplyr::left_join(
          app3 %>% dplyr::select(dplyr::any_of(splot_join_cols)),
          by = "variable_id"
        ) %>%
        dplyr::mutate(
          label_name = dplyr::if_else(
            !is.na(Compound.name) & Compound.name != "",
            Compound.name,
            variable_id
          )
        )

      # 筛选 Top VIP（VIP > threshold 且按 VIP 降序取前 15 个有注释的）
      top_label <- s_df %>%
        dplyr::filter(!is.na(VIP) & VIP > vip_threshold) %>%
        dplyr::arrange(dplyr::desc(VIP)) %>%
        dplyr::slice_head(n = 15)

      if (nrow(top_label) > 0) {
        # 显示化合物名 + 置信度等级（如果可用）
        top_label <- top_label %>%
          dplyr::mutate(
            display_label = dplyr::case_when(
              !is.na(Compound.name) & Compound.name != "" & has_confidence &
                !is.na(confidence_level) ~ {
                  # Shorten confidence level for display
                  cl_short <- gsub("Level (\\d).*", "L\\1", confidence_level)
                  paste0(Compound.name, " [", cl_short, "]")
                },
              !is.na(Compound.name) & Compound.name != "" ~
                paste0(Compound.name, " (", variable_id, ")"),
              TRUE ~ variable_id
            )
          )
      }

      p_splot_ann <- ggplot2::ggplot(
        s_df, ggplot2::aes(x = covariance_p1, y = correlation_pcorr1)
      ) +
        ggplot2::geom_point(
          ggplot2::aes(color = highlight, size = highlight),
          alpha = 0.65
        ) +
        ggplot2::scale_color_manual(
          values = c("FALSE" = "#808180FF", "TRUE" = "#D20A13"),
          labels = c("FALSE" = paste0("VIP ≤ ", vip_threshold),
                     "TRUE"  = paste0("VIP > ", vip_threshold)),
          name   = "Variable Importance"
        ) +
        ggplot2::scale_size_manual(values = c("FALSE" = 1.5, "TRUE" = 3), guide = "none") +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                            color = "grey50", linewidth = 0.5) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                            color = "grey50", linewidth = 0.5)

      # 添加标签
      if (nrow(top_label) > 0) {
        p_splot_ann <- p_splot_ann +
          ggrepel::geom_text_repel(
            data = top_label,
            ggplot2::aes(label = display_label),
            size        = 3,
            max.overlaps = 15,
            box.padding  = 0.4,
            point.padding = 0.3,
            force        = 2,
            segment.color = "grey50",
            min.segment.length = 0.2
          )
      }

      p_splot_ann <- p_splot_ann +
        ggplot2::labs(
          title = paste0("S-Plot (Annotated): ", label),
          x = "Covariance (p[1])",
          y = "Correlation (p(corr)[1])"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          panel.grid      = ggplot2::element_blank(),
          aspect.ratio    = 1,
          plot.title      = ggplot2::element_text(hjust = 0.5, size = 14),
          axis.title      = ggplot2::element_text(size = 13),
          axis.text       = ggplot2::element_text(size = 11),
          legend.position = "right"
        )

      save_plot(file.path(output_dir, paste0(cmp_name, "_splot_annotated")),
                plot = p_splot_ann, width = 10, height = 8)
      cat(sprintf("    S-Plot (标注版) 已保存：标注 %s 个 Top VIP 特征\n",
                  nrow(top_label)))
    }

    # ── 注释版 Score Plot ──
    # Score Plot 是样本级别的，不需要用代谢物名标注
    # 但如果有备注意义，可加上样本名标注，这里跳过
  }

  # ── 合并注释版 VIP 汇总 ──
  if (length(all_vip_annotated) > 0) {
    vip_annot_summary <- dplyr::bind_rows(all_vip_annotated)
    # reorder columns: comparison first, then key annotation columns
    reorder_cols <- c("comparison", "variable_id",
                      "Compound.name", if (has_confidence) "confidence_level" else NULL,
                      "VIP", "p1", "pcorr1")
    reorder_cols <- reorder_cols[!sapply(reorder_cols, is.null)]
    other_cols <- setdiff(colnames(vip_annot_summary), reorder_cols)
    vip_annot_summary <- vip_annot_summary[, c(reorder_cols, other_cols)]
    utils::write.csv(vip_annot_summary,
                     file.path(output_dir, "oplsda_vip_annotated.csv"),
                     row.names = FALSE)
    cat(sprintf("  ✓ 注释版 VIP 汇总表: %s 行 -> %s\n",
                nrow(vip_annot_summary),
                file.path(output_dir, "oplsda_vip_annotated.csv")))

    # 打印置信度分布（如果有）
    if (has_confidence && "confidence_level" %in% colnames(vip_annot_summary)) {
      conf_dist <- table(vip_annot_summary$confidence_level, useNA = "ifany")
      cat("  ✓ VIP 置信度分布:\n")
      for (nm in names(conf_dist)) {
        cat(sprintf("      %s: %d\n", nm, conf_dist[nm]))
      }
    }
  }

  cat("  => OPLS-DA 注释更新完成\n")
}
