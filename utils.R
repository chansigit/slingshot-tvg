# Slingshot Analysis Steps:

# 1. Obtain cell embeddings and cluster labels, Specify starting/ending clusters; 
# 2. Run slingshot algorithm to get the fitted curves and pseudotimes; (`run_slingshot_curves(seu, start.clus)`)
# 3. Draw the identified curves above the cell embedding; (`DimPlot(...) & plot_slingshot_curves_layer(curve_df)`)
# 4. Identify temporally varying genes (TVG) using pseudotime and expression matrix; (`compute_gam_temporal_metrics(expr_mat, pseudotime_vec, genes)`)
# 5. Construct temporal gene expression matrix (TGE) along the curve using TVG; (`compute_temporal_matrix(expr_mat, pseudotime_vec, genes)`)
# 6. Draw a heatmap visualizing TGE.


library(slingshot)
library(dplyr)
library(tibble)
library(purrr)
library(tradeSeq)

library(slingshot)
library(dplyr)
library(tibble)
library(purrr)

run_slingshot_curves <- function(seu,
                                 cluster_col = "celltype3a3",
                                 reduced_name = "umap",
                                 start.clus = NULL,
                                 end.clus = NULL,
                                 omega = TRUE,
                                 omega_scale = 1.2,
                                 stretch = 0,
                                 verbose = TRUE) {

  # 1️⃣ 提取数据
  rd <- seu@reductions[[reduced_name]]@cell.embeddings
  colnames(rd)[1:2] <- c("UMAP1", "UMAP2")
  cl <- seu@meta.data[[cluster_col]]

  # 2️⃣ 运行 slingshot
  pto <- slingshot(
    rd, cl,
    start.clus = start.clus,
    end.clus   = end.clus,
    omega = omega,
    omega_scale = omega_scale,
    stretch = stretch
  )

  sds <- as.SlingshotDataSet(pto)
  n_lineages <- length(slingLineages(sds))
  if (verbose) cat(n_lineages, "lineages identified\n")

  # 3️⃣ 提取曲线与伪时间
  curve_list <- slingCurves(sds)
  curve_df <- map_dfr(seq_along(curve_list), function(i){
    crv <- curve_list[[i]]
    n <- nrow(crv$s)
    data.frame(crv$s[, 1:2]) %>%
      setNames(c("UMAP1", "UMAP2")) %>%
      mutate(
        lineage = paste0("L", i),
        pseudotime = seq(
          min(crv$lambda, na.rm = TRUE),
          max(crv$lambda, na.rm = TRUE),
          length.out = n
        )
      )
  })

  # 返回完整结果
  list(
    pto = pto,
    sds = sds,
    curve_list = curve_list,
    curve_df = curve_df
  )
}





library(ggplot2)
library(dplyr)
plot_slingshot_curves_layer <- function(curve_df,
                                        dim_cols = NULL,
                                        color_by_pseudotime = FALSE,
                                        show_start = TRUE,
                                        linewidth = 1.2,
                                        color = "black",
                                        start_shape = 8,
                                        start_size = 4,
                                        start_color = "gold",
                                        progress_cutoff = 1.0) {
    # ⚙️ 防御式检查
    if (progress_cutoff <= 0 || progress_cutoff > 1) {
        stop("progress_cutoff must be in (0, 1].")
    }

    # 🧹 确保数据是 tibble 且 lineage 为字符
    curve_df <- tibble::as_tibble(curve_df)
    curve_df <- dplyr::mutate(curve_df, lineage = as.character(lineage))

    # ⚙️ 自动识别坐标列
    if (is.null(dim_cols)) {
        dim_candidates <- setdiff(colnames(curve_df), c("lineage", "pseudotime"))
        if (length(dim_candidates) < 2) {
            stop("curve_df must contain at least two numeric coordinate columns.")
        }
        dim_cols <- dim_candidates[1:2]
    }

    colnames(curve_df)[match(dim_cols, colnames(curve_df))] <- c("X", "Y")

    # ✂️ 根据 progress_cutoff 裁剪曲线
    if ("pseudotime" %in% colnames(curve_df) && progress_cutoff < 1) {
        curve_df <- curve_df %>%
            group_by(lineage) %>%
            filter(pseudotime <= quantile(pseudotime, progress_cutoff, na.rm = TRUE)) %>%
            ungroup()
    }

    # 🖊️ 主曲线层
    if (color_by_pseudotime && "pseudotime" %in% colnames(curve_df)) {
        curve_layer <- geom_path(
            data = curve_df,
            aes(x = X, y = Y, group = lineage, color = pseudotime),
            linewidth = linewidth, lineend = "round"
        )
    } else {
        curve_layer <- geom_path(
            data = curve_df,
            aes(x = X, y = Y, group = lineage),
            color = color, linewidth = linewidth, lineend = "round"
        )
    }

    # 🌟 起点层
    if (show_start) {
        curve_starts <- curve_df %>%
            group_by(lineage) %>%
            {if ("pseudotime" %in% colnames(.))
                arrange(., pseudotime, .by_group = TRUE)
             else .} %>%
            slice_head(n = 1) %>%
            ungroup()

        start_layer <- geom_point(
            data = curve_starts,
            aes(x = X, y = Y),
            shape = start_shape,
            size = start_size,
            color = start_color,
            stroke = 1.2
        )

        return(list(curve_layer, start_layer))
    } else {
        return(curve_layer)
    }
}




library(mgcv)
library(dplyr)
library(purrr)

compute_gam_temporal_metrics <- function(expr_mat,
                                         pseudotime_vec,
                                         genes = NULL,
                                         family = gaussian(),
                                         k = 6,
                                         grid_n = 200,
                                         scale_pt = TRUE,
                                         na_action = na.omit) {

    stopifnot(ncol(expr_mat) == length(pseudotime_vec))

    keep_cells <- !is.na(pseudotime_vec)
    if (!any(keep_cells)) {
        stop("No cells with non-NA pseudotime values.")
    }
    X <- expr_mat[, keep_cells, drop = FALSE]
    t <- as.numeric(pseudotime_vec[keep_cells])

    if (!is.null(genes)) {
        genes <- intersect(genes, rownames(X))
        if (length(genes) == 0) {
            stop("None of the input genes found in expression matrix.")
        }
        X <- X[genes, , drop = FALSE]
    }

    if (scale_pt) {
        rng <- range(t, na.rm = TRUE)
        t_for_fit <- if (diff(rng) > 0) (t - rng[1]) / diff(rng) else t * 0
    } else {
        t_for_fit <- t
    }

    t_grid <- seq(min(t_for_fit, na.rm = TRUE),
                  max(t_for_fit, na.rm = TRUE),
                  length.out = grid_n)

    res <- apply(X, 1, function(z_vec) {
        z <- as.numeric(z_vec)
        df_fit <- data.frame(z = z, t = t_for_fit)

        fit <- try(
            mgcv::gam(z ~ s(t, k = k),
                      data = df_fit,
                      family = family,
                      na.action = na_action),
            silent = TRUE
        )
        if (inherits(fit, "try-error")) {
            return(c(pval = NA, dev_expl = NA, edf = NA, F_stat = NA,
                     range = NA, monotonicity = NA, t_peak = NA, max_abs_deriv = NA))
        }

        smry <- summary(fit)
        s_tab <- smry$s.table
        pval <- if (!is.null(s_tab) && nrow(s_tab) >= 1) s_tab[1, "p-value"] else NA
        edf <- if (!is.null(s_tab) && nrow(s_tab) >= 1) s_tab[1, "edf"] else NA
        F_stat <- if (!is.null(s_tab) && nrow(s_tab) >= 1) {
            if ("F" %in% colnames(s_tab)) s_tab[1, "F"]
            else if ("Chi.sq" %in% colnames(s_tab)) s_tab[1, "Chi.sq"]
            else NA
        } else NA
        dev_expl <- suppressWarnings(as.numeric(smry$dev.expl))

        pred <- try(
            predict(fit, newdata = data.frame(t = t_grid), type = "response"),
            silent = TRUE
        )
        if (inherits(pred, "try-error")) {
            return(c(pval = pval, dev_expl = dev_expl, edf = edf, F_stat = F_stat,
                     range = NA, monotonicity = NA, t_peak = NA, max_abs_deriv = NA))
        }
        fhat <- as.numeric(pred)

        dyn_range <- diff(range(fhat, na.rm = TRUE))
        monotonicity <- suppressWarnings(abs(cor(fhat, t_grid, method = "spearman", use = "complete.obs")))
        i_peak <- which.max(fhat)
        t_peak <- if (length(i_peak) > 0) t_grid[i_peak] else NA
        if (length(fhat) >= 2) {
            deriv <- diff(fhat) / diff(t_grid)
            max_abs_deriv <- max(abs(deriv), na.rm = TRUE)
        } else {
            max_abs_deriv <- NA
        }

        c(pval = pval, dev_expl = dev_expl, edf = edf, F_stat = F_stat,
          range = dyn_range, monotonicity = monotonicity,
          t_peak = t_peak, max_abs_deriv = max_abs_deriv)
    })

    res <- as.data.frame(t(res))
    res$gene <- rownames(X)
    res$padj <- p.adjust(res$pval, method = "BH")

    res <- res[, c("gene", "pval", "padj",
                   "dev_expl", "edf", "F_stat",
                   "range", "monotonicity", "t_peak", "max_abs_deriv")]
    rownames(res) <- NULL
    return(res)
}




library(mgcv)
library(dplyr)

# ------------------------------------------------------------------------------
# Builid pseudotime x gene matrix
# ------------------------------------------------------------------------------
compute_temporal_matrix <- function(expr_mat,
                                    pseudotime_vec,
                                    genes,
                                    method = c("bin", "smooth"),
                                    n_bins = 50,
                                    k = 6,
                                    grid_n = 200,
                                    family = gaussian()) {

    method <- match.arg(method)
    stopifnot(ncol(expr_mat) == length(pseudotime_vec))

    # 1️⃣ Input filtering
    genes <- intersect(genes, rownames(expr_mat))
    if (length(genes) == 0) stop("None of the input genes found in expr_mat.")

    keep <- !is.na(pseudotime_vec)
    expr_sub <- expr_mat[genes, keep, drop = FALSE]
    pt <- as.numeric(pseudotime_vec[keep])

    # ------------------------------------------------------------------------------
    # 2️⃣ Binning and averaging
    # ------------------------------------------------------------------------------
    if (method == "bin") {
    
        # 按伪时间切成 n_bins 个区间
        pt_bins <- cut(pt, breaks = n_bins, labels = FALSE)
    
        # 找出每个 bin 内细胞数量 > 0 的有效 bin
        bins_used <- which(tabulate(pt_bins) > 0)
    
        # 计算每个有效 bin 的中心位置（伪时间均值）
        bin_centers <- sapply(bins_used, function(b) mean(pt[pt_bins == b], na.rm = TRUE))
    
        # 计算每个有效 bin 中各基因的平均表达
        expr_bin <- sapply(bins_used, function(b) {
            cells <- which(pt_bins == b)
            rowMeans(expr_sub[genes, cells, drop = FALSE])
        })
    
        # 设置行名和列名
        rownames(expr_bin) <- genes
        colnames(expr_bin) <- sprintf("T%.2f", bin_centers)
    
        # 返回结果矩阵
        return(expr_bin)
    }

    # ------------------------------------------------------------------------------
    # 3️⃣ Smooth prediction (GAM)
    # ------------------------------------------------------------------------------
    if (method == "smooth") {
        t_grid <- seq(min(pt, na.rm = TRUE),
                      max(pt, na.rm = TRUE),
                      length.out = grid_n)

        expr_smooth <- sapply(genes, function(g) {
            z <- as.numeric(expr_sub[g, ])
            df <- data.frame(z = z, t = pt)
            fit <- try(
                mgcv::gam(z ~ s(t, k = k),
                          data = df,
                          family = family),
                silent = TRUE
            )
            if (inherits(fit, "try-error"))
                return(rep(NA, length(t_grid)))
            predict(fit, newdata = data.frame(t = t_grid), type = "response")
        })

        expr_smooth <- t(expr_smooth)
        rownames(expr_smooth) <- genes
        colnames(expr_smooth) <- sprintf("T%.2f", t_grid)
        return(expr_smooth)
    }
}

