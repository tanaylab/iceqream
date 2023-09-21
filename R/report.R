#' Plot scatter plot of observed vs predicted ATAC difference
#'
#' @param traj_model Trajectory model object. Please run \code{\link{regress_trajectory_motifs}} first.
#' @return ggplot2 object
#'
#'
#' @export
plot_prediction_scatter <- function(traj_model) {
    r <- get_obs_pred_df(traj_model)
    plot_df <- r$df
    type_str <- r$type_str

    limits <- c(min(plot_df$observed, plot_df$predicted), max(plot_df$observed, plot_df$predicted))

    p <- ggplot(plot_df, aes(x = observed, y = predicted)) +
        geom_point(alpha = 0.5) +
        geom_abline(intercept = 0, slope = 1, color = "red") +
        labs(x = "Observed ATAC difference", y = "Predicted ATAC difference") +
        xlim(limits) +
        ylim(limits) +
        ggtitle(paste0("R^2 = ", round(summary(lm(predicted ~ observed, data = plot_df))$r.squared, 2), type_str)) +
        theme_classic()

    return(p)
}

get_obs_pred_df <- function(traj_model) {
    type_str <- ""
    if (sum(traj_model@type == "test") == 0) {
        cli::cli_alert_warning("No test intervals found. Please run {.func infer_trajectory_motifs} first. Using train intervals instead.")
        plot_df <- data.frame(
            observed = traj_model@diff_score[traj_model@type == "train"],
            predicted = traj_model@predicted_diff_score[traj_model@type == "train"]
        )
        type_str <- " (train)"
    } else {
        plot_df <- data.frame(
            observed = traj_model@diff_score[traj_model@type == "test"],
            predicted = traj_model@predicted_diff_score[traj_model@type == "test"]
        )
    }

    return(list(df = plot_df, type_str = type_str))
}

#' Plot a boxplot of observed vs predicted ATAC differences
#'
#' This function takes a trajectory model and plots a boxplot of the observed vs predicted ATAC differences.
#' The plot is grouped by the observed ATAC difference, which is divided into n_groups quantiles.
#'
#' @param traj_model A trajectory model object
#' @param n_groups The number of groups to divide the observed ATAC difference into
#'
#' @return A ggplot object representing the boxplot
#'
#' @export
plot_prediction_boxplot <- function(traj_model, n_groups = 5) {
    r <- get_obs_pred_df(traj_model)
    plot_df <- r$df
    type_str <- r$type_str

    plot_df <- plot_df %>%
        mutate(obs_group = cut(observed, quantile(observed, seq(0, 1, length.out = n_groups + 1)), include.lowest = TRUE)) %>%
        gather(key = "type", value = "value", observed, predicted)

    p <- plot_df %>%
        ggplot(aes(x = obs_group, y = value, fill = type)) +
        geom_boxplot() +
        labs(x = "Observed ATAC difference", y = "ATAC difference") +
        scale_fill_manual(name = "", values = c(observed = "blue", predicted = "darkred")) +
        theme_classic() +
        theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))


    return(p)
}

#' Plot scatter of a variable vs ATAC difference
#'
#' @param traj_model Trajectory model object. Please run \code{\link{regress_trajectory_motifs}} first.
#' @param variable Variable to plot.
#' @param point_size Size of the points.
#'
#' @return ggplot2 object
#'
#' @export
plot_variable_vs_response <- function(traj_model, variable, point_size = 0.5) {
    validate_traj_model(traj_model)
    if (!(variable %in% traj_model_all@coefs$variable)) {
        cli_abort("Variable {.val {variable}} not found in the model.")
    }

    if (variable %in% colnames(traj_model@additional_features)) {
        plot_df <- data.frame(
            response = traj_model@diff_score,
            variable = traj_model@additional_features[[variable]]
        )
        p <- ggplot(plot_df, aes(x = variable, y = response)) +
            geom_point(alpha = 0.5) +
            labs(x = variable, y = "ATAC difference") +
            theme_classic()
    } else {
        plot_df <- purrr::map_dfr(c("early", "linear", "late"), ~ data.frame(
            response = traj_model@diff_score,
            variable = traj_model@normalized_energies[, paste0(variable, "_", .x)],
            type = .x
        ))

        p <- ggplot(aes(x = variable, y = response, color = type), data = plot_df) +
            geom_point(size = point_size, alpha = 0.5) +
            scale_color_manual(name = "", values = c(early = "red", linear = "green", late = "blue")) +
            labs(x = variable, y = "ATAC difference") +
            theme_classic()
    }

    return(p)
}

#' Plot a report of trajectory motifs
#'
#' @param traj_model Trajectory model object. Please run \code{regress_trajectory_motifs} first.
#' @param motif_num Number of motifs to plot. If NULL, all motifs will be plotted.
#' @param free_coef_axis Whether to use a free axis for the coefficient plots.
#' @param spatial_freqs Pre-computed spatial frequencies to plot. Use \code{compute_traj_model_spatial_freq} to compute.
#' @param filename Filename to save the plot to. If NULL, the plot will be returned.
#' @param width Width of the plot.
#' @param height Height of the plot.
#' @param dev Device to use for saving the plot. Default: \code{grDevices::pdf}.
#' @param ... Additional arguments to pass to the device.
#'
#' @return ggplot2 object. If filename is not NULL, the plot will be saved to the file and the function will return \code{invisible(NULL)}.
#'
#'
#' @export
plot_motifs_report <- function(traj_model, motif_num = NULL, free_coef_axis = TRUE, spatial_freqs = NULL, filename = NULL, width = NULL, height = NULL, dev = grDevices::pdf, title = NULL, ...) {
    validate_traj_model(traj_model)

    models <- traj_model@motif_models

    if (length(traj_model@features_r2) > 0) {
        sorted_vars <- names(sort(traj_model@features_r2, decreasing = TRUE))
    } else {
        sorted_vars <- traj_model@coefs %>%
            tibble::column_to_rownames("variable") %>%
            as.matrix() %>%
            abs() %>%
            apply(1, max) %>%
            sort(decreasing = TRUE) %>%
            names()
    }


    sorted_vars <- sorted_vars[!(sorted_vars %in% colnames(traj_model@additional_features))]

    if (!is.null(motif_num)) {
        if (motif_num > length(models)) {
            cli_abort("Motif number {.val {motif_num}} is larger than the number of motifs in the model ({.val {length(models)}})")
        }
    } else {
        motif_num <- length(sorted_vars)
    }

    cli_alert_info("Plotting {.val {motif_num}} motifs")
    models <- models[sorted_vars[1:motif_num]]

    if (length(traj_model@features_r2) > 0) {
        spatial_p <- purrr::imap(models, ~ prego::plot_spat_model(.x$spat, title = paste0("R^2=", round(traj_model@features_r2[.y], 6))))
    } else {
        spatial_p <- purrr::imap(models, ~ prego::plot_spat_model(.x$spat))
    }

    motifs_p <- purrr::imap(models, ~ prego::plot_pssm_logo(.x$pssm, title = .y))

    if (free_coef_axis) {
        coef_limits <- NULL
    } else {
        coef_limits <- c(min(as.matrix(traj_model@coefs[, c("early", "linear", "late")])), max(as.matrix(traj_model@coefs[, c("early", "linear", "late")])))
    }
    coefs_p <- purrr::map(names(models), ~ plot_coefs(traj_model, .x, limits = coef_limits, title = ""))

    if (is.null(spatial_freqs)) {
        spatial_freqs <- compute_traj_model_spatial_freq(traj_model, size = 1000, pwm_threshold = 7, top_q = 0.1, bottom_q = 0.1)
    }

    spat_freq_p <- purrr::map(names(models), ~ plot_motif_spatial_freq(spatial_freqs, .x, smooth = 10))

    if ("atac_freq" %in% colnames(spatial_freqs)) {
        atac_spat_freq_p <- purrr::map(names(models), ~ plot_motif_spatial_freq(spatial_freqs, .x, smooth = 10, plot_atac = TRUE))
    } else {
        atac_spat_freq_p <- purrr::map(names(models), ~ ggplot() +
            theme_void())
    }

    if ("k4me3" %in% colnames(spatial_freqs)) {
        k4me3_p <- purrr::map(names(models), ~ plot_epi_spatial_freq(spatial_freqs, .x, "k4me3"))
    } else {
        k4me3_p <- purrr::map(names(models), ~ ggplot() +
            theme_void())
    }

    if ("k27me3" %in% colnames(spatial_freqs)) {
        k27me3_p <- purrr::map(names(models), ~ plot_epi_spatial_freq(spatial_freqs, .x, "k27me3"))
    } else {
        k27me3_p <- purrr::map(names(models), ~ ggplot() +
            theme_void())
    }

    # scatter_p <- purrr::map(names(models), ~ plot_variable_vs_response(traj_model, .x, point_size = 0.001))

    p <- patchwork::wrap_plots(
        A = patchwork::wrap_plots(motifs_p, ncol = 1),
        B = patchwork::wrap_plots(coefs_p, ncol = 1),
        C = patchwork::wrap_plots(spatial_p, ncol = 1),
        D = patchwork::wrap_plots(spat_freq_p, ncol = 1),
        E = patchwork::wrap_plots(atac_spat_freq_p, ncol = 1),
        F = patchwork::wrap_plots(k4me3_p, ncol = 1),
        G = patchwork::wrap_plots(k27me3_p, ncol = 1),
        design = "ABCDEFG",
        widths = c(0.5, 0.1, 0.1, 0.3, 0.3, 0.3, 0.3)
    )

    if (!is.null(title)) {
        p <- p + patchwork::plot_annotation(title = title)
    }

    if (!is.null(filename)) {
        if (is.null(width)) {
            width <- 24
        }
        if (is.null(height)) {
            height <- motif_num * 1.8
        }
        cli_alert_info("Saving plot...")
        dev(filename, width = width, height = height, ...)
        print(p)
        dev.off()
        cli_alert_success("Plot saved to {.file {filename}}")
        invisible(p)
    } else {
        return(p)
    }
}

plot_coefs <- function(traj_model, variable, limits = NULL, title = variable) {
    validate_traj_model(traj_model)
    coef_df <- traj_model@coefs %>% filter(variable == !!variable)

    if (nrow(coef_df) == 0) {
        cli_abort("Coefficient for variable {.val {variable}} not found.")
    }

    coef_df <- coef_df %>%
        gather("type", "value", -variable) %>%
        mutate(type = factor(type, levels = c("late", "linear", "early")))

    p <- ggplot(coef_df, aes(x = type, y = value)) +
        geom_col() +
        labs(x = "", y = "Coefficient") +
        coord_flip() +
        theme_classic() +
        ggtitle(title)

    if (!is.null(limits)) {
        p <- p + ylim(limits)
    }

    return(p)
}

#' Plot trajectory model report
#'
#' @param traj_model A trajectory model object
#' @param dir A directory to save the report files
#' @param k Number of clusters to split the heatmap into
#' @param spatial_freqs A vector of spatial frequencies to plot. Use \code{compute_traj_model_spatial_freq} to compute.
#'
#' @return None.
#'
#' @export
plot_traj_model_report <- function(traj_model, dir, k = 10, spatial_freqs = NULL) {
    validate_traj_model(traj_model)
    e_mat <- traj_model@normalized_energies
    e_mat <- e_mat[, setdiff(colnames(e_mat), colnames(traj_model@additional_features))]
    cm <- tgs_cor(e_mat, pairwise.complete.obs = TRUE)
    cm_no_diag <- cm
    diag(cm_no_diag) <- 0
    hc <- tgs_dist(cm_no_diag) %>% hclust(method = "ward.D2")

    hm <- ComplexHeatmap::Heatmap(cm, name = "features", cluster_rows = hc, cluster_columns = hc, col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), split = k, column_split = k)

    if (dir.exists(dir)) {
        unlink(dir, recursive = TRUE)
    }

    dir.create(dir, showWarnings = FALSE, recursive = TRUE)

    png(file.path(dir, "heatmap.png"), width = 2000, height = 1000)
    if (length(traj_model@features_r2) > 0) {
        ha <- ComplexHeatmap::rowAnnotation(r2 = ComplexHeatmap::anno_points(traj_model@features_r2[rownames(cm)]), width = grid::unit(3, "cm"))
        hm <- ComplexHeatmap::draw(hm + ha, heatmap_legend_side = "left")
    } else {
        hm <- ComplexHeatmap::draw(hm, heatmap_legend_side = "left")
    }

    dev.off()

    clust_df <- ComplexHeatmap::row_order(hm) %>%
        purrr:::imap_dfr(~ tibble(ord = .x, clust = .y)) %>%
        mutate(motif = rownames(cm)[ord])

    plyr::d_ply(clust_df, "clust", function(x) {
        traj_model_clust <- traj_model
        traj_model_clust@motif_models <- traj_model@motif_models[x$motif]
        traj_model_clust@coefs <- traj_model@coefs %>% filter(variable %in% x$motif)
        traj_model_clust@features_r2 <- traj_model@features_r2[x$motif]
        plot_motifs_report(traj_model_clust, filename = file.path(dir, paste0("clust_", x$clust[1], ".pdf")), title = paste0("Cluster ", x$clust[1], " (", nrow(x), " motifs)"), spatial_freqs = spatial_freqs)
    })

    system(glue("ml load python3/3.7.5; /home/aviezerl/tools/pdftools/pdfmerge.py -o {dir}/all.pdf {files}", files = paste0(glue::glue("{dir}/clust_"), 1:k, ".pdf", collapse = " ")))

    invisible(clust_df)
}
