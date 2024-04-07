plot_multi_traj_model_report <- function(multi_traj, filename = NULL, width = NULL, height = NULL, dev = grDevices::pdf, title = NULL, ...) {
    models <- multi_traj@motif_models
    motif_num <- length(models)
    spatial_p <- purrr::imap(models, ~ {
        if (is.null(.x$spat)) {
            return(plot_fake_spat())
        }
        prego::plot_spat_model(.x$spat)
    })

    motifs_p <- purrr::imap(models, ~ prego::plot_pssm_logo(.x$pssm, title = .y))

    partial_resp_list <- purrr::map(multi_traj@models, compute_partial_response)

    pr_p <- purrr::map(names(models), function(motif) {
        purrr::map(names(multi_traj@models), function(model) {
            pr <- partial_resp_list[[model]]
            if (!(motif %in% names(pr))) {
                return(plot_fake_pr() + ggtitle(model))
            }
            plot_e_vs_pr(motif, pr, multi_traj@models[[model]]) + ggtitle(model)
        }) %>%
            purrr::discard(is.null) %>%
            patchwork::wrap_plots(nrow = 1)
    })

    n_traj_models <- length(multi_traj@models)

    p <- patchwork::wrap_plots(
        A = patchwork::wrap_plots(motifs_p, ncol = 1),
        B = patchwork::wrap_plots(spatial_p, ncol = 1),
        C = patchwork::wrap_plots(pr_p, ncol = 1),
        design = "ABC",
        widths = c(3, 0.8, n_traj_models)
    )

    if (!is.null(title)) {
        p <- p + patchwork::plot_annotation(title = title)
    }

    if (!is.null(filename)) {
        if (is.null(width)) {
            width <- 3 + 3 * n_traj_models
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

plot_fake_spat <- function() {
    prego::plot_spat_model(tibble(bin = seq(0, 400, 20), spat_factor = 1 / 21))
}

plot_fake_pr <- function() {
    ggplot() +
        labs(x = "Energy", y = "Partial response") +
        theme_classic() +
        theme(aspect.ratio = 1) +
        xlim(0, 10) +
        ylim(0, 1)
}
