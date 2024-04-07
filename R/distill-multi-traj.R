#' Distills multiple trajectory models
#'
#' This function takes a list of trajectory models and performs a distillation process to create a single set of motifs.
#'
#' @param traj_models A list of trajectory models.
#' @param max_motif_num The maximum number of motifs to be identified. If NULL - the number would be set to the maximum number of motifs in the input models.
#' @param distill_single Logical indicating whether to distill clusters with a single motif.
#' @param use_all_motifs Logical indicating whether to use all motifs in the resulting models. If FALSE, only motifs from clusters which had a motif from the original model are used.
#' @param seed The random seed for reproducibility. Defaults to 60427.
#' @param cluster_report_dir The directory to store cluster reports. If not NULL, a png would be created for each cluster.
#'
#' @return A distilled trajectory model.
#'
#' @inheritParams distill_traj_model
#' @inheritParams filter_traj_model
#' @export
distill_traj_model_multi <- function(traj_models, max_motif_num = NULL, min_diff = 0.1, intra_cor_thresh = 0.6, distill_single = FALSE, use_all_motifs = FALSE, bits_threshold = 1.75, r2_threshold = NULL, seed = 60427, parallel = TRUE, cluster_report_dir = NULL) {
    if (is.null(max_motif_num)) {
        max_motif_num <- max(purrr::map_dbl(traj_models, ~ length(.x@motif_models)))
        cli_alert_info("Setting {.field max_motif_num} to {.val {max_motif_num}}")
    }

    # make sure that all models have the same number of peaks
    if (length(unique(sapply(traj_models, function(x) nrow(x@peak_intervals)))) > 1) {
        cli_abort("All trajectory models must have the same number of peaks.")
    }
    orig_traj_model_names <- names(traj_models)
    names(traj_models) <- paste0("m", 1:length(traj_models))

    # unify energies
    features <- purrr::imap(traj_models, ~ {
        e <- .x@normalized_energies
        e <- e[, colnames(e) %in% names(.x@motif_models)]
        colnames(e) <- paste0(.y, ".", colnames(e))
        e
    }) %>% do.call(cbind, .)

    # unify models
    motif_models <- purrr::imap(traj_models, ~ {
        .x@motif_models
    }) %>% do.call(c, .)

    # unify atac_diff
    atac_diff <- purrr::imap_dfc(traj_models, ~ {
        tibble(!!.y := .x@diff_score)
    }) %>% as.matrix()
    diff_filter <- rowSums(abs(atac_diff) > min_diff) > 0
    cli_alert_info("Filtering out peaks with absolute ATAC-seq difference less than {.val {min_diff}} in at least on of the trajectories. Number of peaks left: {.val {sum(diff_filter)}}")

    # cluster models
    cm <- tgs_cor(features, pairwise.complete.obs = TRUE)
    hc <- hclust(as.dist(1 - cm), method = "complete")

    nclust <- min(ncol(features), max_motif_num)

    clust_map <- data.frame(feat = colnames(features), clust = paste0("c", cutree(hc, k = nclust))) %>%
        mutate(model = stringr::str_extract(feat, "^m\\d+")) %>%
        mutate(motif = gsub("^m\\d+\\.", "", feat))

    # correlate each motif with diff_score of its corresponding model
    motif_cors <- purrr::imap_dfr(traj_models, ~ {
        tgs_cor(.x@normalized_energies, as.matrix(.x@diff_score), pairwise.complete.obs = TRUE) %>%
            enframe("motif", "cor") %>%
            mutate(model = .y)
    })

    clust_map <- clust_map %>%
        left_join(motif_cors, by = c("model", "motif")) %>%
        arrange(clust, desc(abs(cor)))


    if (!is.null(intra_cor_thresh)) {
        avg_intra_cluster_cor <- function(cluster_id, clust_map, corr_matrix) {
            features_in_cluster <- clust_map$feat[clust_map$clust == cluster_id]
            cluster_corr <- corr_matrix[features_in_cluster, features_in_cluster]
            mean(cluster_corr[upper.tri(cluster_corr)], na.rm = TRUE)
        }

        clust_map <- clust_map %>%
            group_by(clust) %>%
            mutate(intra_cor = avg_intra_cluster_cor(clust[1], clust_map, cm)) %>%
            ungroup() %>%
            mutate(intra_cor = ifelse(is.na(intra_cor), 1, intra_cor)) %>%
            add_count(clust)

        to_split <- clust_map %>%
            filter(intra_cor < intra_cor_thresh, n > 1) %>%
            pull(clust) %>%
            unique()

        if (length(to_split) > 0) {
            cli_alert_info("Splitting {.val {length(to_split)}} cluster{?s} with average intra-cluster correlation < {.val {intra_cor_thresh}}")
            clust_map <- clust_map %>%
                group_by(clust) %>%
                mutate(clust_i = ifelse(clust %in% to_split, 1:n(), 1)) %>%
                ungroup() %>%
                tidyr::unite(clust, clust, clust_i, sep = "_") %>%
                mutate(clust = as.integer(as.factor(clust)))
        }
    }

    clust_map <- clust_map %>%
        group_by(clust) %>%
        mutate(clust_name = feat[1]) %>%
        ungroup()

    clust_to_name <- clust_map %>%
        distinct(clust, clust_name) %>%
        deframe()

    sequences <- toupper(misha::gseq.extract(misha.ext::gintervals.normalize(traj_models[[1]]@peak_intervals, traj_models[[1]]@params$peaks_size)))

    if (!is.null(cluster_report_dir)) {
        plot_traj_model_multi_clust(traj_models, clust_map, motif_models, cluster_report_dir)
    }

    prego_distilled <- plyr::dlply(clust_map, "clust", function(x) {
        n_feats <- nrow(x)
        clust_name <- clust_to_name[x$clust[1]]

        if (!distill_single && n_feats == 1) {
            cli_alert("Only one feature in cluster {.val {clust_name}}. Skipping distillation.")
            return(motif_models[[x$feat[1]]])
        }
        cli_alert_info("Running {.field prego} on cluster {.val {clust_name}}, distilling {.val {n_feats}} features")

        clust_models <- names(traj_models)[names(traj_models) %in% x$model]
        partial_y <- purrr::map_dfc(clust_models, ~ {
            model <- traj_models[[.x]]
            feats <- feat_to_variable(model) %>%
                filter(variable %in% x$motif[x$model == .x]) %>%
                pull(feature)

            py <- (model@model_features[, feats, drop = FALSE] %*% coef(model@model, s = model@params$lambda)[feats, , drop = FALSE])[, 1]
            tibble(!!.x := py)
        })

        partial_y <- as.matrix(partial_y)

        (prego_model <- prego::regress_pwm(
            sequences = sequences[diff_filter],
            response = partial_y[diff_filter, ],
            seed = seed,
            match_with_db = FALSE,
            screen_db = FALSE,
            multi_kmers = FALSE,
            spat_num_bins = traj_models[[1]]@params$spat_num_bins,
            spat_bin_size = traj_models[[1]]@params$spat_bin_size,
            kmer_sequence_length = traj_models[[1]]@params$kmer_sequence_length,
            symmetrize_spat = TRUE
        )
        ) %>%
            cli::cli_fmt()

        cli::cli_alert_success("Finished distilling cluster {.val {clust_name}}")

        return(prego::export_regression_model(prego_model))
    }, .parallel = TRUE)

    names(prego_distilled) <- clust_to_name[names(prego_distilled)]

    norm_sequences <- toupper(misha::gseq.extract(misha.ext::gintervals.normalize(traj_models[[1]]@normalization_intervals, traj_models[[1]]@params$peaks_size)))

    cli_alert_info("Infering energies...")
    clust_energies <- infer_energies(sequences, norm_sequences, prego_distilled, traj_models[[1]]@params$min_energy, traj_models[[1]]@params$energy_norm_quantile, traj_models[[1]]@params$norm_energy_max)

    traj_models_new <- purrr::imap(traj_models, ~ {
        nm <- .x
        if (!use_all_motifs) {
            models_to_use <- unique(clust_map$clust_name[clust_map$model == .y])
            clust_energies <- clust_energies[, colnames(clust_energies) %in% models_to_use]
        } else {
            models_to_use <- colnames(clust_energies)
        }
        cli_alert_info("Computing new trajectory model {.val {.y}}. Using {.val {length(models_to_use)}} motifs")
        clust_energies_logist <- create_logist_features(cbind(clust_energies, nm@additional_features))
        atac_diff <- nm@diff_score
        atac_diff_n <- norm01(atac_diff)
        model <- glmnet::glmnet(clust_energies_logist, atac_diff_n, binomial(link = "logit"), alpha = nm@params$alpha, lambda = nm@params$lambda, parallel = parallel, seed = nm@params$seed)

        predicted_diff_score <- logist(glmnet::predict.glmnet(model, newx = clust_energies_logist, type = "link", s = nm@params$lambda))[, 1]
        predicted_diff_score <- norm01(predicted_diff_score)
        predicted_diff_score <- rescale(predicted_diff_score, atac_diff)

        TrajectoryModel(
            model = model,
            motif_models = homogenize_pssm_models(prego_distilled[names(prego_distilled) %in% models_to_use]),
            coefs = get_model_coefs(model),
            normalized_energies = as.matrix(clust_energies),
            model_features = clust_energies_logist,
            type = nm@type,
            normalization_intervals = nm@normalization_intervals,
            additional_features = nm@additional_features,
            diff_score = atac_diff,
            predicted_diff_score = predicted_diff_score,
            initial_prego_models = nm@initial_prego_models,
            peak_intervals = nm@peak_intervals,
            params = nm@params
        )
    })

    names(traj_models_new) <- orig_traj_model_names

    cli_alert_info("Filtering models...")
    traj_models_new <- purrr::imap(traj_models_new, ~ {
        cli_alert("Filtering model {.val {.y}}")
        filter_traj_model(.x, r2_threshold = r2_threshold, bits_threshold = bits_threshold)
    })

    all_model_names <- purrr::map(traj_models_new, ~ names(.x@motif_models)) %>%
        unlist() %>%
        unique()
    prego_distilled <- prego_distilled[names(prego_distilled) %in% all_model_names]

    traj_models_new <- purrr::map(traj_models_new, add_traj_model_stats)
    traj_models <- purrr::map(traj_models, add_traj_model_stats)

    stats <- compute_traj_list_stats(traj_models_new)

    return(
        TrajectoryModelMulti(
            models = traj_models_new,
            motif_models = prego_distilled,
            cluster_map = clust_map,
            stats = stats,
            params = list(
                max_motif_num = max_motif_num,
                min_diff = min_diff,
                intra_cor_thresh = intra_cor_thresh,
                distill_single = distill_single,
                use_all_motifs = use_all_motifs,
                seed = seed
            )
        )
    )
}

#' Infer trajectory motifs for a multi-trajectory model
#'
#' This function infers trajectory motifs for each individual trajectory in a multi-trajectory model.
#'
#' @param traj_multi A TrajectoryModelMulti object.
#'
#' @return A TrajectoryModelMulti object with the result of \code{infer_trajectory_motifs} for each individual trajectory.
#'
#'
#' @inheritParams infer_trajectory_motifs
#' @export
infer_trajectory_motifs_multi <- function(traj_multi, peak_intervals, atac_scores = NULL, bin_start = 1, bin_end = ncol(atac_scores), additional_features = NULL) {
    validate_traj_model_multi(traj_multi)

    if (!is.list(atac_scores)) {
        cli_abort("The {.field atac_scores} must be a list of data frames")
    }

    if (length(atac_scores) != length(traj_multi@models)) {
        cli_abort("Number of elements in the {.field atac_scores} list must be equal to the number of models in the TrajectoryModelMulti object. It is {.val {length(atac_scores)}} while the number of models is {.val {length(traj_multi@models)}}")
    }

    traj_models <- traj_multi@models
    new_models <- purrr::imap(traj_models, ~ {
        infer_trajectory_motifs(.x, peak_intervals, atac_scores = atac_scores[[.y]], bin_start = bin_start, bin_end = bin_end, additional_features = additional_features)
    })
    names(new_models) <- names(traj_models)

    traj_multi@models <- new_models
    traj_multi@models <- purrr::map(traj_multi@models, add_traj_model_stats)
    traj_multi@stats <- compute_traj_list_stats(traj_multi)
    return(traj_multi)
}

compute_traj_list_stats <- function(traj_models) {
    r2_df <- purrr::imap_dfr(traj_models, ~ {
        tibble(
            model = .y,
            r2_train = .x@params$stats$r2_train,
            r2_test = .x@params$stats$r2_test,
            r2_all = .x@params$stats$r2_all
        )
    })

    return(r2_df)
}




plot_traj_model_multi_clust <- function(traj_models, clust_map, prego_distilled, out_dir) {
    cli_alert_info("Plotting cluster motifs to {.val {out_dir}}")
    pssm_db <- map_dfr(names(traj_models), function(m) {
        imap_dfr(traj_models[[m]]@motif_models, ~ .x$pssm %>% mutate(motif = paste0(m, ".", .y)))
    })

    dir.create(out_dir)

    for (i in unique(clust_map$clust)) {
        x <- clust_map %>% filter(clust == i)
        if (nrow(x) > 1) {
            cli::cli_alert("Plotting cluster {.val {x$clust_name[1]}}, intra_cor = {.val {x$intra_cor[1]}}")
            png(file.path(out_dir, paste0(x$clust_name[1], ".png")), width = 500, height = 150 * nrow(x))
            p <- clust_map %>%
                filter(clust == i) %>%
                as.data.frame() %>%
                pull(feat) %>%
                map(~ prego::plot_pssm_logo_dataset(.x, dataset = pssm_db))

            model <- prego_distilled[[x$clust_name[1]]]
            p_all <- prego::plot_pssm_logo(model$pssm) + ggtitle("Distilled")
            p <- c(list(p_all), p)

            p <- p %>%
                patchwork::wrap_plots(ncol = 1)
            p <- p + patchwork::plot_annotation(
                title = glue("Cluster {x$clust_name[1]}, cor = {x$intra_cor[1]}"),
                theme = theme(plot.title = element_text(size = 30))
            )
            print(p)
            dev.off()
        }
    }
}
