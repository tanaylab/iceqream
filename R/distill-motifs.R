#' Distill trajectory model
#'
#' @description This function takes a trajectory model and reduces the number of motifs using 'prego'.
#'
#' @param traj_model A trajectory model object
#' @param max_motif_num The maximum number of motifs to select
#' @param min_diff The minimum difference of scores to use for distillation
#' @param intra_cor_thresh The threshold for the average intra-cluster correlation to split clusters. If NULL, no splitting is done.
#' @param parallel Whether to use parallel processing
#' @param use_non_linear Whether to use non-linear features. If FALSE, a glm model without non-linear features is used.
#'
#' @return A distilled trajectory model object with \code{max_motif_num} models + additional features.
#'
#' @export
distill_traj_model <- function(traj_model, max_motif_num, min_diff = 0.1, intra_cor_thresh = 0.6, use_non_linear = TRUE, parallel = TRUE) {
    if (traj_model_has_test(traj_model)) {
        cli_alert_info("Using only the training set for distillation.")
        traj_model_list <- split_traj_model_to_train_test(traj_model)
        traj_model_train <- traj_model_list$train
        traj_model_test <- traj_model_list$test
        cli_alert("# of train peaks: {.val {nrow(traj_model_train@peak_intervals)}}, # of test peaks: {.val {nrow(traj_model_test@peak_intervals)}}")

        traj_model_d <- distill_traj_model(traj_model = traj_model_train, max_motif_num = max_motif_num, min_diff = min_diff, intra_cor_thresh = intra_cor_thresh, use_non_linear = use_non_linear, parallel = parallel)
        cli_alert_info("Infering test peaks...")
        traj_model_test <- infer_trajectory_motifs(traj_model_d, peak_intervals = traj_model_test@peak_intervals, additional_features = traj_model_test@additional_features, diff_score = traj_model_test@diff_score)

        traj_model <- add_traj_model_stats(traj_model)

        cli_alert_success("R^2 train: {.val {traj_model_test@params$stats$r2_train}} (before: {.val {traj_model@params$stats$r2_train}}), R^2 test: {.val {traj_model_test@params$stats$r2_test}} (before: {.val {traj_model@params$stats$r2_test}}). # of motifs: {.val {length(traj_model_test@motif_models)}} (before: {.val {length(traj_model@motif_models)}})")
        return(traj_model_test)
    }

    pssm_db <- purrr::imap_dfr(traj_model@motif_models, ~ .x$pssm %>% mutate(motif = .y)) %>%
        select(motif, pos, everything())

    params <- traj_model@params
    withr::local_options(list(gmax.data.size = 1e9))
    seqs <- toupper(misha::gseq.extract(misha.ext::gintervals.normalize(traj_model@peak_intervals, params$peaks_size)))
    norm_seqs <- toupper(misha::gseq.extract(misha.ext::gintervals.normalize(traj_model@normalization_intervals, params$peaks_size)))
    atac_diff <- traj_model@diff_score
    atac_diff_n <- norm01(atac_diff)

    diff_filter <- abs(atac_diff) >= min_diff

    if (use_non_linear) {
        glm_model <- traj_model@model
    } else {
        glm_model <- glmnet::glmnet(as.matrix(cbind(traj_model@normalized_energies, traj_model@additional_features)), atac_diff_n, binomial(link = "logit"), alpha = params$alpha, lambda = params$lambda, parallel = parallel, seed = params$seed)
        glm_model <- strip_glmnet(glm_model)
    }


    distilled <- distill_motifs(traj_model@normalized_energies, max_motif_num, glm_model, y = atac_diff_n, diff_filter = diff_filter, seqs = seqs, norm_seqs = norm_seqs, additional_features = traj_model@additional_features, pssm_db = pssm_db, prev_models = traj_model@motif_models, lambda = params$lambda, alpha = params$alpha, energy_norm_quantile = params$energy_norm_quantile, seed = params$seed, spat_num_bins = params$spat_num_bins, spat_bin_size = params$spat_bin_size, kmer_sequence_length = params$kmer_sequence_length, nclust = max_motif_num, n_clust_factor = params$n_clust_factor, distill_single = FALSE, intra_cor_thresh = intra_cor_thresh, use_non_linear = use_non_linear, symmetrize_spat = traj_model@params$symmetrize_spat %||% TRUE)

    clust_energies <- distilled$energies
    clust_energies_logist <- create_logist_features(clust_energies)

    model <- glmnet::glmnet(clust_energies_logist, atac_diff_n, binomial(link = "logit"), alpha = params$alpha, lambda = params$lambda, parallel = parallel, seed = params$seed)
    model <- strip_glmnet(model)

    predicted_diff_score <- logist(glmnet::predict.glmnet(model, newx = clust_energies_logist, type = "link", s = params$lambda))[, 1]
    predicted_diff_score <- norm01(predicted_diff_score)
    predicted_diff_score <- rescale(predicted_diff_score, atac_diff)


    cli_alert_success("Finished running model. Number of non-zero coefficients: {.val {sum(model$beta != 0)}} (out of {.val {ncol(clust_energies_logist)}}). R^2: {.val {cor(predicted_diff_score, atac_diff_n)^2}}")
    params <- traj_model@params
    params$distilled_features <- distilled$features

    # remove additional features from clust_energies
    normalized_energies <- clust_energies[, setdiff(colnames(clust_energies), names(traj_model@additional_features)), drop = FALSE]

    traj_model_distilled <- TrajectoryModel(
        model = model,
        motif_models = homogenize_pssm_models(distilled$motifs),
        coefs = get_model_coefs(model),
        normalized_energies = as.matrix(normalized_energies),
        model_features = clust_energies_logist,
        type = traj_model@type,
        normalization_intervals = traj_model@normalization_intervals,
        additional_features = traj_model@additional_features,
        diff_score = atac_diff,
        predicted_diff_score = predicted_diff_score,
        initial_prego_models = traj_model@initial_prego_models,
        peak_intervals = traj_model@peak_intervals,
        params = params
    )

    return(traj_model_distilled)
}

distill_motifs <- function(features, target_number, glm_model, y, seqs, norm_seqs, diff_filter, additional_features = NULL, pssm_db = prego::all_motif_datasets(), prev_models = list(), lambda = 1e-5, alpha = 1, energy_norm_quantile = 1, norm_energy_max = 10, min_energy = -7, seed = 60427, spat_num_bins = NULL, spat_bin_size = NULL, kmer_sequence_length = NULL, nclust = NULL, n_clust_factor = 1, distill_single = TRUE, intra_cor_thresh = NULL, use_non_linear = FALSE, symmetrize_spat = TRUE) {
    if (is.null(nclust)) {
        nclust <- min(ncol(features), target_number * n_clust_factor)
    }

    cli_alert_info("Clustering {.val {ncol(features)}} features into {.val {nclust}} clusters...")
    features_cm <- tgs_cor(features, pairwise.complete.obs = TRUE)

    hc <- hclust(as.dist(1 - features_cm), method = "complete")
    clust_map <- data.frame(feat = colnames(features), clust = cutree(hc, k = nclust))
    if (use_non_linear) {
        clust_map$beta <- apply(clust_map, 1, function(x) max(abs(glm_model$beta[variable_to_feat(glm_model, x[1]), 1]), na.rm = TRUE))
    } else {
        clust_map$beta <- glm_model$beta[colnames(features), 1]
    }


    if (!is.null(intra_cor_thresh)) {
        avg_intra_cluster_cor <- function(cluster_id, clust_map, corr_matrix) {
            features_in_cluster <- clust_map$feat[clust_map$clust == cluster_id]
            cluster_corr <- corr_matrix[features_in_cluster, features_in_cluster]
            mean(cluster_corr[upper.tri(cluster_corr)], na.rm = TRUE)
        }

        clust_map <- clust_map %>%
            group_by(clust) %>%
            mutate(intra_cor = avg_intra_cluster_cor(clust[1], clust_map, features_cm)) %>%
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
            nclust <- length(unique(clust_map$clust))
            if (target_number < nclust) {
                target_number <- nclust
            }
        }
    }

    cli_alert_info("Choosing top {.val {target_number}} features clusters")
    selected_clust <- clust_map %>%
        group_by(clust) %>%
        summarize(max_beta = max(beta)) %>%
        arrange(desc(max_beta)) %>%
        mutate(new_clust = 1:n()) %>%
        head(target_number)

    clust_map <- clust_map %>%
        inner_join(selected_clust, by = "clust") %>%
        mutate(clust = new_clust) %>%
        select(-new_clust)

    cli_alert_info("Features left: {.val {nrow(clust_map)}}")

    best_clust_map <- clust_map %>%
        arrange(clust, abs(beta)) %>%
        group_by(clust) %>%
        slice(1) %>%
        ungroup() %>%
        as.data.frame()

    cli_alert_info("Learning a model for each motif cluster...")

    best_motifs_prego <- plyr::alply(best_clust_map, 1, function(x) {
        n_feats <- nrow(clust_map %>% filter(clust == x$clust))
        motif <- NULL
        optimize_pwm <- TRUE
        if (!distill_single &&
            n_feats == 1) {
            if (!is.null(prev_models[[x$feat]])) {
                cli_alert("Cluster {.val {x$feat}} has only one feature, skipping...")
                return(prev_models[[x$feat]])
            } else if (nrow(pssm_db %>% filter(motif == x$feat)) > 0) {
                cli_alert_info("Cluster {.val {x$feat}} has only one feature, learning only spatial model...")
                motif <- pssm_db %>%
                    filter(motif == x$feat) %>%
                    select(pos, A, C, G, T)
                optimize_pwm <- FALSE
            } else {
                cli::cli_alert_warning("No current model found for {.val {x$feat}}.")
            }
        } else {
            cli_alert_info("Running {.field prego} on cluster {.val {x$feat}} (fusing {.val {n_feats}} motifs)")
        }

        res <- run_prego_on_clust_residuals(
            model = glm_model,
            feats = features[diff_filter, ],
            clust_motifs = clust_map %>% filter(clust == x$clust) %>% pull(feat),
            sequences = seqs[diff_filter],
            lambda = lambda,
            seed = seed,
            spat_num_bins = spat_num_bins,
            spat_bin_size = spat_bin_size,
            kmer_sequence_length = kmer_sequence_length,
            use_non_linear = use_non_linear,
            motif = motif,
            optimize_pwm = optimize_pwm,
            symmetrize_spat = symmetrize_spat
        )
        cli::cli_alert_success("Finished running {.field prego} on cluster {.val {x$feat}}")
        return(res)
    }, .parallel = getOption("prego.parallel", TRUE))
    names(best_motifs_prego) <- best_clust_map$feat

    cli_alert_info("Infering energies...")
    clust_energies <- infer_energies_new(seqs, norm_seqs, best_motifs_prego, min_energy, energy_norm_quantile, norm_energy_max)

    # add missing features
    missing_features <- setdiff(best_clust_map$feat, colnames(clust_energies))
    if (length(missing_features) > 0) {
        clust_energies <- cbind(clust_energies, features[, missing_features, drop = FALSE])
    }
    clust_energies <- clust_energies[, best_clust_map$feat, drop = FALSE]

    if (!is.null(additional_features) && length(additional_features) > 0 && nrow(additional_features) > 0) {
        additional_features[is.na(additional_features)] <- 0
        clust_energies <- cbind(clust_energies, additional_features)
    }

    distilled_features <- clust_map %>%
        left_join(best_clust_map %>% select(clust, name = feat), by = "clust") %>%
        select(distilled = name, model = feat, beta, max_beta, any_of("intra_cor")) %>%
        arrange(distilled)

    return(list(energies = clust_energies, motifs = best_motifs_prego, features = distilled_features))
}

run_prego_on_clust_residuals <- function(model, feats, clust_motifs, sequences, lambda = 1e-5, seed = 60427, spat_num_bins = NULL, spat_bin_size = NULL, kmer_sequence_length = NULL, use_non_linear = FALSE, motif = NULL, optimize_pwm = TRUE, symmetrize_spat = TRUE) {
    if (use_non_linear) {
        feats <- create_logist_features(feats[, clust_motifs, drop = FALSE])
        partial_y <- (feats %*% coef(model, s = lambda)[colnames(feats), , drop = FALSE])[, 1]
    } else {
        partial_y <- (feats[, clust_motifs, drop = FALSE] %*% coef(model, s = lambda)[clust_motifs, , drop = FALSE])[, 1]
    }

    cli::cli_fmt(prego_model <- prego::regress_pwm(sequences = sequences, response = partial_y, seed = seed, match_with_db = FALSE, screen_db = FALSE, multi_kmers = FALSE, spat_num_bins = spat_num_bins, spat_bin_size = spat_bin_size, kmer_sequence_length = kmer_sequence_length, symmetrize_spat = symmetrize_spat, motif = motif, optimize_pwm = optimize_pwm))

    return(prego::export_regression_model(prego_model))
}
