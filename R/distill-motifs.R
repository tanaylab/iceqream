#' Distill trajectory model
#'
#' @description This function takes a trajectory model and reduces the number of motifs using 'prego'.
#'
#' @param traj_model A trajectory model object
#' @param max_motif_num The maximum number of motifs to select
#' @param min_diff The minimum difference of scores to use for distillation
#' @param parallel Whether to use parallel processing
#'
#' @return A distilled trajectory model object with \code{max_motif_num} models + additional features.
#'
#' @export
distill_traj_model <- function(traj_model, max_motif_num, min_diff = 0.1, parallel = TRUE) {
    if (any(traj_model@type == "test")) {
        cli_abort("Cannot distill a trajectory model with test peaks.")
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

    glm_linear <- glmnet::glmnet(as.matrix(cbind(traj_model@normalized_energies, traj_model@additional_features)), atac_diff_n, binomial(link = "logit"), alpha = params$alpha, lambda = params$lambda, parallel = parallel, seed = params$seed)

    distilled <- distill_motifs(traj_model@normalized_energies, max_motif_num, glm_linear, y = atac_diff_n, diff_filter = diff_filter, seqs = seqs, norm_seqs = norm_seqs, additional_features = traj_model@additional_features, pssm_db = pssm_db, prego_models = traj_model@motif_models, lambda = params$lambda, alpha = params$alpha, energy_norm_quantile = params$energy_norm_quantile, seed = params$seed, spat_num_bins = params$spat_num_bins, spat_bin_size = params$spat_bin_size, kmer_sequence_length = params$kmer_sequence_length, nclust = max_motif_num, n_clust_factor = params$n_clust_factor)

    clust_energies <- distilled$energies
    clust_energies_logist <- create_logist_features(clust_energies)

    model <- glmnet::glmnet(clust_energies_logist, atac_diff_n, binomial(link = "logit"), alpha = params$alpha, lambda = params$lambda, parallel = parallel, seed = params$seed)

    predicted_diff_score <- logist(glmnet::predict.glmnet(model, newx = clust_energies_logist, type = "link", s = params$lambda))[, 1]
    predicted_diff_score <- (predicted_diff_score * max(atac_diff)) + min(atac_diff)


    cli_alert_success("Finished running model. Number of non-zero coefficients: {.val {sum(model$beta != 0)}} (out of {.val {ncol(clust_energies_logist)}}). R^2: {.val {cor(predicted_diff_score, atac_diff_n)^2}}")

    traj_model_distilled <- TrajectoryModel(
        model = model,
        motif_models = homogenize_pssm_models(distilled$motifs),
        coefs = get_model_coefs(model),
        normalized_energies = as.matrix(clust_energies),
        model_features = clust_energies_logist,
        type = traj_model@type,
        additional_features = traj_model@additional_features,
        diff_score = atac_diff,
        predicted_diff_score = predicted_diff_score,
        initial_prego_models = traj_model@initial_prego_models,
        peak_intervals = traj_model@peak_intervals,
        params = traj_model@params
    )

    return(traj_model_distilled)
}

distill_motifs <- function(features, target_number, glm_model, y, seqs, norm_seqs, diff_filter, additional_features = NULL, pssm_db = prego::all_motif_datasets(), prego_models = list(), lambda = 1e-5, alpha = 1, energy_norm_quantile = 1, norm_energy_max = 10, min_energy = -7, seed = 60427, spat_num_bins = NULL, spat_bin_size = NULL, kmer_sequence_length = NULL, nclust = NULL, n_clust_factor = 1) {
    if (is.null(nclust)) {
        nclust <- min(ncol(features), target_number * n_clust_factor)
    }

    cli_alert_info("Clustering {.val {ncol(features)}} features into {.val {nclust}} clusters...")
    features_cm <- tgs_cor(features, pairwise.complete.obs = TRUE)
    hc <- tgs_dist(features_cm) %>% hclust(method = "ward.D2")
    clust_map <- data.frame(feat = colnames(features), clust = cutree(hc, k = nclust), beta = glm_model$beta[colnames(features), 1])

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
    library(glmnet)

    best_motifs_prego <- plyr::alply(best_clust_map, 1, function(x) {
        run_prego_on_clust_residuals(
            x$feat,
            model = glm_model,
            y = y,
            feats = features[diff_filter, ],
            clust_motifs = clust_map %>% filter(clust == x$clust) %>% pull(feat),
            sequences = seqs[diff_filter],
            lambda = lambda,
            pssm_db = pssm_db,
            spat_db = prego_models %>% purrr::map("spat"),
            spat_min = prego_models %>% purrr::map("spat_min"),
            spat_max = prego_models %>% purrr::map("spat_max"),
            seed = seed,
            spat_num_bins = spat_num_bins,
            spat_bin_size = spat_bin_size,
            kmer_sequence_length = kmer_sequence_length
        )
    }, .parallel = TRUE)
    names(best_motifs_prego) <- best_clust_map$feat

    cli_alert_info("Infering energies...")
    clust_energies <- plyr::llply(purrr::discard(best_motifs_prego, is.null), function(x) {
        prego::compute_pwm(seqs, x$pssm, spat = x$spat, spat_min = x$spat_min %||% 1, spat_max = x$spat_max)
    }, .parallel = TRUE)
    names(clust_energies) <- best_clust_map$feat
    clust_energies_raw <- do.call(cbind, clust_energies)

    norm_clust_energies <- plyr::llply(purrr::discard(best_motifs_prego, is.null), function(x) {
        prego::compute_pwm(norm_seqs, x$pssm, spat = x$spat, spat_min = x$spat_min %||% 1, spat_max = x$spat_max)
    }, .parallel = TRUE)
    names(norm_clust_energies) <- best_clust_map$feat
    norm_clust_energies <- do.call(cbind, norm_clust_energies)

    clust_energies <- norm_energy_matrix(clust_energies_raw, norm_clust_energies, min_energy = min_energy, q = energy_norm_quantile, norm_energy_max = norm_energy_max)

    # add missing features
    missing_features <- setdiff(best_clust_map$feat, colnames(clust_energies))
    if (length(missing_features) > 0) {
        clust_energies <- cbind(clust_energies, features[, missing_features, drop = FALSE])
    }
    clust_energies <- clust_energies[, best_clust_map$feat, drop = FALSE]

    if (!is.null(additional_features)) {
        additional_features[is.na(additional_features)] <- 0
        clust_energies <- cbind(clust_energies, additional_features)
    }

    distilled_features <- clust_map %>%
        left_join(best_clust_map %>% select(clust, name = feat), by = "clust") %>%
        select(clust = name, feat, beta, max_beta) %>%
        arrange(clust)

    return(list(energies = clust_energies, motifs = best_motifs_prego, features = distilled_features))
}

run_prego_on_clust_residuals <- function(motif, model, y, feats, clust_motifs, sequences, pssm_db, spat_db = NULL, spat_min = NULL, spat_max = NULL, lambda = 1e-5, seed = 60427, spat_num_bins = NULL, spat_bin_size = NULL, kmer_sequence_length = NULL) {
    cli_alert("Running {.field prego} on cluster {.val {motif}}...")
    pssm <- pssm_db %>%
        filter(motif == !!motif)

    partial_y <- (feats[, clust_motifs, drop = FALSE] %*% coef(model, s = lambda)[clust_motifs, , drop = FALSE])[, 1]

    cli::cli_fmt(prego_model <- prego::regress_pwm(sequences = sequences, response = partial_y, seed = seed, match_with_db = FALSE, screen_db = FALSE, multi_kmers = FALSE, spat_num_bins = spat_num_bins, spat_bin_size = spat_bin_size, kmer_sequence_length = kmer_sequence_length, symmetrize_spat = TRUE))

    cli::cli_alert_success("Finished running {.field prego} on cluster {.val {motif}}")
    return(prego::export_regression_model(prego_model))
}
