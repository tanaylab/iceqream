distill_motifs <- function(features, target_number, glm_model, y, seqs, additional_features = NULL, pssm_db = prego::all_motif_datasets(), prego_models = list(), lambda = 1e-5, alpha = 1, energy_norm_quantile = 1, seed = 60427, spat_num_bins = NULL, spat_bin_size = NULL, kmer_sequence_length = NULL) {
    nclust <- min(ncol(features), target_number * 2)
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
            feats = features,
            clust_motifs = clust_map %>% filter(clust == x$clust) %>% pull(feat),
            sequences = seqs,
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
    clust_energies <- apply(clust_energies_raw, 2, norm_energy, min_energy = -10, q = energy_norm_quantile)

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

    if (nrow(pssm) == 0) { # Motif is part of additional features
        return(NULL)
    }

    partial_y <- (feats[, clust_motifs, drop = FALSE] %*% coef(model, s = lambda)[clust_motifs, , drop = FALSE])[, 1]

    sequences_directed <- direct_sequences(sequences, pssm)

    cli::cli_fmt(prego_model <- prego::regress_pwm(sequences = sequences_directed, response = partial_y, seed = seed, match_with_db = FALSE, screen_db = FALSE, multi_kmers = FALSE, spat_num_bins = spat_num_bins, spat_bin_size = spat_bin_size, kmer_sequence_length = kmer_sequence_length, symmetrize_spat = FALSE))

    cli::cli_alert_success("Finished running {.field prego} on cluster {.val {motif}}")
    return(prego::export_regression_model(prego_model))
}
