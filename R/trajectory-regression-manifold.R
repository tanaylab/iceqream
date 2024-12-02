#' Regress trajectory motifs of a full manifold (multiple trajectories)
#'
#' @param target_traj_motif_num Number of motifs to select for each trajectory. Note that the actual number of motifs can be less or more than this number, depending on the data from other trajectories.
#' @param initial_min_diff minimal ATAC difference for a peak to participate in the initial prego motif inference.
#'
#' @return a TrajectoryModelMulti object
#'
#' @inheritParams regress_trajectory_motifs
#' @inheritDotParams distill_traj_model_multi
#' @export
regress_trajectory_motifs_manifold <- function(peak_intervals,
                                               atac_diff_add_mat,
                                               norm_intervals = NULL,
                                               max_motif_num = 120,
                                               target_traj_motif_num = 30,
                                               n_clust_factor = 1,
                                               motif_energies = NULL,
                                               norm_motif_energies = NULL,
                                               pssm_db = iceqream::motif_db,
                                               additional_features = NULL,
                                               min_tss_distance = 5000,
                                               bin_start = 1,
                                               bin_end = NULL,
                                               min_initial_energy_cor = 0.05,
                                               normalize_energies = TRUE,
                                               energy_norm_quantile = 1,
                                               norm_energy_max = 10,
                                               n_prego_motifs = 0,
                                               traj_prego = NULL,
                                               initial_min_diff = 0.1,
                                               min_diff = initial_min_diff,
                                               prego_sample_for_kmers = TRUE,
                                               prego_sample_fraction = 0.1,
                                               seed = 60427,
                                               feature_selection_beta = 0.003,
                                               lambda = 1e-5,
                                               alpha = 1,
                                               filter_using_r2 = FALSE,
                                               r2_threshold = 0.0005,
                                               parallel = TRUE,
                                               peaks_size = 500,
                                               spat_num_bins = NULL,
                                               spat_bin_size = 2,
                                               kmer_sequence_length = 300,
                                               include_interactions = FALSE,
                                               interaction_threshold = 0.001,
                                               max_motif_interaction_n = NULL,
                                               max_add_interaction_n = NULL,
                                               ...) {
    withr::local_options(list(gmax.data.size = 1e9))


    validate_peak_intervals(peak_intervals)
    validate_additional_features(additional_features, peak_intervals)
    additional_features[is.na(additional_features)] <- 0
    if (is.null(norm_intervals)) {
        norm_intervals <- peak_intervals
    }

    validate_motif_energies(motif_energies, peak_intervals, pssm_db)

    min_energy <- -7

    validate_misha()

    # filter peaks that are too close to TSS
    enhancers_filter <- get_tss_distance_filter(peak_intervals, min_tss_distance)

    peak_intervals_all <- peak_intervals
    peak_intervals <- peak_intervals[enhancers_filter, ]
    atac_diff_add_mat <- atac_diff_add_mat[enhancers_filter, ]
    motif_energies <- motif_energies[enhancers_filter, ]

    if (!is.null(additional_features)) {
        additional_features <- additional_features[enhancers_filter, ]
    }

    cli_alert_info("Number of peaks: {.val {nrow(peak_intervals)}}")

    cli_alert("Extracting sequences...")
    all_seqs <- toupper(misha::gseq.extract(misha.ext::gintervals.normalize(peak_intervals_all, peaks_size)))
    norm_seqs <- toupper(misha::gseq.extract(misha.ext::gintervals.normalize(norm_intervals, peaks_size)))

    motif_energies_all <- motif_energies

    cli_alert_info("Calculating correlations between {.val {ncol(motif_energies)}} motif energies and ATAC differences...")
    atac_diff_add_mat_n <- apply(atac_diff_add_mat, 2, norm01)
    atac_diff_add_mat_n[abs(atac_diff_add_mat_n) < initial_min_diff] <- NA
    cm <- tgs_cor(motif_energies, as.matrix(atac_diff_add_mat_n), pairwise.complete.obs = TRUE)

    cli::cli_alert("Collecting motifs from {.val {ncol(atac_diff_add_mat)}} trajectories")
    traj_model_ls <- purrr::imap(as.data.frame(atac_diff_add_mat), function(atac_diff, traj_nm) {
        cli::cli_h2("Processing trajectory {.val {traj_nm}}")
        # normalize differential accessibility
        atac_diff_n <- norm01(atac_diff)

        diff_filter <- abs(atac_diff) >= initial_min_diff
        diff_filter[is.na(diff_filter)] <- FALSE

        if (is.null(traj_prego) && n_prego_motifs > 0) {
            traj_prego <- learn_traj_prego(peak_intervals, atac_diff,
                n_motifs = n_prego_motifs, min_diff = initial_min_diff,
                sample_for_kmers = prego_sample_for_kmers,
                sample_fraction = prego_sample_fraction,
                energy_norm_quantile = energy_norm_quantile,
                sequences = all_seqs,
                norm_intervals = norm_intervals,
                seed = seed,
                spat_bin_size = spat_bin_size,
                spat_num_bins = spat_num_bins
            )
        }

        if (!is.null(traj_prego)) {
            prego_models <- traj_prego$models
            prego_e <- traj_prego$energies
            prego_pssm <- traj_prego$pssm
            if (!is.null(min_tss_distance)) {
                prego_e <- prego_e[enhancers_filter, ]
            }
            motif_energies <- cbind(motif_energies, prego_e)
            pssm_db <- bind_rows(pssm_db, prego_pssm)
        } else {
            prego_models <- list()
        }

        motifs <- select_motifs_by_correlation(motif_energies, atac_diff, min_initial_energy_cor, target_traj_motif_num, cm = cm[, traj_nm])
        motif_energies <- motif_energies_all[, motifs]

        result <- select_features_by_regression(
            motif_energies = motif_energies,
            atac_diff_n = atac_diff_n,
            additional_features = additional_features,
            feature_selection_beta = feature_selection_beta,
            alpha = alpha,
            lambda = lambda,
            parallel = parallel,
            seed = seed
        )
        features_mat <- result$features_mat
        glm_model2 <- result$glm_model2
        chosen_motifs <- result$chosen_motifs

        chosen_motifs_ls <- purrr::map(chosen_motifs, function(x) {
            list(pssm = pssm_db %>%
                filter(motif == x) %>%
                select(pos, A, C, G, T), spat = NULL)
        })
        names(chosen_motifs_ls) <- chosen_motifs

        clust_energies_logist <- create_logist_features(features_mat)

        TrajectoryModel(
            model = glm_model2,
            motif_models = homogenize_pssm_models(chosen_motifs_ls),
            coefs = get_model_coefs(glm_model2),
            normalized_energies = as.matrix(features_mat),
            model_features = features_mat,
            type = rep("train", length(atac_diff)),
            additional_features = as.data.frame(additional_features),
            diff_score = atac_diff,
            predicted_diff_score = atac_diff,
            initial_prego_models = prego_models,
            peak_intervals = peak_intervals,
            normalization_intervals = norm_intervals,
            interactions = matrix(0, nrow = nrow(features_mat), ncol = 0),
            params = list(
                energy_norm_quantile = energy_norm_quantile,
                norm_energy_max = norm_energy_max,
                min_energy = min_energy,
                alpha = alpha,
                lambda = lambda,
                peaks_size = peaks_size,
                spat_num_bins = spat_num_bins,
                spat_bin_size = spat_bin_size,
                distilled_features = motif_energies,
                n_clust_factor = n_clust_factor,
                include_interactions = include_interactions,
                interaction_threshold = interaction_threshold,
                seed = seed,
                features_type = "linear"
            )
        )
    })
    names(traj_model_ls) <- colnames(atac_diff_add_mat)

    cli::cli_h3("Number of motifs selected for each trajectory:")
    purrr::iwalk(traj_model_ls, ~ cli::cli_bullets("{.val {.y}}: {.val {length(.x@motif_models)}}"))

    n_unique_models <- purrr::map(traj_model_ls, ~ names(.x@motif_models)) %>%
        unlist() %>%
        unique() %>%
        length()
    cli_alert_info("Number of unique motifs selected: {.val {n_unique_models}}")

    cli::cli_h2("Distilling trajectory models...")
    traj_model_multi <- distill_traj_model_multi(traj_model_ls, max_motif_num = max_motif_num, min_diff = min_diff, seed = seed, filter_models = FALSE, ...)

    return(traj_model_multi)
}
