#' Perform motif regression on ATAC trajectories.
#'
#' @description This function performs motif regression on ATAC trajectories. It can work with either ATAC scores on trajectory bins or directly with differential accessibility.
#'
#' The basic inputs for regression are the genomic positions of the peaks, two vectors of ATAC scores (as an \code{atac_scores} matrix) or differential accessibility (\code{atac_diff}), and database energies computed for all the genomics positions.
#' Optionally, additional features such as epigenomic features can be provided.
#'
#' @param peak_intervals A data frame, indicating the genomic positions ('chrom', 'start', 'end') of each peak.
#' @param atac_scores Optional. A numeric matrix, representing mean ATAC score per bin per peak. Rows: peaks, columns: bins. By default iceqream would regress the last column minus the first column. If you want to regress something else, please either change bin_start or bin_end, or provide \code{atac_diff} instead. If \code{normalize_bins} is TRUE, the scores will be normalized to be between 0 and 1.
#' @param atac_diff Optional. A numeric vector representing the differential accessibility between the start and end of the trajectory. Either this or \code{atac_scores} must be provided.
#' @param normalize_bins whether to normalize the ATAC scores to be between 0 and 1. Default: TRUE
#' @param norm_intervals A data frame, indicating the genomic positions ('chrom', 'start', 'end') of peaks used for energy normalization. If NULL, the function will use \code{peak_intervals} for normalization.
#' @param max_motif_num maximum number of motifs to consider. Default: 50
#' @param n_clust_factor factor to divide the number of to keep after clustering. e.g. if n_clust_factor > 1 the number of motifs to keep will be reduced by a factor of n_clust_factor. Default: 1
#' @param motif_energies A numeric matrix, representing the energy of each motif in each peak. If NULL, the function will use \code{pssm_db} to calculate the motif energies. Note that this might take a while.
#' @param norm_motif_energies A numeric matrix, representing the normalized energy of each motif in each interval of \code{norm_intervals}. If NULL, the function will use \code{pssm_db} to calculate the motif energies. Note that this might take a while.
#' @param pssm_db a data frame with PSSMs ('A', 'C', 'G' and 'T' columns), with an additional column 'motif' containing the motif name. All the motifs in \code{motif_energies} (column names) should be present in the 'motif' column. Default: all motifs.
#' @param additional_features A data frame, representing additional genomic features (e.g. CpG content, distance to TSS, etc.) for each peak. Note that NA values would be replaced with 0.
#' @param min_tss_distance distance from Transcription Start Site (TSS) to classify a peak as an enhancer. Default: 5000. If NULL, no filtering will be performed - use this option if your peaks are already filtered. \cr
#' Note that in order to filter peaks that are too close to TSS, the current \code{misha} genome must have an intervals set called \code{intervs.global.tss}.
#' @param bin_start the start of the trajectory. Default: 1
#' @param bin_end the end of the trajectory. Default: the last bin (only used when atac_scores is provided)
#' @param normalize_energies whether to normalize the motif energies. Set this to FALSE if the motif energies are already normalized.
#' @param min_initial_energy_cor minimal correlation between the motif normalized energy and the ATAC difference.
#' @param energy_norm_quantile quantile of the energy used for normalization. Default: 1
#' @param norm_energy_max maximum value of the normalized energy. Default: 10
#' @param n_prego_motifs number of prego motifs (de-novo motifs) to consider.
#' @param traj_prego output of \code{learn_traj_prego}. If provided, no additional prego models would be inferred.
#' @param min_diff minimal ATAC difference for a peak to participate in the initial prego motif inference and in the distillation step (if \code{distill_on_diff} is TRUE).
#' @param distill_on_diff whether to distill motifs based on differential accessibility. If FALSE, all peaks will be used for distillation, if TRUE - only peaks with differential accessibility >= min_diff will be used.
#' @param prego_sample_for_kmers whether to use a sample of the peaks for kmer screening. Default: TRUE
#' @param prego_sample_fraction Fraction of peaks to sample for prego motif inference. A smaller number would be faster but might lead to over-fitting. Default: 0.1
#' @param prego_enrgy_norm_quantile,prego_spat_bin_size,prego_spat_num_bins parameters for prego motif inference. See \code{prego::regress_pwm} for more details.
#' @param seed random seed for reproducibility.
#' @param feature_selection_beta beta parameter used for feature selection.
#' @param filter_using_r2 whether to filter features using R^2.
#' @param r2_threshold minimal R^2 for a feature to be included in the model.
#' @param parallel whether to use parallel processing on glmnet.
#' @param peaks_size size of the peaks to extract sequences from. Default: 500bp
#' @param spat_num_bins number of spatial bins to use.
#' @param spat_bin_size size of each spatial bin.
#' @param kmer_sequence_length length of the kmer sequence to use for kmer screening. By default the full sequence is used.
#' @param include_interactions whether to include interactions between motifs / additional fetures as model features. IQ will create interactions between significant additional features and all motifs, and between significant motifs. Default: FALSE

#' @param max_motif_interaction_n maximum number of motifs to consider for interactions. If NULL, all motifs above the interaction_threshold will be considered. Default: NULL
#' @param max_add_interaction_n maximum number of additional features to consider for interactions. If NULL, all additional features above the interaction_threshold will be considered. Default: NULL
#' @param max_interaction_n maximum number of interactions to consider.
#'
#' @return An instance of `TrajectoryModel` containing model information and results:
#' \itemize{
#'   \item model: The final General Linear Model (GLM) object.
#'   \item motif_models: Named List, PSSM and spatial models for each motif cluster.
#'   \item normalized_energies: Numeric vector, normalized energies of each motif in each peak.
#'   \item additional_features: data frame of the additional features.
#'   \item diff_score: Numeric, normalized score of differential accessibility between 'bin_start' and 'bin_end'.
#'   \item predicted_diff_score: Numeric, predicted differential accessibility score between 'bin_start' and 'bin_end'.
#'   \item initial_prego_models: List, inferred prego models at the initial step of the algorithm.
#'   \item peak_intervals: data frame, indicating the genomic positions ('chrom', 'start', 'end') of each peak used for training.
#' }
#' Print the model to see more details.
#'
#' @inheritParams glmnet::glmnet
#' @inheritParams add_interactions
#' @inheritParams prego::regress_pwm
#' @export
regress_trajectory_motifs <- function(peak_intervals,
                                      atac_scores = NULL,
                                      atac_diff = NULL,
                                      normalize_bins = TRUE,
                                      norm_intervals = NULL,
                                      max_motif_num = 30,
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
                                      min_diff = 0.1,
                                      distill_on_diff = FALSE,
                                      prego_min_diff = min_diff,
                                      prego_sample_for_kmers = TRUE,
                                      prego_sample_fraction = 0.1,
                                      prego_energy_norm_quantile = 1,
                                      prego_spat_bin_size = NULL,
                                      prego_spat_num_bins = NULL,
                                      seed = 60427,
                                      feature_selection_beta = 0.003,
                                      lambda = 1e-5,
                                      alpha = 1,
                                      filter_using_r2 = FALSE,
                                      r2_threshold = 0.0005,
                                      parallel = TRUE,
                                      peaks_size = 500,
                                      spat_num_bins = NULL,
                                      spat_bin_size = NULL,
                                      kmer_sequence_length = 300,
                                      include_interactions = FALSE,
                                      interaction_threshold = 0.001,
                                      max_motif_interaction_n = NULL,
                                      max_add_interaction_n = NULL,
                                      max_interaction_n = NULL,
                                      symmetrize_spat = TRUE) {
    withr::local_options(list(gmax.data.size = 1e9))

    if (is.null(atac_scores) && is.null(atac_diff)) {
        cli_abort("Either 'atac_scores' or 'atac_diff' must be provided.")
    }

    if (!is.null(atac_scores) && !is.null(atac_diff)) {
        cli_abort("Only one of 'atac_scores' or 'atac_diff' should be provided.")
    }

    if (!is.null(atac_scores)) {
        atac_scores <- as.matrix(atac_scores)
        if (is.null(bin_end)) {
            bin_end <- ncol(atac_scores)
        }
        validate_atac_scores(atac_scores, bin_start, bin_end)
        if (normalize_bins) {
            atac_scores[, bin_start] <- norm01(atac_scores[, bin_start])
            atac_scores[, bin_end] <- norm01(atac_scores[, bin_end])
        }
        if (nrow(peak_intervals) != nrow(atac_scores)) {
            cli_abort("Number of rows in {.field {peak_intervals}} must be equal to the number of rows in {.field {atac_scores}}")
        }
        atac_diff <- atac_scores[, bin_end] - atac_scores[, bin_start]
    } else {
        if (length(atac_diff) != nrow(peak_intervals)) {
            cli_abort("Length of {.field {atac_diff}} must be equal to the number of rows in {.field {peak_intervals}}")
        }
        atac_scores <- matrix(0, nrow = length(atac_diff), ncol = 2)
        atac_scores[, 2] <- atac_diff
    }

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
    atac_scores <- atac_scores[enhancers_filter, ]
    atac_diff <- atac_diff[enhancers_filter]
    motif_energies <- motif_energies[enhancers_filter, ]

    if (!is.null(additional_features)) {
        additional_features <- additional_features[enhancers_filter, ]
    }

    cli_alert_info("Number of peaks: {.val {nrow(peak_intervals)}}")

    # normalize differential accessibility
    atac_diff_n <- norm01(atac_diff)

    diff_filter <- abs(atac_diff) >= min_diff
    diff_filter[is.na(diff_filter)] <- FALSE

    cli_alert("Extracting sequences...")
    all_seqs <- prego::intervals_to_seq(peak_intervals_all, peaks_size)
    norm_seqs <- prego::intervals_to_seq(norm_intervals, peaks_size)


    if (is.null(traj_prego) && n_prego_motifs > 0) {
        traj_prego <- learn_traj_prego(peak_intervals, atac_diff,
            n_motifs = n_prego_motifs, min_diff = prego_min_diff,
            sample_for_kmers = prego_sample_for_kmers,
            sample_fraction = prego_sample_fraction, energy_norm_quantile = prego_energy_norm_quantile, sequences = all_seqs, norm_intervals = norm_intervals, seed = seed, spat_bin_size = prego_spat_bin_size, spat_num_bins = prego_spat_num_bins, peaks_size = peaks_size
        )
    }

    if (!is.null(traj_prego)) {
        prego_models <- traj_prego$models
        prego_e <- traj_prego$energies
        prego_pssm <- traj_prego$pssm
        if (!is.null(min_tss_distance) && min_tss_distance > 0) {
            prego_e <- prego_e[enhancers_filter, ]
        }
        motif_energies <- cbind(motif_energies, prego_e)
        pssm_db <- bind_rows(pssm_db, prego_pssm)
    } else {
        prego_models <- list()
    }

    motifs <- select_motifs_by_correlation(motif_energies, atac_diff, min_initial_energy_cor, max_motif_num)

    motif_energies <- motif_energies[, motifs]
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
    features <- result$features_mat
    glm_model2 <- result$glm_model2
    chosen_motifs <- result$chosen_motifs

    if (distill_on_diff) {
        diff_filter <- abs(atac_diff) >= min_diff
        diff_filter[is.na(diff_filter)] <- FALSE
    } else {
        diff_filter <- rep(TRUE, nrow(peak_intervals))
    }

    distilled <- distill_motifs(
        features,
        max_motif_num,
        glm_model2,
        y = atac_diff_n,
        seqs = all_seqs[enhancers_filter],
        norm_seqs = norm_seqs,
        diff_filter,
        additional_features = additional_features,
        pssm_db = pssm_db,
        lambda = lambda,
        alpha = alpha,
        energy_norm_quantile = energy_norm_quantile,
        norm_energy_max = norm_energy_max,
        min_energy = min_energy,
        seed = seed,
        spat_num_bins = spat_num_bins,
        spat_bin_size = spat_bin_size,
        kmer_sequence_length = kmer_sequence_length,
        n_clust_factor = n_clust_factor,
        distill_single = FALSE,
        symmetrize_spat = symmetrize_spat
    )
    clust_energies <- distilled$energies

    if (include_interactions) {
        interactions <- get_significant_interactions(
            energies = clust_energies,
            y = atac_diff_n,
            interaction_threshold = interaction_threshold, additional_features = additional_features,
            max_motif_n = max_motif_interaction_n,
            max_add_n = max_add_interaction_n,
            max_n = max_interaction_n,
            lambda = lambda, alpha = alpha, seed = seed
        )
        clust_energies <- cbind(clust_energies, interactions)
    } else {
        interactions <- matrix(0, nrow = nrow(clust_energies), ncol = 0)
    }

    clust_energies_logist <- create_logist_features(clust_energies)

    cli_alert_info("Running final regression, number of features: {.val {ncol(clust_energies_logist)}}")
    model <- glmnet::glmnet(clust_energies_logist, atac_diff_n, binomial(link = "logit"), alpha = alpha, lambda = lambda, parallel = parallel, seed = seed)
    model <- strip_glmnet(model)

    predicted_diff_score <- logist(glmnet::predict.glmnet(model, newx = clust_energies_logist, type = "link", s = lambda))[, 1]
    predicted_diff_score <- norm01(predicted_diff_score)
    predicted_diff_score <- rescale(predicted_diff_score, atac_diff)

    cli_alert_success("Finished running model. Number of non-zero coefficients: {.val {sum(model$beta != 0)}} (out of {.val {ncol(clust_energies_logist)}}). R^2: {.val {cor(predicted_diff_score, atac_diff_n)^2}}")

    # remove additional features from clust_energies
    normalized_energies <- clust_energies[, setdiff(colnames(clust_energies), colnames(additional_features))]

    if (include_interactions) {
        normalized_energies <- normalized_energies[, setdiff(colnames(normalized_energies), colnames(interactions))]
    }

    traj_model <- TrajectoryModel(
        model = model,
        motif_models = homogenize_pssm_models(distilled$motifs),
        coefs = get_model_coefs(model),
        normalized_energies = as.matrix(normalized_energies),
        model_features = clust_energies_logist,
        type = rep("train", nrow(atac_scores)),
        additional_features = as.data.frame(additional_features),
        diff_score = atac_diff,
        predicted_diff_score = predicted_diff_score,
        initial_prego_models = prego_models,
        peak_intervals = peak_intervals,
        normalization_intervals = norm_intervals,
        interactions = interactions,
        params = list(
            energy_norm_quantile = energy_norm_quantile,
            norm_energy_max = norm_energy_max,
            min_energy = min_energy,
            alpha = alpha,
            lambda = lambda,
            peaks_size = peaks_size,
            spat_num_bins = spat_num_bins,
            spat_bin_size = spat_bin_size,
            distilled_features = distilled$features,
            n_clust_factor = n_clust_factor,
            include_interactions = include_interactions,
            interaction_threshold = interaction_threshold,
            symmetrize_spat = symmetrize_spat,
            seed = seed
        )
    )

    if (filter_using_r2) {
        traj_model <- filter_traj_model(traj_model, r2_threshold = r2_threshold)
    }


    return(traj_model)
}


validate_atac_scores <- function(atac_scores, bin_start, bin_end) {
    # validate bin_start and bin_end
    if (bin_start < 1 || bin_start > ncol(atac_scores)) {
        cli_abort("{.field 'bin_start'} must be between 1 and {.code 'ncol(atac_scores)'}", call = parent.frame(1))
    }
    if (bin_end < 1 || bin_end > ncol(atac_scores)) {
        cli_abort("{.field 'bin_end'} must be between 1 and {.code 'ncol(atac_scores)'}", call = parent.frame(1))
    }
    if (bin_start >= bin_end) {
        cli_abort("{.field 'bin_start'} must be smaller than {.field 'bin_end'}", call = parent.frame(1))
    }
}

validate_peak_intervals <- function(peak_intervals, columns = c("chrom", "start", "end")) {
    # make sure that peak intervals have 'chrom' 'start' 'end' columns
    if (!all(columns %in% colnames(peak_intervals))) {
        cli_abort("{.field 'peak_intervals'} must have '{.val {columns}} columns", call = parent.frame(1))
    }

    peak_sizes <- peak_intervals$end - peak_intervals$start

    # make sure that all peak sizes are equal
    if (length(unique(peak_sizes)) != 1) {
        cli_abort("All peaks must have the same size", call = parent.frame(1))
    }

    # warn if peak sizes are smaller than 100bp
    if (unique(peak_sizes) < 100) {
        cli_warn("Peak sizes are smaller than 100bp", call = parent.frame(1))
    }
}

validate_additional_features <- function(additional_features, peak_intervals) {
    if (!is.null(additional_features)) {
        if (nrow(additional_features) != nrow(peak_intervals)) {
            cli_abort("{.field 'additional_features'} must have the same number of rows as {.field 'peak_intervals'}. intervals number of rows {.val {nrow(peak_intervals)}} != additional_features number of rows {.val {nrow(additional_features)}}}}", call = parent.frame(1))
        }
    }
    if (!is.null(additional_features)) {
        if (any(is.na(additional_features))) {
            cli_warn("NA values in additional features are replaced with 0", call = parent.frame(1))
        }
    }
}


get_model_coefs <- function(model, s = model$lambda) {
    df <- glmnet::coef.glmnet(model, s = s) %>%
        as.matrix() %>%
        as.data.frame() %>%
        tibble::rownames_to_column("variable")
    df <- df %>% filter(variable != "(Intercept)")
    colnames(df)[2] <- "s1"

    df <- df %>%
        mutate(type = sub(".*_", "", variable), variable = sub("_(low-energy|high-energy|higher-energy|sigmoid)$", "", variable)) %>%
        tidyr::spread(type, s1)

    df[is.na(df)] <- 0

    return(df)
}
