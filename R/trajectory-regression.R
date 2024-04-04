#' Perform motif regression on ATAC trajectories.
#'
#' @description This function performs motif regression on ATAC trajectories. Given ATAC scores on trajectory bins, it predicts the differential accessibility between the start and end of the trajectory and returns the motifs that are most likely to be responsible for this differential accessibility.
#'
#' @param atac_scores A numeric matrix, representing mean ATAC score per bin per peak. Rows: peaks, columns: bins.
#' @param peak_intervals A data frame, indicating the genomic positions ('chrom', 'start', 'end') of each peak, with an additional column named "const" indicating whether the peak is constitutive. Optionally, a column named "cluster" can be added with indication of the cluster of each peak.
#' @param norm_intervals A data frame, indicating the genomic positions ('chrom', 'start', 'end') of peaks used for energy normalization. If NULL, the function will use \code{peak_intervals} for normalization.
#' @param motif_energies A numeric matrix, representing the energy of each motif in each peak. If NULL, the function will use \code{pssm_db} to calculate the motif energies. Note that this might take a while.
#' @param norm_motif_energies A numeric matrix, representing the normalized energy of each motif in each interval of \code{norm_intervals}. If NULL, the function will use \code{pssm_db} to calculate the motif energies. Note that this might take a while.
#' @param pssm_db a data frame with PSSMs ('A', 'C', 'G' and 'T' columns), with an additional column 'motif' containing the motif name. All the motifs in \code{motif_energies} (column names) should be present in the 'motif' column. Default: all motifs in the prego package.
#' @param additional_features A data frame, representing additional genomic features (e.g. CpG content, distance to TSS, etc.) for each peak. Note that NA values would be replaced with 0.
#' @param max_motif_num maximum number of motifs to consider. Default: 50
#' @param n_clust_factor factor to divide the number of to keep after clustering. e.g. if n_clust_factor > 1 the number of motifs to keep will be reduced by a factor of n_clust_factor. Default: 1
#' @param min_tss_distance distance from Transcription Start Site (TSS) to classify a peak as an enhancer. Default: 5000. If NULL, no filtering will be performed - use this option if your peaks are already filtered. \cr
#' Note that in order to filter peaks that are too close to TSS, the current \code{misha} genome must have an intervals set called \code{intervs.global.tss}.
#' @param bin_start the start of the trajectory. Default: 1
#' @param bin_end the end of the trajectory. Default: the last bin
#' @param normalize_energies whether to normalize the motif energies. Set this to FALSE if the motif energies are already normalized.
#' @param min_initial_energy_cor minimal correlation between the motif normalized energy and the ATAC difference.
#' @param energy_norm_quantile quantile of the energy used for normalization. Default: 1
#' @param norm_energy_max maximum value of the normalized energy. Default: 10
#' @param n_prego_motifs number of prego motifs to consider.
#' @param traj_prego output of \code{learn_traj_prego}. If provided, no additional prego models would be inferred.
#' @param min_diff minimal ATAC difference for a peak to participate in the initial prego motif inference and in the distillation step (if \code{distill_on_diff} is TRUE).
#' @param distill_on_diff whether to distill motifs based on differential accessibility. If FALSE, all peaks will be used for distillation, if TRUE - only peaks with differential accessibility >= min_diff will be used.
#' @param prego_sample_fraction Fraction of peaks to sample for prego motif inference. A smaller number would be faster but might lead to over-fitting. Default: 0.1
#' @param seed random seed for reproducibility.
#' @param feature_selection_beta beta parameter used for feature selection.
#' @param filter_using_r2 whether to filter features using R^2.
#' @param r2_threshold minimal R^2 for a feature to be included in the model.
#' @param parallel whether to use parallel processing on glmnet.
#' @param peaks_size size of the peaks to extract sequences from. Default: 300bp
#' @param spat_num_bins number of spatial bins to use.
#' @param spat_bin_size size of each spatial bin.
#' @param kmer_sequence_length length of the kmer sequence to use for kmer screening. By default the full sequence is used.
#'
#' @return An instance of `TrajectoryModel` containing:
#' \itemize{
#'   \item{model}{The final General Linear Model (GLM) object.}
#'   \item{motif_models}{Named List, PSSM and spatial models for each motif cluster.}
#'   \item{normalized_energies}{Numeric vector, normalized energies of each motif in each peak.}
#'   \item{additional_features}{data frame of the additional features.}
#'   \item{diff_score}{Numeric, normalized score of differential accessibility between 'bin_start' and 'bin_end'.}
#'   \item{predicted_diff_score}{Numeric, predicted differential accessibility score between 'bin_start' and 'bin_end'.}
#'   \item{initial_prego_models}{List, inferred prego models at the initial step of the algorithm.}
#'   \item{peak_intervals}{data frame, indicating the genomic positions ('chrom', 'start', 'end') of each peak used for training.}
#' }
#'
#' @inheritParams glmnet::glmnet
#' @export
regress_trajectory_motifs <- function(atac_scores,
                                      peak_intervals,
                                      norm_intervals = peak_intervals,
                                      max_motif_num = 50,
                                      n_clust_factor = 1,
                                      motif_energies = NULL,
                                      norm_motif_energies = NULL,
                                      pssm_db = prego::all_motif_datasets(),
                                      additional_features = NULL,
                                      min_tss_distance = 5000,
                                      bin_start = 1,
                                      bin_end = ncol(atac_scores),
                                      min_initial_energy_cor = 0.05,
                                      normalize_energies = TRUE,
                                      energy_norm_quantile = 1,
                                      norm_energy_max = 10,
                                      n_prego_motifs = 4,
                                      traj_prego = NULL,
                                      min_diff = 0.1,
                                      distill_on_diff = FALSE,
                                      prego_sample_fraction = 0.1,
                                      seed = 60427,
                                      feature_selection_beta = 0.003,
                                      lambda = 1e-5,
                                      alpha = 1,
                                      filter_using_r2 = FALSE,
                                      r2_threshold = 0.0005,
                                      parallel = TRUE,
                                      peaks_size = 300,
                                      spat_num_bins = NULL,
                                      spat_bin_size = NULL,
                                      kmer_sequence_length = NULL) {
    withr::local_options(list(gmax.data.size = 1e9))
    atac_scores <- as.matrix(atac_scores)

    validate_atac_scores(atac_scores, bin_start, bin_end)
    validate_peak_intervals(peak_intervals, atac_scores)
    validate_additional_features(additional_features, peak_intervals)

    # Extract features
    motif_energies <- calc_motif_energies(peak_intervals, pssm_db, motif_energies, field_name = "motif_energies", intervals_field_name = "peak_intervals")

    min_energy <- -7

    if (normalize_energies) {
        cli_alert_info("Normalizing motif energies...")
        norm_motif_energies <- calc_motif_energies(norm_intervals, pssm_db, norm_motif_energies, field_name = "norm_motif_energies", intervals_field_name = "norm_intervals")
        motif_energies <- norm_energy_matrix(motif_energies, norm_motif_energies, min_energy = min_energy, q = energy_norm_quantile, norm_energy_max = norm_energy_max)
    }

    # filter peaks that are too close to TSS
    if (!is.null(min_tss_distance)) {
        if (!misha::gintervals.exists("intervs.global.tss")) {
            cli_abort("Please make sure the current genome ({.field {GROOT}}) has an intervals set called {.val intervs.global.tss}")
        }

        tss_dist <- abs(misha::gintervals.neighbors(peak_intervals, "intervs.global.tss", na.if.notfound = TRUE)$dist)
        enhancers_filter <- tss_dist > min_tss_distance
        enhancers_filter[is.na(enhancers_filter)] <- FALSE
        if (sum(!enhancers_filter) > 0) {
            cli_alert_info("{.val {sum(!enhancers_filter)}} peaks were filtered out because they are too close to TSS (<= {.val min_tss_distance}bp)")
        }

        peak_intervals_all <- peak_intervals
        peak_intervals <- peak_intervals[enhancers_filter, ]
        atac_scores <- atac_scores[enhancers_filter, ]
        motif_energies <- motif_energies[enhancers_filter, ]

        if (!is.null(additional_features)) {
            additional_features <- additional_features[enhancers_filter, ]
        }
    } else {
        peak_intervals_all <- peak_intervals
        enhancers_filter <- rep(TRUE, nrow(peak_intervals))
    }

    cli_alert_info("Number of peaks: {.val {nrow(peak_intervals)}}")

    # calculate differential accessibility
    atac_diff <- atac_scores[, bin_end] - atac_scores[, bin_start]
    atac_diff_n <- norm01(atac_diff)

    diff_filter <- abs(atac_diff) >= min_diff
    diff_filter[is.na(diff_filter)] <- FALSE

    cli_alert("Extracting sequences...")
    all_seqs <- toupper(misha::gseq.extract(misha.ext::gintervals.normalize(peak_intervals_all, peaks_size)))
    norm_seqs <- toupper(misha::gseq.extract(misha.ext::gintervals.normalize(norm_intervals, peaks_size)))


    if (is.null(traj_prego) && n_prego_motifs > 0) {
        traj_prego <- learn_traj_prego(peak_intervals, atac_diff, n_motifs = n_prego_motifs, min_diff = min_diff, sample_fraction = prego_sample_fraction, energy_norm_quantile = energy_norm_quantile, sequences = all_seqs, norm_seqeunces = norm_seqs, seed = seed)
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


    cli_alert_info("Calculating correlations between {.val {ncol(motif_energies)}} motif energies and ATAC difference...")
    cm <- tgs_cor(motif_energies, as.matrix(atac_diff), pairwise.complete.obs = TRUE)[, 1]
    motifs <- names(cm[abs(cm) >= min_initial_energy_cor])
    if (length(motifs) == 0) {
        min_initial_energy_cor <- min_initial_energy_cor / 2
        cli::cli_alert_warning("No features with absolute correlation >= {.val {min_initial_energy_cor}}. Trying again with {.val {min_initial_energy_cor}}")
        motifs <- names(cm[abs(cm) >= min_initial_energy_cor])
        if (length(motifs) == 0) {
            n <- min(30, length(motifs))
            motifs <- names(sort(abs(cm), decreasing = TRUE)[1:n])
            cli::cli_alert_warning("No features with absolute correlation >= {.val {min_initial_energy_cor}}. Trying again with top {.val {n}} features")
        }
    }

    cli_alert_info("Selected {.val {length(motifs)}} (out of {.val {ncol(motif_energies)}}) features with absolute correlation >= {.val {min_initial_energy_cor}}")
    motifs <- motifs[!is.na(motifs)]
    motif_energies <- motif_energies[, motifs]

    cli_alert_info("Running first round of regression, # of features: {.val {ncol(motif_energies)}}")
    glm_model1 <- glmnet::glmnet(motif_energies, atac_diff_n, binomial(link = "logit"), alpha = alpha, lambda = lambda, parallel = parallel, seed = seed)

    features <- rownames(glm_model1$beta)[abs(glm_model1$beta[, 1]) >= feature_selection_beta]
    cli_alert_info("Taking {.val {length(features)}} features with beta >= {.val {feature_selection_beta}}")
    if (length(features) == 0) {
        cli::cli_alert_warning("No features with beta >= {.val {feature_selection_beta}}. Using all features.")
        features <- rownames(glm_model1$beta)
    }

    cli_alert_info("Running second round of regression...")
    additional_features[is.na(additional_features)] <- 0
    glm_model2 <- glmnet::glmnet(as.matrix(cbind(motif_energies[, features], additional_features)), atac_diff_n, binomial(link = "logit"), alpha = alpha, lambda = lambda, parallel = parallel, seed = seed)

    chosen_motifs <- rownames(glm_model2$beta)[abs(glm_model2$beta[, 1]) > 0]
    if (length(chosen_motifs) == 0) {
        cli::cli_alert_warning("No features with beta > 0. Using all features.")
        chosen_motifs <- rownames(glm_model2$beta)
    }
    features <- motif_energies[, setdiff(chosen_motifs, colnames(additional_features))]

    if (distill_on_diff) {
        diff_filter <- abs(atac_diff) >= min_diff
        diff_filter[is.na(diff_filter)] <- FALSE
    } else {
        diff_filter <- rep(TRUE, nrow(peak_intervals))
    }
    distilled <- distill_motifs(features, max_motif_num, glm_model2, y = atac_diff_n, seqs = all_seqs[enhancers_filter], norm_seqs = norm_seqs, diff_filter, additional_features = additional_features, pssm_db = pssm_db, lambda = lambda, alpha = alpha, energy_norm_quantile = energy_norm_quantile, norm_energy_max = norm_energy_max, min_energy = min_energy, seed = seed, spat_num_bins = spat_num_bins, spat_bin_size = spat_bin_size, kmer_sequence_length = kmer_sequence_length, n_clust_factor = n_clust_factor, distill_single = FALSE)
    clust_energies <- distilled$energies

    clust_energies_logist <- create_logist_features(clust_energies)

    model <- glmnet::glmnet(clust_energies_logist, atac_diff_n, binomial(link = "logit"), alpha = alpha, lambda = lambda, parallel = parallel, seed = seed)

    predicted_diff_score <- logist(glmnet::predict.glmnet(model, newx = clust_energies_logist, type = "link", s = lambda))[, 1]
    predicted_diff_score <- norm01(predicted_diff_score)
    predicted_diff_score <- rescale(predicted_diff_score, atac_diff)


    cli_alert_success("Finished running model. Number of non-zero coefficients: {.val {sum(model$beta != 0)}} (out of {.val {ncol(clust_energies_logist)}}). R^2: {.val {cor(predicted_diff_score, atac_diff_n)^2}}")

    traj_model <- TrajectoryModel(
        model = model,
        motif_models = homogenize_pssm_models(distilled$motifs),
        coefs = get_model_coefs(model),
        normalized_energies = clust_energies,
        model_features = clust_energies_logist,
        type = rep("train", nrow(atac_scores)),
        additional_features = as.data.frame(additional_features),
        diff_score = atac_diff,
        predicted_diff_score = predicted_diff_score,
        initial_prego_models = prego_models,
        peak_intervals = peak_intervals,
        normalization_intervals = norm_intervals,
        params = list(
            energy_norm_quantile = energy_norm_quantile,
            norm_energy_max = norm_energy_max,
            min_energy = min_energy,
            alpha = alpha,
            lambda = lambda,
            peaks_size = peaks_size,
            spat_num_bins = spat_num_bins,
            spat_bin_size = spat_bin_size,
            n_clust_factor = n_clust_factor,
            seed = seed
        )
    )

    if (filter_using_r2) {
        coefs <- get_model_coefs(model)
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

validate_peak_intervals <- function(peak_intervals, atac_scores, columns = c("chrom", "start", "end", "const")) {
    # make sure that peak intervals have 'chrom' 'start' 'end' and 'const' columns
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

    if (nrow(peak_intervals) != nrow(atac_scores)) {
        cli_abort("{.field 'peak_intervals'} must have the same number of rows as {.field 'atac_scores'}. intervals number of rows {.val {nrow(peak_intervals)}} != atac_scores number of rows {.val {nrow(atac_scores)}}}}", call = parent.frame(1))
    }
}

validate_additional_features <- function(additional_features, peak_intervals) {
    if (!is.null(additional_features)) {
        if (nrow(additional_features) != nrow(peak_intervals)) {
            cli_abort("{.field 'additional_features'} must have the same number of rows as {.field 'peak_intervals'}. intervals number of rows {.val {nrow(peak_intervals)}} != additional_features number of rows {.val {nrow(additional_features)}}}}", call = parent.frame(1))
        }
    }
}


calc_motif_energies <- function(peak_intervals, pssm_db = prego::all_motif_datasets(), motif_energies = NULL, field_name = "motif_energies", intervals_field_name = "peak_intervals") {
    if (is.null(motif_energies)) {
        cli_alert("Computing motif energies (this might take a while)")
        motif_energies <- prego::gextract_pwm(peak_intervals, dataset = pssm_db, prior = 0.01) %>%
            select(-chrom, -start, -end) %>%
            as.matrix()
    }

    if (nrow(motif_energies) != nrow(peak_intervals)) {
        cli_abort("{.field '{field_name}'} must have the same number of rows as {.field '{intervals_field_name}'}. {intervals_field_name} number of rows {.val {nrow(peak_intervals)}} != {field_name} number of rows {.val {nrow(motif_energies)}}}}", call = parent.frame(1))
    }

    missing_motifs <- colnames(motif_energies)[!(colnames(motif_energies) %in% pssm_db$motif)]
    if (length(missing_motifs) > 0) {
        cli_abort("Motifs {.val {missing_motifs}} are missing from the motif database ({.field pssm_db} parameter).", call = parent.frame(1))
    }

    return(motif_energies)
}

get_model_coefs <- function(model) {
    df <- coef(model, s = model$lambda) %>%
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
