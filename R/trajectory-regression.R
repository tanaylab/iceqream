#' TrajectoryModel Class
#'
#' An S4 class to represent a trajectory model.
#'
#' @slot diff_score numeric
#'   A numerical value representing the difference score.
#'
#' @slot predicted_diff_score numeric
#'   A numerical value representing the predicted difference score.
#'
#' @slot model ANY
#'   The model object or any data type that represents the trajectory model.
#'
#' @slot coefs data.frame
#'  A data frame containing the coefficients of the model (not including the intercept). Contains
#'  the coefficients for the 'early', 'linear' and 'late' models
#'
#' @slot normalized_energies matrix
#'   A matrix containing normalized energies.
#'
#' @slot type
#'  A vector the length of the number of peaks, indicating whether each peak is a training ('train') or
#' a prediction peak ('test').
#'
#' @slot motif_models list
#'   A list of models representing different motifs.
#'
#' @slot initial_prego_models list
#'   A list of initial pre-go models.
#'
#' @slot peak_intervals data.frame
#'  A data frame containing the peak intervals.
#'
#' @slot params list
#'  A list of parameters used for training.
#'
#' @slot additional_features data.frame
#'  A data frame containing the additional features.
#'
#'
#'
#' @exportClass TrajectoryModel
TrajectoryModel <- setClass(
    "TrajectoryModel",
    slots = list(
        model = "ANY",
        motif_models = "list",
        coefs = "data.frame",
        normalized_energies = "matrix",
        type = "character",
        diff_score = "numeric",
        predicted_diff_score = "numeric",
        initial_prego_models = "list",
        peak_intervals = "data.frame",
        additional_features = "data.frame",
        params = "list"
    )
)

#' @export
#' @rdname TrajectoryModel-class
setMethod("show", signature = "TrajectoryModel", definition = function(object) {
    cli::cli({
        cli::cli_text("{.cls TrajectoryModel} with {.val {length(object@motif_models)}} motifs and {.val {length(object@additional_features)}} additional features\n")
        cli::cli_text("\n")
        cli::cli_text("Slots include:")
        cli_ul(c("{.field @model}: A GLM model object. Number of non-zero coefficients: {.val {sum(object@model$beta[, 1] != 0)}}"))
        cli_ul(c("{.field @motif_models}: A named list of motif models. Each element contains PSSM and spatial model ({.val {length(object@motif_models)}} models: {.val {names(object@motif_models)}})"))
        cli_ul(c("{.field @additional_features}: A data frame of additional features ({.val {nrow(object@additional_features)}} elements)"))
        cli_ul(c("{.field @coefs}: A data frame of coefficients ({.val {nrow(object@coefs)}} elements)"))
        cli_ul(c("{.field @normalized_energies}: A matrix of normalized energies of the model features (logistic functions of the motif models energies, dimensions: {.val {nrow(object@normalized_energies)}}x{.val {ncol(object@normalized_energies)}})"))
        cli_ul(c("{.field @type}: A vector the length of the number of peaks, indicating whether each peak is a training ('train') or a prediction peak ('test')"))
        cli_ul(c("{.field @diff_score}: A numeric value representing the difference score the model was trained on ({.val {length(object@normalized_energies[,1])}} elements)"))
        if (any(object@type == "test")) {
            cli_ul(c("{.field @predicted_diff_score}: A numeric value representing the predicted difference score (R^2 train: {.val {round(cor(object@diff_score[object@type == 'train'], object@predicted_diff_score[object@type == 'train'])^2, digits = 3)}}; R^2 test: {.val {round(cor(object@diff_score[object@type == 'test'], object@predicted_diff_score[object@type == 'test'])^2, digits = 3)}})"))
        } else {
            cli_ul(c("{.field @predicted_diff_score}: A numeric value representing the predicted difference score (R^2 (train): {.val {round(cor(object@diff_score, object@predicted_diff_score)^2, digits = 3)}})"))
        }
        cli_ul(c("{.field @initial_prego_models}: A list of prego models used in the initial phase of the algorithm ({.val {length(object@initial_prego_models)}} models)"))
        cli_ul(c("{.field @peak_intervals}: A data frame containing the peak intervals ({.val {nrow(object@peak_intervals)}} elements)"))
        cli_ul(c("{.field @params}: A list of parameters used for training (including: {.val {names(object@params)}})"))

        cli::cli_text("\n")
        cli::cli_text("Run {.code predict(object, peak_intervals)} to predict the model on new data.")
    })
})

#' Predict TrajectoryModel on new data
#'
#' Computes the predicted differential accessibility score between the start and end of the trajectory.
#'
#' @param object An instance of `TrajectoryModel`.
#' @param peak_intervals data frame, indicating the genomic positions ('chrom', 'start', 'end') of each peak to predict.
#'
#' @return A numeric vector of predicted differential accessibility scores.
#'
#' @inheritParams regress_trajectory_motifs
#' @export
setMethod("predict", signature = "TrajectoryModel", definition = function(object, peak_intervals, additional_features = NULL) {
    traj_model <- infer_trajectory_motifs(object, peak_intervals, additional_features = additional_features)
    return(traj_model@predicted_diff_score[traj_model@type == "test"])
})

validate_traj_model <- function(object) {
    if (!methods::is(object, "TrajectoryModel")) {
        cli_abort("Please provide an instance of {.cls TrajectoryModel}", call = parent.frame(1))
    }
}


#' Perform motif regression on ATAC trajectories.
#'
#' @description This function performs motif regression on ATAC trajectories. Given ATAC scores on trajectory bins, it predicts the differential accessibility between the start and end of the trajectory and returns the motifs that are most likely to be responsible for this differential accessibility.
#'
#' @param atac_scores A numeric matrix, representing mean ATAC score per bin per peak. Rows: peaks, columns: bins.
#' @param peak_intervals A data frame, indicating the genomic positions ('chrom', 'start', 'end') of each peak, with an additional column named "const" indicating whether the peak is constitutive.
#' @param motif_energies A numeric matrix, representing the energy of each motif in each peak. If NULL, the function will use \code{pssm_db} to calculate the motif energies. Note that this might take a while.
#' @param pssm_db a data frame with PSSMs ('A', 'C', 'G' and 'T' columns), with an additional column 'motif' containing the motif name. All the motifs in \code{motif_energies} (column names) should be present in the 'motif' column. Default: all motifs in the prego package.
#' @param additional_features A data frame, representing additional genomic features (e.g. CpG content, distance to TSS, etc.) for each peak. Note that NA values would be replaced with 0.
#' @param max_motif_num maximum number of motifs to consider. Default: 50
#' @param min_tss_distance distance from Transcription Start Site (TSS) to classify a peak as an enhancer. Default: 5000. If NULL, no filtering will be performed - use this option if your peaks are already filtered. \cr
#' Note that in order to filter peaks that are too close to TSS, the current \code{misha} genome must have an intervals set called \code{intervs.global.tss}.
#' @param bin_start the start of the trajectory. Default: 1
#' @param bin_end the end of the trajectory. Default: the last bin
#' @param normalize_energies whether to normalize the motif energies. Set this to FALSE if the motif energies are already normalized.
#' @param min_initial_energy_cor minimal correlation between the motif normalized energy and the ATAC difference.
#' @param energy_norm_quantile quantile of the energy used for normalization. Default: 1
#' @param n_prego_motifs number of prego motifs to consider.
#' @param traj_prego output of \code{infer_traj_prego}. If provided, no additional prego models would be inferred.
#' @param min_diff minimal ATAC difference for a peak to participate in the initial prego motif inference.
#' @param prego_sample_fraction Fraction of peaks to sample for prego motif inference. A smaller number would be faster but might lead to over-fitting. Default: 0.1
#' @param seed random seed for reproducibility.
#' @param feature_selection_beta beta parameter used for feature selection.
#' @param parallel whether to use parallel processing on glmnet.
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
                                      max_motif_num = 50,
                                      motif_energies = NULL,
                                      pssm_db = prego::all_motif_datasets(),
                                      additional_features = NULL,
                                      min_tss_distance = 5000,
                                      bin_start = 1,
                                      bin_end = ncol(atac_scores),
                                      min_initial_energy_cor = 0.05,
                                      normalize_energies = TRUE,
                                      energy_norm_quantile = 1,
                                      n_prego_motifs = 4,
                                      traj_prego = NULL,
                                      min_diff = 0.2,
                                      prego_sample_fraction = 0.1,
                                      seed = 60427,
                                      feature_selection_beta = 0.003,
                                      lambda = 1e-5,
                                      alpha = 1,
                                      parallel = TRUE) {
    withr::local_options(list(gmax.data.size = 1e9))
    atac_scores <- as.matrix(atac_scores)

    validate_atac_scores(atac_scores, bin_start, bin_end)
    validate_peak_intervals(peak_intervals, atac_scores)
    validate_additional_features(additional_features, peak_intervals)

    # Extract features
    motif_energies <- calc_motif_energies(peak_intervals, pssm_db, motif_energies)

    cli_alert_info("Normalizing motif energies...")
    if (normalize_energies) {
        motif_energies <- apply(motif_energies, 2, norm_energy, min_energy = -10, q = energy_norm_quantile)
    }

    # filter peaks that are too close to TSS
    if (!is.null(min_tss_distance)) {
        if (!misha::gintervals.exists("intervs.global.tss")) {
            cli_abort("Please make sure the current genome ({.field {GROOT}}) has an intervals set called {.val intervs.global.tss}")
        }

        tss_dist <- abs(misha::gintervals.neighbors(peak_intervals, "intervs.global.tss")$dist)
        enhancers_filter <- tss_dist > min_tss_distance
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
        enhancers_filter <- rep(TRUE, nrow(peak_intervals))
    }

    cli_alert_info("Number of peaks: {.val {nrow(peak_intervals)}}")

    # calculate differential accessibility
    atac_diff <- atac_scores[, bin_end] - atac_scores[, bin_start]
    atac_diff_n <- norm01(atac_diff)

    cli_alert("Extracting sequences...")
    all_seqs <- toupper(misha::gseq.extract(peak_intervals_all))


    if (is.null(traj_prego) && n_prego_motifs > 0) {
        traj_prego <- infer_traj_prego(peak_intervals, atac_diff, n_motifs = n_prego_motifs, min_diff = min_diff, sample_fraction = prego_sample_fraction, energy_norm_quantile = energy_norm_quantile, sequences = all_seqs, seed = seed)
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

    cli_alert_info("Selected {.val {length(motifs)}} (out of {.val {ncol(motif_energies)}}) features with absolute correlation >= {.val {min_initial_energy_cor}}")
    motif_energies <- motif_energies[, motifs]

    cli_alert_info("Running first round of regression, # of features: {.val {ncol(motif_energies)}}")
    glm_model1 <- glmnet::glmnet(motif_energies, atac_diff_n, binomial(link = "logit"), alpha = alpha, lambda = lambda, parallel = parallel, seed = seed)

    features <- rownames(glm_model1$beta)[abs(glm_model1$beta[, 1]) >= feature_selection_beta]
    cli_alert_info("Taking {.val {length(features)}} features with beta >= {.val {feature_selection_beta}}")

    cli_alert_info("Running second round of regression...")
    glm_model2 <- glmnet::glmnet(motif_energies[, features], atac_diff_n, binomial(link = "logit"), alpha = alpha, lambda = lambda, parallel = parallel, seed = seed)

    chosen_motifs <- rownames(glm_model2$beta)[abs(glm_model2$beta[, 1]) > 0]
    features <- motif_energies[, chosen_motifs]

    nclust <- min(ncol(features), max_motif_num * 2)
    cli_alert_info("Clustering {.val {ncol(features)}} features into {.val {nclust}} clusters...")
    features_cm <- tgs_cor(features, pairwise.complete.obs = TRUE)
    hc <- tgs_dist(features_cm) %>% hclust(method = "ward.D2")
    clust_map <- data.frame(feat = colnames(features), clust = cutree(hc, k = nclust), beta = glm_model2$beta[chosen_motifs, 1])

    cli_alert_info("Choosing top {.val {max_motif_num}} features clusters")
    selected_clust <- clust_map %>%
        group_by(clust) %>%
        summarize(max_beta = max(beta)) %>%
        arrange(desc(max_beta)) %>%
        mutate(new_clust = 1:n()) %>%
        head(max_motif_num)

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
            model = glm_model2,
            y = atac_diff_n,
            feats = features,
            clust_motifs = clust_map %>% filter(clust == x$clust) %>% pull(feat),
            sequences = all_seqs[enhancers_filter],
            lambda = lambda,
            pssm_db = pssm_db,
            spat_db = prego_models %>% purrr::map("spat"),
            seed = seed
        )
    }, .parallel = TRUE)
    names(best_motifs_prego) <- best_clust_map$feat


    cli_alert_info("Infering energies...")
    clust_energies <- plyr::llply(purrr::discard(best_motifs_prego, is.null), function(x) {
        prego::compute_pwm(all_seqs[enhancers_filter], x$pssm, spat = x$spat, spat_min = x$spat_min %||% 1, spat_max = x$spat_max)
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

    clust_energies_logist <- create_logist_features(clust_energies)

    model <- glmnet::glmnet(clust_energies_logist, atac_diff_n, binomial(link = "logit"), alpha = alpha, lambda = lambda, parallel = parallel, seed = seed)

    predicted_diff_score <- logist(glmnet::predict.glmnet(model, newx = clust_energies_logist, type = "link", s = lambda))[, 1]
    predicted_diff_score <- (predicted_diff_score * max(atac_diff)) + min(atac_diff)

    cli_alert_success("Finished running model. Number of non-zero coefficients: {.val {sum(model$beta != 0)}} (out of {.val {ncol(clust_energies_logist)}}). R^2: {.val {cor(predicted_diff_score, atac_diff_n)^2}}")

    return(TrajectoryModel(
        model = model,
        motif_models = best_motifs_prego,
        coefs = get_model_coefs(model),
        normalized_energies = clust_energies_logist,
        type = rep("train", nrow(atac_scores)),
        additional_features = as.data.frame(additional_features),
        diff_score = atac_diff,
        predicted_diff_score = predicted_diff_score,
        initial_prego_models = prego_models,
        peak_intervals = peak_intervals,
        params = list(
            energy_norm_quantile = energy_norm_quantile
        )
    ))
}


run_prego_on_clust_residuals <- function(motif, model, y, feats, clust_motifs, sequences, pssm_db, spat_db = NULL, lambda = 1e-5, seed = 60427) {
    cli_alert("Running {.field prego} on cluster {.val {motif}}...")
    pssm <- pssm_db %>%
        filter(motif == !!motif)

    if (nrow(pssm) == 0) { # Motif is part of additional features
        return(NULL)
    }

    if (length(clust_motifs) == 1) {
        pssm <- pssm_db %>% filter(motif == !!motif)

        if (!is.null(spat_db)) {
            spat <- spat_db[[motif]]
        } else {
            spat <- NULL
        }
        return(list(pssm = pssm, spat = spat))
    }

    partial_y <- (feats[, clust_motifs, drop = FALSE] %*% coef(model, s = lambda)[clust_motifs, , drop = FALSE])[, 1]
    cli::cli_fmt(prego_model <- prego::regress_pwm(sequences = sequences, response = partial_y, motif = pssm, seed = seed, match_with_db = FALSE, screen_db = FALSE))
    cli::cli_alert_success("Finished running {.field prego} on cluster {.val {motif}}")
    return(prego::export_regression_model(prego_model))
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


calc_motif_energies <- function(peak_intervals, pssm_db = prego::all_motif_datasets(), motif_energies = NULL) {
    if (is.null(motif_energies)) {
        cli_alert("Computing motif energies (this might take a while)")
        motif_energies <- prego::gextract_pwm(peak_intervals, dataset = pssm_db, prior = 0.01) %>%
            select(-chrom, -start, -end) %>%
            as.matrix()
    }

    if (nrow(motif_energies) != nrow(peak_intervals)) {
        cli_abort("{.field 'motif_energies'} must have the same number of rows as {.field 'peak_intervals'}. intervals number of rows {.val {nrow(peak_intervals)}} != motif_energies number of rows {.val {nrow(motif_energies)}}}}", call = parent.frame(1))
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

    df <- df %>%
        mutate(type = sub(".*_", "", variable), variable = sub("_(early|late|linear)$", "", variable)) %>%
        spread(type, s1)

    df[is.na(df)] <- 0

    return(df)
}
