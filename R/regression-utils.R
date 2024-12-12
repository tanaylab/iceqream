validate_misha <- function() {
    # check if misha package is loaded
    if (!requireNamespace("misha", quietly = TRUE)) {
        cli_abort("Please install the misha package to use the 'min_tss_distance' parameter")
    }
    # check that an environment called .misha exists
    if (!exists(".misha", envir = .GlobalEnv)) {
        cli_abort("Please load the misha package {.code library(misha)} and set the genome using {.code misha::gsetroot()}")
    }
}

validate_tss_intervals <- function() {
    # Check for required TSS intervals
    if (!misha::gintervals.exists("intervs.global.tss")) {
        misha_root <- .misha$GROOT
        cli_abort("Please make sure the current genome ({.field {misha_root}}) has an intervals set called {.val intervs.global.tss}")
    }
}

get_tss_distance_filter <- function(peak_intervals, min_tss_distance = NULL) {
    # If no minimum TSS distance specified, return all TRUE
    if (is.null(min_tss_distance)) {
        return(rep(TRUE, nrow(peak_intervals)))
    }

    validate_tss_intervals()

    # Calculate distances to TSS
    tss_dist <- abs(misha::gintervals.neighbors(
        as.data.frame(peak_intervals),
        "intervs.global.tss",
        na.if.notfound = TRUE
    )$dist)

    # Create filter based on TSS distance
    enhancers_filter <- tss_dist > min_tss_distance
    enhancers_filter[is.na(enhancers_filter)] <- FALSE

    # Report filtered peaks
    if (sum(!enhancers_filter) > 0) {
        cli_alert_info("{.val {sum(!enhancers_filter)}} peaks were filtered out because they are too close to TSS (<= {.val {min_tss_distance}}bp)")
    }

    return(enhancers_filter)
}

select_motifs_by_correlation <- function(motif_energies, atac_diff, min_initial_energy_cor, max_motif_num, cm = NULL) {
    if (is.null(cm)) {
        cli_alert_info("Calculating correlations between {.val {ncol(motif_energies)}} motif energies and ATAC difference...")
        cm <- tgs_cor(motif_energies, as.matrix(atac_diff), pairwise.complete.obs = TRUE)[, 1]
    }

    cm <- cm[!is.na(cm)]
    motifs <- names(cm[abs(cm) >= min_initial_energy_cor])

    if (length(motifs) < min(max_motif_num, ncol(motif_energies))) {
        cli::cli_alert_warning("No features with absolute correlation >= {.val {min_initial_energy_cor}}. Trying again with {.val {min_initial_energy_cor/2}}")
        min_initial_energy_cor <- min_initial_energy_cor / 2
        motifs <- names(cm[abs(cm) >= min_initial_energy_cor])
        if (length(motifs) == 0) {
            n <- min(min(max_motif_num, ncol(motif_energies)), length(motifs))
            motifs <- names(sort(abs(cm), decreasing = TRUE)[1:n])
            cli::cli_alert_warning("No features with absolute correlation >= {.val {min_initial_energy_cor}}. Taking top {.val {n}} correlated features")
        }
    }

    cli_alert_info("Selected {.val {length(motifs)}} (out of {.val {ncol(motif_energies)}}) features with absolute correlation >= {.val {min_initial_energy_cor}}")
    motifs <- motifs[!is.na(motifs)]

    return(motifs)
}

select_features_by_regression <- function(motif_energies, atac_diff_n, additional_features,
                                          feature_selection_beta, alpha, lambda,
                                          parallel = TRUE, seed = NULL) {
    cli_alert_info("Running first round of regression, # of features: {.val {ncol(motif_energies)}}")
    glm_model1 <- glmnet::glmnet(motif_energies, atac_diff_n,
        binomial(link = "logit"),
        alpha = alpha,
        lambda = lambda,
        parallel = parallel,
        seed = seed
    )
    glm_model1 <- strip_glmnet(glm_model1)

    features <- rownames(glm_model1$beta)[abs(glm_model1$beta[, 1]) >= feature_selection_beta]
    cli_alert_info("Taking {.val {length(features)}} features with beta >= {.val {feature_selection_beta}}")
    if (length(features) == 0) {
        cli::cli_alert_warning("No features with beta >= {.val {feature_selection_beta}}. Using all features.")
        features <- rownames(glm_model1$beta)
    }

    cli_alert_info("Running second round of regression...")
    additional_features[is.na(additional_features)] <- 0
    glm_model2 <- glmnet::glmnet(
        as.matrix(cbind(motif_energies[, features], additional_features)),
        atac_diff_n,
        binomial(link = "logit"),
        alpha = alpha,
        lambda = lambda,
        parallel = parallel,
        seed = seed
    )
    glm_model2 <- strip_glmnet(glm_model2)

    chosen_motifs <- rownames(glm_model2$beta)[abs(glm_model2$beta[, 1]) > 0]
    if (length(chosen_motifs) == 0) {
        cli::cli_alert_warning("No features with beta > 0. Using all features.")
        chosen_motifs <- rownames(glm_model2$beta)
    }

    final_features <- motif_energies[, setdiff(chosen_motifs, colnames(additional_features))]

    return(list(
        features_mat = final_features,
        glm_model2 = glm_model2,
        chosen_motifs = chosen_motifs
    ))
}
