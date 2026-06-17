has_interactions <- function(traj_model) {
    if (!.hasSlot(traj_model, "interactions")) {
        return(FALSE)
    }
    n_interactions(traj_model) > 0
}

n_interactions <- function(traj_model) {
    ncol(traj_model@interactions)
}

create_specific_terms <- function(energies, terms, scale = 1) {
    term1_matrix <- energies[, terms$term1]
    term2_matrix <- energies[, terms$term2]
    inter <- term1_matrix * term2_matrix
    max_vals <- apply(inter, 2, max, na.rm = TRUE)
    max_vals[max_vals == 0] <- 1
    inter <- t(t(inter) / max_vals)
    inter <- apply(inter, 2, norm01) * scale
    colnames(inter) <- terms$variable
    return(inter)
}

create_features_terms <- function(energies, features, data, scale = 1) {
    # Build each anchor's interaction block as a matrix and cbind them, instead
    # of accumulating a column-wise data.frame via purrr::map_dfc and round-
    # tripping through apply(, 2, norm01). For large motif sets the candidate
    # matrix is nrow(peaks) x O(n^2), so the data.frame dance and the per-column
    # apply() copies dominate. This is the matrix-level equivalent and produces
    # identical values (norm01's constant-column case maps to 0 either way).
    # `energies` may arrive as a data.frame (e.g. cbind(matrix, additional_features));
    # matrixStats needs a matrix. Coercing here is equivalent to the old
    # apply()/t() path, which silently coerced too.
    energies <- as.matrix(energies)
    blocks <- lapply(features, function(f) {
        inter <- energies[, setdiff(colnames(energies), f), drop = FALSE] * data[, f]
        max_vals <- matrixStats::colMaxs(inter, na.rm = TRUE)
        max_vals[max_vals == 0] <- 1
        inter <- sweep(inter, 2, max_vals, `/`)
        colnames(inter) <- paste0(f, ":", colnames(inter))
        inter
    })
    interactions <- do.call(cbind, blocks)

    # Column-wise norm01, vectorized (matches norm01(): constant columns -> 0).
    cmin <- matrixStats::colMins(interactions, na.rm = TRUE)
    rng <- matrixStats::colMaxs(interactions, na.rm = TRUE) - cmin
    rng[rng == 0] <- 1
    interactions <- sweep(sweep(interactions, 2, cmin, `-`), 2, rng, `/`) * scale
    interactions
}


create_interaction_terms <- function(energies, motif_feats = NULL, add_feats = NULL, additional_features = NULL, max_motif_n = NULL, max_add_n = NULL, only_sig_motifs = FALSE, only_sig_add_motifs = TRUE, scale = 1) {
    create_interactions <- function(features, data, max_n) {
        if (is.null(features) || is.null(data)) {
            return(NULL)
        }

        features <- head(features, n = max_n %||% length(features))
        interactions <- create_features_terms(data, features, data, scale = scale)

        interactions
    }

    # When only_sig_add_motifs=FALSE, use all additional feature columns as anchors
    if (!only_sig_add_motifs && !is.null(additional_features)) {
        add_feats <- colnames(additional_features)
    }

    if (length(add_feats) > 0) {
        add_inter <- create_interactions(add_feats, cbind(energies, additional_features), max_add_n)

        if (!is.null(add_inter)) {
            cli::cli_alert_info("Created {.val {ncol(add_inter)}} interactions between additional features and motif features.")
        }
    } else {
        add_inter <- NULL
    }

    if (only_sig_motifs) {
        motif_inter <- create_interactions(motif_feats, energies[, motif_feats], max_motif_n)
    } else {
        motif_inter <- create_interactions(motif_feats, energies, max_motif_n)
    }

    if (!is.null(motif_inter)) {
        cli::cli_alert_info("Created {.val {ncol(motif_inter)}} interactions between motif features.")
    }

    interactions <- cbind(motif_inter, add_inter)
    if (!is.null(interactions)) {
        if (!is.null(rownames(energies))) {
            rownames(interactions) <- rownames(energies)
        }
        cli::cli_alert_info("Created {.val {ncol(interactions)}} interactions in total.")
    }

    return(interactions)
}

get_significant_interactions <- function(
  energies, y, interaction_threshold, max_motif_n = NULL, max_add_n = NULL,
  max_n = NULL,
  additional_features = NULL, lambda = 1e-5, alpha = 1, seed = 60427,
  ignore_feats = c("TT", "CT", "GT", "AT", "TC", "CC", "GC", "AC", "TG", "CG", "GG", "AG", "TA", "CA", "GA", "AA"), idxs = NULL, only_sig_motifs = FALSE, only_sig_add_motifs = TRUE, scale = 1, min_signal_correlation = NULL
) {
    if (is.null(idxs)) {
        idxs <- seq_len(nrow(energies))
    }
    # Interactions are pairwise: with fewer than 2 features there are no pairs to
    # form, and the linear pre-selection glmnet() below would itself error on a
    # single-column design ("x should be a matrix with 2 or more columns"). This
    # guards both callers - add_interactions() and the de-novo regression path.
    if (ncol(energies) < 2) {
        cli::cli_alert_warning("Fewer than 2 features available; no interactions can be formed.")
        return(NULL)
    }
    glm_model_lin <- glmnet::glmnet(as.matrix(energies[idxs, , drop = FALSE]), y[idxs], binomial(link = "logit"), alpha = alpha, lambda = lambda, seed = seed)
    glm_model_lin <- strip_glmnet(glm_model_lin)

    feats_all <- abs(stats::coef(glm_model_lin)[-1])
    names(feats_all) <- rownames(stats::coef(glm_model_lin))[-1]
    sig_feats <- names(feats_all)[feats_all > interaction_threshold]
    sig_feats <- setdiff(sig_feats, ignore_feats)

    if (length(sig_feats) == 0) {
        cli::cli_alert_warning("No significant features to consider for interactions.")
        return(NULL)
    }
    additional_features <- additional_features[, setdiff(colnames(additional_features), ignore_feats)]
    add_feats <- intersect(sig_feats, colnames(additional_features))
    add_feats <- setdiff(add_feats, ignore_feats)
    motif_feats <- setdiff(sig_feats, add_feats)

    cli::cli_alert_info("# of significant features to consider for interactions: {.val {length(sig_feats)}} (out of {.val {ncol(energies)}}) above the threshold of {.val {interaction_threshold}}. Of these, {.val {length(motif_feats)}} are motif features and {.val {length(add_feats)}} are additional features.")

    if (!is.null(additional_features)) {
        # remove the features from energies
        energies <- energies[, setdiff(colnames(energies), colnames(additional_features))]

        # remove the ignored features
        energies <- energies[, setdiff(colnames(energies), ignore_feats)]
    }

    inter <- create_interaction_terms(energies,
        motif_feats = motif_feats, add_feats = add_feats,
        additional_features = additional_features, max_motif_n = max_motif_n, max_add_n = max_add_n, only_sig_motifs = only_sig_motifs, only_sig_add_motifs = only_sig_add_motifs, scale = scale
    )

    if (!is.null(max_n) && ncol(inter) > max_n) {
        cli::cli_alert_info("Selecting top {.val {max_n}} interactions based on correlation with the signal.")
        max_n <- min(max_n, ncol(inter))
        cm <- tgs_cor(inter[idxs, ], as.matrix(y[idxs]))[, 1]
        chosen <- names(sort(abs(cm), decreasing = TRUE)[1:max_n])
        inter <- inter[, chosen]
    }

    if (!is.null(min_signal_correlation) && ncol(inter) > 0) {
        cm_final <- abs(tgs_cor(inter[idxs, ], as.matrix(y[idxs]))[, 1])
        if (!any(is.finite(cm_final))) {
            # All interaction-signal correlations are NA/non-finite (e.g. every
            # candidate column is constant). max(., na.rm = TRUE) would be -Inf
            # and silently drop every interaction; skip the filter instead.
            cli::cli_alert_warning(
                "{.field min_signal_correlation}: all interaction-signal correlations are non-finite; skipping the filter."
            )
        } else {
            threshold_cor <- min_signal_correlation * max(cm_final, na.rm = TRUE)
            keep <- names(cm_final)[is.finite(cm_final) & cm_final > threshold_cor]
            n_dropped <- ncol(inter) - length(keep)
            if (n_dropped > 0) {
                cli::cli_alert_info(
                    "Applied {.field min_signal_correlation} = {.val {min_signal_correlation}}: dropped {.val {n_dropped}} of {.val {ncol(inter)}} interactions below {.val {signif(threshold_cor, 3)}}."
                )
            }
            inter <- inter[, keep, drop = FALSE]
        }
    }

    return(inter)
}

#' Add interactions to a trajectory model
#'
#' This function adds significant interactions to a given trajectory model if they do not already exist.
#' It identifies significant interactions based on the provided threshold and updates the model features
#' with logistic features derived from these interactions. The trajectory model is then re-learned with
#' the new features.
#' @param interaction_threshold threshold for the selecting features to create interactions. IQ learns a linear model on the features and selects the features with coefficients above this threshold. Default: 0.001
#' @param max_motif_n maximum number of motifs to consider for interactions. If NULL, all motifs above the interaction_threshold will be considered. Default: NULL
#' @param max_add_n maximum number of additional features to consider for interactions. If NULL, all additional features above the interaction_threshold will be considered. Default: NULL
#' @param max_n maximum number of interactions to consider. If NULL, all interactions will be considered. If set, the interactions will be selected based on correlation with the signal in the training data. Default: (motif models + additional features) * 10
#' @param interactions A precomputed interaction matrix. If provided, the function will not compute the interactions. Default: NULL
#' @param ignore_feats A character vector of features to ignore when creating interactions. Default: dinucleotides
#' @param force If TRUE, the function will add interactions even if they already exist. Default: FALSE
#' @param logist_interactions Logical indicating whether to transform interactions to logistic functions. Default: FALSE
#' @param only_sig_motifs Logical indicating whether to only consider significant motifs for interactions. Default: FALSE
#' @param only_sig_add_motifs Logical indicating whether to restrict additional-feature interaction anchors to the significant additional features (those whose linear-model coefficient exceeds `interaction_threshold`). If FALSE, all additional feature columns are candidate anchors. Default: TRUE
#' @param alpha The elastic net mixing parameter for glmnet. alpha=1 is the lasso penalty, alpha=0 is the ridge penalty. Default: 1
#' @param interaction_scale_factor A multiplier applied to the normalized interaction matrix before it is joined into the model features. Default: 1 (no scaling). Values >1 make interactions more prominent relative to motif features; values <1 down-weight them.
#' @param min_signal_correlation If non-NULL, after the `max_n` correlation-based top-N cap, drop any interaction column whose absolute correlation with the training `diff_score` is below `min_signal_correlation * max(|cor|)`. For example, `min_signal_correlation = 1/8` keeps only interactions whose |cor| exceeds 1/8 of the best interaction's |cor|, mirroring Akhiad's manual post-hoc filter. Default: NULL (no filter).
#'
#' @inheritParams regress_trajectory_motifs
#' @inheritParams relearn_traj_model
#'
#' @return The updated trajectory model with added interactions.
#' @export
add_interactions <- function(traj_model, interaction_threshold = 0.001, max_motif_n = NULL, max_add_n = NULL, max_n = NULL, lambda = 1e-5, alpha = 1, seed = 60427, interactions = NULL, ignore_feats = c("TT", "CT", "GT", "AT", "TC", "CC", "GC", "AC", "TG", "CG", "GG", "AG", "TA", "CA", "GA", "AA"), force = FALSE, logist_interactions = FALSE, use_cv = FALSE, nfolds = 10, family = "binomial", rescale_pred = TRUE, only_sig_motifs = FALSE, only_sig_add_motifs = TRUE, interaction_scale_factor = 1, min_signal_correlation = NULL) {
    r2_all_before <- cor(traj_model@diff_score, traj_model@predicted_diff_score)^2
    if (traj_model_has_test(traj_model)) {
        r2_train_before <- cor(traj_model@diff_score[traj_model@type == "train"], traj_model@predicted_diff_score[traj_model@type == "train"])^2
        r2_test_before <- cor(traj_model@diff_score[traj_model@type == "test"], traj_model@predicted_diff_score[traj_model@type == "test"])^2
    }

    n_feats <- length(setdiff(colnames(traj_model@additional_features), ignore_feats))
    if (is.null(max_n)) {
        max_n <- (ncol(traj_model@normalized_energies) + n_feats) * 10
        cli::cli_alert("Setting {.field max_n} (maximal number of interactions) to {.val {max_n}}.")
    }

    if (has_interactions(traj_model)) {
        if (force) {
            traj_model <- remove_interactions(traj_model)
        } else {
            cli::cli_alert_warning("Interactions already exist. Set {.field force} to TRUE to overwrite.")
            return(traj_model)
        }
    }

    # Interactions are pairwise, so they need at least 2 candidate features
    # (motifs + non-ignored additional features). A 0/1-feature model has no
    # pairs to form, and the linear pre-selection (glmnet) would itself error on
    # a single-column matrix. Return the model unchanged in that case.
    if (is.null(interactions) && (ncol(traj_model@normalized_energies) + n_feats) < 2) {
        cli::cli_alert_warning("Too few features to form interactions; returning the model unchanged.")
        return(traj_model)
    }

    if (is.null(interactions)) {
        cli::cli_alert("Adding interactions")
        interactions <- get_significant_interactions(
            cbind(traj_model@normalized_energies, traj_model@additional_features), norm01(traj_model@diff_score), interaction_threshold,
            max_motif_n = max_motif_n, max_add_n = max_add_n,
            max_n = max_n,
            additional_features = traj_model@additional_features, lambda = lambda, alpha = alpha, seed = seed, ignore_feats = ignore_feats, idxs = which(traj_model@type == "train"), only_sig_motifs = only_sig_motifs, only_sig_add_motifs = only_sig_add_motifs, scale = interaction_scale_factor, min_signal_correlation = min_signal_correlation
        )
    }

    if (!is.null(interactions)) {
        traj_model@interactions <- interactions
    }

    if (logist_interactions) {
        logist_inter <- create_logist_features(interactions)
        traj_model@model_features <- cbind(traj_model@model_features, logist_inter)
        cli::cli_alert_info("Re-learning the model with the new interactions. Number of features: {.val {ncol(traj_model@model_features)}}, out of which {.val {ncol(traj_model@normalized_energies)}}*4 are motif features, {.val {ncol(traj_model@additional_features)}}*4 are additional features and {.val {ncol(traj_model@interactions)}}*4 are interactions.")
    } else {
        traj_model@model_features <- cbind(traj_model@model_features, interactions)
        cli::cli_alert_info("Re-learning the model with the new interactions. Number of features: {.val {ncol(traj_model@model_features)}}, out of which {.val {ncol(traj_model@normalized_energies)}}*4 are motif features, {.val {ncol(traj_model@additional_features)}}*4 are additional features and {.val {ncol(traj_model@interactions)}} are interactions.")
    }
    traj_model@params$logist_interactions <- logist_interactions


    traj_model <- relearn_traj_model(traj_model, new_energies = FALSE, new_logist = FALSE, use_additional_features = TRUE, use_motifs = TRUE, verbose = FALSE, logist_interactions = logist_interactions, use_cv = use_cv, nfolds = nfolds, lambda = lambda, family = family, rescale_pred = rescale_pred)
    if (traj_model_has_test(traj_model)) {
        r2_train_after <- cor(traj_model@diff_score[traj_model@type == "train"], traj_model@predicted_diff_score[traj_model@type == "train"])^2
        cli::cli_alert_info("R^2 train before: {.val {r2_train_before}}, after: {.val {r2_train_after}}, change: {.val {r2_train_after - r2_train_before}}")
        r2_test_after <- cor(traj_model@diff_score[traj_model@type == "test"], traj_model@predicted_diff_score[traj_model@type == "test"])^2
        cli::cli_alert_info("R^2 test before: {.val {r2_test_before}}, after: {.val {r2_test_after}}, change: {.val {r2_test_after - r2_test_before}}")
    } else {
        r2_after <- cor(traj_model@diff_score, traj_model@predicted_diff_score)^2
        cli::cli_alert_info("R^2 all before: {.val {r2_all_before}}, after: {.val {r2_after}}, change: {.val {r2_after - r2_all_before}}")
    }


    return(traj_model)
}

#' Build base-score / end-score engineered additional features
#'
#' Used by [add_interactions_progressive()] between passes. Relearns the
#' trajectory model twice — once with each end of `atac_scores` as the
#' response — and returns a data frame with three engineered columns
#' suitable for `cbind()` onto `@additional_features`:
#'
#' * `base_pred`: prediction under the bin_start score, scaled to `[0, 10]`.
#' * `end_pred`: prediction under the bin_end score, scaled to `[0, 10]`.
#' * `pred_diff_e_b`: `end_pred - base_pred` scaled to `[0, 10]`.
#'
#' These features mirror the manual pattern in Akhiad's analysis
#' workflow and are what the second pass of
#' [add_interactions_progressive()] anchors interactions against.
#'
#' Requires `atac_scores` with at least two distinct columns. If the
#' inputs don't support it, the caller should fall back to a single
#' [add_interactions()] call.
#'
#' @section Caveat — test-time leakage:
#' Using this function as the `additional_features_builder` of
#' [add_interactions_progressive()] injects `base_pred` / `end_pred` /
#' `pred_diff_e_b` columns that are defined only for the training peaks
#' the `traj_model` was fit on. If you later call
#' [infer_trajectory_motifs()] on test peaks, those columns are imputed
#' to 0 for the test rows, which drops test R^2 severely (measured:
#' −0.27 on the gastrulation vignette). To use this builder correctly
#' you must independently run the base-only / end-only helper models
#' on the test peaks and supply the resulting predictions as
#' `additional_features` at [infer_trajectory_motifs()] time.
#'
#' Because this pitfall is easy to hit, `iq_regression(strategy =
#' "progressive")` errors at call time — it reserves the progressive
#' path for a future release that threads the test-time propagation
#' through automatically. Use this builder only when calling
#' [add_interactions_progressive()] directly with full control over
#' test-time inference.
#'
#' @param traj_model A `TrajectoryModel` already fit with a first pass
#'   of interactions (so the relearns have something meaningful to fit).
#' @param atac_scores A data frame / matrix of per-peak ATAC scores with
#'   one column per bin, rows aligned with `traj_model@peak_intervals`.
#' @param bin_start,bin_end Column index (integer) or name (character) of
#'   the start and end bins in `atac_scores`.
#'
#' @return A data frame with columns `base_pred`, `end_pred`,
#'   `pred_diff_e_b`, row-aligned with `traj_model@peak_intervals`.
#' @export
default_score_split_features <- function(traj_model, atac_scores, bin_start, bin_end) {
    atac_scores <- as.data.frame(atac_scores)
    if (nrow(atac_scores) != nrow(traj_model@peak_intervals)) {
        cli::cli_abort(
            "{.field atac_scores} has {.val {nrow(atac_scores)}} rows but {.field traj_model} has {.val {nrow(traj_model@peak_intervals)}} peaks."
        )
    }
    if (is.character(bin_start)) bin_start <- match(bin_start, colnames(atac_scores))
    if (is.character(bin_end)) bin_end <- match(bin_end, colnames(atac_scores))
    if (is.na(bin_start) || is.na(bin_end) || bin_start == bin_end) {
        cli::cli_abort("{.field bin_start} and {.field bin_end} must resolve to two distinct columns of {.field atac_scores}.")
    }

    base_score <- atac_scores[, bin_start]
    end_score <- atac_scores[, bin_end]

    tm_base <- suppressMessages(
        relearn_traj_model(traj_model, new_energies = FALSE, new_logist = FALSE,
            verbose = FALSE, new_score = base_score, rescale_pred = FALSE)
    )
    tm_end <- suppressMessages(
        relearn_traj_model(traj_model, new_energies = FALSE, new_logist = FALSE,
            verbose = FALSE, new_score = end_score, rescale_pred = FALSE)
    )

    base_pred <- tm_base@predicted_diff_score
    end_pred <- tm_end@predicted_diff_score

    data.frame(
        base_pred = norm01(base_pred) * 10,
        end_pred = norm01(end_pred) * 10,
        pred_diff_e_b = norm01(end_pred - base_pred) * 10,
        row.names = rownames(traj_model@normalized_energies)
    )
}

#' Progressive two-pass interaction selection
#'
#' Implements the two-pass interaction-selection pattern used manually in
#' Akhiad's analysis workflow:
#'
#' 1. First pass: tight `interaction_threshold[1]`, typically with
#'    `only_sig_motifs = TRUE` — a strict selection to seed the model with
#'    high-confidence interactions.
#' 2. (Optional) Call `additional_features_builder(traj_model)` to produce
#'    engineered additional features (e.g. `base_pred`, `end_pred` from
#'    [default_score_split_features()]) and inject them. Relearn once so
#'    the next pass sees an enriched anchor set.
#' 3. Second pass: looser `interaction_threshold[2]`, typically with
#'    `only_sig_motifs = FALSE` and `force = TRUE` — broader selection
#'    that can pick interactions anchored on the newly-injected features.
#'
#' Without a feature builder, the two passes are equivalent to a single
#' [add_interactions()] call at the final (loosest) threshold, so the
#' primary use case is with `default_score_split_features` or a custom
#' builder. For data without multi-bin `atac_scores`, prefer a single
#' [add_interactions()] call with `interaction_threshold = thresholds[1]`.
#'
#' @param traj_model A `TrajectoryModel`.
#' @param thresholds A numeric vector of `interaction_threshold` values,
#'   one per pass. Default `c(0.01, 0.0005)` mirrors Akhiad's pattern.
#' @param only_sig_motifs A logical vector matched to `thresholds`.
#'   Default `c(TRUE, FALSE)`.
#' @param only_sig_add_motifs A logical vector matched to `thresholds`.
#'   Default `c(TRUE, TRUE)`.
#' @param additional_features_builder NULL or a function taking the
#'   traj_model after pass 1 and returning a data frame to cbind onto
#'   `@additional_features`. Use [default_score_split_features()] as a
#'   ready-made builder when `atac_scores` is available.
#' @param interaction_scale_factor,min_signal_correlation Forwarded to
#'   each [add_interactions()] call.
#' @param seed Integer seed forwarded to [add_interactions()].
#' @param ... Additional arguments forwarded to [add_interactions()]
#'   (e.g. `max_motif_n`, `max_add_n`, `max_n`, `logist_interactions`).
#'
#' @return The updated trajectory model with interactions and any
#'   engineered additional features from the builder.
#' @export
add_interactions_progressive <- function(
    traj_model,
    thresholds = c(0.01, 0.0005),
    only_sig_motifs = c(TRUE, FALSE),
    only_sig_add_motifs = c(TRUE, TRUE),
    additional_features_builder = NULL,
    interaction_scale_factor = 1,
    min_signal_correlation = NULL,
    seed = 60427,
    ...
) {
    if (length(only_sig_motifs) == 1) {
        only_sig_motifs <- rep(only_sig_motifs, length(thresholds))
    }
    if (length(only_sig_add_motifs) == 1) {
        only_sig_add_motifs <- rep(only_sig_add_motifs, length(thresholds))
    }
    if (length(thresholds) != length(only_sig_motifs) ||
        length(thresholds) != length(only_sig_add_motifs)) {
        cli::cli_abort(
            "{.field thresholds}, {.field only_sig_motifs}, and {.field only_sig_add_motifs} must have the same length."
        )
    }
    if (length(thresholds) < 1) {
        cli::cli_abort("{.field thresholds} must have at least one value.")
    }

    # Guard against ... collisions with the args we set explicitly on each
    # inner add_interactions() call. These two are the only args forwarded
    # via ... that we override per-pass — the others (only_sig_motifs,
    # interaction_scale_factor, seed, etc.) are formal args of this
    # function, so R matches them to the formal slot before ... can
    # capture them. interaction_threshold and force are NOT formal args
    # of progressive, so a user could pass them via ... and hit an
    # ambiguous "matched by multiple actual arguments" error.
    reserved_via_dots <- c("interaction_threshold", "force")
    extra_args <- list(...)
    clash <- intersect(names(extra_args), reserved_via_dots)
    if (length(clash) > 0) {
        cli::cli_abort(c(
            "{.arg {clash}} cannot be passed to {.fn add_interactions_progressive}.",
            "i" = "{.arg interaction_threshold} is controlled per-pass via the {.arg thresholds} argument. {.arg force} is set internally ({.val FALSE} on pass 1, {.val TRUE} on passes 2+)."
        ))
    }

    # Single-threshold: identical to a single add_interactions call.
    if (length(thresholds) == 1) {
        return(add_interactions(
            traj_model,
            interaction_threshold = thresholds[1],
            only_sig_motifs = only_sig_motifs[1],
            only_sig_add_motifs = only_sig_add_motifs[1],
            interaction_scale_factor = interaction_scale_factor,
            min_signal_correlation = min_signal_correlation,
            seed = seed,
            ...
        ))
    }

    cli::cli_alert("Progressive interactions: pass 1 of {.val {length(thresholds)}} at threshold {.val {thresholds[1]}} ({.field only_sig_motifs} = {.val {only_sig_motifs[1]}})")
    traj_model <- add_interactions(
        traj_model,
        interaction_threshold = thresholds[1],
        only_sig_motifs = only_sig_motifs[1],
        only_sig_add_motifs = only_sig_add_motifs[1],
        interaction_scale_factor = interaction_scale_factor,
        min_signal_correlation = min_signal_correlation,
        seed = seed,
        force = FALSE,
        ...
    )

    if (!is.null(additional_features_builder)) {
        cli::cli_alert("Building engineered additional features between passes")
        new_feats <- additional_features_builder(traj_model)
        if (!is.null(new_feats) && nrow(new_feats) > 0 && ncol(new_feats) > 0) {
            existing <- traj_model@additional_features
            if (nrow(existing) == 0) {
                existing <- data.frame(row.names = rownames(traj_model@normalized_energies))
            }
            # Drop any already-present columns so we don't double-register.
            dup_cols <- intersect(colnames(new_feats), colnames(existing))
            if (length(dup_cols) > 0) {
                existing <- existing[, setdiff(colnames(existing), dup_cols), drop = FALSE]
            }
            traj_model@additional_features <- cbind(existing, new_feats)
            # Relearn so model_features reflects the new additional features.
            # logist_dinucs = TRUE matches how infer_trajectory_motifs expands
            # features at inference time (create_logist_features on every
            # additional-feature column, including dinucleotides). Without
            # this, the between-pass relearn leaves raw "TT"/"CT"/... columns
            # in @model_features while inference emits "TT_low-energy"/..., and
            # infer_trajectory_motifs aborts with a Missing-columns error.
            traj_model <- suppressMessages(
                relearn_traj_model(
                    traj_model,
                    new_energies = FALSE,
                    new_logist = TRUE,
                    logist_dinucs = TRUE,
                    logist_interactions = isTRUE(traj_model@params$logist_interactions),
                    verbose = FALSE
                )
            )
        }
    }

    for (k in seq.int(2, length(thresholds))) {
        cli::cli_alert("Progressive interactions: pass {.val {k}} of {.val {length(thresholds)}} at threshold {.val {thresholds[k]}} ({.field only_sig_motifs} = {.val {only_sig_motifs[k]}})")
        traj_model <- add_interactions(
            traj_model,
            interaction_threshold = thresholds[k],
            only_sig_motifs = only_sig_motifs[k],
            only_sig_add_motifs = only_sig_add_motifs[k],
            interaction_scale_factor = interaction_scale_factor,
            min_signal_correlation = min_signal_correlation,
            seed = seed,
            force = TRUE,
            ...
        )
    }

    traj_model
}

remove_interactions <- function(traj_model) {
    if (has_interactions(traj_model)) {
        if (is.null(traj_model@params$logist_interactions) || !traj_model@params$logist_interactions) {
            inter_terms <- colnames(traj_model@interactions)
        } else {
            inter_terms <- purrr::map(c("low-energy", "high-energy", "higher-energy", "sigmoid"), ~ {
                paste0(colnames(traj_model@interactions), "_", .x)
            }) %>% do.call(c, .)
        }
        traj_model@model_features <- traj_model@model_features[, !(colnames(traj_model@model_features) %in% inter_terms)]
        traj_model@interactions <- matrix(nrow = 0, ncol = 0)
        traj_model@params$logist_interactions <- FALSE
    }
    return(traj_model)
}

