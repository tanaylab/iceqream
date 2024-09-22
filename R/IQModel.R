#' IQmodel class
#'
#' This class represents an IQ model including IQ features and PBMs,
#' along with necessary model parameters.
#'
#' @slot iq_features A list of IQFeature objects.
#' @slot pbms A list of PBM objects.
#' @slot intercept The intercept term from the model.
#' @slot func The function used for computing energies (e.g., "logSumExp" or "max").
#' @slot lambda The regularization parameter.
#' @slot alpha The elastic net mixing parameter.
#' @slot min_pred Minimum value of training predictions before normalization.
#' @slot max_pred Maximum value of training predictions before normalization.
#' @slot norm_factors Normalization factors for scaling predictions.
#'
#' @exportClass IQmodel
IQmodel <- setClass(
    "IQmodel",
    slots = list(
        iq_features = "list",
        pbms = "list",
        intercept = "numeric",
        func = "character",
        lambda = "numeric",
        alpha = "numeric",
        min_pred = "numeric",
        max_pred = "numeric",
        norm_factors = "numeric"
    )
)

#' Show method for IQmodel
#'
#' This method defines how an IQmodel object should be displayed.
#'
#' @param object An IQmodel object
#'
#' @rdname IQmodel-class
#' @exportMethod show
setMethod("show", signature = "IQmodel", function(object) {
    cli::cli({
        cli::cli_text("An {.cls IQmodel} object")
        n_groups <- sum(purrr::map_lgl(object@iq_features, ~ inherits(.x, "IQFeatureGroup")))
        if (n_groups > 0) {
            cli::cli_text("Contains {.val {length(object@pbms)}} PBMs ({.code @pbms}), and {.val {length(object@iq_features)}} IQ features ({.code @iq_features}), including {.val {n_groups}} feature groups")
        } else {
            cli::cli_text("Contains {.val {length(object@pbms)}} PBMs ({.code @pbms}), and {.val {length(object@iq_features)}} IQ features ({.code @iq_features})")
        }

        cli::cli_text("Intercept: {.val {object@intercept}} ({.code @intercept})")
        cli::cli_text("Energy computation function: {.val {object@func}} ({.code @func})")
        cli::cli_text("Regularization parameter: {.val {object@lambda}} ({.code @lambda})")
        cli::cli_text("Elastic net mixing parameter: {.val {object@alpha}} ({.code @alpha})")
        cli::cli_text("Min prediction value: {.val {object@min_pred}} ({.code @min_pred})")
        cli::cli_text("Max prediction value: {.val {object@max_pred}} ({.code @max_pred})")
        cli::cli_text("Normalization factors: {.val {object@norm_factors}} ({.code @norm_factors})")

        cli::cli_text("\n")
        if (length(object@iq_features) > 0) {
            cli::cli_text("Run {.code predict(object, intervals = new_intervals, new_data = new_data)} to make predictions using this model.")
        } else {
            cli::cli_text("Run {.code predict(object, intervals = new_intervals)} to make predictions using this model.")
        }
    })
})

#' Create an IQmodel from a trajectory model
#'
#' This function creates an IQmodel object from a trajectory model.
#'
#' @param traj_model A trajectory model object
#' @param func The function to use for computing energies (default is "logSumExp")
#'
#' @return An IQmodel object
#'
#' @export
create_iq_model <- function(traj_model, func = "logSumExp") {
    traj_model <- split_traj_model_to_train_test(traj_model)$train
    iq_features <- traj_model_to_iq_feature_list(traj_model)
    pbms <- traj_model_to_pbm_list(traj_model, func = func)
    intercept <- coef(traj_model@model, s = traj_model@params$lambda)["(Intercept)", ]
    lambda <- traj_model@params$lambda
    alpha <- traj_model@params$alpha

    # Calculate min and max prediction values
    feats <- traj_model@model_features[, colnames(traj_model@model_features), drop = FALSE]
    pred_train <- logist(glmnet::predict.glmnet(traj_model@model, newx = feats, type = "link", s = lambda))[, 1]
    min_val <- min(pred_train)
    max_val <- max(pred_train)

    # Normalize predictions
    normalized_pred <- (pred_train - min_val) / (max_val - min_val)

    # Calculate min and max of scaled predictions
    norm_factors <- min(traj_model@diff_score, na.rm = TRUE)
    norm_factors <- c(norm_factors, max(traj_model@diff_score - norm_factors, na.rm = TRUE))

    IQmodel(
        iq_features = iq_features,
        pbms = pbms,
        intercept = intercept,
        func = func,
        lambda = lambda,
        alpha = alpha,
        min_pred = min_val,
        max_pred = max_val,
        norm_factors = norm_factors
    )
}


#' Validate input for IQmodel prediction
#'
#' @param intervals Optional genomic intervals
#' @param sequences Optional DNA sequences
#' @param pbm_responses Optional pre-computed PBM responses
#'
#' @return A logical indicating whether to use intervals
#' @noRd
validate_predict_input <- function(intervals, sequences, pbm_responses) {
    if (is.null(intervals) && is.null(sequences) && is.null(pbm_responses)) {
        cli::cli_abort("Either {.field intervals}, {.field sequences}, or {.field pbm_responses} must be provided.")
    }

    if (sum(!is.null(intervals), !is.null(sequences), !is.null(pbm_responses)) > 1) {
        cli::cli_abort("Only one of {.field intervals}, {.field sequences}, or {.field pbm_responses} should be provided.")
    }

    return(is.null(sequences) && !is.null(intervals))
}

#' Compute PBM responses
#'
#' @param pbms List of PBM objects
#' @param sequences DNA sequences
#' @param func Energy computation function
#' @param pbm_responses Pre-computed PBM responses
#'
#' @return Computed or provided PBM responses
#' @noRd
compute_pbm_responses <- function(pbms, sequences, func, pbm_responses) {
    if (is.null(pbm_responses)) {
        pbm_responses <- pbm_list.compute(pbms, sequences, response = TRUE, func = func)
    }
    return(pbm_responses)
}

#' Get all feature names from IQ features
#'
#' @param iq_features List of IQ features
#'
#' @return Vector of all feature names
#' @noRd
get_all_feature_names <- function(iq_features) {
    unlist(lapply(iq_features, function(feat) {
        if (inherits(feat, "IQFeatureGroup")) {
            return(names(feat@features))
        } else {
            return(feat@name)
        }
    }))
}

#' Identify missing features
#'
#' @param all_feature_names Vector of all feature names
#' @param new_data Data frame of provided features
#'
#' @return Vector of missing feature names
#' @noRd
identify_missing_features <- function(all_feature_names, new_data) {
    setdiff(all_feature_names, colnames(new_data))
}

#' Create list of features to compute
#'
#' @param iq_features List of IQ features
#' @param missing_features Vector of missing feature names
#'
#' @return List of features to compute
#' @noRd
create_features_to_compute <- function(iq_features, missing_features) {
    features_to_compute <- list()
    for (i in 1:length(iq_features)) {
        feat <- iq_features[[i]]
        if (inherits(feat, "IQFeatureGroup")) {
            missing_group_features <- intersect(names(feat@features), missing_features)
            if (length(missing_group_features) > 0) {
                features_to_compute[[names(iq_features)[i]]] <- feat
            }
        } else if (feat@name %in% missing_features) {
            features_to_compute[[feat@name]] <- feat
        }
    }
    return(features_to_compute)
}

#' Compute missing features
#'
#' @param features_to_compute List of features to compute
#' @param use_intervals Logical indicating whether to use intervals
#' @param intervals Genomic intervals
#' @param sequences DNA sequences
#'
#' @return Data frame of computed features
#' @noRd
compute_missing_features <- function(features_to_compute, use_intervals, intervals, sequences) {
    if (use_intervals) {
        iq_feature_list.compute(features_to_compute, intervals = intervals)
    } else {
        iq_feature_list.compute(features_to_compute, sequences = sequences)
    }
}

#' Compute IQ responses
#'
#' @param iq_features List of IQ features
#' @param new_data Data frame of feature values
#' @param use_intervals Logical indicating whether to use intervals
#' @param intervals Genomic intervals
#' @param sequences DNA sequences
#'
#' @return Data frame of IQ responses
#' @noRd
compute_iq_responses <- function(iq_features, new_data, use_intervals, intervals, sequences) {
    if (is.null(new_data)) {
        new_data <- data.frame()
    }

    all_feature_names <- get_all_feature_names(iq_features)
    missing_features <- identify_missing_features(all_feature_names, new_data)

    if (length(missing_features) > 0) {
        cli::cli_alert("Computing missing IQ features: {.val {missing_features}}")
        features_to_compute <- create_features_to_compute(iq_features, missing_features)

        computed_features <- compute_missing_features(features_to_compute, use_intervals, intervals, sequences)
        new_data <- cbind(new_data, computed_features)
    }

    iq_feature_list.compute_response(iq_features, new_data)
}

#' @noRd
handle_missing_iq_responses <- function(iq_responses, all_feature_names) {
    missing_features <- setdiff(all_feature_names, colnames(iq_responses))
    if (length(missing_features) > 0) {
        cli::cli_warn("Some features are still missing in the IQ responses. They will be set to 0: {.val {missing_features}}")
    }
    return(iq_responses)
}

#' Predict method for IQmodel
#'
#' This method makes predictions using an IQmodel object on new data.
#'
#' @param object An IQmodel object
#' @param new_data Optional data frame or matrix of new data for IQ features
#' @param intervals Optional genomic intervals
#' @param sequences Optional DNA sequences
#' @param pbm_responses Optional pre-computed PBM responses
#' @param iq_responses Optional pre-computed IQ responses
#'
#' @return A vector of normalized and rescaled predictions
#'
#' @exportMethod predict
setMethod("predict", signature = "IQmodel", function(object, new_data = NULL, intervals = NULL, sequences = NULL, pbm_responses = NULL, iq_responses = NULL) {
    use_intervals <- validate_predict_input(intervals, sequences, pbm_responses)

    if (is.null(sequences) && !is.null(intervals)) {
        sequences <- prego::intervals_to_seq(intervals)
    }

    pbm_responses <- compute_pbm_responses(object@pbms, sequences, object@func, pbm_responses)
    all_features <- pbm_responses

    if (length(object@iq_features) > 0) {
        if (is.null(iq_responses)) {
            iq_responses <- compute_iq_responses(object@iq_features, new_data, use_intervals, intervals, sequences)
        }

        all_feature_names <- get_all_feature_names(object@iq_features)
        iq_responses <- handle_missing_iq_responses(iq_responses, all_feature_names)

        if (nrow(iq_responses) != nrow(pbm_responses)) {
            cli::cli_abort("The number of rows in IQ responses ({.val {nrow(iq_responses)}}) should match the number of rows in PBM responses ({.val {nrow(pbm_responses)}}).")
        }

        all_features <- cbind(all_features, iq_responses)
    }

    pred <- rowSums(all_features, na.rm = TRUE) + object@intercept
    pred <- logist(pred)
    pred <- (pred - object@min_pred) / (object@max_pred - object@min_pred)
    pred <- pred * object@norm_factors[2] + object@norm_factors[1]

    return(pred)
})
