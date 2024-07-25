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
#' @export
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
#' @export
setMethod("show", signature = "IQmodel", function(object) {
    cli::cli({
        cli::cli_text("An {.cls IQmodel} object")
        cli::cli_text("Contains {.val {length(object@pbms)}} PBMs ({.code @pbms}), and {.val {length(object@iq_features)}} IQ features ({.code @iq_features})")
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
#' @export
setMethod("predict", signature = "IQmodel", function(object, new_data = NULL, intervals = NULL, sequences = NULL, pbm_responses = NULL, iq_responses = NULL) {
    # Validate input
    if (is.null(intervals) && is.null(sequences) && is.null(pbm_responses)) {
        cli::cli_abort("Either {.field intervals}, {.field sequences}, or {.field pbm_responses} must be provided.")
    }

    if (sum(!is.null(intervals), !is.null(sequences), !is.null(pbm_responses)) > 1) {
        cli::cli_abort("Only one of {.field intervals}, {.field sequences}, or {.field pbm_responses} should be provided.")
    }

    # Compute or use provided PBM responses
    if (is.null(pbm_responses)) {
        if (is.null(sequences) && !is.null(intervals)) {
            sequences <- prego::intervals_to_seq(intervals)
        }
        pbm_responses <- pbm_list.compute(object@pbms, sequences, response = TRUE, func = object@func)
    }

    # Initialize all_features with PBM responses
    all_features <- pbm_responses

    # Compute or use provided IQ feature responses
    if (length(object@iq_features) > 0) {
        if (is.null(iq_responses)) {
            if (is.null(new_data)) {
                cli::cli_abort("{.field new_data} must be provided when the model includes IQ features and {.field iq_responses} is not provided.")
            }
            # make sure new data has all the needed features
            if (!all(names(new_data) %in% names(object@iq_features))) {
                cli::cli_warn("Some features are missing in the new data. They will be set to 0: {.val {setdiff(names(object@iq_features), names(new_data))}")
            }

            # make sure the number of rows in new data matches the number of rows in pbm_responses
            if (nrow(new_data) != nrow(pbm_responses)) {
                cli::cli_abort("The number of rows in {.field new_data} should match the number of rows in {.field pbm_responses}.")
            }

            iq_responses <- iq_feature_list.compute(object@iq_features, new_data)
        }

        # make sure iq_responses has all the needed features
        if (!all(names(iq_responses) %in% names(object@iq_features))) {
            cli::cli_warn("Some features are missing in the IQ responses. They will be set to 0: {.val {setdiff(names(object@iq_features), names(iq_responses))}")
        }

        # make sure the number of rows in iq_responses matches the number of rows in pbm_responses
        if (nrow(iq_responses) != nrow(pbm_responses)) {
            cli::cli_abort("The number of rows in {.field iq_responses} should match the number of rows in {.field pbm_responses}.")
        }

        all_features <- cbind(all_features, iq_responses)
    }

    # Make predictions
    if (is.matrix(all_features)) {
        pred <- rowSums(all_features, na.rm = TRUE) + object@intercept
    } else {
        pred <- all_features + object@intercept
    }

    # Apply logistic function
    pred <- logist(pred)

    # Normalize predictions
    pred <- (pred - object@min_pred) / (object@max_pred - object@min_pred)

    # Rescale predictions to match the original scaling
    pred <- pred * object@norm_factors[2] + object@norm_factors[1]

    return(pred)
})
