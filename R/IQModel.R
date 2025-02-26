#' IQmodel class
#'
#' This class represents an IQ model including IQ features and PBMs,
#' along with necessary model parameters.
#'
#' @slot features A list of IQFeature objects.
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
        features = "list",
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
        n_groups <- sum(purrr::map_lgl(object@features, ~ inherits(.x, "IQFeatureGroup")))
        if (n_groups > 0) {
            cli::cli_text("Contains {.val {length(object@pbms)}} PBMs ({.code @pbms}), and {.val {length(object@features)}} IQ features ({.code @features}), including {.val {n_groups}} feature groups")
        } else {
            cli::cli_text("Contains {.val {length(object@pbms)}} PBMs ({.code @pbms}), and {.val {length(object@features)}} IQ features ({.code @features})")
        }

        cli::cli_text("Intercept: {.val {object@intercept}} ({.code @intercept})")
        cli::cli_text("Energy computation function: {.val {object@func}} ({.code @func})")
        cli::cli_text("Regularization parameter: {.val {object@lambda}} ({.code @lambda})")
        cli::cli_text("Elastic net mixing parameter: {.val {object@alpha}} ({.code @alpha})")
        cli::cli_text("Min prediction value: {.val {object@min_pred}} ({.code @min_pred})")
        cli::cli_text("Max prediction value: {.val {object@max_pred}} ({.code @max_pred})")
        cli::cli_text("Normalization factors: {.val {object@norm_factors}} ({.code @norm_factors})")

        cli::cli_text("\n")
        if (length(object@features) > 0) {
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
    if (has_interactions(traj_model)) {
        cli::cli_abort("Creating an IQ model from a trajectory model with interactions is not implemented yet.")
    }

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
        features = iq_features,
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
    if (is.null(new_data)) {
        return(all_feature_names)
    }
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
    all_feature_names <- get_all_feature_names(iq_features)
    missing_features <- identify_missing_features(all_feature_names, new_data)

    if (length(missing_features) > 0) {
        cli::cli_alert("Computing missing IQ features: {.val {missing_features}}")
        features_to_compute <- create_features_to_compute(iq_features, missing_features)

        computed_features <- compute_missing_features(features_to_compute, use_intervals, intervals, sequences)
        if (!is.null(computed_features)) {
            if (is.null(new_data)) {
                new_data <- computed_features
            } else {
                new_data <- cbind(new_data, computed_features)
            }
        }
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
#' @param rescale Logical indicating whether to rescale predictions
#'
#' @return A vector of normalized and rescaled predictions
#'
#' @exportMethod predict
setMethod("predict", signature = "IQmodel", function(object, new_data = NULL, intervals = NULL, sequences = NULL, pbm_responses = NULL, iq_responses = NULL, rescale = TRUE) {
    use_intervals <- validate_predict_input(intervals, sequences, pbm_responses)

    if (is.null(sequences) && !is.null(intervals)) {
        sequences <- prego::intervals_to_seq(intervals)
    }

    pbm_responses <- compute_pbm_responses(object@pbms, sequences, object@func, pbm_responses)
    all_features <- pbm_responses

    if (length(object@features) > 0) {
        if (is.null(iq_responses)) {
            iq_responses <- compute_iq_responses(object@features, new_data, use_intervals, intervals, sequences)
        }

        all_feature_names <- get_all_feature_names(object@features)
        iq_responses <- handle_missing_iq_responses(iq_responses, all_feature_names)

        if (nrow(iq_responses) != nrow(pbm_responses)) {
            cli::cli_abort("The number of rows in IQ responses ({.val {nrow(iq_responses)}}) should match the number of rows in PBM responses ({.val {nrow(pbm_responses)}}).")
        }

        all_features <- cbind(all_features, iq_responses)
    }

    pred <- rowSums(all_features, na.rm = TRUE) + object@intercept
    pred <- logist(pred)
    if (rescale) {
        pred <- (pred - object@min_pred) / (object@max_pred - object@min_pred)
        pred <- pred * object@norm_factors[2] + object@norm_factors[1]
    }


    return(pred)
})

#' Compute all features required by an IQModel for sequences or intervals
#'
#' This function computes all features required by a given IQModel for prediction,
#' based on the model's configuration. It can return either raw features or model responses.
#'
#' @param model An IQModel object that defines which features to compute.
#' @param sequences Optional vector of DNA sequences. Either sequences or intervals must be provided.
#' @param intervals Optional data frame of genomic intervals with columns 'chrom', 'start', and 'end'.
#' @param size Optional size for intervals. If not provided, will be determined from the model.
#' @param return_responses Logical indicating whether to return responses (TRUE) or raw features (FALSE).
#'        Default is FALSE to return responses.
#' @param return_separate Logical indicating whether to return separate data frames for PBM and IQ features.
#'        Default is FALSE.
#'
#' @return If return_separate is FALSE, a single data frame containing all computed features/responses.
#'         If return_separate is TRUE, a list with components for PBMs and IQ features.
#'
#' @export
iq_model.compute_features <- function(model, sequences = NULL, intervals = NULL, size = NULL,
                                      return_responses = FALSE, return_separate = FALSE) {
    # Validate input parameters
    if (is.null(sequences) && is.null(intervals)) {
        cli::cli_abort("Either sequences or intervals must be provided.")
    }

    # Prepare sequences from intervals if needed
    if (is.null(sequences) && !is.null(intervals)) {
        # Determine size if not provided
        if (is.null(size)) {
            # Try to get size from IQ features
            if (!is.null(model@features) && length(model@features) > 0) {
                for (feat in model@features) {
                    if ("size" %in% slotNames(feat)) {
                        size <- feat@size
                        cli::cli_alert_info("Using size {.val {size}} from IQ features")
                        break
                    }
                }
            }

            # If not found in IQ features, try PBMs
            if (is.null(size) && !is.null(model@pbms) && length(model@pbms) > 0) {
                size <- model@pbms[[1]]@size
                cli::cli_alert_info("Using size {.val {size}} from PBMs")
            }

            # If still not found, throw an error
            if (is.null(size)) {
                cli::cli_abort("Size could not be determined from model. Please provide the size parameter.")
            }
        }

        # Convert intervals to sequences with the specified size
        cli::cli_alert_info("Converting intervals to sequences with size {.val {size}}")
        sequences <- prego::intervals_to_seq(intervals, size = size)
    }

    # Process PBMs
    if (!is.null(model@pbms) && length(model@pbms) > 0) {
        cli::cli_alert_info("Computing PBM features")

        if (return_responses) {
            # Compute PBM responses directly
            pbm_results <- pbm_list.compute(model@pbms, sequences, response = TRUE, func = model@func)
        } else {
            # Compute raw PBM features (binding scores) without aggregating with the response function
            pbm_results <- pbm_list.compute(model@pbms, sequences, response = FALSE)
        }
    } else {
        pbm_results <- NULL
    }

    # Process IQ features
    if (!is.null(model@features) && length(model@features) > 0) {
        # First compute the raw feature values
        if (!is.null(intervals)) {
            cli::cli_alert_info("Computing IQ features from intervals")
            iq_features <- iq_feature_list.compute(model@features, intervals = intervals)
        } else {
            cli::cli_alert_info("Computing IQ features from sequences")
            iq_features <- iq_feature_list.compute(model@features, sequences = sequences)
        }

        if (return_responses && !is.null(iq_features)) {
            # Compute responses based on the model coefficients
            cli::cli_alert_info("Computing IQ responses")
            iq_results <- iq_feature_list.compute_response(model@features, iq_features)
        } else {
            # Use the raw features
            iq_results <- iq_features
        }
    } else {
        iq_results <- NULL
    }

    # Return results based on format preference
    if (return_separate) {
        if (return_responses) {
            return(list(
                pbm_responses = pbm_results,
                iq_responses = iq_results
            ))
        } else {
            return(list(
                pbm_features = pbm_results,
                iq_features = iq_results
            ))
        }
    } else {
        # Combine all results into a single data frame
        combined_results <- NULL

        if (!is.null(pbm_results)) {
            combined_results <- pbm_results
        }

        if (!is.null(iq_results)) {
            if (is.null(combined_results)) {
                combined_results <- iq_results
            } else {
                # Check for dimension compatibility
                if (nrow(pbm_results) == nrow(iq_results)) {
                    combined_results <- cbind(combined_results, iq_results)
                } else {
                    cli::cli_warn("PBM data ({.val {nrow(pbm_results)}} rows) and IQ data ({.val {nrow(iq_results)}} rows) have different dimensions. Returning them separately.")
                    if (return_responses) {
                        return(list(
                            pbm_responses = pbm_results,
                            iq_responses = iq_results
                        ))
                    } else {
                        return(list(
                            pbm_features = pbm_results,
                            iq_features = iq_results
                        ))
                    }
                }
            }
        }

        # If intervals had rownames, add them back
        if (!is.null(intervals) && !is.null(rownames(intervals)) && !is.null(combined_results)) {
            if (nrow(combined_results) == nrow(intervals)) {
                rownames(combined_results) <- rownames(intervals)
            }
        }

        return(as.data.frame(combined_results))
    }
}
