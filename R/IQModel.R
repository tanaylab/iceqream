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
    # Split model into train and test components
    traj_model <- split_traj_model_to_train_test(traj_model)$train

    # Extract features and PBMs
    iq_features <- traj_model_to_iq_feature_list(traj_model)
    pbms <- traj_model_to_pbm_list(traj_model, func = func)

    # Extract model parameters
    intercept <- glmnet::coef.glmnet(traj_model@model, s = traj_model@params$lambda)["(Intercept)", ]
    lambda <- traj_model@params$lambda
    alpha <- traj_model@params$alpha

    # Add interactions if present
    if (has_interactions(traj_model)) {
        cli::cli_alert_info("Adding interactions to IQ model")
        iq_interactions <- traj_model_interactions_to_iq_feature_list(traj_model)
        if (length(iq_interactions) > 0) {
            cli::cli_alert_success("Added {.val {length(iq_interactions)}} interaction features")
            iq_features <- c(iq_features, iq_interactions)
        } else {
            cli::cli_alert_warning("No valid interactions found to add")
        }
    }

    # Calculate min and max prediction values
    feats <- traj_model@model_features[, colnames(traj_model@model_features), drop = FALSE]
    pred_train <- logist(glmnet::predict.glmnet(traj_model@model, newx = feats, type = "link", s = lambda))[, 1]
    min_val <- min(pred_train)
    max_val <- max(pred_train)

    # Normalize predictions
    normalized_pred <- (pred_train - min_val) / (max_val - min_val)

    # Calculate normalization factors
    norm_factors <- min(traj_model@diff_score, na.rm = TRUE)
    norm_factors <- c(norm_factors, max(traj_model@diff_score - norm_factors, na.rm = TRUE))

    # Create and return the IQmodel object
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
        } else if (inherits(feat, "IQInteraction")) {
            # For interactions, we need both the interaction name and its terms
            return(c(feat@name, feat@term1, feat@term2))
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
    # Separate regular features from interaction features
    interaction_features <- list()
    regular_features <- list()

    for (i in 1:length(iq_features)) {
        feat <- iq_features[[i]]
        feat_name <- names(iq_features)[i]

        if (inherits(feat, "IQFeatureGroup")) {
            missing_group_features <- intersect(names(feat@features), missing_features)
            if (length(missing_group_features) > 0) {
                regular_features[[feat_name]] <- feat
            }
        } else if (inherits(feat, "IQInteraction")) {
            # For interactions, check if the interaction itself is missing
            # or if any of its terms are missing
            if (feat@name %in% missing_features ||
                feat@term1 %in% missing_features ||
                feat@term2 %in% missing_features) {
                interaction_features[[feat_name]] <- feat
            }
        } else if (feat@name %in% missing_features) {
            regular_features[[feat_name]] <- feat
        }
    }

    # Return a list with separated features for ordered computation
    return(list(
        regular = regular_features,
        interaction = interaction_features
    ))
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
compute_iq_responses <- function(iq_features, new_data, use_intervals, intervals, sequences, pbm_data = NULL) {
    # Separate interaction features from regular features
    regular_features <- iq_features[!sapply(iq_features, function(x) inherits(x, "IQInteraction"))]
    interaction_features <- iq_features[sapply(iq_features, function(x) inherits(x, "IQInteraction"))]

    # Calculate required features and identify missing ones
    all_feature_names <- get_all_feature_names(iq_features)
    missing_features <- identify_missing_features(all_feature_names, new_data)

    # Create a complete data frame that includes PBM data if available
    complete_data <- new_data
    if (!is.null(pbm_data)) {
        if (is.null(complete_data)) {
            complete_data <- pbm_data
        } else if (nrow(pbm_data) == nrow(complete_data)) {
            # Only add columns that don't already exist
            pbm_cols_to_add <- setdiff(colnames(pbm_data), colnames(complete_data))
            if (length(pbm_cols_to_add) > 0) {
                complete_data <- cbind(complete_data, pbm_data[, pbm_cols_to_add, drop = FALSE])
            }
        } else {
            cli::cli_warn("PBM data ({.val {nrow(pbm_data)}} rows) and feature data ({.val {nrow(complete_data)}} rows) have different dimensions. PBM data will not be included in interaction computations.")
        }
    }

    # Update missing features based on complete_data
    missing_features <- identify_missing_features(all_feature_names, complete_data)

    if (length(missing_features) > 0) {
        features_to_compute <- create_features_to_compute(iq_features, missing_features)

        # First compute regular features
        if (length(features_to_compute$regular) > 0) {
            cli::cli_alert_info("Computing regular features...")
            computed_features <- NULL
            if (use_intervals) {
                computed_features <- iq_feature_list.compute(features_to_compute$regular, intervals = intervals)
            } else {
                computed_features <- iq_feature_list.compute(features_to_compute$regular, sequences = sequences)
            }

            if (!is.null(computed_features)) {
                # Add computed features to both new_data and complete_data
                if (is.null(new_data)) {
                    new_data <- computed_features
                } else {
                    new_data <- cbind(new_data, computed_features[, setdiff(colnames(computed_features), colnames(new_data)), drop = FALSE])
                }

                if (is.null(complete_data)) {
                    complete_data <- computed_features
                } else {
                    complete_data <- cbind(complete_data, computed_features[, setdiff(colnames(computed_features), colnames(complete_data)), drop = FALSE])
                }
            }
        }

        # Then compute interaction features if needed
        if (length(features_to_compute$interaction) > 0) {
            cli::cli_alert_info("Computing interaction features...")
            # Use complete_data to ensure PBM features are included
            interaction_values <- compute_interaction_values(features_to_compute$interaction, complete_data)

            if (!is.null(interaction_values)) {
                # Add interaction values to both new_data and complete_data
                if (is.null(new_data)) {
                    new_data <- interaction_values
                } else {
                    new_data <- cbind(new_data, interaction_values[, setdiff(colnames(interaction_values), colnames(new_data)), drop = FALSE])
                }

                if (is.null(complete_data)) {
                    complete_data <- interaction_values
                } else {
                    complete_data <- cbind(complete_data, interaction_values[, setdiff(colnames(interaction_values), colnames(complete_data)), drop = FALSE])
                }
            }
        }
    }

    # Compute responses for regular features using complete_data
    if (length(regular_features) > 0) {
        reg_responses <- iq_feature_list.compute_response(regular_features, complete_data)
    } else {
        reg_responses <- NULL
    }

    # Compute responses for interaction features
    if (length(interaction_features) > 0) {
        # Use complete_data to extract interaction values
        int_values <- complete_data[, intersect(names(interaction_features), colnames(complete_data)), drop = FALSE]
        int_responses <- compute_interaction_responses(interaction_features, int_values)
    } else {
        int_responses <- NULL
    }

    # Combine responses
    if (!is.null(reg_responses) && !is.null(int_responses)) {
        return(cbind(reg_responses, int_responses))
    } else if (!is.null(reg_responses)) {
        return(reg_responses)
    } else if (!is.null(int_responses)) {
        return(int_responses)
    } else {
        return(NULL)
    }
}


#' Predict method for IQmodel
#'
#' This method makes predictions using an IQmodel object on new data.
#'
#' @param object An IQmodel object
#' @param new_data Optional data frame or matrix of pre-computed features
#' @param intervals Optional genomic intervals
#' @param sequences Optional DNA sequences
#' @param rescale Logical indicating whether to rescale predictions to the original scale
#'
#' @return A vector of normalized and rescaled predictions
#'
#' @exportMethod predict
setMethod("predict", signature = "IQmodel", function(object, new_data = NULL, intervals = NULL,
                                                     sequences = NULL, rescale = TRUE) {
    # Get feature responses directly
    responses <- iq_model.compute_features(
        object,
        new_data = new_data,
        sequences = sequences,
        intervals = intervals,
        return_responses = TRUE
    )

    # Compute final prediction
    pred <- rowSums(responses, na.rm = TRUE) + object@intercept
    pred <- logist(pred)

    # Rescale if requested
    if (rescale) {
        pred <- (pred - object@min_pred) / (object@max_pred - object@min_pred)
        pred <- pred * object@norm_factors[2] + object@norm_factors[1]
    }

    return(pred)
})

#' Compute all features or responses required by an IQModel
#'
#' This function computes features or feature responses for an IQModel based on
#' the provided inputs.
#'
#' @param model An IQModel object
#' @param sequences Optional vector of DNA sequences
#' @param intervals Optional data frame of genomic intervals
#' @param new_data Optional data frame of pre-computed features
#' @param return_responses Logical indicating whether to return responses (TRUE) or raw features (FALSE)
#' @param return_separate Logical indicating whether to return separate data frames
#' @param normalize_energies Logical indicating whether to normalize PBM energies
#'
#' @return A data frame of computed features/responses or list of separate components
#'
#' @export
iq_model.compute_features <- function(model, sequences = NULL, intervals = NULL,
                                      new_data = NULL, return_responses = FALSE,
                                      return_separate = FALSE, normalize_energies = TRUE) {
    # Validate inputs
    if (is.null(sequences) && is.null(intervals) && is.null(new_data)) {
        cli::cli_abort("Either sequences, intervals, or new_data must be provided.")
    }

    # Prepare sequences from intervals if needed
    if (is.null(sequences) && !is.null(intervals)) {
        size <- get_size_from_model(model)
        sequences <- prego::intervals_to_seq(intervals, size = size)
    }

    # Initialize containers
    features_list <- list()

    # STEP 1: Compute all raw features

    # 1a. Process PBMs if present (always compute raw features, not responses)
    if (length(model@pbms) > 0 && (!is.null(sequences) || !is.null(new_data))) {
        pbm_results <- NULL

        if (!is.null(new_data)) {
            # Check which PBMs are already in new_data
            existing_pbms <- intersect(names(model@pbms), colnames(new_data))
            missing_pbms <- setdiff(names(model@pbms), colnames(new_data))

            # Extract existing PBM results from new_data
            if (length(existing_pbms) > 0) {
                pbm_results <- new_data[, existing_pbms, drop = FALSE]
            }

            # Only compute missing PBMs if we have sequences
            if (length(missing_pbms) > 0 && !is.null(sequences)) {
                cli::cli_alert("Computing missing PBMs: {.val {missing_pbms}}")
                pbms_to_compute <- model@pbms[missing_pbms]

                # Compute missing PBMs
                computed_pbms <- pbm_list.compute(
                    pbms_to_compute,
                    sequences,
                    response = FALSE, # Always compute raw features first
                    func = model@func,
                    normalize_energies = normalize_energies
                )

                # Combine existing and newly computed PBM results
                if (!is.null(pbm_results)) {
                    pbm_results <- cbind(pbm_results, computed_pbms)
                } else {
                    pbm_results <- computed_pbms
                }
            }
        } else if (!is.null(sequences)) {
            # Compute all PBMs from scratch when no new_data provided
            pbm_results <- pbm_list.compute(
                model@pbms,
                sequences,
                response = FALSE, # Always compute raw features first
                func = model@func,
                normalize_energies = normalize_energies
            )
        }

        features_list$pbm <- pbm_results
    }

    # 1b. Process regular features
    reg_features <- model@features[!sapply(model@features, inherits, "IQInteraction")]
    if (length(reg_features) > 0) {
        # Get features from new_data or compute as needed
        if (!is.null(new_data)) {
            # Use existing data and compute only what's missing
            iq_results <- new_data
            missing_feats <- get_missing_features(reg_features, new_data)

            if (length(missing_feats) > 0 && (!is.null(sequences) || !is.null(intervals))) {
                feats_to_compute <- filter_features_by_names(reg_features, missing_feats)
                cli::cli_alert("Computing missing features: {.val {missing_feats}}")

                if (length(feats_to_compute) > 0) {
                    # Compute missing features
                    computed <- iq_feature_list.compute(feats_to_compute, sequences = sequences, intervals = intervals)

                    if (!is.null(computed)) {
                        new_cols <- setdiff(colnames(computed), colnames(iq_results))
                        if (length(new_cols) > 0) {
                            iq_results <- cbind(iq_results, computed[, new_cols, drop = FALSE])
                        }
                    }
                }
            }
        } else if (!is.null(sequences) || !is.null(intervals)) {
            # Compute all features from scratch
            iq_results <- iq_feature_list.compute(reg_features, sequences = sequences, intervals = intervals)
        }

        features_list$iq <- iq_results
    }

    # 1c. Process interaction features
    int_features <- sapply(model@features, inherits, "IQInteraction")
    if (length(int_features) > 0) {
        int_features <- model@features[int_features]
    }

    if (length(int_features) > 0) {
        # Combine all available feature data
        combined_features <- combine_feature_data(features_list, new_data)

        if (!is.null(combined_features)) {
            # Compute interactions
            int_values <- iq_model.compute_interactions(model, combined_features)
            features_list$interaction <- int_values
        }
    }

    # Combine all raw features
    combined_features <- combine_feature_data(features_list, NULL)

    # Add row names from intervals if applicable
    if (!is.null(intervals) && !is.null(rownames(intervals)) &&
        !is.null(combined_features) && nrow(combined_features) == nrow(intervals)) {
        rownames(combined_features) <- rownames(intervals)
    }

    # STEP 2: Apply response calculations if requested
    if (return_responses && !is.null(combined_features)) {
        responses <- iq_feature_list.compute_response(c(model@pbms, model@features), combined_features)
        rownames(responses) <- rownames(combined_features)
        if (return_separate) {
            return(list(
                features = as.data.frame(combined_features),
                responses = as.data.frame(responses)
            ))
        }
        return(responses)
    }
    return(as.data.frame(combined_features))
}

#' Get size parameter from model
#' @noRd
get_size_from_model <- function(model) {
    # Try features first
    for (feat in model@features) {
        if ("size" %in% slotNames(feat)) {
            return(feat@size)
        }
    }
    # Then try PBMs
    if (length(model@pbms) > 0) {
        return(model@pbms[[1]]@size)
    }
    cli::cli_abort("Could not determine sequence size from model.")
}

#' Get list of features missing from provided data
#' @noRd
get_missing_features <- function(features, data) {
    needed_features <- unlist(lapply(features, function(feat) {
        if (inherits(feat, "IQFeatureGroup")) {
            return(names(feat@features))
        } else {
            return(feat@name)
        }
    }))
    setdiff(needed_features, colnames(data))
}

#' Filter features list to only include those with specified names
#' @noRd
filter_features_by_names <- function(features, names_to_include) {
    result <- list()
    for (i in seq_along(features)) {
        feat <- features[[i]]
        if (inherits(feat, "IQFeatureGroup")) {
            missing_in_group <- intersect(names(feat@features), names_to_include)
            if (length(missing_in_group) > 0) {
                result[[length(result) + 1]] <- feat
            }
        } else if (feat@name %in% names_to_include) {
            result[[length(result) + 1]] <- feat
        }
    }
    return(result)
}

#' Combine feature data from multiple sources
#' @noRd
combine_feature_data <- function(result_list, new_data = NULL) {
    combined <- NULL

    # First add from result_list
    for (key in names(result_list)) {
        if (is.null(combined)) {
            combined <- result_list[[key]]
        } else if (!is.null(result_list[[key]])) {
            new_cols <- setdiff(colnames(result_list[[key]]), colnames(combined))
            if (length(new_cols) > 0) {
                combined <- cbind(combined, result_list[[key]][, new_cols, drop = FALSE])
            }
        }
    }

    # Then add from new_data if provided
    if (!is.null(new_data)) {
        if (is.null(combined)) {
            combined <- new_data
        } else {
            new_cols <- setdiff(colnames(new_data), colnames(combined))
            if (length(new_cols) > 0) {
                combined <- cbind(combined, new_data[, new_cols, drop = FALSE])
            }
        }
    }

    return(combined)
}
