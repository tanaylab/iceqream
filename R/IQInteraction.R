#' IQInteraction class
#'
#' This class represents an interaction feature for IQ models, inheriting from IQFeature.
#'
#' @slot name The name of the IQ interaction feature (character).
#' @slot coefs The coefficients of the IQ interaction feature (numeric).
#' @slot term1 The name of the first term in the interaction (character).
#' @slot term2 The name of the second term in the interaction (character).
#' @slot term1_type The type of the first term (motif, additional, etc.).
#' @slot term2_type The type of the second term (motif, additional, etc.).
#'
#' @exportClass IQInteraction
IQInteraction <- setClass(
    "IQInteraction",
    slots = list(
        term1 = "character",
        term2 = "character",
        term1_type = "character",
        term2_type = "character"
    ),
    contains = "IQFeature"
)

#' Show method for IQInteraction
#'
#' This method defines how an IQInteraction object should be displayed.
#'
#' @param object An IQInteraction object
#'
#' @rdname IQInteraction-class
#' @exportMethod show
setMethod("show", signature = "IQInteraction", definition = function(object) {
    cli::cli({
        cli::cli_text("An {.cls IQInteraction} object named {.val {object@name}}")
        cli::cli_text("Interaction between {.val {object@term1}} ({.val {object@term1_type}}) and {.val {object@term2}} ({.val {object@term2_type}})")
        if (length(object@coefs) > 0) {
            cli::cli_text("Contains {.val {length(object@coefs)}} coefficients: {.val {names(object@coefs)}}")
        } else {
            cli::cli_text("Contains no coefficients")
        }
    })
})

#' Check if an IQmodel has interaction features
#'
#' This function checks if an IQmodel contains any interaction features.
#'
#' @param model An IQmodel object
#' @return Logical indicating whether the model has interactions
#'
#' @export
iq_model.has_interactions <- function(model) {
    any(sapply(model@features, function(x) inherits(x, "IQInteraction")))
}

#' Get interaction features from an IQmodel
#'
#' This function extracts all interaction features from an IQmodel.
#'
#' @param model An IQmodel object
#' @return A list of IQInteraction objects
#'
#' @export
iq_model.get_interactions <- function(model) {
    model@features[sapply(model@features, function(x) inherits(x, "IQInteraction"))]
}

#' Get regular features from an IQmodel
#'
#' This function extracts all non-interaction features from an IQmodel.
#'
#' @param model An IQmodel object
#' @return A list of IQFeature objects (excluding interactions)
#'
#' @export
iq_model.get_regular_features <- function(model) {
    model@features[!sapply(model@features, function(x) inherits(x, "IQInteraction"))]
}

#' Get all feature names from IQ features with categorization by type
#'
#' @param iq_features List of IQ features
#' @param include_interaction_terms Logical indicating whether to include interaction terms
#'
#' @return A list with named vectors for different feature types
#' @export
iq_model.feature_names_by_type <- function(iq_features, include_interaction_terms = FALSE) {
    result <- list(
        regular = character(),
        interaction = character(),
        interaction_terms = character()
    )

    for (i in seq_along(iq_features)) {
        feat <- iq_features[[i]]
        feat_name <- names(iq_features)[i]

        if (inherits(feat, "IQFeatureGroup")) {
            result$regular <- c(result$regular, names(feat@features))
        } else if (inherits(feat, "IQInteraction")) {
            result$interaction <- c(result$interaction, feat@name)
            if (include_interaction_terms) {
                result$interaction_terms <- c(result$interaction_terms, c(feat@term1, feat@term2))
            }
        } else {
            result$regular <- c(result$regular, feat@name)
        }
    }

    # Remove duplicates
    result$regular <- unique(result$regular)
    result$interaction <- unique(result$interaction)
    result$interaction_terms <- unique(result$interaction_terms)

    return(result)
}

#' Compute all interaction features
#'
#' This function computes all interaction features for a given dataset.
#'
#' @param model An IQmodel object
#' @param feature_data A data frame containing feature values needed for interactions
#' @return A data frame containing computed interaction values
#'
#' @export
iq_model.compute_interactions <- function(model, feature_data) {
    # Get interaction features
    interaction_features <- iq_model.get_interactions(model)

    if (length(interaction_features) == 0) {
        return(NULL)
    }

    # Compute interaction values
    result <- matrix(0, nrow = nrow(feature_data), ncol = length(interaction_features))
    colnames(result) <- names(interaction_features)

    # Process each interaction
    for (i in seq_along(interaction_features)) {
        feat <- interaction_features[[i]]
        feat_name <- names(interaction_features)[i]

        # Check if the required terms are available
        if (feat@term1 %in% colnames(feature_data) && feat@term2 %in% colnames(feature_data)) {
            # Compute interaction
            result[, i] <- iq_interaction.compute(
                feat,
                feature_data[, feat@term1],
                feature_data[, feat@term2]
            )
        } else {
            # Issue a warning for missing terms
            missing_terms <- setdiff(c(feat@term1, feat@term2), colnames(feature_data))
            cli::cli_warn("Cannot compute interaction {.val {feat_name}}: missing terms {.val {missing_terms}}")

            # Set to NA
            result[, i] <- NA
        }
    }

    return(as.data.frame(result))
}

#' Compute an IQInteraction
#'
#' This function computes the value of an IQ interaction feature given the values of its terms.
#'
#' @param iq An IQInteraction object.
#' @param term1_values A vector of values for the first term.
#' @param term2_values A vector of values for the second term.
#'
#' @return A vector of computed and normalized interaction values.
#'
#' @export
iq_interaction.compute <- function(iq, term1_values, term2_values) {
    # Compute interaction by multiplying term values
    inter <- term1_values * term2_values

    # Normalize the interaction
    inter <- inter / max(inter, na.rm = TRUE)
    inter <- norm01(inter) * 10

    return(inter)
}

#' Convert trajectory model interactions to IQ feature list
#'
#' This function takes a trajectory model with interactions and converts them to a list of IQ interaction features.
#'
#' @param traj_model The trajectory model object to convert.
#' @return A list of IQ interaction features, where each feature contains the interaction name and its corresponding coefficients.
#'
#' @export
traj_model_interactions_to_iq_feature_list <- function(traj_model) {
    if (!has_interactions(traj_model)) {
        return(list())
    }

    # Extract feature-to-variable mapping with type information
    f2v <- feat_to_variable(traj_model, add_types = TRUE)

    # Get only interaction features
    f2v_inter <- f2v %>%
        filter(type == "interaction")

    # Create a list to store IQInteraction objects
    iq_interactions <- list()

    # Process each interaction variable
    for (i in 1:nrow(f2v_inter)) {
        var_row <- f2v_inter[i, ]
        var_name <- var_row$variable

        # Extract variables for this interaction
        var_features <- f2v %>%
            filter(variable == var_name) %>%
            pull(feature)

        # Extract coefficients for this interaction
        coefs <- glmnet::coef.glmnet(traj_model@model, s = traj_model@params$lambda)[var_features, , drop = TRUE]
        names(coefs) <- gsub(paste0("^", var_name, "_"), "", names(coefs))

        # Get term information
        term1 <- var_row$term1
        term2 <- var_row$term2
        term1_type <- var_row$term1_type
        term2_type <- var_row$term2_type

        # Create an IQInteraction object for this interaction
        iq_interactions[[var_name]] <- IQInteraction(
            name = var_name,
            coefs = coefs,
            term1 = term1,
            term2 = term2,
            term1_type = term1_type,
            term2_type = term2_type
        )
    }

    return(iq_interactions)
}

#' Compute interaction values
#'
#' @param interaction_features List of IQInteraction features
#' @param feature_data Data frame containing feature values
#'
#' @return Data frame of computed interaction values
#' @noRd
compute_interaction_values <- function(interaction_features, feature_data) {
    if (length(interaction_features) == 0) {
        return(NULL)
    }

    result <- matrix(0, nrow = nrow(feature_data), ncol = length(interaction_features))
    colnames(result) <- names(interaction_features)

    for (i in seq_along(interaction_features)) {
        feat <- interaction_features[[i]]
        feat_name <- names(interaction_features)[i]

        # Check if both terms are available in the data
        if (feat@term1 %in% colnames(feature_data) && feat@term2 %in% colnames(feature_data)) {
            # Compute the interaction
            result[, i] <- iq_interaction.compute(
                feat,
                feature_data[, feat@term1],
                feature_data[, feat@term2]
            )
        } else {
            cli::cli_warn("Cannot compute interaction {.val {feat_name}}: missing terms {.val {setdiff(c(feat@term1, feat@term2), colnames(feature_data))}}")
            result[, i] <- NA
        }
    }

    return(as.data.frame(result))
}

#' Compute interaction responses
#'
#' @param interaction_features List of IQInteraction features
#' @param interaction_values Data frame of interaction values
#'
#' @return Data frame of computed interaction responses
#' @noRd
compute_interaction_responses <- function(interaction_features, interaction_values) {
    if (length(interaction_features) == 0 || is.null(interaction_values)) {
        return(NULL)
    }

    result <- matrix(0, nrow = nrow(interaction_values), ncol = length(interaction_features))
    colnames(result) <- names(interaction_features)

    for (i in seq_along(interaction_features)) {
        feat <- interaction_features[[i]]
        feat_name <- names(interaction_features)[i]

        if (feat_name %in% colnames(interaction_values)) {
            # Apply the logistic transformation
            logist_e <- create_logist_features(as.matrix(interaction_values[, feat_name, drop = FALSE]))
            result[, i] <- (logist_e %*% feat@coefs)[, 1]
        } else {
            result[, i] <- NA
        }
    }

    return(as.data.frame(result))
}
