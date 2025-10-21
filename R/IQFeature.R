#' IQFeature class
#'
#' This class represents an IQ feature
#'
#' @slot name The name of the IQ feature (character).
#' @slot coefs The coefficients of the IQ feature (numeric).
#'
#' @exportClass IQFeature
IQFeature <- setClass(
    "IQFeature",
    slots = list(
        name = "character",
        coefs = "numeric"
    )
)

#' Show method for IQFeature
#'
#' This method defines how an IQFeature object should be displayed.
#'
#' @param object An IQFeature object
#'
#' @rdname IQFeature-class
#' @exportMethod show
setMethod("show", signature = "IQFeature", definition = function(object) {
    cli::cli({
        cli::cli_text("An {.cls IQFeature} object named {.val {object@name}}")
        if (length(object@coefs) > 0) {
            cli::cli_text("Contains {.val {length(object@coefs)}} coefficients: {.val {names(object@coefs)}}")
        } else {
            cli::cli_text("Contains no coefficients")
        }
    })
})

#' Convert trajectory model to IQ feature list
#'
#' This function takes a trajectory model object and converts its additional features into a list of IQ features.
#'
#' @param traj_model The trajectory model object to convert.
#' @return A list of IQ features, where each feature contains the feature name and its corresponding coefficients.
#'
#' @export
traj_model_to_iq_feature_list <- function(traj_model) {
    f2v <- feat_to_variable(traj_model)
    all_feats <- colnames(traj_model@additional_features)

    dinucs <- c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")
    dinuc_feats <- all_feats[all_feats %in% dinucs]

    # Check if gc_content is among the additional features
    has_gc_content <- "gc_content" %in% all_feats

    # Check if cg_cont is among the additional features
    has_cg_cont <- "cg_cont" %in% all_feats

    # Initialize the IQ features list
    iq_features <- list()

    if (length(dinuc_feats) > 0) {
        # Create dinucleotide feature group, including GC content if available
        iq_features$dinucleotides <- create_dinuc_feature_group(traj_model,
            dinucleotides = dinuc_feats,
            include_gc = has_gc_content
        )

        # Remove dinucleotide features and GC content (if included) from all_feats
        feats_to_remove <- c(dinuc_feats)
        if (has_gc_content) feats_to_remove <- c(feats_to_remove, "gc_content")
        all_feats <- all_feats[!(all_feats %in% feats_to_remove)]
    }

    # Handle the cg_cont feature if present
    if (has_cg_cont) {
        # Create the CG content feature
        iq_features$cg_cont <- create_cg_content_feature(traj_model)

        # Remove cg_cont from all_feats to avoid processing it twice
        all_feats <- all_feats[all_feats != "cg_cont"]
    }

    # Process all remaining features
    other_features <- purrr::map(all_feats, function(name) {
        variables <- f2v %>%
            filter(variable == name) %>%
            pull(feature)
        coefs <- glmnet::coef.glmnet(traj_model@model, s = traj_model@params$lambda)[variables, , drop = TRUE]
        names(coefs) <- gsub(paste0("^", name, "_"), "", names(coefs))
        IQFeature(name = name, coefs = coefs)
    })
    names(other_features) <- all_feats

    iq_features <- c(iq_features, other_features)

    return(iq_features)
}

#' Compute the IQ feature response
#'
#' This function computes the IQ feature response for a given set of values.
#'
#' @param iq An IQFeature object.
#' @param values A vector of values to compute the IQ feature.
#'
#' @return A vector of computed IQ feature responses.
#'
#' @export
iq_feature.compute_response <- function(iq, values) {
    if (length(iq@coefs) == 1) {
        return(values * iq@coefs)
    } else {
        logist_e <- create_logist_features(as.matrix(values))
        return((logist_e %*% iq@coefs)[, 1])
    }
}


#' Compute IQ feature list response
#'
#' This function computes the IQ response for a given list of IQ features and a matrix of values.
#'
#' @param iq_list A list of IQFeature objects
#' @param mat The input matrix
#'
#' @return A matrix containing the computed IQ feature responses.
#'
#' @export
iq_feature_list.compute_response <- function(iq_list, mat) {
    compute <- function(.x, .y) {
        if (!(.y %in% colnames(mat))) {
            return(rep(NA, nrow(mat)))
        }
        iq_feature.compute_response(.x, mat[, .y])
    }

    resp <- purrr::imap(iq_list, ~ {
        if (inherits(.x, "IQFeatureGroup")) {
            return(purrr::imap(.x@features, compute))
        }
        compute(.x, .y)
    })

    result <- list()
    for (i in 1:length(resp)) {
        if (is.list(resp[[i]])) {
            result <- c(result, resp[[i]])
        } else {
            result <- c(result, list(resp[[i]]))
            names(result)[length(result)] <- names(resp)[i]
        }
    }
    result <- do.call(cbind, result)

    return(result)
}

#' Compute features for a list of IQFeature and IQFeatureGroup objects
#'
#' This function computes features for a given list of IQFeature and IQFeatureGroup objects,
#' using either provided sequences or intervals.
#'
#' @param feature_list A list of IQFeature and/or IQFeatureGroup objects.
#' @param sequences A vector of DNA sequences. Optional if intervals are provided.
#' @param intervals A data frame of genomic intervals with columns 'chrom', 'start', and 'end'. Optional if sequences are provided.
#'
#' @return A data frame containing the computed features for all objects in the feature_list.
#'
#' @export
iq_feature_list.compute <- function(feature_list, sequences = NULL, intervals = NULL) {
    if (is.null(sequences) && is.null(intervals)) {
        cli::cli_abort("Either sequences or intervals must be provided.")
    }

    # Initialize an empty list to store results
    results <- list()

    # Iterate through the feature_list
    for (feature_obj in feature_list) {
        if (inherits(feature_obj, "IQFeatureGroup")) {
            # Compute features for IQFeatureGroup
            group_results <- iq_feature_group.compute(feature_obj, sequences, intervals)
            results[[length(results) + 1]] <- group_results
        } else if (inherits(feature_obj, "IQSeqFeature")) {
            # Compute features for individual IQFeature or IQSeqFeature
            if (is.null(sequences)) {
                sequences <- prego::intervals_to_seq(gintervals.normalize(intervals, feature_obj@size))
            }
            feature_result <- iq_seq_feature.compute(feature_obj, sequences)
            results[[length(results) + 1]] <- data.frame(feature_result)
            colnames(results[[length(results)]]) <- feature_obj@name
        } else if (inherits(feature_obj, "IQFeature")) {
            # check if the feature has a compute function
            if ("compute_func" %in% slotNames(feature_obj)) {
                feature_result <- feature_obj@compute_func(intervals)
                results[[length(results) + 1]] <- data.frame(feature_result)
                colnames(results[[length(results)]]) <- feature_obj@name
            } else {
                cli::cli_warn("Skipping IQFeature object without a compute function: {.val {feature_obj@name}}. Note that the model accuracy may be affected.")
            }
        } else {
            cli::cli_warn("Skipping unknown object type in feature_list: {.type {class(feature_obj)}}")
        }
    }

    # Combine all results into a single data frame
    combined_results <- do.call(cbind, results)

    if (is.null(combined_results)) {
        return(NULL)
    }

    # make sure that the number of rows in the results match the number of intervals
    if (!is.null(intervals)) {
        if (nrow(combined_results) != nrow(intervals)) {
            cli::cli_abort("The number of rows in the results of the IQ features ({.val {nrow(combined_results)}}) did not match the number of intervals ({.val {nrow(intervals)}}).")
        }
    }

    # if intervals had rownames then add them back
    if (!is.null(intervals)) {
        rownames(combined_results) <- rownames(intervals)
    }

    return(as.data.frame(combined_results))
}
