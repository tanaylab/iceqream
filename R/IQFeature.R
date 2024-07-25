#' IQFeature class
#'
#' This class represents an IQ feature
#'
#' @slot name The name of the IQ feature (character).
#' @slot coefs The coefficients of the IQ feature (numeric).
#'
#' @export
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
#' @export
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
    iq_features <- purrr::map(colnames(traj_model@additional_features), function(name) {
        variables <- f2v %>%
            filter(variable == name) %>%
            pull(feature)
        coefs <- coef(traj_model@model, s = traj_model@params$lambda)[variables, , drop = TRUE]
        names(coefs) <- gsub(paste0("^", name, "_"), "", names(coefs))
        IQFeature(name = name, coefs = coefs)
    })
    names(iq_features) <- colnames(traj_model@additional_features)

    return(iq_features)
}

#' Compute the IQ feature
#'
#' This function computes the IQ feature response for a given set of values.
#'
#' @param iq An IQFeature object.
#' @param values A vector of values to compute the IQ feature.
#'
#' @return A vector of computed IQ feature responses.
#'
#' @export
iq_feature.compute <- function(iq, values) {
    logist_e <- create_logist_features(as.matrix(values))
    return((logist_e %*% iq@coefs)[, 1])
}


#' Compute IQ feature list
#'
#' This function computes the IQ response for a given list of IQ features and a matrix of values.
#'
#' @param iq_list A list of IQFeature objects
#' @param mat The input matrix
#'
#' @return A matrix containing the computed IQ feature responses.
#'
#' @export
iq_feature_list.compute <- function(iq_list, mat) {
    resp <- purrr::imap(iq_list, ~ {
        if (!(.y %in% colnames(mat))) {
            return(rep(NA, nrow(mat)))
        }
        iq_feature.compute(.x, mat[, .y])
    })
    return(do.call(cbind, resp))
}
