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
#' @slot model_features
#'  A matrix containing the features used for training the model.
#'
#' @slot normalized_energies matrix
#'   A matrix containing normalized energies. If additional variables were used, they are also included.
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
#' @slot features_r2 numeric
#'  A numeric vector of R^2 values for each feature.
#'
#' @slot normalization_intervals data.frame
#'  A data frame containing the intervals used for energy normalization.
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
        model_features = "matrix",
        normalized_energies = "matrix",
        type = "character",
        diff_score = "numeric",
        predicted_diff_score = "numeric",
        initial_prego_models = "list",
        peak_intervals = "data.frame",
        normalization_intervals = "data.frame",
        additional_features = "data.frame",
        features_r2 = "numeric",
        params = "list"
    )
)

#' @param object An instance of `TrajectoryModel`.
#' @rdname TrajectoryModel-class
setMethod("show", signature = "TrajectoryModel", definition = function(object) {
    cli::cli({
        cli::cli_text("{.cls TrajectoryModel} with {.val {length(object@motif_models)}} motifs and {.val {length(object@additional_features)}} additional features\n")
        cli::cli_text("\n")
        cli::cli_text("Slots include:")
        cli_ul(c("{.field @model}: A GLM model object. Number of non-zero coefficients: {.val {sum(object@model$beta[, 1] != 0)}}"))
        cli_ul(c("{.field @motif_models}: A named list of motif models. Each element contains PSSM and spatial model ({.val {length(object@motif_models)}} models: {.val {names(object@motif_models)}})"))
        cli_ul(c("{.field @additional_features}: A data frame of additional features ({.val {ncol(object@additional_features)}} features)"))
        cli_ul(c("{.field @coefs}: A data frame of coefficients ({.val {nrow(object@coefs)}} elements)"))
        cli_ul(c("{.field @model_features}: A matrix of the model features (logistic functions of the motif models energies, dimensions: {.val {nrow(object@model_features)}}x{.val {ncol(object@model_features)}})"))
        cli_ul(c("{.field @normalized_energies}: A matrix of normalized energies of the model features (dimensions: {.val {nrow(object@normalized_energies)}}x{.val {ncol(object@normalized_energies)}})"))
        cli_ul(c("{.field @type}: A vector the length of the number of peaks, indicating whether each peak is a training ('train') or a prediction peak ('test')"))
        cli_ul(c("{.field @diff_score}: A numeric value representing the difference score the model was trained on ({.val {length(object@normalized_energies[,1])}} elements)"))
        cli_ul(c("{.field @predicted_diff_score}: A numeric value representing the predicted difference score"))
        cli_ul(c("{.field @initial_prego_models}: A list of prego models used in the initial phase of the algorithm ({.val {length(object@initial_prego_models)}} models)"))
        cli_ul(c("{.field @peak_intervals}: A data frame containing the peak intervals ({.val {nrow(object@peak_intervals)}} elements)"))
        if ("normalization_intervals" %in% slotNames(object)) { # here for backwards compatibility
            cli_ul(c("{.field @normalization_intervals}: A data frame containing the intervals used for energy normalization ({.val {nrow(object@normalization_intervals)}} elements)"))
        }
        if (length(object@features_r2) > 0) {
            cli_ul(c("{.field @features_r2}: A numeric vector of the added R^2 values for each feature ({.val {length(object@features_r2)}} elements)"))
        } else {
            cli_ul(c("{.field @features_r2}: Model was not filtered, no R^2 values available."))
        }
        cli_ul(c("{.field @params}: A list of parameters used for training (including: {.val {names(object@params)}})"))

        cli::cli_text("\n")
        if (any(object@type == "train")) {
            cli::cli_text("R^2 train: {.val {round(cor(object@diff_score[object@type == 'train'], object@predicted_diff_score[object@type == 'train'], use = 'pairwise.complete.obs')^2, digits = 3)}}")
        }

        if (!any(object@type == "test")) {
            cli::cli_text("\n")
            cli::cli_text("Run {.code predict(object, peak_intervals)} to predict the model on new data.")
            cli::cli_text("Run {.code infer_trajectory_motifs(object, peak_intervals)} to create an object that includes the predicted peaks.")
        } else {
            cli::cli_text("R^2 test: {.val {round(cor(object@diff_score[object@type == 'test'], object@predicted_diff_score[object@type == 'test'], use = 'pairwise.complete.obs')^2, digits = 3)}}")
        }
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
setMethod("predict", signature = "TrajectoryModel", definition = function(object, peak_intervals, atac_scores = NULL, bin_start = 1, bin_end = ncol(atac_scores), additional_features = NULL) {
    traj_model <- infer_trajectory_motifs(object, peak_intervals, additional_features = additional_features)
    return(traj_model@predicted_diff_score[traj_model@type == "test"])
})

validate_traj_model <- function(object) {
    if (!methods::is(object, "TrajectoryModel")) {
        cli_abort("Please provide an instance of {.cls TrajectoryModel}", call = parent.frame(1))
    }
}
