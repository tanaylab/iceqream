#' TrajectoryModelMulti Class
#'
#' This class represents a multi-trajectory model, which is a collection of individual trajectory models.
#'
#' @slot models A list of individual trajectory models.
#' @slot models A list of individual trajectory models with full features.
#' @slot motif_models A list of motif models corresponding to the individual trajectory models.
#' @slot cluster_map a data frame that maps original motif to its distilled cluster.
#' @slot stats A data frame with r2 statistics for each trajectory model after distillation.
#' @slot params additional parameters for the multi-trajectory distillation.
#'
#' @export
TrajectoryModelMulti <- setClass(
    "TrajectoryModelMulti",
    slots = list(
        models = "list",
        models_full = "list",
        motif_models = "list",
        cluster_map = "data.frame",
        stats = "data.frame",
        params = "list"
    )
)

#' @param object An instance of `TrajectoryModelMulti`.
#' @rdname TrajectoryModelMulti-class
setMethod(
    "show", "TrajectoryModelMulti",
    function(object) {
        cli::cli({
            cli::cli_text("{.cls TrajectoryModelMulti} with {.val {length(object@models)}} individual models\n")
            cli::cli_text("\n")
            cli::cli_text("Slots include:")
            cli_ul(c("{.field @models}: A list of individual trajectory models ({.val {length(object@models)}} models). models: {.val {names(object@models)}}."))
            cli_ul(c("{.field @models_full}: A list of individual trajectory models with full features ({.val {length(object@models_full)}} models: {.val {names(object@models_full)}})."))
            cli_ul(c("{.field @motif_models}: A named list of motif models. Each element contains PSSM and spatial model ({.val {length(object@motif_models)}} models: {.val {names(object@motif_models)}})."))
            cli_ul(c("{.field @cluster_map}: A data frame that maps original motif to its distilled cluster."))
            cli_ul(c("{.field @stats}: A data frame with r2 statistics for each trajectory model after distillation."))
            cli_ul(c("{.field @params}: A list of parameters used for multi-trajectory distillation (including: {.val {names(object@params)}})."))
        })
    }
)

validate_traj_model_multi <- function(object) {
    if (!methods::is(object, "TrajectoryModelMulti")) {
        cli_abort("Please provide an instance of {.cls TrajectoryModelMulti}", call = parent.frame(1))
    }
}
