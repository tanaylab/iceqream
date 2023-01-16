#' EFlow object
#'
#'
#' @description Eflow object is ...
#'
#'
#' @slot id an identifier for the object, e.g. "gastrulation".
#' @slot description description of the object, e.g. "Gastrulation flow"
#'
#' @exportClass Eflow
Eflow <- setClass(
    "Eflow",
    slots = c(
        id = "character",
        description = "character",
        mc_metadata = "data.frame",
        peaks_metadata = "data.frame",
        mc_flow_graph = "matrix",
        matched_mcs_nms = "data.frame",
        mc_rna = "matrix",
        mc_atac = "matrix",
        mc_atac_s = "matrix",
        mc_atac_norm = "matrix",
        pwm_metadata = "data.frame",
        pwms = "character",
        peaks_energy = "matrix",
        pwm_activity = "matrix"
    )
)

setMethod(
    "initialize",
    signature = "Eflow",
    definition = function(.Object, id = NULL, description = NULL, mc_metadata = NULL, peaks_metadata = NULL, mc_flow_graph = NULL, matched_mcs_nms=NULL, mc_rna = NULL, mc_atac = NULL, mc_atac_s = NULL, mc_atac_norm = NULL, pwm_metadata = NULL, pwms = NULL, peaks_energy = NULL, pwm_activity = NULL) {
        .Object <- make_eflow_object(.Object, id, description, mc_metadata, peaks_metadata, mc_flow_graph, matched_mcs_nms, mc_rna , mc_atac, mc_atac_s, mc_atac_norm, pwm_metadata, pwms, peaks_energy, pwm_activity = NULL)
        return(.Object)
    }
)

make_eflow_object <- function(.Object, id = NULL, description = NULL, mc_metadata = NULL, peaks_metadata = NULL, mc_flow_graph = NULL, matched_mcs_nms=NULL, mc_rna = NULL, mc_atac = NULL, mc_atac_s = NULL,  mc_atac_norm = NULL, pwm_metadata = NULL, pwms = NULL, peaks_energy = NULL, pwm_activity = NULL) {
    # TODO: this is just a skeleton - need to validate the input of these things!
    if (!is.null(id)) {
        .Object@id <- id
    }
    if (!is.null(description)) {
        .Object@description <- description
    }
    if (!is.null(mc_metadata)) {
        .Object@mc_metadata <- mc_metadata
    }
    if (!is.null(peaks_metadata)) {
        .Object@peaks_metadata <- peaks_metadata
    }
    if (!is.null(mc_flow_graph)) {
        .Object@mc_flow_graph <- mc_flow_graph
    }
    if (!is.null(matched_mcs_nms)) {
      .Object@matched_mcs_nms <- matched_mcs_nms
    }
    if (!is.null(mc_rna)) {
        .Object@mc_rna <- mc_rna
    }
    if (!is.null(mc_atac)) {
        .Object@mc_atac <- mc_atac
    }
    if (!is.null(mc_atac_s)) {
        .Object@mc_atac_s <- mc_atac_s
    }
    if (!is.null(mc_atac_norm)) {
      .Object@mc_atac_norm <- mc_atac_norm
    }
    if (!is.null(pwm_metadata)) {
        .Object@pwm_metadata <- pwm_metadata
    }
    if (!is.null(pwms)) {
        .Object@pwms <- pwms
    }
    if (!is.null(peaks_energy)) {
        .Object@peaks_energy <- peaks_energy
    }
    if (!is.null(pwm_activity)) {
        .Object@pwm_activity <- pwm_activity
    }
    return(.Object)
}

#' Initialize an EFlow object
#'
#' @param mc2_file Path to metacell2 (metacells python package) Anndata file
#' @param md A data frame with a column named "metacell" with the metacell IDs, and additional columns with metadata, or a file
#' which contains such data frame.
#' @param id an identifier for the object, e.g. "gastrulation".
#' @param description description of the object, e.g. "Gastrulation flow"
#'
#' @return an Eflow object
#'
#' @export
init <- function(mc2_file, md, id = NULL, description = NULL) {
    # TODO

    return(new("Eflow", id = id, description = description))
}
