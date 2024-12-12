#' Homogenize a list of PSSM models
#'
#' This function adjusts each PSSM model in the input list so that the GC content is not higher
#' than the AT content. If a model's GC content is higher than its AT content, the function applies
#' a reverse complement to the model using the prego pssm_rc function.
#'
#' @param models A list of PSSM models. Each model should be a list with a pssm element, which
#'               is a data frame containing columns 'A', 'C', 'G', 'T', and 'pos'.
#' @return A list of homogenized PSSM models.
#'
#' @examples
#' # Create simulated data
#' pssm1 <- data.frame(
#'     pos = 1:4,
#'     A = c(0.1, 0.2, 0.3, 0.4),
#'     C = c(0.3, 0.3, 0.2, 0.1),
#'     G = c(0.3, 0.3, 0.3, 0.3),
#'     T = c(0.3, 0.2, 0.2, 0.2)
#' )
#' pssm2 <- data.frame(
#'     pos = 1:4,
#'     A = c(0.1, 0.2, 0.3, 0.4),
#'     C = c(0.1, 0.1, 0.1, 0.1),
#'     G = c(0.2, 0.2, 0.2, 0.2),
#'     T = c(0.6, 0.5, 0.4, 0.3)
#' )
#'
#' models <- list(list(pssm = pssm1), list(pssm = pssm2))
#'
#' # Homogenize the models
#' homogenized_models <- homogenize_pssm_models(models)
#'
#' @export
homogenize_pssm_models <- function(models) {
    models <- models %>%
        purrr::map(homogenize_model)

    return(models)
}

homogenize_model <- function(model) {
    # if A is higher than T reverse complement the model
    if (sum(model$pssm$T) > sum(model$pssm$A)) {
        model$pssm <- prego::pssm_rc(model$pssm)
    }
    return(model)
}



compute_directed_pwm <- function(
    sequences, pssm, spat = NULL, spat_min = 1, spat_max = NULL,
    prior = 0.01, func = "logSumExp") {
    sequences <- direct_sequences(sequences, pssm)
    return(prego::compute_pwm(sequences, pssm,
        spat = spat, spat_min = spat_min, spat_max = spat_max,
        bidirect = TRUE, prior = prior, func = func
    ))
}
