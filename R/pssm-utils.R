#' Homogenize a list of PSSM models
#'
#' This function adjusts each PSSM model in the input list so that the GC content is not higher
#' than the AT content. If a model's GC content is higher than its AT content, the function applies
#' a reverse complement to the model using the pssm_rc function.
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
    homogenize_model <- function(model) {
        # if GC is higher than AT, we need to reverse complement the model
        gc <- sum(model$pssm$C + model$pssm$G)
        at <- sum(model$pssm$A + model$pssm$T)
        if (gc > at) {
            model$pssm <- pssm_rc(model$pssm)
        }
        return(model)
    }

    models <- models %>%
        purrr::map(homogenize_model)

    return(models)
}

#' Reverse complement a PSSM
#'
#' @param pssm A PSSM. Data frame with columns 'A', 'C', 'G', 'T' and 'pos'.
#' @return A PSSM with the same format, but reverse complemented.
#'
#' @examples
#' # Create simulated PSSM data frame
#' pssm <- data.frame(
#'     pos = 1:4,
#'     A = c(0.1, 0.2, 0.3, 0.1),
#'     C = c(0.1, 0.3, 0.2, 0.1),
#'     G = c(0.1, 0.3, 0.3, 0.7),
#'     T = c(0.7, 0.2, 0.2, 0.1)
#' )
#'
#' # Reverse complement the PSSM
#' rc_pssm <- pssm_rc(pssm)
#'
#' @export
pssm_rc <- function(pssm) {
    pssm <- pssm %>%
        mutate(tmp_A = A, tmp_C = C, tmp_G = G, tmp_T = T) %>%
        mutate(A = tmp_T, T = tmp_A, C = tmp_G, G = tmp_C) %>%
        select(-starts_with("tmp_")) %>%
        arrange(desc(pos)) %>%
        mutate(pos = n() - pos + 1)
    return(pssm)
}
