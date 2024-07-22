#' PBM Class
#'
#' This class represents a PBM (Physical Binding Model) object.
#'
#' @slot name A character string representing the name of the PBM.
#' @slot pssm A matrix representing the position-specific scoring matrix.
#' @slot max_energy A numeric value representing the maximum energy for value normalization
#' @slot min_energy A numeric value representing the minimal energy for value normalization
#' @slot energy_range A numeric value representing the energy range for value normalization
#' @slot spat A data frame representing the spatial distribution (optional)
#' @slot spat_min A numeric value representing the spatial minimum (optional).
#' @slot spat_max A numeric value representing the spatial maximum (optional)
#' @slot seq_length A numeric value representing the length of the sequence for spatial distribution (optional)
#' @slot coefs A set of coefficients learned by IceQream trajectory model (optional)
#'
#' @export
PBM <- setClass(
    "PBM",
    slots = list(
        name = "character",
        pssm = "matrix",
        max_energy = "numeric",
        min_energy = "numeric",
        energy_range = "numeric",
        spat = "data.frame",
        spat_min = "numeric",
        spat_max = "numeric",
        seq_length = "numeric",
        coefs = "numeric"
    )
)

setMethod("show", signature = "PBM", definition = function(object) {
    cli::cli({
        cli::cli_text("a {.cls PBM} object named {.val {object@name}} with {.val {nrow(object@pssm)}} positions ({.code @pssm})")
        cli::cli_text("Energy normalization max = {.val {object@max_energy}} ({.code @max_energy})")
        if (nrow(object@spat) > 0) {
            cli::cli_text("Spatial distribution with {.val {nrow(object@spat)}} spatial factors, from position {.val {object@spat_min}} to {.val {object@spat_max}} ({.val {object@seq_length}} bp) ({.code @spat})")
        }
        if (length(object@coefs) > 0) {
            cli::cli_text("Includes a model with {.val {length(object@coefs)}} coefficients ({.val {names(object@coefs)}}) ({.code @coefs})")
        }
    })
})

#' Validate a PBM object
#'
#' This function checks if the given object is a valid PBM object and contains all necessary components.
#'
#' @param pbm A PBM object to validate
#'
#' @return This function doesn't return a value. It raises an error if the PBM object is invalid.
#'
#' @noRd
validate_pbm <- function(pbm) {
    if (!methods::is(pbm, "PBM")) {
        cli::cli_abort("The object is not a {.cls PBM} object.")
    }

    if (length(pbm@name) == 0) {
        cli::cli_abort("The PBM object does not contain a name.")
    }

    if (nrow(pbm@pssm) == 0) {
        cli::cli_abort("The PBM object does not contain a PSSM.")
    }

    if (length(pbm@max_energy) == 0) {
        cli::cli_abort("The PBM object does not contain the maximum energy.")
    }

    if (length(pbm@min_energy) == 0) {
        cli::cli_abort("The PBM object does not contain the minimum energy.")
    }

    if (length(pbm@energy_range) != 2) {
        cli::cli_abort("The PBM object does not contain the energy range.")
    }

    if (nrow(pbm@spat) > 0) {
        spat <- pbm@spat
        if (!is.data.frame(spat)) {
            cli_abort("The {.field spat} argument should be a data frame")
        }
        if (!all(c("bin", "spat_factor") %in% colnames(spat))) {
            cli_abort("The {.field spat} data frame should have columns {.val bin} and {.val spat_factor}")
        }
        if (!is.numeric(spat$bin) || !is.numeric(spat$spat_factor)) {
            cli_abort("The {.field spat} data frame should have columns {.val bin} and {.val spat_factor} of type numeric")
        }
        binsize <- unique(diff(spat$bin))
        if (length(binsize) > 1) {
            cli_abort("The bins in {.field spat} should be of equal size")
        }

        if (pbm@spat_min >= pbm@spat_max) {
            cli::cli_abort("The spatial minimum must be less than the spatial maximum.")
        }
        if (pbm@seq_length <= 0) {
            cli::cli_abort("The sequence length must be greater than 0.")
        }
    }
}

#' Convert a trajectory model to a list of PBM objects
#'
#' This function takes a trajectory model and converts it into a list of PBM objects.
#'
#' @param traj_model A trajectory model object
#' @param func The function to use for computing the energies "logSumExp" or "max". Default is "logSumExp"
#' @param normalization_intervals The normalization intervals to use for computing the energies. Default is the normalization intervals of the trajectory model.
#' @param bits_threshold The threshold for trimming the PSSM (default is NULL, no trimming)
#'
#' @return A list of PBM objects
#'
#' @export
traj_model_to_pbm_list <- function(
    traj_model, func = "logSumExp",
    normalization_intervals = traj_model@normalization_intervals,
    bits_threshold = NULL) {
    f2v <- feat_to_variable(traj_model)
    norm_sequences <- prego::intervals_to_seq(misha.ext::gintervals.normalize(normalization_intervals, traj_model@params$peaks_size))
    cli_alert_info("Computing motif energies for {.val {length(traj_model@motif_models)}} motifs on {.val {nrow(normalization_intervals)}} normalization intervals")

    if (!is.null(bits_threshold)) {
        traj_model@motif_models <- purrr::map(traj_model@motif_models, ~ {
            .x$pssm <- prego::trim_pssm(.x$pssm, bits_threshold)
            .x
        })
    }

    norm_energies <- plyr::llply(traj_model@motif_models, function(x) {
        prego::compute_pwm(norm_sequences, x$pssm, spat = x$spat, spat_min = x$spat_min %||% 1, spat_max = x$spat_max, func = func)
    }, .parallel = TRUE)
    names(norm_energies) <- names(traj_model@motif_models)
    norm_energies <- do.call(cbind, norm_energies)

    norm_energies <- norm_energies / log(2)

    max_e <- log2(matrixStats::colQuantiles(2^norm_energies, probs = traj_model@params$energy_norm_quantile, na.rm = TRUE))
    norm_energies <- pmax(pmin(sweep(norm_energies, 2, max_e, `-`), 0), traj_model@params$min_energy)
    range_data <- matrixStats::colRanges(norm_energies, na.rm = TRUE)


    pbm_list <- purrr::imap(traj_model@motif_models, function(model, name) {
        pssm <- model$pssm %>%
            select(A, C, G, T) %>%
            as.matrix()
        spat <- model$spat
        spat_min <- model$spat_min
        spat_max <- model$spat_max
        seq_length <- model$seq_length
        max_energy <- max_e[name]
        min_energy <- traj_model@params$min_energy
        energy_range <- range_data[name, ]
        variables <- f2v %>%
            filter(variable == name) %>%
            pull(feature)
        coefs <- coef(traj_model@model, s = traj_model@params$lambda)[variables, , drop = TRUE]
        names(coefs) <- gsub(paste0("^", name, "_"), "", names(coefs))
        pbm <- PBM(name = name, pssm = pssm, max_energy = max_energy, min_energy = min_energy, energy_range = energy_range, spat = spat, spat_min = spat_min, spat_max = spat_max, seq_length = seq_length, coefs = coefs)
        validate_pbm(pbm)
        pbm
    })
    return(pbm_list)
}

#' Normalize energies for a PBM
#'
#' This function normalizes the energies for a given PBM.
#'
#' @param pbm A PBM object
#' @param energies A vector of energy values to normalize
#' @param max_energy The maximum energy value (default is pbm@max_energy)
#' @param energy_range The range of energy values (default is pbm@energy_range)
#' @param norm_energy_max The maximum value for the normalized energy (default is 10)
#'
#' @return A vector of normalized energy values
#'
#' @export
pbm.normalize_energies <- function(pbm, energies, max_energy = pbm@max_energy, min_energy = pbm@min_energy, energy_range = pbm@energy_range, norm_energy_max = 10) {
    x <- energies / log(2)
    y <- pmax(pmin(x - max_energy, 0), min_energy)
    range <- energy_range[2] - energy_range[1]
    y <- (y - energy_range[1]) / range * norm_energy_max
    y <- pmax(y, 0)
    y <- pmin(y, norm_energy_max)

    return(y)
}

#' Compute energies / response for a PBM on given sequences
#'
#'
#' @param pbm  A PBM object
#' @param sequences A set of sequences on which to compute the energies.
#' @param response A logical flag indicating whether to compute the response. Default is FALSE.
#' @return A data frame containing the computed energies.
#'
#' @inheritParams prego::compute_pwm
#' @export
pbm.compute <- function(pbm, sequences, response = FALSE, func = "logSumExp") {
    pssm <- pbm@pssm %>%
        as.data.frame() %>%
        mutate(pos = row_number() - 1) %>%
        select(pos, everything())
    energies <- prego::compute_pwm(sequences, pssm, spat = pbm@spat, spat_min = pbm@spat_min %||% 1, spat_max = pbm@spat_max, func = func)
    energies <- pbm.normalize_energies(pbm, energies)

    if (response) {
        logist_e <- create_logist_features(as.matrix(energies))
        energies <- (logist_e %*% pbm@coefs)[, 1]
    }

    return(energies)
}

#' Extract PBM scores for genomic intervals
#'
#' This function extracts PBM scores of a given PBM for genomic intervals and optionally computes the response.
#'
#' @param pbm A PBM object
#' @param intervals Genomic intervals to extract
#' @param response Logical, whether to compute the response (default is FALSE)
#'
#' @return The intervals, with an additional column containing the PBM scores / response
#'
#' @inheritParams pbm.compute
#' @export
pbm.gextract <- function(pbm, intervals, response = FALSE, func = "logSumExp") {
    sequences <- prego::intervals_to_seq(intervals)
    energies <- intervals
    energies[, pbm@name] <- pbm.compute(pbm, sequences, response, func = func)
    return(energies)
}

#' Compute energies / response for a list of PBMs on given sequences
#'
#' This function computes the energies for a list of PBMs on given sequences.
#'
#' @param pbm_list A list of PBM objects
#'
#' @inheritParams pbm.compute
#'
#' @return A matrix containing the computed energies for each PBM
#'
#' @export
pbm_list.compute <- function(pbm_list, sequences, response = FALSE, func = "logSumExp") {
    cli::cli_alert_info("Computing energies for {.val {length(pbm_list)}} PBMs on {.val {length(sequences)}} sequences")
    energies <- plyr::llply(pbm_list, function(pbm) {
        pbm.compute(pbm, sequences, response, func = func)
    }, .parallel = TRUE)
    energies <- do.call(cbind, energies)

    return(energies)
}

#' Compute energy for a list of pbm lists (multiple trajectories)
#'
#'
#' @param multi_traj A list of PBM lists
#'
#' @return A matrix containing the computed energies for each PBM
#'
#' @inheritParams pbm_list.compute
#' @export
pbm_list.multi_traj.compute_energy <- function(multi_traj, sequences, func = "logSumExp") {
    pbm_list <- purrr::flatten(multi_traj)
    pbm_list <- pbm_list[unique(names(pbm_list))]
    pbm_list.compute(pbm_list, sequences, response = FALSE, func = func)
}

#' Extract energy for a list of pbm lists (multiple trajectories)
#'
#' @return A matrix containing the computed energies for each PBM
#'
#' @inheritParams pbm.gextract
#' @inheritParams pbm_list.multi_traj.compute_energy
#'
#' @export
pbm_list.multi_traj.gextract_energy <- function(multi_traj, intervals, func = "logSumExp") {
    sequences <- prego::intervals_to_seq(intervals)
    pbm_list.multi_traj.compute_energy(multi_traj, sequences, func)
}


#' Extract PBM scores / response for a list of PBMs
#'
#' This function extracts  PBM scores for genomic intervals and optionally computes the response for a list of PBMs.
#'
#' @param pbm_list A list of PBM objects
#'
#' @inheritParams pbm.gextract
#'
#' @return The intervals, with an additional column containing the PBM scores / response for each PBM
#'
#' @export
pbm_list.gextract <- function(pbm_list, intervals, response = FALSE, func = "logSumExp") {
    sequences <- prego::intervals_to_seq(intervals)
    energies <- pbm_list.compute(pbm_list, sequences, response, func)
    energies <- cbind(intervals, energies)

    return(energies)
}

#' Normalize energies locally for a PBM
#'
#' This function performs local normalization of energies for a given PBM.
#'
#' @param pbm A PBM object
#' @param energies A matrix of local energy values to normalize
#'
#' @return A matrix of locally normalized energy values
#'
#' @export
pbm.normalize_energies_local <- function(pbm, energies) {
    pbm.normalize_energies(pbm, energies, pbm@max_energy)
}

#' Compute local PBM scores / response
#'
#' This function computes local energies for a given PBM and set of sequences.
#'
#' @param pbm A PBM object
#' @param sequences A vector of DNA sequences
#' @param response Logical, whether to compute the response (default is FALSE)
#' @param scaling_q Scaling quantile for normalization (default is NULL)
#'
#' @return A matrix of computed local energies
#'
#' @export
pbm.compute_local <- function(pbm, sequences, response = FALSE, scaling_q = NULL) {
    pssm <- pbm@pssm %>%
        as.data.frame() %>%
        mutate(pos = row_number() - 1) %>%
        select(pos, everything())
    energies <- prego::compute_local_pwm(sequences, pssm)

    # energies <- pbm.normalize_energies(pbm, energies)
    energies <- pbm.normalize_energies_local(pbm, energies)

    if (!is.null(scaling_q)) {
        energies <- norm0q(energies, scaling_q) * 10
    }

    if (response) {
        empty_cols <- colSums(!is.na(energies)) == 0
        resp <- apply(energies[, !empty_cols, drop = FALSE], 2, function(x) (create_logist_features(as.matrix(x)) %*% pbm@coefs)[, 1])
        energies[, !empty_cols] <- resp
    }

    return(energies)
}

#' Extract local PBM scores / response
#'
#' This function extracts local PBM scores for genomic intervals and optionally computes the response.
#'
#' @param intervals Genomic intervals to extract
#'
#' @return A matrix of extracted local energies / response
#'
#' @inheritParams pbm.compute_local
#'
#' @export
pbm.gextract_local <- function(pbm, intervals, response = FALSE, scaling_q = 0.995) {
    sequences <- prego::intervals_to_seq(intervals)

    pbm.compute_local(pbm, sequences, response, scaling_q)
}

#' Extract local PBM scores / response for a list of PBMs given sequences
#'
#' This function extracts local PBM scores for sequences and optionally computes the response for a list of PBMs.
#'
#' @param pbm_list A list of PBM objects
#'
#' @inheritParams pbm.compute_local
#'
#' @return A list of matrices containing extracted local energies for each PBM
#'
#' @export
pbm_list.compute_local <- function(pbm_list, sequences, response = FALSE, scaling_q = NULL) {
    energies <- plyr::llply(pbm_list, function(pbm) {
        pbm.compute_local(pbm, sequences, response, scaling_q)
    }, .parallel = TRUE)

    return(energies)
}

#' Extract local PBM scores / response for a list of PBMs
#'
#' This function extracts local PBM scores for genomic intervals and optionally computes the response for a list of PBMs.
#'
#' @param pbm_list A list of PBM objects
#'
#' @inheritParams pbm.gextract_local
#'
#' @return A list of matrices containing extracted local energies for each PBM
#'
#' @export
pbm_list.gextract_local <- function(pbm_list, intervals, response = FALSE, scaling_q = NULL) {
    sequences <- prego::intervals_to_seq(intervals)
    pbm_list.compute_local(pbm_list, sequences, response, scaling_q)
}


#' Trim a PBM
#'
#' This function trims a PBM object by removing positions with low information content.
#'
#' @param pbm A PBM object
#' @param bits_thresh The threshold for trimming positions (default is 0.1)
#'
#' @return A trimmed PBM object
#'
#' @export
pbm.trim_pssm <- function(pbm, bits_thresh = 0.1) {
    pbm@pssm <- as.matrix(prego::trim_pssm(pbm@pssm, bits_thresh) %>% select(-pos))
    return(pbm)
}

#' Trim a list of PBMs
#'
#' @param pbm_list A list of PBM objects
#'
#' @inheritParams pbm.trim_pssm
#'
#' @return A list of trimmed PBM objects
#'
#' @export
pbm_list.trim_pssm <- function(pbm_list, bits_thresh = 0.1) {
    purrr::map(pbm_list, ~ pbm.trim_pssm(.x, bits_thresh))
}

tn5bias <- function(sequences = NULL, intervals = NULL) {
    if (is.null(sequences) && is.null(intervals)) {
        cli::cli_abort("Either sequences or intervals must be provided")
    }

    if (is.null(sequences) && !is.null(intervals)) {
        sequences <- prego::intervals_to_seq(intervals)
    }

    dn <- prego::calc_sequences_dinuc_dist(sequences)
    bias <- log2(dn$GC / dn$AT)

    return(bias)
}
