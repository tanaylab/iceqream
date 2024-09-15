compute_motif_energies <- function(peak_intervals, db, normalization_intervals = peak_intervals, prior = 0.01, normalize = TRUE, energy_norm_quantile = 1, norm_energy_max = 10, min_energy = -7) {
    cli_alert("Computing motif energies (this might take a while). You can set the number of cores to use with {.code prego::set_parallel()}")
    motif_energies <- prego::gextract_pwm(peak_intervals %>% select(chrom, start, end), dataset = db, prior = prior) %>%
        select(-chrom, -start, -end) %>%
        as.matrix()

    if (normalize) {
        missing_intervals <- normalization_intervals %>% anti_join(peak_intervals, by = c("chrom", "start", "end"))
        if (nrow(missing_intervals) > 0) {
            cli::cli_alert_info("Extracting {.val {nrow(missing_intervals)}} intervals to normalize the motif energies")
            missing_energies <- prego::gextract_pwm(missing_intervals %>% select(chrom, start, end), dataset = db, prior = prior) %>%
                select(-chrom, -start, -end) %>%
                as.matrix()
            norm_energies <- rbind(motif_energies, missing_energies)
        } else {
            norm_energies <- motif_energies
        }

        cli_alert_info("Normalizing {.val {nrow(motif_energies)}} motif energies using {.val {nrow(norm_energies)}} intervals")
        motif_energies <- norm_energy_matrix(motif_energies, norm_energies, min_energy = min_energy, q = energy_norm_quantile, norm_energy_max = norm_energy_max)
    }


    return(motif_energies)
}

validate_motif_energies <- function(motif_energies, peak_intervals, pssm_db = motif_db) {
    if (nrow(motif_energies) != nrow(peak_intervals)) {
        cli_abort("Motif energies must have the same number of rows as peak_intervals. peak_intervals number of rows {.val {nrow(peak_intervals)}} != motif_energies number of rows {.val {nrow(motif_energies)}}", call = parent.frame(1))
    }

    missing_motifs <- colnames(motif_energies)[!(colnames(motif_energies) %in% pssm_db$motif)]
    if (length(missing_motifs) > 0) {
        cli_abort("Motifs {.val {missing_motifs}} are missing from the motif database ({.field pssm_db} parameter).", call = parent.frame(1))
    }
}
