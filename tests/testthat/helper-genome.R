# Shared, memoized fixture for the genome-backed tests (test-genome-*.R).
#
# Gated exactly like test-e2e.R: callers are skipped on CRAN and when the misha
# genome or the cached vignette data is unavailable. Override the genome with
# ICEQREAM_GENOME_ROOT (defaults to mm10). The two base trajectory models are
# built once per test-run and shared across every genome test file.

genome_root_path <- function() {
    Sys.getenv("ICEQREAM_GENOME_ROOT", unset = "/home/aviezerl/mm10")
}

# Returns a list with:
#   base, base2   - trajectory models (1->4 and 1->2) with no additional features
#   base_af       - a smaller model trained WITH 3 additional features (for the
#                   inference additional-feature alignment tests)
#   new_peaks     - 40 held-out peaks disjoint from the training set
#   new_af        - the matching additional features for new_peaks
# Calls skip_*() when prerequisites are missing. Built once per test run.
genome_traj_models <- local({
    cache <- NULL
    function() {
        if (!is.null(cache)) {
            return(cache)
        }
        skip_on_cran()
        skip_if_not_installed("misha")
        skip_if_not_installed("prego")
        skip_if_not_installed("misha.ext")
        groot <- genome_root_path()
        skip_if(!dir.exists(groot), paste("misha genome not found at", groot))
        cdir <- file.path(tools::R_user_dir("iceqream", "cache"), "gastrulation-example")
        needed <- file.path(cdir, c(
            "peak_intervals.rds", "atac_scores.rds", "motif_energies.rds", "additional_features.rds"
        ))
        skip_if(!all(file.exists(needed)), "cached vignette data not available")

        misha::gsetroot(groot)
        pk <- readr::read_rds(needed[1])
        at <- as.matrix(readr::read_rds(needed[2]))
        me <- readr::read_rds(needed[3])
        af_all <- readr::read_rds(needed[4])

        set.seed(1)
        ri <- sort(sample(nrow(pk), 600))
        ci <- sort(sample(ncol(me), 80))
        pk2 <- pk[ri, c("chrom", "start", "end", "peak_name")]
        at2 <- as.matrix(at[ri, ])
        me2 <- me[ri, ci]
        rownames(me2) <- pk2$peak_name

        mk <- function(bin_start, bin_end) {
            suppressWarnings(suppressMessages(regress_trajectory_motifs(
                peak_intervals = pk2, atac_scores = at2, motif_energies = me2,
                additional_features = NULL, norm_intervals = pk2,
                min_tss_distance = NULL, max_motif_num = 10,
                bin_start = bin_start, bin_end = bin_end,
                n_prego_motifs = 0, peaks_size = 500, seed = 60427
            )))
        }

        # Smaller model trained WITH additional features.
        af <- af_all[ri, 1:3]
        rownames(af) <- pk2$peak_name
        base_af <- suppressWarnings(suppressMessages(regress_trajectory_motifs(
            peak_intervals = pk2, atac_scores = at2, motif_energies = me2[, 1:40],
            additional_features = af, norm_intervals = pk2,
            min_tss_distance = NULL, max_motif_num = 8,
            n_prego_motifs = 0, peaks_size = 500, seed = 60427
        )))

        # Held-out peaks disjoint from the training set, with matching features.
        newi <- sort(sample(setdiff(seq_len(nrow(pk)), ri), 40))
        new_peaks <- pk[newi, c("chrom", "start", "end", "peak_name")]
        new_af <- af_all[newi, 1:3]
        rownames(new_af) <- new_peaks$peak_name

        cache <<- list(
            base = mk(1, 4), base2 = mk(1, 2), base_af = base_af,
            new_peaks = new_peaks, new_af = new_af
        )
        cache
    }
})
