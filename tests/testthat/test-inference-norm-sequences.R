# Regression guard for the normalization-sequence misalignment in
# extract_traj_model_sequences(): `intervalID` indexes the DEDUPED interval set,
# so the sequence vector it subscripts has to be extracted from that same set.
# Extracting the full (duplicated) set instead handed every norm interval at or
# after the first duplicate somebody else's sequence.

test_that("extract_traj_model_sequences aligns norm sequences when norm overlaps peaks", {
    skip_on_cran()
    skip_if_not(dir.exists(genome_root_path()), "misha genome unavailable")
    misha::gsetroot(genome_root_path())

    chrom <- misha::gintervals.all()$chrom[1]
    ivs <- data.frame(
        chrom = chrom,
        start = seq(3000000, by = 5000, length.out = 8),
        stringsAsFactors = FALSE
    )
    ivs$end <- ivs$start + 500

    peaks <- ivs[1:5, , drop = FALSE]
    # norm set deliberately shares intervals 3:5 with the peaks, then adds 3 new
    norms <- ivs[3:8, , drop = FALSE]

    tm <- methods::new("TrajectoryModel",
        peak_intervals = peaks,
        normalization_intervals = norms,
        diff_score = rep(0, nrow(peaks)),
        params = list(peaks_size = 500)
    )
    got <- suppressMessages(extract_traj_model_sequences(tm, peaks))

    expect_equal(unname(got$sequences), unname(prego::intervals_to_seq(peaks, 500)))
    expect_equal(unname(got$norm_sequences), unname(prego::intervals_to_seq(norms, 500)))
})
