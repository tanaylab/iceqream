library(iceqream)
gsetroot("example_data/mm10")

normalization_intervals <- readr::read_tsv("example_data/gastrulation-example/gastrulation_intervals.tsv", show_col_types = FALSE)
motif_energies_raw <- compute_motif_energies(normalization_intervals, normalize = FALSE)

quants <- seq(0.98, 1, by = 0.001)
mq <- matrixStats::colQuantiles(motif_energies_raw, probs = quants)
colnames(mq) <- quants

mouse_db_quantiles <- mq
usethis::use_data(mouse_db_quantiles, overwrite = TRUE)
