library(iceqream)

clust_map_scenic <- motif_db %>%
    filter(!is.na(cluster)) %>%
    distinct(motif, cluster) %>%
    select(motif, clust = cluster)

set.seed(60427)
chosen_motifs_scenic <- clust_map_scenic %>%
    sample_n(n()) %>%
    mutate(bad = grepl("Unknown", motif) | grepl("SeqBias", motif) | grepl("bias", motif) | grepl("Hotspot", motif) | grepl("Satellite", motif) | grepl("MA00", motif)) %>%
    arrange(clust, bad) %>%
    group_by(clust) %>%
    slice(1) %>%
    pull(motif)
head(chosen_motifs_scenic)

motif_db_scenic_clusters <- motif_db %>% filter(motif %in% chosen_motifs_scenic)

usethis::use_data(motif_db_scenic_clusters, overwrite = TRUE)
