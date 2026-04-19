# Shared fixtures for tests. Loaded automatically by testthat via "helper-*.R".

# Larger traj_model with real motif-motif interaction structure in the response,
# used as a regression baseline for add_interactions.
create_interaction_traj_model <- function(n_peaks = 200, n_motifs = 10, train_frac = 0.7, seed = 1) {
    set.seed(seed)
    n_train <- round(n_peaks * train_frac)

    motif_names <- paste0("motif", seq_len(n_motifs))
    norm_energies <- matrix(runif(n_peaks * n_motifs, 0, 10), nrow = n_peaks, ncol = n_motifs)
    colnames(norm_energies) <- motif_names
    rownames(norm_energies) <- paste0("peak", seq_len(n_peaks))

    # Response has both linear contributions (first 3 motifs) and a real
    # motif-motif interaction (motif1 * motif2), so add_interactions has
    # something to find.
    linear_signal <- 0.4 * norm_energies[, 1] - 0.3 * norm_energies[, 2] + 0.25 * norm_energies[, 3]
    interaction_signal <- 0.05 * norm_energies[, 1] * norm_energies[, 2]
    y <- norm01(linear_signal + interaction_signal + rnorm(n_peaks, sd = 0.3))

    model_features <- create_logist_features(norm_energies)
    type_vec <- c(rep("train", n_train), rep("test", n_peaks - n_train))

    X_train <- model_features[type_vec == "train", , drop = FALSE]
    y_train <- y[type_vec == "train"]

    lambda_val <- 1e-5
    alpha_val <- 1
    glm_model <- suppressWarnings(glmnet::glmnet(
        x = X_train, y = y_train,
        family = binomial(link = "logit"),
        alpha = alpha_val, lambda = lambda_val, seed = seed
    ))
    pred <- logist(glmnet::predict.glmnet(glm_model, newx = model_features, type = "link", s = lambda_val))[, 1]

    motif_models <- lapply(motif_names, function(nm) {
        list(
            pssm = data.frame(pos = 0:4,
                A = c(0.7, 0.1, 0.1, 0.1, 0.25),
                C = c(0.1, 0.7, 0.1, 0.1, 0.25),
                G = c(0.1, 0.1, 0.7, 0.1, 0.25),
                T = c(0.1, 0.1, 0.1, 0.7, 0.25)),
            spat = data.frame(bin = numeric(0), spat_factor = numeric(0)),
            spat_min = numeric(0), spat_max = numeric(0), seq_length = numeric(0)
        )
    })
    names(motif_models) <- motif_names

    peak_intervals <- data.frame(
        chrom = rep("chr1", n_peaks),
        start = seq(1, by = 500, length.out = n_peaks),
        end = seq(500, by = 500, length.out = n_peaks),
        peak_name = paste0("peak", seq_len(n_peaks)),
        stringsAsFactors = FALSE
    )

    TrajectoryModel(
        model = glm_model,
        motif_models = motif_models,
        coefs = iceqream:::get_model_coefs(glm_model, s = lambda_val),
        normalized_energies = norm_energies,
        model_features = model_features,
        type = type_vec,
        diff_score = y,
        predicted_diff_score = pred,
        initial_prego_models = list(),
        peak_intervals = peak_intervals,
        normalization_intervals = data.frame(),
        additional_features = data.frame(row.names = paste0("peak", seq_len(n_peaks))),
        features_r2 = numeric(0),
        interactions = matrix(nrow = n_peaks, ncol = 0),
        params = list(
            lambda = lambda_val, alpha = alpha_val, seed = seed, peaks_size = 500
        )
    )
}

create_mock_traj_model <- function(n_peaks = 50, n_motifs = 2, train_frac = 0.8) {
    set.seed(42)
    n_train <- round(n_peaks * train_frac)
    n_test <- n_peaks - n_train

    motif_names <- paste0("motif", seq_len(n_motifs))

    norm_energies <- matrix(runif(n_peaks * n_motifs, 0, 10), nrow = n_peaks, ncol = n_motifs)
    colnames(norm_energies) <- motif_names
    rownames(norm_energies) <- paste0("peak", seq_len(n_peaks))

    model_features <- create_logist_features(norm_energies)

    y <- norm01(rowSums(norm_energies) + rnorm(n_peaks, sd = 0.5))

    type_vec <- c(rep("train", n_train), rep("test", n_test))

    X_train <- model_features[type_vec == "train", , drop = FALSE]
    y_train <- y[type_vec == "train"]

    lambda_val <- 0.01
    alpha_val <- 0.5

    glm_model <- suppressWarnings(glmnet::glmnet(
        x = X_train,
        y = y_train,
        family = binomial(link = "logit"),
        alpha = alpha_val,
        lambda = lambda_val
    ))

    pred <- logist(glmnet::predict.glmnet(glm_model, newx = model_features, type = "link", s = lambda_val))[, 1]
    pred <- norm01(pred)
    pred <- rescale(pred, y)

    motif_models <- lapply(motif_names, function(nm) {
        pssm <- matrix(
            c(
                0.7, 0.1, 0.1, 0.1,
                0.1, 0.7, 0.1, 0.1,
                0.1, 0.1, 0.7, 0.1,
                0.1, 0.1, 0.1, 0.7,
                0.25, 0.25, 0.25, 0.25
            ),
            nrow = 5, ncol = 4, byrow = TRUE,
            dimnames = list(NULL, c("A", "C", "G", "T"))
        )
        list(
            pssm = as.data.frame(pssm) %>% dplyr::mutate(pos = 0:4) %>% dplyr::select(pos, A, C, G, T),
            spat = data.frame(bin = numeric(0), spat_factor = numeric(0)),
            spat_min = numeric(0),
            spat_max = numeric(0),
            seq_length = numeric(0)
        )
    })
    names(motif_models) <- motif_names

    peak_intervals <- data.frame(
        chrom = rep("chr1", n_peaks),
        start = seq(1, by = 500, length.out = n_peaks),
        end = seq(500, by = 500, length.out = n_peaks),
        peak_name = paste0("peak", seq_len(n_peaks)),
        stringsAsFactors = FALSE
    )

    coefs_df <- iceqream:::get_model_coefs(glm_model, s = lambda_val)

    TrajectoryModel(
        model = glm_model,
        motif_models = motif_models,
        coefs = coefs_df,
        normalized_energies = norm_energies,
        model_features = model_features,
        type = type_vec,
        diff_score = y,
        predicted_diff_score = pred,
        initial_prego_models = list(),
        peak_intervals = peak_intervals,
        normalization_intervals = data.frame(),
        additional_features = data.frame(row.names = paste0("peak", seq_len(n_peaks))),
        features_r2 = numeric(0),
        interactions = matrix(nrow = n_peaks, ncol = 0),
        params = list(
            lambda = lambda_val,
            alpha = alpha_val,
            seed = 60427,
            peaks_size = 500
        )
    )
}
