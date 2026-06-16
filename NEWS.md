# iceqream 0.0.8

## Breaking changes

* `regress_trajectory_motifs()` and `regress_trajectory_motifs_manifold()` no
  longer accept `normalize_energies` / `norm_motif_energies`. They were ignored
  (energy normalization happens in `compute_motif_energies()`); remove them from
  your calls.

## Bug fixes

* `filter_traj_model()` no longer crashes on small or aggressively-filtered
  models: it is a no-op on models with fewer than 2 motifs, never removes every
  motif (it keeps at least one so the model stays fittable), and its bits-based
  removal is now applied even when no feature also fails the R^2 threshold
  (previously it was silently skipped in that case).
* `add_interactions()` now rescales `predicted_diff_score` to the trajectory's
  accessibility-difference range (`rescale_pred = TRUE`), matching
  `regress_trajectory_motifs()`. Previously the default left predictions on the
  raw `[0, 1]` logistic scale, so adding interactions silently changed the units
  of `predicted_diff_score` (reported R^2 was unaffected, as it is
  correlation-based).
* `pbm_list.gextract()`, `pbm.gextract()`, and
  `pbm_list.multi_traj.gextract_energy()` now extract sequences at the model's
  trained `size` rather than the input intervals' width. Previously, passing
  intervals whose width differed from the training `peaks_size` (e.g. 300bp
  peaks for a model trained at 500bp) produced incorrect energies in the
  documented inference path.
* `merge_trajectory_motifs()` now removes the merged-away motifs from the model
  entirely (previously they were dropped from the model features but left behind
  in `@motif_models` and `@normalized_energies`).
* `regress_trajectory_motifs()` and `regress_trajectory_motifs_manifold()` no
  longer error when `additional_features` is omitted (`NULL`).
* `add_interactions(min_signal_correlation = ...)` no longer drops all
  interactions when none of them correlate with the signal.
* `compute_motif_energies(db_quantiles = ...)` is now documented as using a
  fixed-range normalization that can differ slightly from the default path.

## Improvements

* Faster model reports and `add_interactions()` on large motif sets; results are
  unchanged.

## Breaking changes

* Removed three exported functions that had no internal callers, no tests,
  and were never referenced from any vignette or public workflow:
  `add_motif_models_to_traj()`, `adjust_energies()`, `adjust_motif_seq_lengths()`.
* `iq_regression(include_interactions = TRUE)` now uses Akhiad-inspired
  tight single-pass defaults instead of the previous loose defaults:
    * `interaction_threshold`: `0.01` (was `0.001`).
    * `interaction_only_sig_motifs`: `TRUE` (new argument, previously
      hard-coded to the `add_interactions()` default of `FALSE`).
  For exact numeric parity with 0.0.6 / paper / pycqream Phase 2 results,
  pass `interaction_threshold = 0.001, interaction_only_sig_motifs = FALSE`
  explicitly. Benchmarked on the gastrulation vignette data
  (5000 peaks, 21862-motif db) under two configurations:
    * At the default `max_n_interactions` cap (both configs bounded at
      260 interactions): R^2 train/test 0.531/0.351 (0.0.7) vs
      0.529/0.352 (0.0.6). Parity within measurement noise — expected,
      because the top-N correlation selection dominates when both
      configurations hit the same cap.
    * With `max_n_interactions = 5000` (unbounded, probes the selection
      criteria directly): 0.0.7 produces **492 interactions** at
      R^2 test 0.367; 0.0.6 produces 517 interactions at R^2 test 0.360.
      The tight default is more parsimonious (25 fewer interactions)
      and generalizes slightly better (+0.007 test R^2) — this is the
      evidence that motivates the default change.
* `iq_regression(strategy = "progressive")` now errors at call time
  instead of warning. The default progressive builder
  (`default_score_split_features()`) causes silent test-R^2 collapse
  (~0.27 on gastrulation) because its helper-model predictions are
  defined only for training peaks and imputed to 0 at test. Errors are
  louder than warnings. **The `strategy = "progressive"` path inside
  `iq_regression` is currently disabled — it will error unconditionally.**
  Power users who want the two-pass workflow should call
  `add_interactions_progressive()` directly on a model that already spans
  train + test peaks and supply `additional_features` (with `base_pred` /
  `end_pred`) for both splits at inference time.

## New features

* Exported `add_interactions_progressive()` — N-pass interaction selection
  with a pluggable between-pass feature-injection hook. Mirrors the manual
  workflow used in Akhiad's analysis notebooks.
* Exported `default_score_split_features()` — helper that builds
  `base_pred`, `end_pred`, `pred_diff_e_b` engineered additional features
  from start- and end-bin ATAC scores by relearning the trajectory model
  twice. Using this as an `additional_features_builder` requires
  propagating the helper-model predictions to test peaks yourself — see
  `?default_score_split_features`. Not safe to use inside `iq_regression`
  until test-time propagation lands (see the breaking-change note above
  — `iq_regression(strategy = "progressive")` errors unconditionally).
* `iq_regression()` gains a `strategy = c("single", "progressive")`
  argument. `"single"` is the default. `"progressive"` is **currently
  disabled** (errors at call time — see the breaking-change note above);
  reserved for a future release once test-time propagation of the
  engineered features lands. Power users who want the two-pass workflow
  today should call `add_interactions_progressive()` directly on a
  traj_model that already spans train + test peaks and supply the
  engineered additional features explicitly at inference time.
* `add_interactions()` gains `interaction_scale_factor` (multiplier on
  the normalized interaction matrix) and `min_signal_correlation`
  (post-filter: drop interactions with |cor(col, diff_score[train])| <
  threshold × max(|cor|)). `min_signal_correlation = 1/8` mirrors the
  manual post-filter in Akhiad's notebook.

## Bug fixes

* Fix: `remove_interactions()` silently failed to strip the logist-expanded
  interaction columns from `@model_features` when `logist_interactions = TRUE`,
  leaving orphan columns behind after a subsequent `add_interactions(force = TRUE)`.
* Fix: `infer_trajectory_motifs()` produced a
  "not valid for @'interactions'" S4 validity error on any
  interaction-augmented model because `rbind(matrix, data.frame)` was
  assigned back to the matrix-only slot. Latent until now — no prior test
  exercised `iq_regression(include_interactions = TRUE)` through
  inference on the final model.
* Fix: `escape_vars()` used a regex that made the escaped dot optional
  (`\\.?`), causing name-prefix matches to incorrectly hit variables that
  shared a prefix with the intended target (e.g. `motif.low-energy` also
  matched `motiflow-energy`).
* Fix: `cli_warn` in `infer_trajectory_motifs` printed the full additional-
  features data frame instead of just the feature names.
* Fix: `plot_traj_model_clusters_report()` would silently `unlink(dir, recursive = TRUE)`
  any user-supplied path — now refuses to delete `/`, `$HOME`, `getwd()`, or
  any path shallower than 3 components.
* Fix: `regress_trajectory_motifs()` and `regress_trajectory_motifs_manifold()`
  (and therefore `iq_regression()`) now clamp `kmer_sequence_length` to
  `peaks_size` when it would exceed the extracted sequence length. Previously
  the mismatch (default `kmer_sequence_length = 300`, user `peaks_size < 300`)
  surfaced deep inside `distill_motifs()`'s parallel `plyr::llply` loop as the
  opaque `task 1 failed - "kmer_sequence_length cannot be greater than the
  length"` originating from `prego::regress_pwm`.

## Improvements

* `relearn_traj_model()` with `use_cv = TRUE` no longer refits `glmnet` after
  `cv.glmnet`; the path fit stored inside `cv_model$glmnet.fit` is reused.
  Saves one full fit per CV-enabled relearn.
* Finished the partial `tidyr::gather` → `pivot_longer` and `plyr::*` →
  `purrr::*` migrations started in 0.0.6. `plyr::llply(..., .parallel = TRUE)`
  is retained in the hot paths where parallelism is active
  (`distill-motifs.R`, `distill-multi-traj.R`, `inference.R`, `PBM.R`).

## Deferred

* Integration with misha's `feat/glm-pred` branch
  (`glm_batch_quantiles()`, `glm_extract_features()`, `glm_pred.create()`
  virtual tracks): deferred until the branch merges to misha master.
  Expected wins: 5–10× on per-motif quantile normalization, faster
  feature-matrix construction, and post-fit genome-wide scoring without
  per-position R round-trip. See `.a5c/runs/deep-code-review-2026-04-19.md`
  §3.
* `cv.glmnet` parallelization via `doParallel`: deferred due to the
  OpenMP-oversubscription / deadlock risk when glmnet's internal
  parallelism stacks with R-level `foreach`. Safe rollout requires
  forcing `OMP_NUM_THREADS = 1` inside each worker and exposing it as
  opt-in via a package option. See `.a5c/runs/deep-code-review-2026-04-19.md`
  §5.2.
* Test-time `base_pred` / `end_pred` propagation so that
  `strategy = "progressive"` with `default_score_split_features` doesn't
  overfit on test. Would require inferring the base-only / end-only
  helper models on test peaks inside `infer_trajectory_motifs`.

# iceqream 0.0.6

## Bug fixes

* Fix: `prob1_thresh` parameter in `preprocess_data` was silently ignored due to being passed as NULL instead of the computed value.
* Fix: `select_motifs_by_correlation` fallback path always returned zero motifs.
* Fix: `distill_traj_model` selected the weakest motif per cluster instead of the strongest (wrong sort order).
* Fix: `regress_trajectory_motifs_manifold` modified the wrong object (`mm` instead of `mm_new`), discarding motif model updates.
* Fix: `iq_regression` applied redundant TSS distance filtering twice (once in outer and inner call).
* Fix: divide-by-zero in `create_specific_terms` and `create_features_terms` when interaction normalization factors are zero.

## Improvements

* New exported functions: `norm_energy_dataset` and `strip_traj_model`.
* Migrated deprecated ggplot2 APIs: `size` aesthetic in line geoms replaced with `linewidth`; numeric `legend.position` replaced with `legend.position.inside`.
* Migrated `tidyr::gather` calls to `tidyr::pivot_longer`.
* Virtual tracks created during `preprocess_data`, `normalize_regional`, and plotting functions are now cleaned up automatically via `on.exit()`.
* Graphics devices opened during report generation are now closed reliably via `on.exit()`.
* S4 validity checks added for `IQFeature`, `PBM`, and `TrajectoryModel` classes.
* Fixed typo in function name: `create_specifc_terms` renamed to `create_specific_terms`.

## Tests

* Added 127 validation tests covering utility functions, S4 validity, refactored helpers, and API migration correctness (total: 575 tests).

# iceqream 0.0.5

## Bug fixes

* Fix: `plot_normalization_scatters` no longer crashes when `anchor_cell_type` is NULL.
* Fix: `preprocess_data` now correctly excludes the marginal track from cell type names, preventing size mismatch errors.
* Fix: `names(new_models_full)` was incorrectly assigned in `infer_trajectory_motifs_multi`, leaving `new_models_full` unnamed.
* Fix: S4 copy-on-modify bug in `regress_trajectory_motifs_manifold` where loop modifications to model params were silently discarded.
* Fix: undefined variable `symmetrize_spat` in `iq_regression` when learning prego motifs de-novo.
* Fix: `add_type` partial argument match corrected to `add_types` in `infer_trajectory_motifs` and `get_model_coefs`.

## Improvements

* `preprocess_data`: renamed `peaks` parameter to `peak_intervals` for API consistency.
* `regress_trajectory_motifs`: `bin_start` and `bin_end` now accept column names (character) in addition to integer indices.
* `validate_atac_scores`: relaxed validation to allow `bin_start == bin_end` only when they refer to different columns (removed `bin_start < bin_end` constraint).
* Improved documentation linking `preprocess_data()$atac_norm_prob` to the `atac_scores` parameter in `iq_regression` and `regress_trajectory_motifs`.
* Added input validation for `norm_energy_matrix` (matrix coercion, column name alignment).
* Added input validation for interaction terms in `infer_trajectory_motifs`.

## Tests

* Added 423 unit tests covering energy utilities, S4 classes (PBM, TrajectoryModel, IQFeature family), validation functions, and preprocessing helpers.

# iceqream 0.0.4

* Fix: `iq_regression` now filters also test set peaks that are too close to TSS.

# iceqream 0.0.3

* Paper version.
