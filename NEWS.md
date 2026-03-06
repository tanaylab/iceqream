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
