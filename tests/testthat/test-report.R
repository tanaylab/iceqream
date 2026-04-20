test_that("plot_traj_model_clusters_report refuses to recursively delete dangerous dirs", {
    skip_if_not_installed("ComplexHeatmap")
    tm <- create_mock_traj_model()

    expect_error(
        plot_traj_model_clusters_report(tm, dir = "/"),
        "Refusing to recursively delete"
    )
    expect_error(
        plot_traj_model_clusters_report(tm, dir = normalizePath("~", mustWork = FALSE)),
        "Refusing to recursively delete"
    )
})
