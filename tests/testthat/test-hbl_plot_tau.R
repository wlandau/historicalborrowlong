test_that("hbl_plot_tau()", {
  skip_on_cran()
  set.seed(0)
  data <- hbl_sim_independent(n_continuous = 2)$data
  tmp <- utils::capture.output(
    suppressWarnings(
      mcmc <- hbl_mcmc_hierarchical(
        data,
        chains = 1,
        warmup = 10,
        iter = 20,
        seed = 0
      )
    )
  )
  out <- hbl_plot_tau(mcmc)
  expect_s3_class(out, "ggplot")
})
