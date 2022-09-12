test_that("hbl_plot_borrow()", {
  skip_on_cran()
  set.seed(0)
  data <- hbl_sim_independent(
    n_study = 2,
    n_group = 2,
    n_patient = 5,
    n_rep = 3
  )$data
  tmp <- utils::capture.output(
    suppressWarnings(
      mcmc_borrow <- hbl_mcmc_hierarchical(
        data,
        chains = 1,
        warmup = 10,
        iter = 20,
        seed = 0
      )
    )
  )
  tmp <- utils::capture.output(
    suppressWarnings(
      mcmc_pool <- hbl_mcmc_pool(
        data,
        chains = 1,
        warmup = 10,
        iter = 20,
        seed = 0
      )
    )
  )
  tmp <- utils::capture.output(
    suppressWarnings(
      mcmc_independent <- hbl_mcmc_independent(
        data,
        chains = 1,
        warmup = 10,
        iter = 20,
        seed = 0
      )
    )
  )
  borrow <- hbl_summary(mcmc_borrow, data)
  pool <- hbl_summary(mcmc_pool, data)
  independent <- hbl_summary(mcmc_independent, data)
  out <- hbl_plot_borrow(
    borrow = borrow,
    pool = pool,
    independent = independent
  )
  expect_s3_class(out, "ggplot")
})
