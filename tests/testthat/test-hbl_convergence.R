test_that("hbl_convergence()", {
  skip_on_cran()
  set.seed(0)
  data <- hbl_sim_pool(
    n_study = 2,
    n_group = 2,
    n_patient = 5,
    n_rep = 3
  )$data
  tmp <- utils::capture.output(
    suppressWarnings(
      mcmc <- hbl_mcmc_pool(
        data,
        chains = 1,
        warmup = 10,
        iter = 20,
        seed = 0
      )
    )
  )
  out <- hbl_convergence(mcmc)
  expect_equal(dim(out), c(1, 3))
  expect_equal(
    sort(colnames(out)),
    sort(c("max_rhat", "min_ess_bulk", "min_ess_tail"))
  )
})
