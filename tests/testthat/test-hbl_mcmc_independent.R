test_that("hbl_mcmc_independent() with betas + ar1 + diagonal", {
  set.seed(0)
  data <- hbl_sim_independent(
    n_study = 2,
    n_group = 3,
    n_patient = 5,
    n_rep = 4,
    n_continuous = 1,
    n_binary = 0,
    s_alpha = 1,
    s_delta = 1,
    s_beta = 1,
    s_sigma = 1,
    constraint = FALSE
  )$data
  tmp <- utils::capture.output(
    suppressWarnings(
      out <- hbl_mcmc_independent(
        data,
        chains = 2,
        warmup = 10,
        iter = 20,
        seed = 0,
        constraint = FALSE,
        covariance_current = "ar1",
        covariance_historical = "diagonal"
      )
    )
  )
  lapply(out, function(x) true(is.numeric(x) && all(is.finite(x))))
  sigma_grid <- tidyr::expand_grid(
    study = seq_len(2),
    rep = seq_len(4)
  )
  exp <- c(
    ".chain",
    ".draw",
    ".iteration",
    "lp__",
    sprintf("alpha[%s]", seq_len(8)),
    sprintf("delta[%s]", seq_len(8)),
    sprintf("beta[%s]", seq_len(2)),
    sprintf(
      "sigma[%s,%s]",
      sigma_grid$study,
      sigma_grid$rep
    ),
    "rho_current[1]"
  )
  expect_equal(sort(colnames(out)), sort(exp))
})
