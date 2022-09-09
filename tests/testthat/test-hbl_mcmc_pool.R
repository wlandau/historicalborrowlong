test_that("hbl_mcmc_pool() unstructured + ar1", {
  set.seed(0)
  data <- hbl_sim_pool(
    n_study = 2,
    n_group = 3,
    n_patient = 5,
    n_rep = 4,
    n_continuous = 0,
    n_binary = 0,
    s_alpha = 1,
    s_delta = 1,
    s_beta = 1,
    s_sigma = 1,
    constraint = FALSE
  )$data
  tmp <- utils::capture.output(
    suppressWarnings(
      out <- hbl_mcmc_pool(
        data,
        chains = 2,
        warmup = 10,
        iter = 20,
        seed = 0,
        constraint = FALSE,
        covariance_current = "unstructured",
        covariance_historical = "ar1"
      )
    )
  )
  lapply(out, function(x) true(is.numeric(x) && all(is.finite(x))))
  sigma_grid <- tidyr::expand_grid(
    study = seq_len(2),
    rep = seq_len(4)
  )
  lambda_current_grid <- tidyr::expand_grid(
    study = seq_len(1),
    rep1 = seq_len(4),
    rep2 = seq_len(4)
  )
  lambda_current_grid <- dplyr::filter(
    lambda_current_grid,
    rep1 >= rep2,
    rep1 + rep2 > 2
  )
  exp <- c(
    ".chain",
    ".draw",
    ".iteration",
    "lp__",
    sprintf("alpha[%s]", seq_len(4)),
    sprintf("delta[%s]", seq_len(8)),
    sprintf(
      "sigma[%s,%s]",
      sigma_grid$study,
      sigma_grid$rep
    ),
    sprintf(
      "lambda_current[%s,%s,%s]",
      lambda_current_grid$study,
      lambda_current_grid$rep1,
      lambda_current_grid$rep2
    ),
    "rho_historical[1]"
  )
  expect_equal(sort(colnames(out)), sort(exp))
})
