test_that("hbl_mcmc_sge() + diagonal + unstructured", {
  skip_on_cran()
  set.seed(0)
  data <- hbl_sim_hierarchical(
    n_study = 2,
    n_group = 3,
    n_patient = 5,
    n_rep = 4,
    n_continuous = 0,
    n_binary = 0,
    s_mu = 1,
    s_tau = 1,
    s_delta = 1,
    s_beta = 1,
    s_sigma = 1,
    constraint = FALSE
  )$data
  tmp <- suppressMessages(
    utils::capture.output(
      suppressWarnings(
        all <- hbl_mcmc_sge(
          data,
          chains = 2,
          cores = 1,
          warmup = 10,
          iter = 20,
          seed = 0,
          constraint = FALSE,
          covariance_current = "diagonal",
          covariance_historical = "unstructured",
          scheduler = "local"
        )
      )
    )
  )
  # hierarchical
  out <- all$hierarchical
  expect_equal(attr(out, "model"), "hierarchical")
  lapply(out, function(x) true(is.numeric(x) && all(is.finite(x))))
  sigma_grid <- tidyr::expand_grid(
    study = seq_len(2),
    rep = seq_len(4)
  )
  lambda_historical_grid <- tidyr::expand_grid(
    study = seq_len(1),
    rep1 = seq_len(4),
    rep2 = seq_len(4)
  )
  lambda_historical_grid <- dplyr::filter(
    lambda_historical_grid,
    rep1 >= rep2,
    rep1 + rep2 > 2
  )
  exp <- c(
    ".chain",
    ".draw",
    ".iteration",
    "lp__",
    sprintf("alpha[%s]", seq_len(8)),
    sprintf("delta[%s]", seq_len(8)),
    sprintf(
      "sigma[%s,%s]",
      sigma_grid$study,
      sigma_grid$rep
    ),
    sprintf(
      "lambda_historical[%s,%s,%s]",
      lambda_historical_grid$study,
      lambda_historical_grid$rep1,
      lambda_historical_grid$rep2
    ),
    sprintf("mu[%s]", seq_len(4)),
    sprintf("tau[%s]", seq_len(4))
  )
  expect_equal(sort(colnames(out)), sort(exp))
  # independent
  out <- all$independent
  expect_equal(attr(out, "model"), "independent")
  lapply(out, function(x) true(is.numeric(x) && all(is.finite(x))))
  sigma_grid <- tidyr::expand_grid(
    study = seq_len(2),
    rep = seq_len(4)
  )
  lambda_historical_grid <- tidyr::expand_grid(
    study = seq_len(1),
    rep1 = seq_len(4),
    rep2 = seq_len(4)
  )
  lambda_historical_grid <- dplyr::filter(
    lambda_historical_grid,
    rep1 >= rep2,
    rep1 + rep2 > 2
  )
  exp <- c(
    ".chain",
    ".draw",
    ".iteration",
    "lp__",
    sprintf("alpha[%s]", seq_len(8)),
    sprintf("delta[%s]", seq_len(8)),
    sprintf(
      "sigma[%s,%s]",
      sigma_grid$study,
      sigma_grid$rep
    ),
    sprintf(
      "lambda_historical[%s,%s,%s]",
      lambda_historical_grid$study,
      lambda_historical_grid$rep1,
      lambda_historical_grid$rep2
    )
  )
  expect_equal(sort(colnames(out)), sort(exp))
  # pool
  out <- all$pool
  expect_equal(attr(out, "model"), "pool")
  lapply(out, function(x) true(is.numeric(x) && all(is.finite(x))))
  sigma_grid <- tidyr::expand_grid(
    study = seq_len(2),
    rep = seq_len(4)
  )
  lambda_historical_grid <- tidyr::expand_grid(
    study = seq_len(1),
    rep1 = seq_len(4),
    rep2 = seq_len(4)
  )
  lambda_historical_grid <- dplyr::filter(
    lambda_historical_grid,
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
      "lambda_historical[%s,%s,%s]",
      lambda_historical_grid$study,
      lambda_historical_grid$rep1,
      lambda_historical_grid$rep2
    )
  )
  expect_equal(sort(colnames(out)), sort(exp))
})

test_that("hbl_mcmc_sge_model() bad model", {
  skip_on_cran()
  expect_error(
    hbl_mcmc_sge_model("nope", list()),
    class = "hbl_error"
  )
})
