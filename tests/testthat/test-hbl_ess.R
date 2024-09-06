test_that("hbl_ess()", {
  set.seed(0L)
  data <- hbl_sim_independent(n_study = 3L)$data
  data$response[1L] <- NA_real_
  data$group <- sprintf("group%s", data$group)
  data$study <- sprintf("study%s", data$study)
  tmp <- utils::capture.output(
    suppressWarnings(
      pool <- hbl_mcmc_pool(
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
      hierarchical <- hbl_mcmc_hierarchical(
        data,
        chains = 1,
        warmup = 10,
        iter = 20,
        seed = 0
      )
    )
  )
  out <- hbl_ess(
    mcmc_pool = pool,
    mcmc_hierarchical = hierarchical,
    data = data
  )
  v0 <- c(
    mean(
      (
        (pool$`sigma[1,1]` ^ (-2)) +
          (pool$`sigma[2,1]` ^ (-2)) +
          (pool$`sigma[3,1]` ^ (-2))
      ) ^ (-1)
    ),
    mean(
      (
        (pool$`sigma[1,2]` ^ (-2)) +
          (pool$`sigma[2,2]` ^ (-2)) +
          (pool$`sigma[3,2]` ^ (-2))
      ) ^ (-1)
    ),
    mean(
      (
        (pool$`sigma[1,3]` ^ (-2)) +
          (pool$`sigma[2,3]` ^ (-2)) +
          (pool$`sigma[3,3]` ^ (-2))
      ) ^ (-1)
    ),
    mean(
      (
        (pool$`sigma[1,4]` ^ (-2)) +
          (pool$`sigma[2,4]` ^ (-2)) +
          (pool$`sigma[3,4]` ^ (-2))
      ) ^ (-1)
    )
  )
  expect_equal(out$v0, v0)
  exp <- c(
    mean(hierarchical[["tau[1]"]]^2 + var(hierarchical[["mu[1]"]])),
    mean(hierarchical[["tau[2]"]]^2 + var(hierarchical[["mu[2]"]])),
    mean(hierarchical[["tau[3]"]]^2 + var(hierarchical[["mu[3]"]])),
    mean(hierarchical[["tau[4]"]]^2 + var(hierarchical[["mu[4]"]]))
  )
  expect_equal(order(out$v_tau), order(exp))
  expect_equal(out$n, c(199L, rep(200L, 3L)))
  expect_equal(out$weight, out$v0 / out$v_tau)
  expect_equal(out$ess, out$n * out$weight)
})
