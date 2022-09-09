test_that("stan_covariance()", {
  expect_equal(stan_covariance("unstructured"), 1L)
  expect_equal(stan_covariance("ar1"), 2L)
  expect_equal(stan_covariance("diagonal"), 3L)
  expect_error(stan_covariance("nope"), class = "hbl_error")
})

test_that("ar1_cholesky()", {
  skip("long test")
  rstan::expose_stan_functions(stanmodels$model)
  rho <- -0.9
  factor <- ar1_cholesky(n = 3, rho = rho)
  out <- factor %*% t(factor)
  expect_equal(diag(out), rep(1, 3))
  expect_equal(out[lower.tri(out)], out[upper.tri(out)])
  expect_equal(out[1, 2], rho)
  expect_equal(out[1, 3], rho ^ 2)
  expect_equal(out[2, 3], rho ^ 1)
})
