test_that("sim pool data", {
  set.seed(0)
  out <- hbl_sim_pool(
    n_study = 2,
    n_group = 3,
    n_patient = 5,
    n_rep = 2,
    n_continuous = 2,
    n_binary = 2
  )$data
  expect_equal(dim(out), c(40, 13))
  expect_equal(out$study, rep(c(1, 2, 2, 2), each = 10))
  expect_equal(out$group, rep(c(1, 1, 2, 3), each = 10))
  expect_equal(out$patient, rep(seq_len(20), each = 2))
  expect_equal(out$rep, rep(seq_len(2), times = 20))
  expect_true(is.numeric(out$response))
  expect_false(anyNA(out$response))
  cols <- c(
    "covariate_study1_continuous1",
    "covariate_study1_continuous2",
    "covariate_study1_binary1",
    "covariate_study1_binary2"
  )
  for (col in cols) {
    expect_true(any(abs(out[[col]][seq_len(10)]) > 0))
    expect_equal(out[[col]][seq(11, 40)], rep(0, 30))
  }
  cols <- c(
    "covariate_study2_continuous1",
    "covariate_study2_continuous2",
    "covariate_study2_binary1",
    "covariate_study2_binary2"
  )
  for (col in cols) {
    expect_true(any(abs(out[[col]][seq(11, 40)]) > 0))
    expect_equal(out[[col]][seq(1, 10)], rep(0, 10))
  }
  cols <- c(
    "covariate_study1_binary1",
    "covariate_study1_binary2",
    "covariate_study2_binary1",
    "covariate_study2_binary2"
  )
  for (col in cols) {
    expect_equal(sort(unique(out[[col]])), sort(unique(c(0, 1))))
  }
  cols <- c(
    "covariate_study1_continuous1",
    "covariate_study1_continuous2",
    "covariate_study2_continuous1",
    "covariate_study2_continuous2"
  )
  for (col in cols) {
    expect_false(identical(sort(unique(out[[col]])), sort(unique(c(0, 1)))))
  }
})

test_that("sim pool parameters with constraint", {
  set.seed(0)
  out <- hbl_sim_pool(
    n_study = 2,
    n_group = 3,
    n_patient = 5,
    n_rep = 3,
    n_continuous = 2,
    n_binary = 2,
    covariance_current = "unstructured",
    covariance_historical = "unstructured",
    constraint = TRUE
  )$parameters
  lapply(out, function(x) expect_true(is.numeric(x)))
  expect_equal(length(out$alpha), 3)
  expect_equal(length(out$delta), 4)
  expect_equal(length(out$beta), 8)
  expect_equal(dim(out$sigma), c(2, 3))
  expect_equal(dim(out$lambda_current), c(1, 3, 3))
  expect_equal(dim(out$lambda_historical), c(1, 3, 3))
  expect_null(out$rho_current)
  expect_null(out$rho_historical)
  expect_null(out$mu)
  expect_null(out$tau)
})

test_that("sim pool parameters without constraint", {
  set.seed(0)
  out <- hbl_sim_pool(
    n_study = 2,
    n_group = 3,
    n_patient = 5,
    n_rep = 3,
    n_continuous = 2,
    n_binary = 2,
    covariance_current = "unstructured",
    covariance_historical = "unstructured",
    constraint = FALSE
  )$parameters
  lapply(out, function(x) expect_true(is.numeric(x)))
  expect_equal(length(out$alpha), 3)
  expect_equal(length(out$delta), 6)
  expect_equal(length(out$beta), 8)
  expect_equal(dim(out$sigma), c(2, 3))
  expect_equal(dim(out$lambda_current), c(1, 3, 3))
  expect_equal(dim(out$lambda_historical), c(1, 3, 3))
  expect_null(out$rho_current)
  expect_null(out$rho_historical)
  expect_null(out$mu)
  expect_null(out$tau)
})

test_that("sim pool unstructured + ar1", {
  set.seed(0)
  out <- hbl_sim_pool(
    n_study = 2,
    n_group = 3,
    n_patient = 5,
    n_rep = 3,
    n_continuous = 2,
    n_binary = 2,
    constraint = FALSE,
    covariance_current = "unstructured",
    covariance_historical = "ar1"
  )$parameters
  lapply(out, function(x) expect_true(is.numeric(x)))
  expect_equal(length(out$alpha), 3)
  expect_equal(length(out$delta), 6)
  expect_equal(length(out$beta), 8)
  expect_equal(dim(out$sigma), c(2, 3))
  expect_equal(dim(out$lambda_current), c(1, 3, 3))
  expect_null(out$lambda_historical)
  expect_null(out$rho_current)
  expect_equal(length(out$rho_historical), 1)
  expect_null(out$mu)
  expect_null(out$tau)
})

test_that("sim pool ar1 + diagonal", {
  set.seed(0)
  out <- hbl_sim_pool(
    n_study = 2,
    n_group = 3,
    n_patient = 5,
    n_rep = 3,
    n_continuous = 2,
    n_binary = 2,
    constraint = FALSE,
    covariance_current = "ar1",
    covariance_historical = "diagonal"
  )$parameters
  lapply(out, function(x) expect_true(is.numeric(x)))
  expect_equal(length(out$alpha), 3)
  expect_equal(length(out$delta), 6)
  expect_equal(length(out$beta), 8)
  expect_equal(dim(out$sigma), c(2, 3))
  expect_null(out$lambda_current)
  expect_null(out$lambda_historical)
  expect_equal(length(out$rho_current), 1)
  expect_null(out$rho_historical)
  expect_null(out$mu)
  expect_null(out$tau)
})

test_that("sim pool matrices", {
  set.seed(0)
  out <- hbl_sim_pool(
    n_study = 2,
    n_group = 3,
    n_patient = 5,
    n_continuous = 2,
    n_binary = 2,
    constraint = FALSE
  )$matrices
  x_alpha <- out$x_alpha
  exp <- do.call(rbind, replicate(10, diag(4), simplify = FALSE))
  exp <- rbind(exp, matrix(0, nrow = 40, ncol = 4))
  dimnames(x_alpha) <- NULL
  expect_equal(x_alpha, exp)
  x_delta <- do.call(rbind, replicate(5, diag(4), simplify = FALSE))
  x_delta <- as.matrix(Matrix::bdiag(list(x_delta, x_delta)))
  x_delta <- rbind(matrix(0, nrow = 40, ncol = 8), x_delta)
  dimnames(out$x_delta) <- NULL
  expect_equal(out$x_delta, x_delta)
  out <- tibble::as_tibble(out$x_beta)
  cols <- c(
    "study1_covariate_study1_continuous1",
    "study1_covariate_study1_continuous2",
    "study1_covariate_study1_binary1",
    "study1_covariate_study1_binary2"
  )
  for (col in cols) {
    expect_true(any(abs(out[[col]][seq_len(20)]) > 0))
    expect_equal(out[[col]][seq(21, 80)], rep(0, 60))
  }
  cols <- c(
    "study2_covariate_study2_continuous1",
    "study2_covariate_study2_continuous2",
    "study2_covariate_study2_binary1",
    "study2_covariate_study2_binary2"
  )
  for (col in cols) {
    expect_true(any(abs(out[[col]][seq(21, 80)]) > 0))
    expect_equal(out[[col]][seq(1, 20)], rep(0, 20))
  }
  cols <- c(
    "study1_covariate_study1_binary1",
    "study1_covariate_study1_binary2",
    "study2_covariate_study2_binary1",
    "study2_covariate_study2_binary2"
  )
  for (col in cols) {
    expect_equal(length(unique(out[[col]])), 3)
  }
  cols <- c(
    "study1_covariate_study1_continuous1",
    "study1_covariate_study1_continuous2",
    "study2_covariate_study2_continuous1",
    "study2_covariate_study2_continuous2"
  )
  for (col in cols) {
    expect_gt(length(unique(out[[col]])), 3)
  }
})
