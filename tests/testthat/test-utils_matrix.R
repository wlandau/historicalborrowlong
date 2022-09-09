test_that("get_x_alpha(constraint = FALSE)", {
  data <- tibble::tibble(
    study = rep(c(1, 2), each = 8),
    group = rep(c(1, 2, 1, 2), each = 4),
    rep = rep(seq_len(4), times = 4)
  )
  x <- get_x_alpha(data, constraint = FALSE)
  expect_equal(dim(x), c(16, 8))
  zero <- matrix(0, nrow = 4, ncol = 4)
  block <- rbind(diag(4), zero)
  exp <- as.matrix(Matrix::bdiag(block, block))
  dimnames(x) <- NULL
  dimnames(exp) <- NULL
  expect_equal(x, exp)
})

test_that("get_x_alpha(constraint = TRUE)", {
  data <- tibble::tibble(
    study = rep(c(1, 2), each = 8),
    group = rep(c(1, 2, 1, 2), each = 4),
    rep = rep(seq_len(4), times = 4)
  )
  x <- get_x_alpha(data, constraint = TRUE)
  expect_equal(dim(x), c(16, 8))
  zero <- matrix(0, nrow = 4, ncol = 4)
  block <- rbind(diag(4), zero)
  exp <- as.matrix(Matrix::bdiag(block, block))
  exp[5, 1] <- 1
  exp[13, 5] <- 1
  dimnames(x) <- NULL
  dimnames(exp) <- NULL
  expect_equal(x, exp)
})

test_that("get_x_delta(constraint = FALSE)", {
  data <- tibble::tibble(
    study = c(1, 1, 1, 2, 2, 3, 3),
    group = c(1, 1, 2, 1, 2, 1, 3)
  )
  data <- tidyr::expand_grid(data, rep = seq_len(2))
  x <- get_x_delta(data, constraint = FALSE)
  expect_equal(dim(x), c(14, 6))
  expect_equal(x[, 1, drop = TRUE], as.integer(seq_len(14) == 5))
  expect_equal(x[, 2, drop = TRUE], as.integer(seq_len(14) == 6))
  expect_equal(x[, 3, drop = TRUE], as.integer(seq_len(14) == 9))
  expect_equal(x[, 4, drop = TRUE], as.integer(seq_len(14) == 10))
  expect_equal(x[, 5, drop = TRUE], as.integer(seq_len(14) == 13))
  expect_equal(x[, 6, drop = TRUE], as.integer(seq_len(14) == 14))
})

test_that("get_x_delta(constraint = TRUE)", {
  data <- tibble::tibble(
    study = c(1, 1, 1, 2, 2, 3, 3),
    group = c(1, 1, 2, 1, 2, 1, 3)
  )
  data <- tidyr::expand_grid(data, rep = seq_len(2))
  x <- get_x_delta(data, constraint = TRUE)
  expect_equal(dim(x), c(14, 3))
  expect_equal(x[, 1, drop = TRUE], as.integer(seq_len(14) == 6))
  expect_equal(x[, 2, drop = TRUE], as.integer(seq_len(14) == 10))
  expect_equal(x[, 3, drop = TRUE], as.integer(seq_len(14) == 14))
})

test_that("get_x_beta() null case", {
  data <- tibble::tibble(
    study = c(1, 1, 1, 2, 2, 3, 3),
    group = c(1, 1, 2, 1, 2, 1, 3)
  )
  data <- tidyr::expand_grid(data, rep = seq_len(2))
  x <- get_x_alpha(data, constraint = FALSE)
  x_delta <- get_x_delta(data, constraint = FALSE)
  x <- get_x_beta(data, x_alpha, x_beta)
  expect_equal(x, matrix(0, nrow = nrow(data), ncol = 0))
})

test_that("get_x_beta() non-null but repeating", {
  data <- tibble::tibble(
    study = c(1, 1, 1, 1, 2, 2, 2, 2),
    group = c(1, 1, 2, 2, 1, 1, 1, 1),
    covariate_a = seq_len(8),
    covariate_b = c(rep(0, 7), 1),
    nope = seq_len(8)
  )
  data <- expand_grid(data, rep = seq_len(2))
  x_alpha <- get_x_alpha(data, constraint = FALSE)
  x_delta <- get_x_delta(data, constraint = TRUE)
  x <- get_x_beta(data, x_alpha, x_delta)
  expect_equal(dim(x), c(16, 3))
  expect_equal(
    sort(colnames(x)),
    sort(c("study1_covariate_a", "study2_covariate_a", "study2_covariate_b"))
  )
  expect_equal(
    x[, 1, drop = TRUE],
    c(as.numeric(scale(data$covariate_a[seq_len(8)])), rep(0, 8))
  )
  expect_equal(
    x[, 2, drop = TRUE],
    c(rep(0, 8), as.numeric(scale(data$covariate_a[seq(9, 16)])))
  )
  expect_equal(
    x[, 3, drop = TRUE],
    c(rep(0, 8), as.numeric(scale(rep(c(0, 0, 0, 1), each = 2))))
  )
})

test_that("get_x_beta() non-null non-repeating", {
  data <- tibble::tibble(
    study = c(1, 1, 1, 1, 2, 2, 2, 2),
    group = c(1, 1, 2, 2, 1, 1, 1, 1),
    covariate_a = c(seq_len(4), c(1, 2, 1, 2)),
    covariate_b = c(0, 0, 1, 0, 1, 0, 0, 0)
  )
  data <- expand_grid(data, rep = seq_len(2))
  x_alpha <- get_x_alpha(data, constraint = FALSE)
  x_delta <- get_x_delta(data, constraint = TRUE)
  x <- get_x_beta(data, x_alpha, x_delta)
  expect_equal(dim(x), c(16, 4))
  expect_equal(
    sort(colnames(x)),
    sort(
      c(
        "study1_covariate_a",
        "study1_covariate_b",
        "study2_covariate_a",
        "study2_covariate_b"
      )
    )
  )
  expect_equal(
    x[, 1, drop = TRUE],
    c(as.numeric(scale(data$covariate_a[seq_len(8)])), rep(0, 8))
  )
  expect_equal(
    x[, 2, drop = TRUE],
    c(as.numeric(scale(rep(c(0, 0, 1, 0), each = 2))), rep(0, 8))
  )
  expect_equal(
    x[, 3, drop = TRUE],
    c(rep(0, 8), as.numeric(scale(rep(c(1, 2, 1, 2), each = 2))))
  )
  expect_equal(
    x[, 4, drop = TRUE],
    c(rep(0, 8), as.numeric(scale(rep(c(1, 0, 0, 0), each = 2))))
  )
})

test_that("get_x_beta() with all of a factor level missing", {
  data <- tibble::tibble(
    response = c(NA_real_, rep(0, 7)),
    study = c(1, 1, 1, 1, 2, 2, 2, 2),
    group = c(1, 1, 2, 2, 1, 1, 1, 1),
    covariate_level = c(1, 0, 0, 0, 1, 0, 0, 0)
  )
  data <- expand_grid(data, rep = seq_len(2))
  x_alpha <- get_x_alpha(data, constraint = FALSE)
  x_delta <- get_x_delta(data, constraint = TRUE)
  # Without accounting for the NA in the response,
  # x would have 2 columns instead of 1.
  x <- get_x_beta(data, x_alpha, x_delta)
  expect_equal(dim(x), c(16, 1))
  expect_equal(
    sort(colnames(x)),
    sort("study2_covariate_level")
  )
  expect_equal(
    x[, 1, drop = TRUE],
    c(rep(0, 8), as.numeric(scale(rep(c(1.5, -0.5, -0.5, -0.5), each = 2))))
  )
})
