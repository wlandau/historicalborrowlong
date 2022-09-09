test_that("hb_sim_lambda(n_matrix = 1)", {
  set.seed(0)
  out <- hbl_sim_lambda(n_matrix = 1, n_rep = 4, s_lambda = 1)
  expect_equal(dim(out), c(1, 4, 4))
  out <- out[1,, .drop = TRUE] # nolint
  expect_true(all(abs(out[lower.tri(out)]) > 0))
  expect_true(all(abs(out[upper.tri(out)]) == 0))
  cor <- out %*% t(out)
  expect_equal(cor, t(cor))
  expect_equal(diag(cor), rep(1, nrow(cor)))
})

test_that("hb_sim_lambda(n_matrix = 3)", {
  set.seed(0)
  array <- hbl_sim_lambda(n_matrix = 3, n_rep = 4, s_lambda = 1)
  expect_equal(dim(array), c(3, 4, 4))
  for (index in seq_len(3)) {
    out <- array[index,, .drop = TRUE] # nolint
    expect_equal(dim(out), c(4, 4))
    expect_true(all(abs(out[lower.tri(out)]) > 0))
    expect_true(all(abs(out[upper.tri(out)]) == 0))
    cor <- out %*% t(out)
    expect_equal(cor, t(cor))
    expect_equal(diag(cor), rep(1, nrow(cor)))
  }
})

test_that("ar1_correlation()", {
  rho <- -0.9
  out <- ar1_correlation(n = 3, rho = rho)
  expect_equal(diag(out), rep(1, 3))
  expect_equal(out[lower.tri(out)], out[upper.tri(out)])
  expect_equal(out[1, 2], rho)
  expect_equal(out[1, 3], rho ^ 2)
  expect_equal(out[2, 3], rho ^ 1)
})
