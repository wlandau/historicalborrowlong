test_that("hbl_metrics()", {
  borrow <- tibble::tibble(
    group = c(1, 1),
    rep = c(1, 2),
    rep_label = c("x", "y"),
    response_mean = c(1.5, 1.52),
    response_variance = c(1.75, 1.752)
  )
  pool <- tibble::tibble(
    group = c(1, 1),
    rep = c(1, 2),
    rep_label = c("x", "y"),
    response_mean = c(1, 1.2),
    response_variance = c(1.3, 1.32)
  )
  independent <- tibble::tibble(
    group = c(1, 1),
    rep = c(1, 2),
    rep_label = c("x", "y"),
    response_mean = c(2, 2.2),
    response_variance = c(1.8, 1.82)
  )
  out <- hbl_metrics(
    borrow = borrow,
    pool = pool,
    independent = independent
  )
  expect_equal(
    out,
    tibble::tibble(
      rep = c(1, 2),
      rep_label = c("x", "y"),
      mean_shift_ratio = c(0.5, 0.68),
      variance_shift_ratio = c(0.1, 0.136)
    )
  )
})
