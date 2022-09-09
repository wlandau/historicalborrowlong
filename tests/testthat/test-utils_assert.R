test_that("true()", {
  expect_silent(true(TRUE))
  expect_error(true(FALSE), class = "hbl_error")
  expect_silent(true(c(2, 3), . > 1, . > 0))
  expect_error(true(2, . < 1), class = "hbl_error")
})

test_that("hbl_warn_identifiable()", {
  x_alpha <- diag(4)
  response <- c(1, 2, 3, 4)
  expect_silent(
    hbl_warn_identifiable(
      response = response,
      x_alpha = x_alpha,
      x_delta = NULL,
      x_beta = NULL
    )
  )
  expect_warning(
    hbl_warn_identifiable(
      response = response,
      x_alpha = cbind(x_alpha, x_alpha),
      x_delta = NULL,
      x_beta = NULL
    ),
    class = "hbl_warn"
  )
  response[3] <- NA_real_
  expect_warning(
    hbl_warn_identifiable(
      response = response,
      x_alpha = x_alpha,
      x_delta = NULL,
      x_beta = NULL
    ),
    class = "hbl_warn"
  )
})
