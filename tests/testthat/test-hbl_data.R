test_that("as_index_min() on characters", {
  x <- c("d", "d", "z", "a", "b", "b")
  out <- as_index_min(x, min = "b")
  expect_true(is.integer(out))
  expect_equal(out, c(3L, 3L, 4L, 2L, 1L, 1L))
})

test_that("as_index_min() on numerics", {
  x <- c(3.7, 3.7, 9999, 1.2, 2.2, 2.2)
  out <- as_index_min(x, min = 2.2)
  expect_true(is.integer(out))
  expect_equal(out, c(3L, 3L, 4L, 2L, 1L, 1L))
})

test_that("as_index_max() on characters", {
  x <- c("d", "d", "z", "a", "b", "b")
  out <- as_index_max(x, max = "b")
  expect_true(is.integer(out))
  expect_equal(out, c(2L, 2L, 3L, 1L, 4L, 4L))
})

test_that("as_index_max() on numerics", {
  x <- c(3.7, 3.7, 9999, 1.2, 2.2, 2.2)
  out <- as_index_max(x, max = 2.2)
  expect_true(is.integer(out))
  expect_equal(out, c(2L, 2L, 3L, 1L, 4L, 4L))
})

test_that("hbl_data() repeat patient ID", {
  data <- tibble::tibble(
    trial = rep(c("study_current", "study_historical"), each = 4),
    arm = rep(c("zzz", "treatment", rep("zzz", 2)), each = 2),
    subject = paste0("patient_", seq_len(8)),
    outcome = c(rnorm(n = 7), NA_real_),
    block1 = seq_len(8),
    block2 = seq_len(8),
    block3 = seq_len(8)
  )
  data <- tidyr::expand_grid(data, visit = c("visit1", "visit2"))
  data$subject[13] <- data$subject[15]
  data$subject[14] <- data$subject[15]
  expect_error(
    hbl_data(
      data = data,
      response = "outcome",
      study = "trial",
      study_reference = "study_current",
      group = "arm",
      group_reference = "zzz",
      patient = "subject",
      rep = "visit",
      rep_reference = "visit1",
      covariates = c("block1", "block3")
    ),
    class = "hbl_error"
  )
})

test_that("hbl_data()", {
  data <- tibble::tibble(
    trial = rep(c("study_current", "study_historical"), each = 4),
    arm = rep(c("zzz", "treatment", rep("zzz", 2)), each = 2),
    subject = paste0("patient_", seq_len(8)),
    outcome = c(rnorm(n = 7), NA_real_),
    block1 = seq_len(8),
    block2 = seq_len(8),
    block3 = seq_len(8)
  )
  data <- tidyr::expand_grid(data, visit = c("visit1", "visit2"))
  data$id <- paste(data$subject, data$visit)
  out <- hbl_data(
    data = data,
    response = "outcome",
    study = "trial",
    study_reference = "study_current",
    group = "arm",
    group_reference = "zzz",
    patient = "subject",
    rep = "visit",
    rep_reference = "visit1",
    covariates = c("block1", "block3")
  )
  expect_equal(dim(out), c(16, 11))
  expect_equal(
    sort(colnames(out)),
    sort(
      c(
        "response",
        "study",
        "study_label",
        "group",
        "group_label",
        "patient",
        "patient_label",
        "rep",
        "rep_label",
        "covariate_block1",
        "covariate_block3"
      )
    )
  )
  out$id <- paste(out$patient_label, out$rep_label)
  exp <- data[match(data$id, out$id), ]
  expect_equal(exp$outcome, out$response)
  expect_equal(out$covariate_block1, exp$block1)
  expect_equal(out$covariate_block3, exp$block3)
  expect_equal(exp$arm, c(rep("zzz", 12), rep("treatment", 4)))
  expect_equal(out$group, c(rep(1, 12), rep(2, 4)))
  expect_equal(out$patient_label, exp$subject)
  expect_equal(out$study, rep(c(1, 2), each = 8))
  expect_equal(out$rep, rep(seq_len(2), times = 8))
  expect_equal(out$rep_label, rep(c("visit1", "visit2"), times = 8))
  expect_equal(exp$trial, rep(c("study_historical", "study_current"), each = 8))
})

test_that("hbl_data() completes the grid", {
  set.seed(0)
  data <- tibble::tibble(
    trial = rep(c("study_current", "study_historical"), each = 4),
    arm = rep(c("zzz", "treatment", rep("zzz", 2)), each = 2),
    subject = paste0("patient_", seq_len(8)),
    block1 = seq_len(8),
    block2 = seq_len(8),
    block3 = seq_len(8)
  )
  data <- tidyr::expand_grid(data, visit = c("visit1", "visit2"))
  data$outcome <- rnorm(n = 16)
  data_full <- data
  data_full$id <- paste(data_full$subject, data_full$visit)
  data <- data[-15,, drop = FALSE] # nolint
  out <- hbl_data(
    data = data,
    response = "outcome",
    study = "trial",
    study_reference = "study_current",
    group = "arm",
    group_reference = "zzz",
    patient = "subject",
    rep = "visit",
    rep_reference = "visit1",
    covariates = c("block1", "block3")
  )
  expect_equal(dim(out), c(16, 11))
  expect_equal(
    sort(colnames(out)),
    sort(
      c(
        "response",
        "study",
        "study_label",
        "group",
        "group_label",
        "patient",
        "patient_label",
        "rep",
        "rep_label",
        "covariate_block1",
        "covariate_block3"
      )
    )
  )
  out$id <- paste(out$patient_label, out$rep_label)
  exp <- data_full[match(data_full$id, out$id), ]
  exp$outcome[exp$id == "patient_8 visit1"] <- NA_real_
  expect_equal(exp$outcome, out$response)
  expect_equal(out$covariate_block1, exp$block1)
  expect_equal(out$covariate_block3, exp$block3)
  expect_equal(exp$arm, c(rep("zzz", 12), rep("treatment", 4)))
  expect_equal(out$group, c(rep(1, 12), rep(2, 4)))
  expect_equal(out$patient_label, exp$subject)
  expect_equal(out$study, rep(c(1, 2), each = 8))
  expect_equal(out$rep, rep(seq_len(2), times = 8))
  expect_equal(out$rep_label, rep(c("visit1", "visit2"), times = 8))
  expect_equal(exp$trial, rep(c("study_historical", "study_current"), each = 8))
})

test_that("hbl_data() on different reference levels", {
  data <- tibble::tibble(
    trial = rep(c("study_current", "study_historical"), each = 4),
    arm = rep(c("zzz", "treatment", rep("zzz", 2)), each = 2),
    subject = paste0("patient_", seq_len(8)),
    outcome = c(rnorm(n = 7), NA_real_),
    block1 = seq_len(8),
    block2 = seq_len(8),
    block3 = seq_len(8)
  )
  data <- tidyr::expand_grid(data, visit = c("visit1", "visit2"))
  data$id <- paste(data$subject, data$visit)
  out <- hbl_data(
    data = data,
    response = "outcome",
    study = "trial",
    study_reference = "study_historical",
    group = "arm",
    group_reference = "treatment",
    patient = "subject",
    rep = "visit",
    rep_reference = "visit2",
    covariates = "block2"
  )
  expect_equal(dim(out), c(16, 10))
  expect_equal(
    sort(colnames(out)),
    sort(
      c(
        "response",
        "study",
        "study_label",
        "group",
        "group_label",
        "patient",
        "patient_label",
        "rep",
        "rep_label",
        "covariate_block2"
      )
    )
  )
  out$id <- paste(out$patient_label, out$rep_label)
  exp <- data[match(data$id, out$id), ]
  expect_equal(exp$outcome, out$response)
  expect_equal(out$covariate_block2, exp$block2)
  expect_equal(exp$arm, c(rep("treatment", 4), rep("zzz", 12)))
  expect_equal(out$group, c(rep(1, 4), rep(2, 12)))
  expect_equal(out$patient_label, exp$subject)
  expect_equal(out$study, rep(c(1, 2), each = 8))
  expect_equal(out$rep, rep(c(1, 2), times = 8))
  expect_equal(out$rep_label, rep(c("visit2", "visit1"), times = 8))
  expect_equal(exp$trial, rep(c("study_current", "study_historical"), each = 8))
})

test_that("hbl_data_enforce_baseline_covariates()", {
  data <- tibble::tibble(
    trial = rep(c("study_current", "study_historical"), each = 4),
    arm = rep(c("zzz", "treatment", rep("zzz", 2)), each = 2),
    subject = paste0("patient_", seq_len(8)),
    outcome = c(rnorm(n = 7), NA_real_),
    block1 = seq_len(8),
    block2 = seq_len(8),
    block3 = seq_len(8)
  )
  data <- tidyr::expand_grid(data, visit = c("visit1", "visit2"))
  data$id <- paste(data$subject, data$visit)
  out <- hbl_data(
    data = data,
    response = "outcome",
    study = "trial",
    study_reference = "study_current",
    group = "arm",
    group_reference = "zzz",
    patient = "subject",
    rep = "visit",
    rep_reference = "visit1",
    covariates = c("block1", "block3")
  )
  expect_silent(hbl_data_enforce_baseline_covariates(out))
  expect_equal(out$covariate_block1[1], 5L)
  expect_equal(out$covariate_block1[2], 5L)
  out$covariate_block1[2] <- -9999L
  expect_equal(out$covariate_block1[1], 5L)
  expect_equal(out$covariate_block1[2], -9999L)
  expect_warning(
    out2 <- hbl_data_enforce_baseline_covariates(out),
    class = "hbl_warn"
  )
  expect_equal(out$covariate_block1[1], 5L)
  expect_equal(out$covariate_block1[2], -9999L)
  expect_equal(out2$covariate_block1[1], 5L)
  expect_equal(out2$covariate_block1[2], 5L)
})
