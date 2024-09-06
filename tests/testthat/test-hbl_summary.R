test_that("hbl_summary() pool with raw response type", {
  skip_on_cran()
  set.seed(0)
  data <- hbl_sim_pool(
    n_study = 2,
    n_group = 2,
    n_patient = 5,
    n_rep = 3,
    constraint = FALSE
  )$data
  data$group <- sprintf("group%s", data$group)
  data$rep <- sprintf("rep%s", data$rep)
  tmp <- utils::capture.output(
    suppressWarnings(
      mcmc <- hbl_mcmc_pool(
        data,
        chains = 1,
        warmup = 10,
        iter = 20,
        seed = 0,
        constraint = FALSE
      )
    )
  )
  out <- hbl_summary(
    mcmc,
    data,
    eoi = c(0, 1),
    direction = c(">", "<"),
    constraint = FALSE
  )
  expect_equal(out$group_label, rep(c("group1", "group2"), each = 3))
  expect_equal(out$group, rep(c(1, 2), each = 3))
  expect_equal(out$rep_label, rep(c("rep1", "rep2", "rep3"), times = 2))
  expect_equal(out$rep, rep(seq_len(3), times = 2))
  expect_equal(dim(out), c(6, 49))
  cols <- c(
    "group", "group_label",
    "rep", "rep_label",
    "data_mean", "data_lower", "data_upper",
    "data_n", "data_N", "data_sd",
    "data_n_study_1", "data_n_study_2",
    "data_N_study_1", "data_N_study_2",
    "response_mean", "response_sd", "response_variance", "response_lower",
    "response_upper", "response_mean_mcse", "response_sd_mcse",
    "response_lower_mcse",  "response_upper_mcse",
    "change_mean", "change_lower",
    "change_upper",
    "change_percent_mean", "change_percent_lower",
    "change_percent_upper",
    "change_percent_mean_mcse", "change_percent_lower_mcse",
    "change_percent_upper_mcse",
    "change_mean_mcse",
    "change_lower_mcse", "change_upper_mcse",
    "diff_mean", "diff_lower",
    "diff_upper",  "diff_mean_mcse", "diff_lower_mcse", "diff_upper_mcse",
    "P(diff > 0)",  "P(diff < 1)", "effect_mean", "effect_lower",
    "effect_upper",  "effect_mean_mcse", "effect_lower_mcse",
    "effect_upper_mcse"
  )
  expect_equal(sort(cols), sort(colnames(out)))
})

test_that("hbl_summary() independent with change response type", {
  set.seed(0)
  data <- hbl_sim_independent(
    n_study = 2,
    n_group = 2,
    n_patient = 5,
    n_rep = 3,
    constraint = FALSE
  )$data
  data$group <- sprintf("group%s", data$group)
  data$rep <- sprintf("rep%s", data$rep)
  tmp <- utils::capture.output(
    suppressWarnings(
      mcmc <- hbl_mcmc_independent(
        data,
        chains = 1,
        warmup = 10,
        iter = 20,
        seed = 0,
        constraint = FALSE
      )
    )
  )
  out <- hbl_summary(
    mcmc,
    data,
    eoi = c(0, 1),
    direction = c(">", "<"),
    constraint = FALSE,
    response_type = "change"
  )
  expect_equal(out$group_label, rep(c("group1", "group2"), each = 3))
  expect_equal(out$group, rep(c(1, 2), each = 3))
  expect_equal(out$rep_label, rep(c("rep1", "rep2", "rep3"), times = 2))
  expect_equal(out$rep, rep(seq_len(3), times = 2))
  expect_equal(dim(out), c(6, 37))
  cols <- c(
    "group", "group_label",
    "rep", "rep_label",
    "data_mean", "data_lower", "data_upper",
    "data_n", "data_N", "data_sd",
    "data_n_study_1", "data_n_study_2",
    "data_N_study_1", "data_N_study_2",
    "response_mean", "response_sd", "response_variance", "response_lower",
    "response_upper", "response_mean_mcse", "response_sd_mcse",
    "response_lower_mcse",  "response_upper_mcse",
    "diff_mean", "diff_lower",
    "diff_upper",  "diff_mean_mcse", "diff_lower_mcse", "diff_upper_mcse",
    "P(diff > 0)",  "P(diff < 1)", "effect_mean", "effect_lower",
    "effect_upper",  "effect_mean_mcse", "effect_lower_mcse",
    "effect_upper_mcse"
  )
  expect_equal(sort(cols), sort(colnames(out)))
})

test_that("hbl_summary() data counts", {
  set.seed(0)
  data <- hbl_sim_independent(
    n_study = 2,
    n_group = 2,
    n_patient = 5,
    n_rep = 3,
    constraint = FALSE
  )$data
  data$group <- sprintf("group%s", data$group)
  data$rep <- sprintf("rep%s", data$rep)
  tmp <- utils::capture.output(
    suppressWarnings(
      mcmc <- hbl_mcmc_independent(
        data,
        chains = 1,
        warmup = 10,
        iter = 20,
        seed = 0,
        constraint = FALSE
      )
    )
  )
  data$response <- NA
  data$response[1] <- 1
  out <- hbl_summary(
    mcmc,
    data,
    eoi = c(0, 1),
    direction = c(">", "<"),
    constraint = FALSE,
    response_type = "change"
  )
  expect_equal(out$data_n, c(1, rep(0, 5)))
  expect_equal(out$data_N, rep(c(10, 5), each = 3))
  expect_equal(out$data_n_study_1, c(1, rep(0, 5)))
  expect_equal(out$data_n_study_2, rep(0, 6))
  expect_equal(out$data_N_study_1, rep(c(5, 0), each = 3))
  expect_equal(out$data_N_study_2, rep(5, 6))
})

test_that("hbl_summary() hierarchical with raw response type", {
  skip_on_cran()
  set.seed(0)
  data <- hbl_sim_hierarchical(
    n_study = 2,
    n_group = 2,
    n_patient = 5,
    n_rep = 3,
    constraint = FALSE
  )$data
  data$group <- sprintf("group%s", data$group)
  data$rep <- sprintf("rep%s", data$rep)
  tmp <- utils::capture.output(
    suppressWarnings(
      mcmc <- hbl_mcmc_hierarchical(
        data,
        chains = 1,
        warmup = 10,
        iter = 20,
        seed = 0,
        constraint = FALSE
      )
    )
  )
  out <- hbl_summary(
    mcmc,
    data,
    eoi = c(0, 1),
    direction = c(">", "<"),
    constraint = FALSE
  )
  expect_equal(out$group_label, rep(c("group1", "group2"), each = 3))
  expect_equal(out$group, rep(c(1, 2), each = 3))
  expect_equal(out$rep_label, rep(c("rep1", "rep2", "rep3"), times = 2))
  expect_equal(out$rep, rep(seq_len(3), times = 2))
  expect_equal(dim(out), c(6, 52))
  cols <- c(
    "group", "group_label",
    "rep", "rep_label",
    "data_mean", "data_lower", "data_upper",
    "data_n", "data_N", "data_sd",
    "data_n_study_1", "data_n_study_2",
    "data_N_study_1", "data_N_study_2",
    "response_mean", "response_sd", "response_variance", "response_lower",
    "response_upper", "response_mean_mcse", "response_sd_mcse",
    "response_lower_mcse",  "response_upper_mcse",
    "change_mean", "change_lower",
    "change_upper",
    "change_percent_mean", "change_percent_lower",
    "change_percent_upper",
    "change_percent_mean_mcse", "change_percent_lower_mcse",
    "change_percent_upper_mcse",
    "change_mean_mcse",
    "change_lower_mcse", "change_upper_mcse",
    "diff_mean", "diff_lower",
    "diff_upper",  "diff_mean_mcse", "diff_lower_mcse", "diff_upper_mcse",
    "P(diff > 0)",  "P(diff < 1)", "effect_mean", "effect_lower",
    "effect_upper",  "effect_mean_mcse", "effect_lower_mcse",
    "effect_upper_mcse",
    "precision_ratio", "precision_ratio_lower", "precision_ratio_upper"
  )
  expect_equal(sort(cols), sort(colnames(out)))
})

test_that("hbl_summary() pool mock mcmc", {
  set.seed(0)
  data <- hbl_sim_pool(
    n_study = 2,
    n_group = 2,
    n_patient = 5,
    n_rep = 3,
    n_continuous = 1,
    constraint = FALSE
  )$data
  tmp <- utils::capture.output(
    suppressWarnings(
      mcmc <- hbl_mcmc_pool(
        data,
        chains = 2,
        warmup = 100,
        iter = 200,
        seed = 0,
        constraint = FALSE
      )
    )
  )
  mcmc <- head(mcmc, 6)
  for (col in colnames(mcmc)) {
    mcmc[[col]] <- seq_len(6) / 10
  }
  out <- hbl_summary(
    mcmc,
    data,
    eoi = c(0, 1),
    direction = c(">", "<"),
    constraint = FALSE
  )
  expect_equal(out$group, rep(seq_len(2), each = 3))
  data_current <- dplyr::filter(data, study == 2)
  x_alpha <- get_x_alpha_pool(data, constraint = FALSE)
  x_delta <- get_x_delta(data, constraint = FALSE)
  grid <- tidyr::expand_grid(
    group = seq_len(2),
    rep = seq_len(3),
    sample = seq_len(6)
  )
  grid <- dplyr::group_by(grid, group, rep, sample)
  grid <- dplyr::group_modify(grid, ~{
    name_alpha <- sprintf("alpha[%s]", .y$rep)
    name_delta <- sprintf("delta[%s]", .y$rep)
    index <- data$study == 2 & data$group == .y$group & data$rep == .y$rep
    if (.y$group == 1) {
      y <- mcmc[[name_alpha]][.y$sample]
    } else {
      y <- mcmc[[name_delta]][.y$sample]
    }
    tibble::tibble(value = y)
  })
  response <- dplyr::ungroup(grid)
  response <- dplyr::group_by(response, group, rep, sample)
  response <- dplyr::summarize(response, value = mean(value), .groups = "drop")
  response <- dplyr::arrange(response, group, rep, sample)
  for (i in seq_len(6)) {
    j <- data$group == out$group[i] & data$rep == out$rep[i]
    group <- out$group[i]
    rep <- out$rep[i]
    k <- response$group == group & response$rep == rep
    s <- j & data$study == max(data$study)
    expect_equal(out$data_mean[i], mean(data$response[s]))
    expect_equal(out$data_sd[i], sd(data$response[s]))
    q_ <- stats::qnorm(p = 0.975)
    m_ <- mean(data$response[s], na.rm = TRUE)
    s_ <- sd(data$response[s], na.rm = TRUE)
    n_ <- sum(!is.na(data$response[s]))
    expect_equal(out$data_lower[i], m_ - q_ * s_ / sqrt(n_))
    expect_equal(out$data_upper[i], m_ + q_ * s_ / sqrt(n_))
    expect_equal(
      out$response_mean[i],
      mean(response$value[k])
    )
    expect_equal(
      out$response_lower[i],
      quantile(response$value[k], 0.025)
    )
    expect_equal(
      out$response_upper[i],
      quantile(response$value[k], 0.975)
    )
    expect_equal(
      out$response_sd[i],
      sd(response$value[k])
    )
    expect_equal(
      out$response_variance[i],
      var(response$value[k])
    )
  }
  value <- response$value
  baseline <- value[c(rep(seq_len(6), times = 3), rep(seq(19, 24), times = 3))]
  change <- value - baseline
  for (i in c(2, 3, 5, 6)) {
    j <- data$group == out$group[i] & data$rep == out$rep[i]
    group <- out$group[i]
    rep <- out$rep[i]
    k <- response$group == group & response$rep == rep
    expect_equal(
      out$change_mean[i],
      mean(change[k])
    )
    expect_equal(
      out$change_lower[i],
      quantile(change[k], 0.025)
    )
    expect_equal(
      out$change_upper[i],
      quantile(change[k], 0.975)
    )
  }
  change_percent <- (value - baseline) / baseline
  for (i in c(2, 3, 5, 6)) {
    j <- data$group == out$group[i] & data$rep == out$rep[i]
    group <- out$group[i]
    rep <- out$rep[i]
    k <- response$group == group & response$rep == rep
    expect_equal(
      out$change_percent_mean[i],
      mean(change_percent[k])
    )
    expect_equal(
      out$change_percent_lower[i],
      quantile(change_percent[k], 0.025)
    )
    expect_equal(
      out$change_percent_upper[i],
      quantile(change_percent[k], 0.975)
    )
  }
  control <- change[rep(seq_len(18), times = 2)]
  diff <- change - control
  for (i in c(5, 6)) {
    j <- data$group == out$group[i] & data$rep == out$rep[i]
    group <- out$group[i]
    rep <- out$rep[i]
    k <- response$group == group & response$rep == rep
    expect_equal(
      out$diff_mean[i],
      mean(diff[k])
    )
    expect_equal(
      out$diff_lower[i],
      quantile(diff[k], 0.025)
    )
    expect_equal(
      out$diff_upper[i],
      quantile(diff[k], 0.975)
    )
  }
  expect_equal(out[["P(diff > 0)"]], c(rep(NA_real_, 4), 0, 0))
  sigma <- c(rep(Inf, 24), mcmc[["sigma[2,2]"]], mcmc[["sigma[2,3]"]])
  effect <- diff / sigma
  for (i in c(5, 6)) {
    j <- data$group == out$group[i] & data$rep == out$rep[i]
    group <- out$group[i]
    rep <- out$rep[i]
    k <- response$group == group & response$rep == rep
    expect_equal(
      out$effect_mean[i],
      mean(effect[k])
    )
    expect_equal(
      out$effect_lower[i],
      quantile(effect[k], 0.025)
    )
    expect_equal(
      out$effect_upper[i],
      quantile(effect[k], 0.975)
    )
  }
})

test_that("hbl_summary() precision ratio", {
  skip_on_cran()
  set.seed(0)
  data <- hbl_sim_hierarchical(
    n_study = 2,
    n_group = 2,
    n_patient = 5,
    n_rep = 3,
    n_continuous = 1,
    constraint = FALSE
  )$data
  tmp <- utils::capture.output(
    suppressWarnings(
      mcmc <- hbl_mcmc_hierarchical(
        data,
        chains = 1,
        warmup = 10,
        iter = 20,
        seed = 0,
        constraint = FALSE
      )
    )
  )
  mcmc <- head(mcmc, 6)
  for (col in colnames(mcmc)) {
    mcmc[[col]] <- mcmc[[col]] + seq_len(6) / 10
  }
  out <- hbl_summary(
    mcmc,
    data,
    eoi = c(0, 1),
    direction = c(">", "<"),
    constraint = FALSE
  )
  for (rep in seq_len(3)) {
    sigma <- mcmc[[sprintf("sigma[2,%s]", rep)]]
    tau <- mcmc[[sprintf("tau[%s]", rep)]]
    n <- 10
    samples <- (1 / tau ^ 2) / ((1 / tau ^ 2) + 1 / (sigma ^ 2 / n))
    expect_equal(
      out$precision_ratio[rep],
      mean(samples)
    )
    expect_equal(
      out$precision_ratio_lower[rep],
      quantile(samples, 0.025)
    )
    expect_equal(
      out$precision_ratio_upper[rep],
      quantile(samples, 0.975)
    )
  }
})
