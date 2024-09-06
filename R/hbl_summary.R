#' @title Model summary
#' @export
#' @family summary
#' @description Summarize a fitted model in a table.
#' @details The `hb_summary()` function post-processes the results from
#'   the model. It estimates marginal means of the response,
#'   treatment effect, and other quantities of interest.
#' @return A tidy data frame with one row per group (e.g. treatment arm)
#'   and the columns in the following list. Unless otherwise specified,
#'   the quantities are calculated at the group-by-rep level.
#'   Some are calculated for the current (non-historical) study only,
#'   while others pertain to the combined dataset which includes
#'   all historical studies.
#'   * `group`: group index.
#'   * `group_label`: original group label in the data.
#'   * `rep`: rep index.
#'   * `rep_label`: original rep label in the data.
#'   * `data_mean`: observed mean of the response specific to the current
#'     study.
#'   * `data_sd`: observed standard deviation of the response
#'     specific to the current study.
#'   * `data_lower`: lower bound of a simple frequentist 95% confidence
#'     interval of the observed data mean specific to the current study.
#'   * `data_upper`: upper bound of a simple frequentist 95% confidence
#'     interval of the observed data mean specific to the current study.
#'   * `data_n`: number of non-missing observations in the combined dataset
#'     (all studies).
#'   * `data_N`: total number of observations (missing and non-missing)
#'     in the combined dataset (all studies).
#'   * `data_n_study_*`: number of non-missing observations in each study.
#'     The suffixes of these column names are integer study indexes.
#'     Call `dplyr::distinct(hbl_data(your_data), study, study_label)`
#'     to see which study labels correspond to these integer indexes.
#'   * `data_N_study_*`: total number of observations
#'     (missing and non-missing) within each study.
#'     The suffixes of these column names are integer study indexes.
#'     Call `dplyr::distinct(hbl_data(your_data), study, study_label)`
#'     to see which study labels correspond to these integer indexes.
#'   * `response_mean`: Estimated posterior mean of the response
#'     from the model. (Here, the response variable in the data
#'     should be a change from baseline outcome.)
#'     Specific to the current study.
#'   * `response_sd`: Estimated posterior standard deviation of the mean
#'     response from the model.
#'     Specific to the current study.
#'   * `response_variance`: Estimated posterior variance of the mean
#'     response from the model.
#'      Specific to the current study.
#'   * `response_lower`: Lower bound of a 95% posterior interval on the mean
#'     response from the model.
#'      Specific to the current study.
#'   * `response_upper`: Upper bound of a 95% posterior interval on the mean
#'     response from the model.
#'      Specific to the current study.
#'   * `response_mean_mcse`: Monte Carlo standard error of `response_mean`.
#'   * `response_sd_mcse`: Monte Carlo standard error of `response_sd`.
#'   * `response_lower_mcse`: Monte Carlo standard error of `response_lower`.
#'   * `response_upper_mcse`: Monte Carlo standard error of `response_upper`.
#'   * `change_*`: same as the `response_*` columns, but for change
#'     from baseline instead of the response. Not included if `response_type`
#'     is `"change"` because in that case the response is already
#'     change from baseline.
#'   * `change_percent_*`: same as the `change_*` columns, but for the
#'     *percent* change from baseline (from 0% to 100%).
#'     Not included if `response_type`
#'     is `"change"` because in that case the response is already
#'     change from baseline.
#'      Specific to the current study.
#'   * `diff_*`: same as the `response_*` columns, but for treatment effect.
#'   * `P(diff > EOI)`, `P(diff < EOI)`: CSF probabilities on the
#'     treatment effect specified with the `eoi` and `direction`
#'     arguments.
#'      Specific to the current study.
#'   * `effect_mean`: same as the `response_*` columns, but for the effect size
#'     (diff / residual standard deviation).
#'      Specific to the current study.
#'   * `precision_ratio*`: same as the `response_*` columns,
#'     but for the precision ratio, which compares within-study variance
#'     to among-study variance. Only returned for the hierarchical model.
#'     Specific to the current study.
#' @inheritParams hbl_mcmc_pool
#' @param mcmc A wide data frame of posterior samples returned by
#'   [hbl_mcmc_hierarchical()] or similar MCMC function.
#' @param response_type Character of length 1: `"raw"` if the response
#'   column in the data is the raw response, `"change"` if the response
#'   columns is change from baseline. In the latter case, the `change_*`
#'   columns in the output table are omitted because the response
#'   is already a change from baseline. Must be one of `"raw"` or `"change"`.
#' @param eoi Numeric of length at least 1,
#'   vector of effects of interest (EOIs) for critical success factors (CSFs).
#' @param direction Character of length `length(eoi)` indicating how
#'   to compare the treatment effect to each EOI. `">"` means
#'   Prob(treatment effect > EOI), and `"<"` means
#'   Prob(treatment effect < EOI). All elements of `direction`
#'   must be either `">"` or `"<"`.
#' @examples
#' if (!identical(Sys.getenv("HBL_TEST", unset = ""), "")) {
#' set.seed(0)
#' data <- hbl_sim_pool(
#'   n_study = 2,
#'   n_group = 2,
#'   n_patient = 5,
#'   n_rep = 3
#' )$data
#' tmp <- utils::capture.output(
#'   suppressWarnings(
#'     mcmc <- hbl_mcmc_hierarchical(
#'       data,
#'       chains = 1,
#'       warmup = 10,
#'       iter = 20,
#'       seed = 0
#'     )
#'   )
#' )
#' hbl_summary(mcmc, data)
#' }
hbl_summary <- function(
  mcmc,
  data,
  response = "response",
  response_type = "raw",
  study = "study",
  study_reference = max(data[[study]]),
  group = "group",
  group_reference = min(data[[group]]),
  patient = "patient",
  rep = "rep",
  rep_reference = min(data[[rep]]),
  covariates = grep("^covariate", colnames(data), value = TRUE),
  constraint = FALSE,
  eoi = 0,
  direction = "<"
) {
  true(constraint, is.logical(.), length(.) == 1, !anyNA(.))
  true(mcmc, is.data.frame(.), !is.null(colnames(.)))
  true(eoi, is.numeric(.), is.finite(.))
  true(all(direction %in% c(">", "<")))
  true(length(eoi) == length(direction))
  true(length(eoi) > 0)
  true(response_type, is.character(.), length(.) == 1, !anyNA(.))
  true(response_type %in% c("raw", "change"))
  data <- hbl_data(
    data = data,
    response = response,
    study = study,
    study_reference = study_reference,
    group = group,
    group_reference = group_reference,
    patient = patient,
    rep = rep,
    rep_reference = rep_reference,
    covariates = covariates
  )
  n_rep <- length(unique(data$rep))
  x_alpha <- if_any(
    sum(grepl("^alpha", colnames(mcmc))) == n_rep,
    get_x_alpha_pool(data, constraint = constraint),
    get_x_alpha(data, constraint = constraint)
  )
  x_delta <- get_x_delta(data, constraint = constraint)
  x_beta <- get_x_beta(data = data, x_alpha = x_alpha, x_delta = x_delta)
  samples_response <- get_samples_response(
    mcmc = mcmc,
    data = data,
    x_alpha = x_alpha,
    x_delta = x_delta
  )
  samples_change <- if_any(
    identical(response_type, "raw"),
    get_samples_change(samples_response),
    samples_response
  )
  samples_diff <- get_samples_diff(samples_change)
  samples_sigma <- get_samples_sigma(mcmc)
  samples_effect <- get_samples_effect(samples_diff, samples_sigma)
  table_data <- get_table_data(data)
  table_data_n_study <- get_table_data_n_study(data)
  table_data_N_study <- get_table_data_N_study(data)
  table_data_current <- get_table_data_current(data)
  table_response <- get_table_response(samples_response)
  if (identical(response_type, "raw")) {
    table_change <- get_table_change(samples_change)
  }
  table_diff <- get_table_diff(samples_diff)
  table_eoi <- get_table_eoi(samples_diff, eoi, direction)
  table_effect <- get_table_effect(samples_effect)
  out <- tidyr::expand_grid(
    group = sort(unique(samples_response$group)),
    rep = sort(unique(samples_response$rep))
  )
  by <- c("group", "rep")
  out <- dplyr::left_join(out, y = table_data, by = by)
  out <- dplyr::left_join(out, y = table_data_n_study, by = by)
  out <- dplyr::left_join(out, y = table_data_N_study, by = by)
  out <- dplyr::left_join(out, y = table_data_current, by = by)
  out <- dplyr::left_join(out, y = table_response, by = by)
  if (identical(response_type, "raw")) {
    out <- dplyr::left_join(out, y = table_change, by = by)
  }
  out <- dplyr::left_join(out, y = table_diff, by = by)
  out <- dplyr::left_join(out, y = table_eoi, by = by)
  out <- dplyr::left_join(out, y = table_effect, by = by)
  if (any(grepl("^tau", colnames(mcmc)))) {
    samples_tau <- get_samples_tau(mcmc)
    samples_precision_ratio <- get_samples_precision_ratio(
      samples_sigma = samples_sigma,
      samples_tau = samples_tau,
      data = data
    )
    table_precision_ratio <- get_table_precision_ratio(samples_precision_ratio)
    out <- dplyr::left_join(out, y = table_precision_ratio, by = by)
  }
  groups <- dplyr::distinct(data, group, group_label, rep, rep_label)
  out <- dplyr::left_join(x = out, y = groups, by = by)
  dplyr::select(
    out,
    group,
    group_label,
    rep,
    rep_label,
    tidyselect::everything()
  )
}

get_samples_response <- function(mcmc, data, x_alpha, x_delta) {
  index_max <- data$study == max(data$study)
  data <- data[index_max,, drop = FALSE] # nolint
  x_alpha <- x_alpha[index_max,, drop = FALSE] # nolint
  x_delta <- x_delta[index_max,, drop = FALSE] # nolint
  alpha <- t(as.matrix(mcmc[, grepl("^alpha", colnames(mcmc)), drop = FALSE]))
  delta <- t(as.matrix(mcmc[, grepl("^delta", colnames(mcmc)), drop = FALSE]))
  beta <- t(as.matrix(mcmc[, grepl("^beta", colnames(mcmc)), drop = FALSE]))
  gc()
  x_alpha <- Matrix::Matrix(x_alpha, sparse = TRUE)
  x_delta <- Matrix::Matrix(x_delta, sparse = TRUE)
  gc()
  fitted <- x_alpha %*% alpha
  rm(x_alpha)
  rm(alpha)
  gc()
  fitted <- fitted + x_delta %*% delta
  rm(x_delta)
  rm(delta)
  gc()
  groups <- tibble::tibble(
    study = data$study,
    group = data$group,
    rep = data$rep
  )
  groups$index <- paste(groups$study, groups$group, groups$rep)
  groups$index <- ordered(groups$index, levels = unique(groups$index))
  unique_groups <- groups[!duplicated(groups$index), ]
  out <- apply(fitted, 2, function(sample) tapply(sample, groups$index, mean))
  rm(fitted)
  gc()
  colnames(out) <- paste0("sample", seq_len(ncol(out)))
  out <- tibble::as_tibble(out)
  out$study <- unique_groups$study
  out$group <- unique_groups$group
  out$rep <- unique_groups$rep
  tidyr::pivot_longer(
    data = out,
    cols = tidyselect::starts_with("sample"),
    names_to = "sample"
  )
}

get_samples_change <- function(samples_response) {
  baseline <- dplyr::filter(samples_response, rep == min(rep))
  baseline <- dplyr::rename(baseline, value_baseline = value)
  baseline$rep <- NULL
  post <- dplyr::filter(samples_response, rep > min(rep))
  out <- dplyr::left_join(
    x = post,
    y = baseline,
    by = c("study", "group", "sample")
  )
  out$value_percent <- 100 * (out$value - out$value_baseline) /
    out$value_baseline
  out$value <- out$value - out$value_baseline
  out$value_baseline <- NULL
  out
}

get_samples_diff <- function(samples_change) {
  control <- dplyr::filter(samples_change, group == min(group))
  control <- dplyr::rename(control, value_control = value)
  control$group <- NULL
  treatment <- dplyr::filter(samples_change, group > min(group))
  out <- dplyr::left_join(
    x = treatment,
    y = control,
    by = c("study", "rep", "sample")
  )
  out$value <- out$value - out$value_control
  out$value_control <- NULL
  out
}

get_samples_sigma <- function(mcmc) {
  out <- dplyr::select(mcmc, tidyselect::starts_with("sigma"))
  out$sample <- paste0("sample", seq_len(nrow(out)))
  out <- tidyr::pivot_longer(out, tidyselect::starts_with("sigma"))
  out$study <- as.integer(gsub("sigma\\[|,.*$", "", out$name))
  out$rep <- as.integer(gsub(".*,|\\]", "", out$name))
  out <- dplyr::filter(out, study == max(study))
  out$name <- NULL
  out
}

get_samples_tau <- function(mcmc) {
  out <- dplyr::select(mcmc, tidyselect::starts_with("tau"))
  out$sample <- paste0("sample", seq_len(nrow(out)))
  out <- tidyr::pivot_longer(out, tidyselect::starts_with("tau"))
  out$rep <- as.integer(gsub("tau\\[|\\]", "", out$name))
  out$name <- NULL
  out
}

get_samples_effect <- function(samples_diff, samples_sigma) {
  samples_diff$value_diff <- samples_diff$value
  samples_diff$value <- NULL
  samples_sigma$value_sigma <- samples_sigma$value
  samples_sigma$value <- NULL
  out <- dplyr::left_join(
    x = samples_diff,
    y = samples_sigma,
    by = c("sample", "rep")
  )
  out$value <- out$value_diff / out$value_sigma
  out$value_diff <- NULL
  out$value_sigma <- NULL
  out
}

get_samples_precision_ratio <- function(samples_sigma, samples_tau, data) {
  n <- dplyr::count(
    filter(data, study == max(study), !is.na(response)),
    rep
  )
  samples_sigma <- dplyr::rename(samples_sigma, sigma = value)
  samples_tau <- dplyr::rename(samples_tau, tau = value)
  out <- dplyr::left_join(
    x = samples_sigma,
    y = samples_tau,
    by = c("rep", "sample")
  )
  out <- dplyr::left_join(x = out, y = n, by = "rep")
  out$value <- (1 / out$tau ^ 2) /
    ((1 / out$tau ^ 2) + 1 / (out$sigma ^ 2 / out$n))
  out
}

get_table_data <- function(data) {
  dplyr::summarize(
    dplyr::group_by(data, group, rep),
    data_n = sum(!is.na(response)),
    data_N = dplyr::n(),
    .groups = "drop"
  )
}

get_table_data_current <- function(data) {
  q <- stats::qnorm(p = 0.975)
  out <- dplyr::summarize(
    dplyr::group_by(filter(data, study == max(study)), group, rep),
    data_mean = mean(response, na.rm = TRUE),
    data_sd = sd(response, na.rm = TRUE),
    n = sum(!is.na(response)),
    data_lower = data_mean - q * data_sd / sqrt(n),
    data_upper = data_mean + q * data_sd / sqrt(n),
    .groups = "drop"
  )
  out$n <- NULL
  out
}

get_table_data_n_study <- function(data) {
  out <- dplyr::summarize(
    dplyr::group_by(data, study, group, rep),
    n = sum(!is.na(response)),
    .groups = "drop"
  )
  tidyr::pivot_wider(
    out,
    names_from = "study",
    values_from = "n",
    names_prefix = "data_n_study_",
    values_fill = 0
  )
}

get_table_data_N_study <- function(data) {
  out <- dplyr::summarize(
    dplyr::group_by(data, study, group, rep),
    n = dplyr::n(),
    .groups = "drop"
  )
  tidyr::pivot_wider(
    out,
    names_from = "study",
    values_from = "n",
    names_prefix = "data_N_study_",
    values_fill = 0
  )
}

get_table_response <- function(samples_response) {
  dplyr::summarize(
    dplyr::group_by(samples_response, group, rep),
    response_mean = mean(value),
    response_variance = stats::var(value),
    response_sd = stats::sd(value),
    response_lower = quantile(value, 0.025),
    response_upper = quantile(value, 0.975),
    response_mean_mcse = posterior::mcse_mean(value),
    response_sd_mcse = posterior::mcse_sd(value),
    response_lower_mcse = posterior::mcse_quantile(value, 0.025),
    response_upper_mcse = posterior::mcse_quantile(value, 0.975),
    .groups = "drop"
  )
}

get_table_change <- function(samples_change) {
  dplyr::summarize(
    dplyr::group_by(samples_change, group, rep),
    change_mean = mean(value),
    change_lower = quantile(value, 0.025),
    change_upper = quantile(value, 0.975),
    change_mean_mcse = posterior::mcse_mean(value),
    change_lower_mcse = posterior::mcse_quantile(value, 0.025),
    change_upper_mcse = posterior::mcse_quantile(value, 0.975),
    change_percent_mean = mean(value_percent),
    change_percent_lower = quantile(value_percent, 0.025),
    change_percent_upper = quantile(value_percent, 0.975),
    change_percent_mean_mcse = posterior::mcse_mean(value_percent),
    change_percent_lower_mcse = posterior::mcse_quantile(value_percent, 0.025),
    change_percent_upper_mcse = posterior::mcse_quantile(value_percent, 0.975),
    .groups = "drop"
  )
}

get_table_diff <- function(samples_diff) {
  dplyr::summarize(
    dplyr::group_by(samples_diff, group, rep),
    diff_mean = mean(value),
    diff_lower = quantile(value, 0.025),
    diff_upper = quantile(value, 0.975),
    diff_mean_mcse = posterior::mcse_mean(value),
    diff_lower_mcse = posterior::mcse_quantile(value, 0.025),
    diff_upper_mcse = posterior::mcse_quantile(value, 0.975),
    .groups = "drop"
  )
}

get_table_eoi <- function(samples_diff, eoi, direction) {
  out <- list()
  for (index in seq_along(eoi)) {
    section <- dplyr::summarize(
      dplyr::group_by(samples_diff, group, rep),
      value = if_any(
        direction[index] == ">",
        mean(value > eoi[index]),
        mean(value < eoi[index])
      )
    )
    name <- sprintf("P(diff %s %s)", direction[index], eoi[index])
    section[[name]] <- section$value
    section$value <- NULL
    out[[index]] <- section
  }
  Reduce(
    f = function(x, y) dplyr::left_join(x, y, by = c("group", "rep")),
    x = out
  )
}

get_table_effect <- function(samples_effect) {
  dplyr::summarize(
    dplyr::group_by(samples_effect, group, rep),
    effect_mean = mean(value),
    effect_lower = quantile(value, 0.025),
    effect_upper = quantile(value, 0.975),
    effect_mean_mcse = posterior::mcse_mean(value),
    effect_lower_mcse = posterior::mcse_quantile(value, 0.025),
    effect_upper_mcse = posterior::mcse_quantile(value, 0.975),
    .groups = "drop"
  )
}

get_table_precision_ratio <- function(samples_precision_ratio) {
  dplyr::summarize(
    dplyr::group_by(samples_precision_ratio, rep),
    group = 1,
    precision_ratio = mean(value),
    precision_ratio_lower = quantile(value, 0.025),
    precision_ratio_upper = quantile(value, 0.975)
  )
}
