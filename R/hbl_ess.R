#' @title Effective sample size (ESS)
#' @export
#' @family summary
#' @description Quantify borrowing with effective sample size (ESS)
#'   as cited and explained in the methods vignette at
#'   <https://wlandau.github.io/historicalborrowlong/articles/methods.html>.
#' @return A data frame with one row per discrete time point ("rep")
#'   and the following columns:
#'   * `v0`: posterior predictive variance of the control group mean of a
#'     hypothetical new study given the pooled model.
#'     Calculated as the mean over MCMC samples of `1 / sum(sigma_i ^ 2)`,
#'     where each `sigma_i` is the residual standard deviation of
#'     study `i` estimated from the pooled model.
#'   * `v_tau`: posterior predictive variance of a hypothetical
#'     new control group mean under the hierarchical model.
#'     Calculated by averaging over predictive draws,
#'     where each predictive draw is from
#'     `rnorm(n = 1, mean = mu_, sd = tau_)` and `mu_` and `tau_` are the
#'     `mu` and `tau` components of an MCMC sample.
#'   * `n`: number of non-missing historical control patients.
#'   * `weight`: strength of borrowing as a ratio of variances: `v0 / v_tau`.
#'   * `ess`: strength of borrowing as a prior effective sample size:
#'      `n v0 / v_tau`, where `n` is the number of non-missing historical
#'      control patients.
#' @inheritParams hbl_data
#' @param mcmc_pool A fitted model from [hbl_mcmc_pool()].
#' @param mcmc_hierarchical A fitted model from [hbl_mcmc_hierarchical()].
#' @examples
#'   set.seed(0)
#'   data <- hbl_sim_independent(n_continuous = 2)$data
#'   data$group <- sprintf("group%s", data$group)
#'   data$study <- sprintf("study%s", data$study)
#'   data$rep <- sprintf("rep%s", data$rep)
#'   tmp <- utils::capture.output(
#'     suppressWarnings(
#'       pool <- hbl_mcmc_pool(
#'         data,
#'         chains = 1,
#'         warmup = 10,
#'         iter = 20,
#'         seed = 0
#'       )
#'     )
#'   )
#'   tmp <- utils::capture.output(
#'     suppressWarnings(
#'       hierarchical <- hbl_mcmc_hierarchical(
#'         data,
#'         chains = 1,
#'         warmup = 10,
#'         iter = 20,
#'         seed = 0
#'       )
#'     )
#'   )
#'   hbl_ess(
#'     mcmc_pool = pool,
#'     mcmc_hierarchical = hierarchical,
#'     data = data
#'   )
hbl_ess <- function(
  mcmc_pool,
  mcmc_hierarchical,
  data,
  response = "response",
  study = "study",
  study_reference = max(data[[study]]),
  group = "group",
  group_reference = min(data[[group]]),
  patient = "patient",
  rep = "rep",
  rep_reference = min(data[[rep]])
) {
  true(
    mcmc_hierarchical,
    is.data.frame(.),
    !is.null(colnames(.)),
    all(c("mu[1]", "tau[1]") %in% colnames(.))
  )
  true(mcmc_pool, is.data.frame(.), !is.null(colnames(.)))
  true(data, is.data.frame(.), !is.null(colnames(.)))
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
    covariates = character(0L)
  )
  v0 <- hbl_ess_v0(data, mcmc_pool)
  v_tau <- hbl_ess_v_tau(data, mcmc_hierarchical)
  weight <- v0 / v_tau
  historical_controls <- dplyr::filter(
    data,
    data$group == min(data$group),
    data$study < max(data$study)
  )
  grouped <- dplyr::group_by(historical_controls, rep, rep_label)
  out <- dplyr::summarize(grouped, n = sum(!is.na(response)))
  out$v0 <- v0
  out$v_tau <- v_tau
  out$weight <- weight
  out$ess <- out$n * weight
  out
}

hbl_ess_v0 <- function(data, mcmc_pool) {
  sigma <- dplyr::select(mcmc_pool, tidyselect::starts_with("sigma["))
  precision <- as.data.frame(lapply(sigma, function(x) x ^ (-2)))
  out <- numeric(0L)
  studies <- sort(unique(data$study))
  for (rep in sort(unique(data$rep))) {
    names <- sprintf("sigma[%s,%s]", studies, rep)
    true(all(names %in% colnames(sigma)))
    precision_rep <- precision[, colnames(sigma) %in% names]
    out[rep] <- mean(1 / rowSums(precision_rep))
  }
  out
}

hbl_ess_v_tau <- function(data, mcmc_hierarchical) {
  mu <- dplyr::select(mcmc_hierarchical, tidyselect::starts_with("mu["))
  tau <- dplyr::select(mcmc_hierarchical, tidyselect::starts_with("tau["))
  out <- numeric(0L)
  for (rep in sort(unique(data$rep))) {
    name_mu <- sprintf("mu[%s]", rep)
    name_tau <- sprintf("tau[%s]", rep)
    true(all(c(name_mu, name_tau) %in% colnames(mcmc_hierarchical)))
    out[rep] <- stats::var(
      stats::rnorm(
        n = nrow(mcmc_hierarchical),
        mean = mcmc_hierarchical[[name_mu]],
        sd = mcmc_hierarchical[[name_tau]]
      )
    )
  }
  out
}
