#' @title Non-longitudinal hierarchical simulations.
#' @export
#' @family simulate
#' @description Simulate from the non-longitudinal hierarchical model.
#' @return A list with the following elements:
#'   * `data`: tidy long-form dataset with the patient-level data.
#'     one row per patient per rep and indicator columns for the study,
#'     group (e.g. treatment arm), patient ID, and rep. The `response`
#'     columns is the patient response. The other columns are
#'     baseline covariates. The control group is the one with
#'     the `group` column equal to 1, and the current study (non-historical)
#'     is the one with the maximum value of the `study` column.
#'     Only the current study has any non-control-group patients,
#'     the historical studies have only the control group.
#'   * `parameters`: named list of model parameter values.
#'     See the model specification vignette for details.
#'   * `matrices`: A named list of model matrices.
#'     See the model specification vignette for details.
#' @inheritParams hbl_sim_pool
#' @param s_mu Numeric of length 1,
#'   prior standard deviation of `mu`.
#' @param s_tau Numeric of length 1,
#'   Upper bound on `tau`.
#' @param d_tau Positive numeric of length 1. Degrees of freedom of the
#'   Student t prior of `tau` if `prior_tau` is `"half_t"`.
#' @param s_tau Non-negative numeric of length 1.
#'   If `prior_tau` is `"half_t"`, then `s_tau` is the scale parameter of
#'   the Student t prior of `tau` and analogous to the `sigma` parameter of
#'   the Student-t parameterization given at
#'   <https://mc-stan.org/docs/functions-reference/unbounded_continuous_distributions.html>. # nolint
#'   If `prior_tau` is `"uniform"`, then `s_tau` is the upper bound of `tau`.
#'   Upper bound on `tau` if `prior_tau` is `"uniform"`.
#' @param prior_tau Character string, family of the prior of `tau`.
#'   If `prior_tau` equals `"uniform"`, then the prior on `tau` is
#'   a uniform prior with lower bound 0 and upper bound `s_tau`.
#'   If `prior_tau` equals `"half_t"`, then the prior on `tau` is a
#'   half Student-t prior with center 0, lower bound 0, scale parameter
#'   `s_tau`, and degrees of freedom `d_tau`. The scale parameter `s_tau`
#'   is analogous to the `sigma` parameter of
#'   the Student-t parameterization given at
#'   <https://mc-stan.org/docs/functions-reference/unbounded_continuous_distributions.html>. # nolint
#' @param mu Numeric of length `n_rep`,
#'   mean of the control group means `alpha` for each rep.
#' @param tau Numeric of length `n_rep`,
#'   standard deviation of the control group means `alpha` for each rep.
#' @examples
#' hbl_sim_hierarchical(n_continuous = 1)$data
hbl_sim_hierarchical <- function(
  n_study = 5,
  n_group = 3,
  n_patient = 100,
  n_rep = 4,
  n_continuous = 0,
  n_binary = 0,
  constraint = FALSE,
  s_delta = 1,
  s_beta = 1,
  s_sigma = 1,
  s_lambda = 1,
  s_mu = 1,
  s_tau = 1,
  d_tau = 4,
  prior_tau = "half_t",
  covariance_current = "unstructured",
  covariance_historical = "unstructured",
  alpha = NULL,
  delta = stats::rnorm(
    n = (n_group - 1) * (n_rep - as.integer(constraint)),
    mean = 0,
    sd = s_delta
  ),
  beta = stats::rnorm(
    n = n_study * (n_continuous + n_binary),
    mean = 0,
    sd = s_delta
  ),
  sigma = stats::runif(n = n_study * n_rep, min = 0, max = s_sigma),
  mu = stats::rnorm(n = n_rep, mean = 0, sd = s_mu),
  tau = NULL,
  rho_current = stats::runif(n = 1, min = -1, max = 1),
  rho_historical = stats::runif(n = n_study - 1, min = -1, max = 1)
) {
  true(n_study, length(.) == 1, is.finite(.), is.numeric(.), . > 0)
  true(n_group, length(.) == 1, is.finite(.), is.numeric(.), . > 0)
  true(n_patient, length(.) == 1, is.finite(.), is.numeric(.), . > 0)
  true(n_rep, length(.) == 1, is.finite(.), is.numeric(.), . > 0)
  true(n_continuous, length(.) == 1, is.finite(.), is.numeric(.), . >= 0)
  true(n_binary, length(.) == 1, is.finite(.), is.numeric(.), . >= 0)
  true(n_study, length(.) == 1, is.finite(.), is.numeric(.), . > 0)
  true(constraint, is.logical(.), length(.) == 1L, !anyNA(.))
  true(s_delta, length(.) == 1, is.finite(.), is.numeric(.), . > 0)
  true(s_beta, length(.) == 1, is.finite(.), is.numeric(.), . > 0)
  true(s_sigma, length(.) == 1, is.finite(.), is.numeric(.), . > 0)
  true(s_mu, length(.) == 1, is.finite(.), is.numeric(.))
  true(s_tau, length(.) == 1, is.finite(.), is.numeric(.), . > 0)
  true(d_tau, length(.) == 1, is.finite(.), is.numeric(.), . > 0)
  if (is.null(tau)) {
    if (identical(as.character(prior_tau), "half_t")) {
      tau <- abs(stats::rt(n = n_rep, df = d_tau)) * s_tau
    } else if (identical(as.character(prior_tau), "uniform")) {
      tau <- stats::runif(n = n_rep, min = 0, max = s_tau)
    }
  }
  true(tau, is.numeric(.), is.finite(.), length(.) == n_rep, . > 0)
  if (is.null(alpha)) {
    alpha <- stats::rnorm(
      n = n_study * n_rep,
      mean = rep(mu, times = n_study),
      sd = rep(tau, times = n_study)
    )
  }
  true(alpha, is.finite(.), length(.) == n_study * n_rep)
  true(
    prior_tau,
    is.character(.),
    length(.) == 1L,
    !anyNA(.),
    as.character(.) %in% c("half_t", "uniform")
  )
  true(
    covariance_current,
    is.character(.),
    !anyNA(.),
    length(.) == 1L,
    . %in% c("unstructured", "ar1", "diagonal")
  )
  true(
    covariance_historical,
    is.character(.),
    !anyNA(.),
    length(.) == 1L,
    . %in% c("unstructured", "ar1", "diagonal")
  )
  true(
    delta,
    is.finite(.),
    is.numeric(.),
    length(.) == (n_group - 1) * (n_rep - as.integer(constraint))
  )
  true(beta, (all(is.finite(.)) || !length(beta)), is.numeric(.))
  true(length(beta) == n_study * (n_continuous + n_binary))
  true(sigma, is.finite(.), is.numeric(.), length(.) == n_study * n_rep)
  true(rho_historical, all(is.finite(.)) || !length(.), is.numeric(.))
  true(length(rho_historical) == n_study - 1)
  true(mu, is.numeric(.), is.finite(.), length(.) == n_rep)
  data <- hbl_sim_grid(n_study, n_group, n_patient, n_rep)
  x_alpha <- get_x_alpha(data, constraint = constraint)
  x_delta <- get_x_delta(data, constraint = constraint)
  covariates <- hbl_sim_x_beta(
    data = data,
    x_alpha = x_alpha,
    x_delta = x_delta,
    n_continuous = n_continuous,
    n_binary = n_binary
  )
  data <- dplyr::bind_cols(data, tibble::as_tibble(covariates))
  x_beta <- get_x_beta(data = data, x_alpha = x_alpha, x_delta = x_delta)
  sigma <- stats::runif(n = n_study * n_rep, min = 0, max = s_sigma)
  sigma <- matrix(sigma, nrow = n_study)
  lambda_current <- hbl_sim_lambda(
    n_matrix = 1,
    n_rep = n_rep,
    s_lambda = s_lambda
  )
  lambda_historical <- hbl_sim_lambda(
    n_matrix = n_study - 1,
    n_rep = n_rep,
    s_lambda = s_lambda
  )
  data$response <- hbl_sim_response(
    data = data,
    covariance_current = covariance_current,
    covariance_historical = covariance_historical,
    x_alpha = x_alpha,
    x_delta = x_delta,
    x_beta = x_beta,
    alpha = alpha,
    delta = delta,
    beta = beta,
    sigma = sigma,
    rho_current = rho_current,
    rho_historical = rho_historical,
    lambda_current = lambda_current,
    lambda_historical = lambda_historical
  )
  hbl_warn_identifiable(
    response = data$response,
    x_alpha = x_alpha,
    x_delta = x_delta,
    x_beta = x_beta
  )
  parameters <- list(
    alpha = alpha,
    delta = delta,
    beta = beta,
    sigma = sigma,
    mu = mu,
    tau = tau
  )
  if (covariance_current == "unstructured") {
    parameters$lambda_current <- lambda_current
  } else if (covariance_current == "ar1") {
    parameters$rho_current <- rho_current
  }
  if (covariance_historical == "unstructured") {
    parameters$lambda_historical <- lambda_historical
  } else if (covariance_historical == "ar1") {
    parameters$rho_historical <- rho_historical
  }
  matrices <- list(
    x_alpha = x_alpha,
    x_delta = x_delta,
    x_beta = x_beta
  )
  list(
    data = data,
    parameters = parameters,
    matrices = matrices
  )
}
