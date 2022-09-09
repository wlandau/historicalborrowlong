#' @title Longitudinal pooled MCMC
#' @export
#' @family mcmc
#' @description Run the longitudinal pooled model with MCMC.
#' @return A tidy data frame of parameter samples from the
#'   posterior distribution. Columns `.chain`, `.iteration`,
#'   and `.draw` have the meanings documented in the
#'   `posterior` package.
#' @inheritSection hbl_data Data processing
#' @inheritParams rstan::stan
#' @inheritParams hbl_data
#' @inheritParams hbl_sim_pool
#' @param data Tidy data frame with one row per patient per rep,
#'   indicator columns for the response variable,
#'   study, group, patient, rep,
#'   and covariates. All columns must be atomic vectors
#'   (e.g. not lists).
#' @param ... Additional named arguments of `rstan::sampling()`.
#'   See the documentation of `rstan::sampling()` for details.
#' @examples
#' if (!identical(Sys.getenv("HBL_TEST", unset = ""), "")) {
#' set.seed(0)
#' data <- hbl_sim_pool(
#'   n_study = 3,
#'   n_group = 2,
#'   n_patient = 5,
#'   n_rep = 3
#' )$data
#' tmp <- utils::capture.output(
#'   suppressWarnings(
#'     mcmc <- hbl_mcmc_pool(
#'       data,
#'       chains = 1,
#'       warmup = 10,
#'       iter = 20,
#'       seed = 0
#'     )
#'   )
#' )
#' mcmc
#' }
hbl_mcmc_pool <- function(
  data,
  response = "response",
  study = "study",
  study_reference = max(data[[study]]),
  group = "group",
  group_reference = min(data[[group]]),
  patient = "patient",
  rep = "rep",
  rep_reference = min(data[[rep]]),
  covariates = grep("^covariate", colnames(data), value = TRUE),
  constraint = FALSE,
  s_alpha = 30,
  s_delta = 30,
  s_beta = 30,
  s_sigma = 30,
  s_lambda = 1,
  covariance_current = "unstructured",
  covariance_historical = "unstructured",
  control = list(max_treedepth = 17, adapt_delta = 0.99),
  ...
) {
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
  true(constraint, is.logical(.), length(.) == 1L, !anyNA(.))
  true(s_alpha, length(.) == 1, is.finite(.), is.numeric(.), . > 0)
  true(s_delta, length(.) == 1, is.finite(.), is.numeric(.), . > 0)
  true(s_beta, length(.) == 1, is.finite(.), is.numeric(.), . > 0)
  true(s_sigma, length(.) == 1, is.finite(.), is.numeric(.), . > 0)
  true(s_lambda, length(.) == 1, is.finite(.), is.numeric(.), . > 0)
  true(is.list(control))
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
  x_alpha <- get_x_alpha_pool(data, constraint = constraint)
  x_delta <- get_x_delta(data, constraint = constraint)
  x_beta <- get_x_beta(data = data, x_alpha = x_alpha, x_delta = x_delta)
  alpha_rep_index <- rep(1L, length(unique(data$rep)))
  args <- list(control = control, ...)
  stan_mcmc(
    model_type = model_pool,
    data = data,
    n_mu = 0L,
    n_tau = 0L,
    x_alpha = x_alpha,
    x_delta = x_delta,
    x_beta = x_beta,
    s_alpha = s_alpha,
    s_delta = s_delta,
    s_beta = s_beta,
    s_sigma = s_sigma,
    s_lambda = s_lambda,
    s_mu = 30,
    s_tau = 30,
    alpha_rep_index = alpha_rep_index,
    covariance_current = covariance_current,
    covariance_historical = covariance_historical,
    args = args
  )
}
