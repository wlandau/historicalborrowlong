#' @title Longitudinal pooled simulations.
#' @export
#' @family simulate
#' @description Simulate from the longitudinal pooled model.
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
#' @param n_study Number of studies to simulate.
#' @param n_group Number of groups (e.g. study arms) to simulate per study.
#' @param n_patient Number of patients to simulate per study per group.
#' @param n_rep Number of repeated measures (time points) per patient.
#' @param n_continuous Number of continuous covariates to simulate
#'   (all from independent standard normal distributions).
#' @param n_binary Number of binary covariates to simulate
#'   (all from independent Bernoulli distributions with p = 0.5).
#' @param constraint Logical of length 1, whether to pool all study arms
#'   at baseline (first rep). Appropriate when the response is the raw
#'   response (as opposed to change from baseline) and the first rep
#'   (i.e. time point) is prior to treatment.
#' @param s_alpha Numeric of length 1, prior standard deviation
#'   of the study-specific control group mean parameters `alpha`.
#' @param s_delta Numeric of length 1, prior standard deviation
#'   of the study-by-group effect parameters `delta`.
#' @param s_beta Numeric of length 1, prior standard deviation
#'   of the fixed effects `beta`.
#' @param s_sigma Numeric of length 1, prior upper bound
#'   of the residual standard deviations.
#' @param s_lambda shape parameter of the LKJ priors
#'   on the unstructured correlation matrices.
#' @param covariance_current Character of length 1,
#'   covariance structure of the current study.
#'   Possible values are `"unstructured"` for fully parameterized
#'   covariance matrices, `"ar1"` for AR(1) covariance matrices,
#'   and `"diagonal"` for residuals independent across time within
#'   each patient. In MCMC (e.g. [hbl_mcmc_hierarchical()]),
#'   the covariance structure affects computational speed.
#'   Unstructured covariance is slower than AR(1), and AR(1) is slower
#'   than diagonal. This is particularly true for `covariance_historical`
#'   if there are many historical studies in the data.
#' @param covariance_historical Same as `covariance_current`,
#'   but for the covariance structure of each separate historical study.
#'   Each historical study has its own separate covariance matrix.
#' @param alpha Numeric vector of length `n_rep` for the pooled
#'   and model and length `n_study * n_rep` for the
#'   independent and hierarchical models.
#'   `alpha` is the vector of control group mean parameters.
#'   `alpha` enters the model by multiplying with
#'   `$matrices$x_alpha` (see the return value).
#'   The control group in the data is the one with the
#'   `group` column equal to 1.
#' @param delta Numeric vector of length
#'   `(n_group - 1) * (n_rep - as.integer(constraint))`
#'   of treatment effect parameters.
#'   `delta` enters the model by multiplying with
#'   `$matrices$x_delta` (see the return value).
#'   The control (non-treatment) group in the data is the one with the
#'   `group` column equal to 1.
#' @param beta Numeric vector of `n_study * (n_continuous + n_binary)`
#'   fixed effect parameters. Within each study,
#'   the first `n_continuous` betas
#'   are for the continuous covariates, and the rest are for
#'   the binary covariates. All the `beta`s for one study
#'   appear before all the `beta`s for the next study,
#'   and studies are arranged in increasing order of
#'   the sorted unique values in `$data$study` in the output.
#'   `betas` enters the model by multiplying with
#'   `$matrices$x_alpha` (see the return value).
#' @param sigma Numeric vector of `n_study * n_rep`
#'   residual standard deviation parameters for each study
#'   and rep. The elements are sorted with all the standard deviations
#'   of study 1 first (all the reps), then all the reps of study 2, etc.
#' @param rho_current Numeric of length 1 between -1 and 1,
#'    AR(1) residual correlation parameter for the current study.
#' @param rho_historical Numeric of length `n_study - 1` between -1 and 1,
#'    AR(1) residual correlation parameters for the historical studies.
#' @examples
#' hbl_sim_pool(n_continuous = 1)$data
hbl_sim_pool <- function(
  n_study = 5,
  n_group = 3,
  n_patient = 100,
  n_rep = 4,
  n_continuous = 0,
  n_binary = 0,
  constraint = FALSE,
  s_alpha = 1,
  s_delta = 1,
  s_beta = 1,
  s_sigma = 1,
  s_lambda = 1,
  covariance_current = "unstructured",
  covariance_historical = "unstructured",
  alpha = stats::rnorm(n = n_rep, mean = 0, sd = s_alpha),
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
  rho_current = stats::runif(n = 1, min = -1, max = 1),
  rho_historical = stats::runif(n = n_study - 1, min = -1, max = 1)
) {
  true(n_study, length(.) == 1, is.finite(.), is.numeric(.), . > 0)
  true(n_group, length(.) == 1, is.finite(.), is.numeric(.), . > 0)
  true(n_patient, length(.) == 1, is.finite(.), is.numeric(.), . > 0)
  true(n_rep, length(.) == 1, is.finite(.), is.numeric(.), . > 0)
  true(n_continuous, length(.) == 1, is.finite(.), is.numeric(.), . >= 0)
  true(n_binary, length(.) == 1, is.finite(.), is.numeric(.), . >= 0)
  true(n_study, length(.) == 1, is.finite(.), is.numeric(.), (. > 0))
  true(constraint, is.logical(.), length(.) == 1L, !anyNA(.))
  true(s_alpha, length(.) == 1, is.finite(.), is.numeric(.), (. > 0))
  true(s_delta, length(.) == 1, is.finite(.), is.numeric(.), (. > 0))
  true(s_beta, length(.) == 1, is.finite(.), is.numeric(.), (. > 0))
  true(s_sigma, length(.) == 1, is.finite(.), is.numeric(.), (. > 0))
  true(s_lambda, length(.) == 1, is.finite(.), is.numeric(.), (. > 0))
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
  true(alpha, length(.) == n_rep, is.finite(.))
  true(
    delta,
    is.finite(.),
    is.numeric(.),
    length(.) == (n_group - 1) * (n_rep - as.integer(constraint))
  )
  true(beta, all(is.finite(.)) || !length(.), is.numeric(.))
  true(length(beta) == n_study * (n_continuous + n_binary))
  true(sigma, is.finite(.), is.numeric(.), length(.) == n_study * n_rep)
  true(rho_current, all(is.finite(.)) || !length(.), is.numeric(.))
  true(length(rho_current) == 1)
  true(rho_historical, all(is.finite(.)) || !length(.), is.numeric(.))
  true(length(rho_historical) == n_study - 1)
  data <- hbl_sim_grid(n_study, n_group, n_patient, n_rep)
  x_alpha <- get_x_alpha_pool(data, constraint = constraint)
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
    sigma = sigma
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
