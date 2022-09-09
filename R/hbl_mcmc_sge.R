#' @title Run all MCMCs on a Sun Grid Engine (SGE) cluster.
#' @export
#' @family mcmc
#' @description Run all MCMCs on a Sun Grid Engine (SGE) cluster.
#'   Different models run in different jobs, and different chains run on
#'   different cores.
#' @return A list of tidy data frames of parameter samples from the
#'   posterior distribution.
#'   Columns `.chain`, `.iteration`,
#'   and `.draw` have the meanings documented in the
#'   `posterior` package.
#' @inheritSection hbl_data Data processing
#' @inheritParams hbl_mcmc_hierarchical
#' @inheritParams hbl_mcmc_independent
#' @inheritParams hbl_mcmc_pool
#' @inheritParams rstan::stan
#' @param data Tidy data frame with one row per patient per rep,
#'   indicator columns for the response variable,
#'   study, group, patient, rep,
#'   and covariates. All columns must be atomic vectors
#'   (e.g. not lists).
#' @param scheduler Either `"sge"` or `"local"`, high-performance computing
#'   scheduler / resource manager to use. Choose `"sge"` for serious use cases
#'   with a Sun Grid Engine (SGE) cluster. Otherwise, to run models
#'   sequentially on the current node, choose `"local"`.
#' @param log Character of length 1, path to a directory (with a trailing `/`)
#'   or a single file path. The SGE log files go here. Only works if
#'   `scheduler` is `"sge"`.
#' @examples
#' if (identical(Sys.getenv("HBL_SGE"), "true")) {
#' if (!identical(Sys.getenv("HBL_TEST", unset = ""), "")) {
#' set.seed(0)
#' data <- hbl_sim_hierarchical(
#'   n_study = 2,
#'   n_group = 2,
#'   n_patient = 5,
#'   n_rep = 3
#' )$data
#' tmp <- utils::capture.output(
#'   suppressWarnings(
#'     mcmc <- hbl_mcmc_sge(
#'       data,
#'       chains = 2,
#'       warmup = 10,
#'       iter = 20,
#'       seed = 0,
#'       scheduler = "local" # change to "sge" for serious runs
#'     )
#'   )
#' )
#' mcmc
#' }
#' }
hbl_mcmc_sge <- function(
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
  s_mu = 30,
  s_tau = 30,
  covariance_current = "unstructured",
  covariance_historical = "unstructured",
  control = list(max_treedepth = 17, adapt_delta = 0.99),
  log = "/dev/null",
  scheduler = "sge",
  chains = 1,
  cores = chains,
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
  true(s_mu, length(.) == 1, is.finite(.), is.numeric(.), . > 0)
  true(s_tau, length(.) == 1, is.finite(.), is.numeric(.), . > 0)
  true(is.list(control))
  true(log, is.character(.), length(.) == 1, !anyNA(.), nzchar(.))
  true(scheduler, is.character(.), length(.) == 1, !anyNA(.), nzchar(.))
  true(scheduler %in% c("sge", "local"))
  args <- list(
    data = data,
    response = response,
    study = study,
    study_reference = study_reference,
    group = group,
    group_reference = group_reference,
    patient = patient,
    rep = rep,
    rep_reference = rep_reference,
    covariates = covariates,
    constraint = constraint,
    s_alpha = s_alpha,
    s_delta = s_delta,
    s_beta = s_beta,
    s_sigma = s_sigma,
    s_lambda = s_lambda,
    s_mu = s_mu,
    s_tau = s_tau,
    covariance_current = covariance_current,
    covariance_historical = covariance_historical,
    control = control,
    chains = chains,
    cores = cores,
    ...
  )
  template <- system.file(
    "clustermq.tmpl",
    package = "historicalborrowlong",
    mustWork = TRUE
  )
  withr::local_options(
    .new = list(
      clustermq.scheduler = scheduler,
      clustermq.template = template
    )
  )
  models <- c(
    "hierarchical",
    "independent",
    "pool"
  )
  out <- clustermq::Q(
    fun = hbl_mcmc_sge_model,
    model = models,
    const = list(args = args),
    n_jobs = length(models),
    template = list(
      log = log,
      slots = args$cores %|||% 1,
      r_version = as.character(getRversion())
    )
  )
  names(out) <- models
  out
}

hbl_mcmc_sge_model <- function(model, args) {
  if (identical(model, "hierarchical")) {
    args$s_alpha <- NULL
    out <- do.call(
      what = historicalborrowlong::hbl_mcmc_hierarchical,
      args = args
    )
  } else if (identical(model, "pool")) {
    args$s_mu <- NULL
    args$s_tau <- NULL
    out <- do.call(
      what = historicalborrowlong::hbl_mcmc_pool,
      args = args
    )
  } else if (identical(model, "independent")) {
    args$s_mu <- NULL
    args$s_tau <- NULL
    out <- do.call(
      what = historicalborrowlong::hbl_mcmc_independent,
      args = args
    )
  } else {
    hbl_error(paste("unsupported model:", model))
  }
  attr(out, "model") <- model
  out
}
