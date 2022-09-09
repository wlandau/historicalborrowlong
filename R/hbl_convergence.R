#' @title Check convergence diagnostics
#' @export
#' @family mcmc
#' @description Check the convergence diagnostics on a model.
#' @return A data frame of summarized convergence diagnostics.
#'   `max_rhat` is the maximum univariate Gelman/Rubin potential scale
#'   reduction factor over all the parameters of the model,
#'   `min_ess_bulk` is the minimum bulk effective sample size over the
#'   parameters, and `min_ess_tail` is the minimum tail effective
#'   sample size. `max_rhat` should be below 1.01, and the ESS metrics
#'   should both be above 100 times the number of MCMC chains. If
#'   any of these conditions are not true, the MCMC did not converge,
#'   and it is recommended to try running the model for more saved
#'   iterations (and if `max_rhat` is high, possibly more warmup
#'   iterations).
#' @inheritParams hbl_summary
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
#'     mcmc <- hbl_mcmc_pool(
#'       data,
#'       chains = 1,
#'       warmup = 10,
#'       iter = 20,
#'       seed = 0
#'     )
#'   )
#' )
#' hbl_convergence(mcmc)
#' }
hbl_convergence <- function(mcmc) {
  out <- posterior::summarize_draws(
    mcmc,
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )
  tibble::tibble(
    max_rhat = max(out$rhat, na.rm = TRUE),
    min_ess_bulk = min(out$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(out$ess_tail, na.rm = TRUE)
  )
}
