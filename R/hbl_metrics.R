#' @title Borrowing metrics
#' @export
#' @family summary
#' @description Calculate historical borrowing metrics using
#'   summary output from a fitted borrowing model and
#'   analogous summaries from the benchmark models.
#' @return A data frame with borrowing metrics.
#' @param borrow A data frame returned by [hbl_summary()]
#'   for the hierarchical model.
#' @param pool A data frame returned by [hbl_summary()]
#'   for the pooled model.
#' @param independent A data frame returned by [hbl_summary()]
#'   for the independent model.
#' @examples
#' if (!identical(Sys.getenv("HBL_TEST", unset = ""), "")) {
#' set.seed(0)
#' data <- hbl_sim_independent(
#'   n_study = 2,
#'   n_group = 2,
#'   n_patient = 5,
#'   n_rep = 3
#' )$data
#' tmp <- utils::capture.output(
#'   suppressWarnings(
#'     mcmc_borrow <- hbl_mcmc_hierarchical(
#'       data,
#'       chains = 1,
#'       warmup = 10,
#'       iter = 20,
#'       seed = 0
#'     )
#'   )
#' )
#' tmp <- utils::capture.output(
#'   suppressWarnings(
#'     mcmc_pool <- hbl_mcmc_pool(
#'       data,
#'       chains = 1,
#'       warmup = 10,
#'       iter = 20,
#'       seed = 0
#'     )
#'   )
#' )
#' tmp <- utils::capture.output(
#'   suppressWarnings(
#'     mcmc_independent <- hbl_mcmc_independent(
#'       data,
#'       chains = 1,
#'       warmup = 10,
#'       iter = 20,
#'       seed = 0
#'     )
#'   )
#' )
#' borrow <- hbl_summary(mcmc_borrow, data)
#' pool <- hbl_summary(mcmc_pool, data)
#' independent <- hbl_summary(mcmc_independent, data)
#' hbl_metrics(
#'   borrow = borrow,
#'   pool = pool,
#'   independent = independent
#' )
#' }
hbl_metrics <- function(borrow, pool, independent) {
  true(is.data.frame(borrow))
  true(is.data.frame(pool))
  true(is.data.frame(independent))
  true(nrow(borrow) == nrow(pool))
  true(nrow(borrow) == nrow(independent))
  true("group" %in% colnames(borrow))
  true("rep" %in% colnames(borrow))
  true(all(borrow$group == pool$group))
  true(all(borrow$group == independent$group))
  true(all(borrow$rep == pool$rep))
  true(all(borrow$rep == independent$rep))
  for (name in c("response_mean", "response_variance")) {
    true(name %in% colnames(borrow))
    true(name %in% colnames(pool))
    true(name %in% colnames(independent))
  }
  borrow <- dplyr::filter(borrow, group == min(group))
  pool <- dplyr::filter(pool, group == min(group))
  independent <- dplyr::filter(independent, group == min(group))
  mean_shift_ratio <- (borrow$response_mean - independent$response_mean) /
    (pool$response_mean - independent$response_mean)
  variance_shift_ratio <-
    (borrow$response_variance - independent$response_variance) /
    (pool$response_variance - independent$response_variance)
  tibble::tibble(
    rep = borrow$rep,
    rep_label = borrow$rep_label,
    mean_shift_ratio = mean_shift_ratio,
    variance_shift_ratio = variance_shift_ratio
  )
}
