#' @title Plot the groups of the hierarchical model
#'   and its benchmark models.
#' @export
#' @family plot
#' @description Plot the groups against one another
#'   for a hierarchical model.
#'   and the independent and pooled benchmark models.
#' @return A `ggplot` object
#' @inheritParams hbl_plot_borrow
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
#' hbl_plot_group(
#'   borrow = borrow,
#'   pool = pool,
#'   independent = independent
#' )
#' }
hbl_plot_group <- function(
  borrow,
  pool,
  independent,
  outcome = c("response", "change", "diff")
) {
  outcome <- match.arg(outcome)
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
  for (name in paste0(outcome, c("_mean", "_lower", "_upper"))) {
    true(name %in% colnames(borrow))
    true(name %in% colnames(pool))
    true(name %in% colnames(independent))
  }
  borrow <- dplyr::select(
    borrow,
    group,
    group_label,
    rep,
    rep_label,
    tidyselect::starts_with(outcome)
  )
  pool <- dplyr::select(
    pool,
    group,
    group_label,
    rep,
    rep_label,
    tidyselect::starts_with(outcome)
  )
  independent <- dplyr::select(
    independent,
    group,
    group_label,
    rep,
    rep_label,
    tidyselect::starts_with(outcome)
  )
  out <- dplyr::bind_rows(
    `1-independent` = independent,
    `2-borrow` = borrow,
    `3-pool` = pool,
    .id = "Model"
  )
  out$Group <- as.character(out$group_label)
  out$Rep <- as.character(out$rep_label)
  out <- out[!is.na(out[[paste0(outcome, "_mean")]]),, drop = FALSE] # nolint
  args_point <- list(
    x = as.symbol("Rep"),
    y = as.symbol(paste0(outcome, "_mean")),
    color = as.symbol("Group")
  )
  args_errorbar <- list(
    x = as.symbol("Rep"),
    ymin = as.symbol(paste0(outcome, "_lower")),
    ymax = as.symbol(paste0(outcome, "_upper")),
    color = as.symbol("Group")
  )
  ggplot2::ggplot(out) +
    ggplot2::geom_point(
      do.call(what = ggplot2::aes, args = args_point),
      position = ggplot2::position_dodge(width = 0.5)
    ) +
    ggplot2::geom_errorbar(
      do.call(what = ggplot2::aes, args = args_errorbar),
      position = ggplot2::position_dodge(width = 0.5)
    ) +
    ggplot2::facet_wrap(~Model) +
    ggplot2::ylab(paste("Posterior", outcome)) +
    ggplot2::theme_gray(20) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
