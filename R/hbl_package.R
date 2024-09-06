#' historicalborrowlong: Bayesian longitudinal historical borrowing models
#'   for clinical studies.
#' @description Bayesian longitudinal historical borrowing models for
#'   clinical studies.
#' @name historicalborrowlong-package
#' @family help
#' @useDynLib historicalborrowlong, .registration = TRUE
#' @importFrom clustermq Q
#' @importFrom dplyr across arrange bind_cols bind_rows distinct filter group_by
#'   group_modify left_join mutate n rename select summarize ungroup
#' @importFrom ggplot2 aes geom_errorbar facet_wrap geom_point ggplot
#'   position_dodge theme_gray xlab ylab
#' @importFrom MASS mvrnorm
#' @importFrom Matrix bdiag
#' @importFrom methods callNextMethod
#' @importFrom posterior as_draws_df mcse_mean mcse_sd mcse_quantile
#' @importFrom Rcpp compileAttributes
#' @importFrom RcppParallel CxxFlags
#' @importFrom rlang abort
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom stats as.formula model.matrix qnorm quantile rnorm runif
#'   sd update var
#' @importFrom tibble as_tibble
#' @importFrom tidyr complete expand_grid fill nesting pivot_longer pivot_wider
#' @importFrom tidyselect any_of everything starts_with
#' @importFrom trialr rlkjcorr
#' @importFrom utils capture.output globalVariables
#' @importFrom withr local_options
#' @importFrom zoo na.locf
NULL

utils::globalVariables(
  c(
    ".",
    "..density..",
    "Group",
    "group",
    "group_label",
    "data_mean",
    "data_sd",
    "data_lower",
    "data_upper",
    "Model",
    "patient",
    "patient_label",
    "Rep",
    "rep",
    "rep_label",
    "response",
    "response_mean",
    "response_lower",
    "response_upper",
    "study",
    "study_index",
    "study_label",
    "tau",
    "Value",
    "value",
    "value_percent"
  )
)
