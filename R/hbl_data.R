#' @title Standardize data
#' @export
#' @family data
#' @description Standardize a tidy input dataset.
#' @details Users do not normally need to call this function.
#'   It mainly serves exposes the indexing behavior of
#'   studies and group levels to aid in interpreting
#'   summary tables.
#' @return A standardized tidy data frame with one row per patient
#'   and the following columns:
#'   * `response`: continuous response/outcome variable. (Should be
#'     change from baseline of an outcome of interest.)
#'   * `study_label`: human-readable label of the study.
#'   * `study`: integer study index with the max index equal to the
#'     current study (at `study_reference`).
#'   * `group_label`: human-readable group label (e.g. treatment arm name).
#'   * `group`: integer group index with an index of 1 equal to the control
#'     group (at `group_reference`).
#'   * `patient_label`: original patient ID.
#'   * `patient`: integer patient index.
#'   * `rep_label`: original rep ID (e.g. time point or patient visit).
#'   * `rep`: integer rep index.
#'   * `covariate_*`: baseline covariate columns.
#' @section Data processing:
#'   Before running the MCMC, dataset is pre-processed.
#'   This includes expanding the rows of the data so every rep
#'   of every patient gets an explicit row. So if your
#'   original data has irregular rep IDs, e.g. unscheduled
#'   visits in a clinical trial that few patients attend,
#'   please remove them before the analysis. Only the most
#'   common rep IDs should be added.
#'
#'   After expanding the rows, the function fills in missing
#'   values for every column except the response. That includes
#'   covariates. Missing covariate values are filled in,
#'   first with last observation carried forward,
#'   then with last observation carried backward.
#'   If there are still missing values after this process,
#'   the program throws an informative error.
#' @param data A tidy data frame or `tibble` with the data.
#' @param response Character of length 1,
#'   name of the column in `data` with the response/outcome variable.
#'   `data[[response]]` must be a continuous variable,
#'   and it *should* be the change from baseline of a
#'   clinical endpoint of interest, as opposed to just
#'   the raw response. Treatment differences
#'   are computed directly from this scale, please supply
#'   change from baseline unless you are absolutely certain
#'   that treatment differences computed directly from
#'   this quantity are clinically meaningful.
#' @param study Character of length 1,
#'   name of the column in `data` with the study ID.
#' @param study_reference Atomic of length 1,
#'   element of the `study` column that indicates
#'   the current study.
#'   (The other studies are historical studies.)
#' @param group Character of length 1,
#'   name of the column in `data` with the group ID.
#' @param group_reference Atomic of length 1,
#'   element of the `group` column that indicates
#'   the control group.
#'   (The other groups may be treatment groups.)
#' @param patient Character of length 1,
#'   name of the column in `data` with the patient ID.
#' @param rep Character of length 1,
#'   name of the column in `data` with the rep ID.
#' @param rep_reference Atomic of length 1,
#'   element of the `rep` column that indicates
#'   baseline, i.e. the first rep chronologically.
#'   (The other reps may be post-baseline study visits or time points.)
#' @param covariates Character vector of column names
#'   in `data` with the columns with baseline covariates.
#'   These can be continuous, categorical, or binary.
#'   Regardless, `historicalborrowlong` derives the appropriate
#'   model matrix.
#'
#'   Each baseline covariate column must truly be a *baseline* covariate:
#'   elements must be equal for all time points within each patient
#'   (after the steps in the "Data processing" section).
#'   In other words, covariates must not be time-varying.
#'
#'   A large number of covariates, or a large number of levels in a
#'   categorical covariate, can severely slow down the computation.
#'   Please consider carefully if you really need to include
#'   such complicated baseline covariates.
#' @examples
#' set.seed(0)
#' data <- hbl_sim_independent(n_continuous = 1, n_study = 2)$data
#' data <- dplyr::select(
#'   data,
#'   study,
#'   group,
#'   rep,
#'   patient,
#'   response,
#'   tidyselect::everything()
#' )
#' data <- dplyr::rename(
#'   data,
#'   change = response,
#'   trial = study,
#'   arm = group,
#'   subject = patient,
#'   visit = rep,
#'   cov1 = covariate_study1_continuous1,
#'   cov2 = covariate_study2_continuous1
#' )
#' data$trial <- paste0("trial", data$trial)
#' data$arm <- paste0("arm", data$arm)
#' data$subject <- paste0("subject", data$subject)
#' data$visit <- paste0("visit", data$visit)
#' hbl_data(
#'   data = data,
#'   response = "change",
#'   study = "trial",
#'   study_reference = "trial1",
#'   group = "arm",
#'   group_reference = "arm1",
#'   patient = "subject",
#'   rep = "visit",
#'   rep_reference = "visit1",
#'   covariates = c("cov1", "cov2")
#' )
hbl_data <- function(
  data,
  response,
  study,
  study_reference,
  group,
  group_reference,
  patient,
  rep,
  rep_reference,
  covariates
) {
  true(is.data.frame(data))
  true(response, !anyNA(.), is.character(.), length(.) == 1)
  true(study, !anyNA(.), is.character(.), length(.) == 1)
  true(group, !anyNA(.), is.character(.), length(.) == 1)
  true(patient, !anyNA(.), is.character(.), length(.) == 1)
  true(rep, !anyNA(.), is.character(.), length(.) == 1)
  true(covariates, !anyNA(.), is.character(.))
  true(study_reference, !anyNA(.), length(.) == 1)
  true(group_reference, !anyNA(.), length(.) == 1)
  true(rep_reference, !anyNA(.), length(.) == 1)
  true(dim(data), length(.) == 2, . > 0)
  true(response %in% colnames(data))
  true(study %in% colnames(data))
  true(group %in% colnames(data))
  true(patient %in% colnames(data))
  true(rep %in% colnames(data))
  true(all(covariates %in% colnames(data)))
  true(data[[response]], is.atomic(.), is.vector(.))
  true(data[[study]], is.atomic(.), is.vector(.), !anyNA(.))
  true(data[[group]], is.atomic(.), is.vector(.), !anyNA(.))
  true(data[[patient]], is.atomic(.), is.vector(.), !anyNA(.))
  true(data[[rep]], is.atomic(.), is.vector(.), !anyNA(.))
  true(study_reference %in% data[[study]])
  true(group_reference %in% data[[group]])
  true(rep_reference %in% data[[rep]])
  true(is.numeric(data[[response]]))
  for (covariate in covariates) {
    true(data[[covariate]], is.atomic(.), is.vector(.), !anyNA(.))
  }
  out <- tibble::tibble(
    response = data[[response]],
    study_label = data[[study]],
    group_label = data[[group]],
    patient_label = data[[patient]],
    rep_label = data[[rep]],
    study = as_index_max(data[[study]], max = study_reference),
    group = as_index_min(data[[group]], min = group_reference),
    patient = as_index_min(data[[patient]]),
    rep = as_index_min(data[[rep]], min = rep_reference),
  )
  hbl_data_assert_unique_patients(out)
  for (x in covariates) {
    name <- if_any(grepl("^covariate", x), x, paste0("covariate_", x))
    out[[name]] <- data[[x]]
  }
  out <- tidyr::complete(
    out,
    tidyr::nesting(patient, patient_label),
    tidyr::nesting(rep, rep_label)
  )
  exclude <- c("response", "patient", "patient_label", "rep", "rep_label")
  out <- dplyr::arrange(out, patient, rep)
  for (col in setdiff(colnames(out), exclude)) {
    out[[col]] <- hbl_data_fill_column(out[[col]], out$patient)
    # Impossible to test but keeping anyway for defensiveness:
    # nocov start
    if (anyNA(out[[col]])) {
      hbl_error(
        sprintf("Missing values in column %s could not be filled in.", col)
      )
    }
    # nocov end
  }
  out <- hbl_data_enforce_baseline_covariates(out)
  dplyr::arrange(out, study, group, patient, rep)
}

as_index_min <- function(x, min = base::min(x)) {
  levels <- c(min, setdiff(sort(unique(x)), min))
  out <- as.integer(ordered(x, levels = levels))
  out - min(out) + 1L
}

as_index_max <- function(x, max = base::max(x)) {
  levels <- c(setdiff(sort(unique(x)), max), max)
  out <- as.integer(ordered(x, levels = levels))
  out - min(out) + 1L
}

hbl_data_fill_column <- function(x, index) {
  out <- tapply(
    X = x,
    INDEX = index,
    FUN = hbl_locf
  )
  unlist(out, use.names = FALSE)
}

hbl_locf <- function(x) {
  x <- zoo::na.locf(x, fromLast = FALSE, na.rm = FALSE)
  x <- zoo::na.locf(x, fromLast = TRUE, na.rm = FALSE)
  x
}

hbl_data_assert_unique_patients <- function(data) {
  dplyr::group_walk(
    dplyr::group_by(data, patient),
    ~ true(
      !anyDuplicated(.x$rep),
      message = paste(
        "Patient IDs must be unique.",
        "No duplicated rep ID for the same patient ID."
      )
    )
  )
}

hbl_data_enforce_baseline_covariates <- function(data) {
  out <- dplyr::group_modify(dplyr::group_by(data, study, patient), ~ {
    for (x in grep("^covariate_", colnames(data), value = TRUE)) {
      if (length(unique(.x[[x]])) != 1L) {
        hbl_warn(
          paste(
            "covariate", x, "of patient", .x$patient_label[1],
            "is time-varying.",
            "All covariates must be *baseline* covariates, not time-varying.",
            "All elements of", x, "and other covariates",
            "must be the same within each patient. The model assumes this fact",
            "and only uses the baseline level of each covariate."
          )
        )
        .x[[x]] <- .x[[x]][.x$rep == min(.x$rep)]
      }
    }
    .x
  })
  dplyr::ungroup(out)
}
