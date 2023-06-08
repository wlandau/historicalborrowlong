stan_mcmc <- function(
    model_type,
    data,
    n_mu,
    n_tau,
    x_alpha,
    x_delta,
    x_beta,
    s_alpha,
    s_delta,
    s_beta,
    s_sigma,
    s_lambda,
    s_mu,
    s_tau,
    alpha_rep_index,
    covariance_current,
    covariance_historical,
    args
) {
  response <- data$response
  hbl_warn_identifiable(
    response = response,
    x_alpha = x_alpha,
    x_delta = x_delta,
    x_beta = x_beta
  )
  which_missing <- which(is.na(response))
  missing <- rep(0L, length(response))
  missing[which_missing] <- 1L
  response[which_missing] <- -99999
  count_missing <- cumsum(missing)
  n_missing <- sum(missing)
  study_index <- vapply(
    X = data$study,
    FUN = match,
    table = sort(unique(data$study)),
    FUN.VALUE = integer(1)
  )
  n_study <- length(unique(data$study))
  n_patient <- length(unique(data$patient))
  n_rep <- length(unique(data$rep))
  n_lambda_current <- as.integer(identical(covariance_current, "unstructured"))
  n_lambda_historical <- (n_study - 1L) * as.integer(
    identical(covariance_historical, "unstructured")
  )
  n_rho_current <- as.integer(identical(covariance_current, "ar1"))
  n_rho_historical <- (n_study - 1L) * as.integer(
    identical(covariance_historical, "ar1")
  )
  patients <- dplyr::distinct(data, study, patient)
  patients <- dplyr::arrange(patients, study)
  # Use index vectors instead of multiplying x_alpha or x_delta.
  x_alpha_sweep <- sweep(x_alpha, 2, seq_len(ncol(x_alpha)), "*")
  alpha_data_index <- rowSums(x_alpha_sweep) + 1L
  x_delta_sweep <- sweep(x_delta, 2, seq_len(ncol(x_delta)), "*")
  delta_data_index <- rowSums(x_delta_sweep) + 1L
  # Multiply x_beta block by block and use the transpose for efficiency.
  rownames(x_beta) <- study_index
  x_beta <- t(x_beta[seq_len(n_patient) * n_rep, ])
  col <- as.integer(colnames(x_beta))
  row <- as.integer(gsub("^study|_.*", "", rownames(x_beta)))
  studies <- sort(unique(row))
  n_study_x_beta <- length(studies)
  x_beta_col_n <- as.integer(table(col))[studies]
  x_beta_row_n <- as.integer(table(row))
  x_beta_col_index <- (cumsum(x_beta_col_n) - x_beta_col_n + 1L)[studies]
  x_beta_row_index <- cumsum(x_beta_row_n) - x_beta_row_n + 1L
  study_patient_match <- dplyr::distinct(data, patient, study)
  study_patient_match <- dplyr::arrange(study_patient_match, patient)
  data_stan <- list(
    model_type = model_type,
    n_alpha = ncol(x_alpha),
    n_mu = n_mu,
    n_tau = n_tau,
    n_delta = ncol(x_delta),
    n_beta = nrow(x_beta),
    n_observe = nrow(data),
    n_missing = n_missing,
    n_patient = n_patient,
    n_rep = n_rep,
    n_study = n_study,
    n_study_x_beta = n_study_x_beta,
    n_lambda_current = n_lambda_current,
    n_lambda_historical = n_lambda_historical,
    n_rho_current = n_rho_current,
    n_rho_historical = n_rho_historical,
    n_patient_study = as.integer(table(patients$study)),
    index_patient_study = match(seq_len(n_study), patients$study),
    index_patient = rep(seq_len(n_patient), each = n_rep),
    s_alpha = s_alpha,
    s_delta = s_delta,
    s_beta = s_beta,
    s_sigma = s_sigma,
    s_lambda = s_lambda,
    s_mu = s_mu,
    s_tau = s_tau,
    missing = missing,
    alpha_rep_index = alpha_rep_index,
    alpha_data_index = alpha_data_index,
    delta_data_index = delta_data_index,
    x_beta_col_index = x_beta_col_index,
    x_beta_row_index = x_beta_row_index,
    x_beta_col_n = x_beta_col_n,
    x_beta_row_n = x_beta_row_n,
    count_missing = count_missing,
    study_index = study_index,
    study_patient = study_patient_match$study,
    covariance_current = stan_covariance(covariance_current),
    covariance_historical = stan_covariance(covariance_historical),
    y = response,
    x_alpha = x_alpha,
    x_delta = x_delta,
    x_beta = x_beta,
    covariance_unstructured = 1L,
    covariance_ar1 = 2L,
    covariance_diagonal = 3L
  )
  args$object <- eval(parse(text = "stanmodels$historicalborrowlong"))
  args$data <- data_stan
  pars <- c(
    "alpha",
    "delta",
    "beta",
    "sigma",
    "lambda_current",
    "lambda_historical",
    "rho_current",
    "rho_historical",
    "mu",
    "tau",
    "lp__"
  )
  args$pars <- pars
  fit <- do.call(what = rstan::sampling, args = args)
  out <- tibble::as_tibble(as.data.frame(fit))
  n <- args$iter - args$warmup
  out$.chain <- rep(seq_len(args$chains), each = n)
  out$.iteration <- rep(seq_len(n), times = args$chains)
  out$.draw <- seq_len(nrow(out))
  out <- dplyr::select(
    out,
    tidyselect::starts_with(pars), tidyselect::starts_with(".")
  )
  out <- dplyr::select(
    out,
    -tidyselect::contains("_latent"),
    -tidyselect::contains("_raw")
  )
  stan_prune_lambda(out)
}

stan_prune_lambda <- function(mcmc) {
  columns_lambda <- grep(
    "^lambda_current|^lambda_historical",
    colnames(mcmc),
    value = TRUE
  )
  columns_not_lambda <- setdiff(colnames(mcmc), columns_lambda)
  root <- gsub("^[^,]*,|\\]", "", columns_lambda)
  row <- as.integer(gsub(",.*", "", root))
  column <- as.integer(gsub(".*,", "", root))
  columns_lambda_keep <- columns_lambda[row >= column & (row + column > 2)]
  mcmc[, c(columns_not_lambda, columns_lambda_keep)]
}

stan_covariance <- function(covariance) {
  switch(
    covariance,
    unstructured = 1L,
    ar1 = 2L,
    diagonal = 3L,
    hbl_error(paste("unsupported covariance choice:", covariance))
  )
}

model_pool <- 1L
model_independent <- 2L
model_hierarchical <- 3L
