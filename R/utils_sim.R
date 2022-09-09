hbl_sim_grid <- function(n_study, n_group, n_patient, n_rep) {
  out <- tidyr::expand_grid(
    study = seq_len(n_study),
    group = seq_len(n_group),
    patient = seq_len(n_patient)
  )
  out <- dplyr::filter(out, study == max(study) | group == min(group))
  out$patient <- seq_len(nrow(out))
  out <- tidyr::expand_grid(out, rep = seq_len(n_rep))
  out
}

hbl_sim_x_beta <- function(
  n_continuous,
  n_binary,
  data,
  x_alpha,
  x_delta
) {
  data_full <- data
  x_alpha_full <- x_alpha
  x_delta_full <- x_delta
  n_rep <- length(unique(data_full$rep))
  index <- data_full$rep == min(data_full$rep)
  data <- data_full[index,, drop = FALSE] # nolint
  x_alpha <- x_alpha[index,, drop = FALSE] # nolint
  x_delta <- x_delta[index,, drop = FALSE] # nolint
  out <- NULL
  try_x_beta_full <- 0
  while (
    is.null(out) || !is_full_rank(cbind(x_alpha_full, x_delta_full, out))
  ) {
    try_x_beta_full <- try_x_beta_full + 1
    true(try_x_beta_full < 1000)
    out <- list()
    for (study in sort(unique(data$study))) {
      out[[study]] <- hbl_sim_x_beta_study(
        n_continuous = n_continuous,
        n_binary = n_binary,
        data = data[data$study == study,, drop = FALSE], # nolint
        x_alpha = x_alpha[data$study == study,, drop = FALSE], # nolint
        x_delta = x_delta[data$study == study,, drop = FALSE] # nolint
      )
    }
    out <- as.matrix(Matrix::bdiag(out))
    studies <- rep(sort(unique(data$study)), each = n_continuous + n_binary)
    covariates <- c(
      sprintf("continuous%s", seq_len(n_continuous)),
      sprintf("binary%s", seq_len(n_binary))
    )
    covariates <- rep(covariates, times = length(unique(studies)))
    colnames(out) <- sprintf(
      "covariate_study%s_%s",
      studies,
      covariates
    )
    index <- rep(seq_len(nrow(data)), each = n_rep)
    out <- out[index,, drop = FALSE] # nolint
  }
  out
}

hbl_sim_x_beta_study <- function(
  n_continuous,
  n_binary,
  data,
  x_alpha,
  x_delta
) {
  out <- NULL
  try_x_beta <- 0
  x_alpha <- drop_zero_columns(x_alpha)
  x_delta <- drop_zero_columns(x_delta)
  while (is.null(out) || !is_full_rank(cbind(x_alpha, x_delta, out))) {
    try_x_beta <- try_x_beta + 1
    true(try_x_beta < 1000)
    out <- matrix(numeric(0), nrow = nrow(data), ncol = 0)
    for (index in seq_len(n_continuous)) {
      x <- stats::rnorm(n = nrow(data), mean = 0, sd = 1)
      out <- cbind(out, x)
    }
    for (index in seq_len(n_binary)) {
      x <- sample(c(0L, 1L), size = nrow(data), replace = TRUE)
      out <- cbind(out, x)
    }
  }
  out
}

hbl_sim_response <- function(
  data,
  covariance_current,
  covariance_historical,
  x_alpha,
  x_delta,
  x_beta,
  alpha,
  delta,
  beta,
  sigma,
  rho_current,
  rho_historical,
  lambda_current,
  lambda_historical
) {
  n_study <- max(data$study)
  n_patient <- max(data$patient)
  n_rep <- max(data$rep)
  true(all(data$patient == rep(seq_len(max(data$patient)), each = n_rep)))
  true(all(data$rep == rep(seq_len(n_rep), times = n_patient)))
  means <- x_alpha %*% alpha +
    x_delta %*% delta +
    x_beta %*% beta
  means <- as.numeric(means)
  residuals <- numeric(0)
  correlations <- hbl_get_correlations(
    n_study = n_study,
    covariance_current = covariance_current,
    covariance_historical = covariance_historical,
    n_rep = n_rep,
    rho_current = rho_current,
    rho_historical = rho_historical,
    lambda_current = lambda_current,
    lambda_historical = lambda_historical
  )
  for (patient in seq_len(n_patient)) {
    study <- unique(data$study[data$patient == patient])
    true(length(study) == 1)
    sigma_patient <- diag(sigma[study,, drop = TRUE]) # nolint
    correlation_patient <- correlations[[study]]
    covariance <- sigma_patient %*% correlation_patient %*% sigma_patient
    epsilon <- MASS::mvrnorm(n = 1L, mu = rep(0, n_rep), Sigma = covariance)
    epsilon <- as.numeric(epsilon)
    residuals <- c(residuals, epsilon)
  }
  means + residuals
}

hbl_get_correlations <- function(
  n_study,
  covariance_current,
  covariance_historical,
  n_rep,
  rho_current,
  rho_historical,
  lambda_current,
  lambda_historical
) {
  out <- list()
  for (study in seq_len(n_study)) {
    if (study == n_study) {
      covariance <- covariance_current
      rho <- rho_current
      lambda <- lambda_current[1,, .drop = TRUE] # nolint
    } else {
      covariance <- covariance_historical
      rho <- rho_historical[study]
      lambda <- lambda_historical[study,, .drop = TRUE] # nolint
    }
    if (covariance == "unstructured") {
      out[[study]] <- lambda %*% t(lambda)
    } else if (covariance == "ar1") {
      out[[study]] <- ar1_correlation(n = n_rep, rho = rho)
    } else if (covariance == "diagonal") {
      out[[study]] <- diag(n_rep)
    }
  }
  out
}

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4455603/
# section 5.1: Generating random AR(1) structures
ar1_correlation <- function(rho, n) {
  scale <- sqrt(1 - (rho * rho))
  out <- matrix(0, nrow = n, ncol = n)
  out[1, 1] <- 1
  for (i in seq(2, n)) {
    out[i, 1] <- rho ^ (i - 1)
  }
  for (i in seq(2, n)) {
    for (j in seq(2, i)) {
      out[i, j] <- scale * rho ^ (i - j)
    }
  }
  out %*% t(out)
}

hbl_sim_lambda <- function(n_matrix, n_rep, s_lambda) {
  out <- trialr::rlkjcorr(n = n_matrix, K = n_rep, eta = s_lambda)
  if (length(dim(out)) == 2L) {
    return(array(t(chol(out)), dim = c(1, n_rep, n_rep)))
  }
  for (index in seq_len(n_matrix)) {
    out[index,, .drop = TRUE] <- # nolint
      t(chol(out[index,, .drop = TRUE])) # nolint
  }
  out
}
