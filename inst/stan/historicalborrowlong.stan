functions {
/**
  * Compute the normal kernel of multiple observations from a
  *   multivariate normal density
  *   up to a proportionality constant on the log scale:
  *   `-0.5 * epsilon' * (LL')^{-1} * e) - 0.5 * log(det(LL'))`
  *   for each column e (iid observation) of matrix epsilon.
  *   Here, L is the lower Cholesky factor of the covariance matrix,
  *   e = y - E(y), y is the uncentered observation, and det() is the
  *   determinant of a matrix (det(LL') is just the product of the
  *   diagonals of L because L is triangular.)
  *   The expression `epsilon' * (LL')^{-1} * epsilon` factors into
  *   the inner product of `L^{-1} * epsilon` with itself. Below,
  *   `mdivide_left_tri_low(cholesky, epsilon)` computes `L^{-1} * epsilon`,
  *   and then sum(columns_dot_self()) computes the inner product.
  * @return Real scalar with the density evaluated at epsilon.
  * @param cholesky Lower-triangular square matrix with one row per
  *   dimension of the multivariate normal. Cholesky factor of the
  *   covariance matrix.
  * @param epsilon Matrix with one row per dimension of the multivariate normal
  *   and one column per data record.
  */
  row_vector log_normal_kernel(matrix cholesky, matrix epsilon) {
    return -0.5 * columns_dot_self(mdivide_left_tri_low(cholesky, epsilon))
      -sum(log(diagonal(cholesky)));
  }
/** 
  * Generate a Cholesky factor of an AR(1) correlation matrix.
  * @references
  *   Madar (2016). "Direct formulation to Cholesky decomposition
  *     of a general nonsingular correlation matrix."
  *     Section 5.1: Generating random AR(1) structures.
  *     Statistics and Probability Letters.
  *     https://doi.org/10.1016%2Fj.spl.2015.03.014
  *     https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4455603/
  * @return An n by n lower-triangular Cholesky factor of an AR(1)
  *   correlation matrix with correlation rho.
  * @param rho Real scalar, correlation between adjacent discrete time points.
  * @param n Positive integer scalar, number of rows of the returned matrix.
  */
  matrix ar1_cholesky(real rho, int n) {
    real scale = sqrt(1 - (rho * rho));
    matrix[n, n] out = rep_matrix(0, n, n);
    out[1, 1] = 1.0;
    for (i in 2:n) {
      out[i, 1] = rho ^ (i - 1);
    }
    for (i in 2:n) {
      for (j in 2:i) {
        out[i, j] = scale * rho ^ (i - j);
      }
    }
    return out;
  }
}
data {
  int<lower=1,upper=3> model_type;
  int<lower=0> n_alpha;
  int<lower=0> n_mu;
  int<lower=0> n_tau;
  int<lower=0> n_delta;
  int<lower=0> n_beta;
  int<lower=0> n_observe;
  int<lower=0> n_missing;
  int<lower=0> n_patient;
  int<lower=0> n_rep;
  int<lower=0> n_study;
  int<lower=0> n_study_x_beta;
  int<lower=0> n_lambda_current;
  int<lower=0> n_lambda_historical;
  int<lower=0> n_rho_current;
  int<lower=0> n_rho_historical;
  array[n_study] int<lower=0> n_patient_study;
  array[n_study] int<lower=0> index_patient_study;
  array[n_observe] int<lower=0> index_patient;
  real<lower=0> s_alpha;
  real<lower=0> s_mu;
  real<lower=0> s_tau;
  real<lower=0> s_beta;
  real<lower=0> s_delta;
  real<lower=0> s_sigma;
  real<lower=0> s_lambda;
  array[n_observe] int<lower=0> missing;
  array[n_observe] int<lower=0> count_missing;
  array[n_observe] int<lower=0> study_index;
  array[n_alpha] int<lower=0> alpha_rep_index;
  array[n_observe] int<lower=0> alpha_data_index;
  array[n_observe] int<lower=0> delta_data_index;
  array[n_study_x_beta] int<lower=0> x_beta_col_index;
  array[n_study_x_beta] int<lower=0> x_beta_row_index;
  array[n_study_x_beta] int<lower=0> x_beta_col_n;
  array[n_study_x_beta] int<lower=0> x_beta_row_n;
  array[n_patient] int<lower=0,upper=n_patient> study_patient;
  vector[n_observe] y;
  matrix[n_observe, n_alpha] x_alpha;
  matrix[n_observe, n_delta] x_delta;
  matrix[n_beta, n_patient] x_beta;
  int<lower=1,upper=3> covariance_current;
  int<lower=1,upper=3> covariance_historical;
  int covariance_unstructured;
  int covariance_ar1;
  int covariance_diagonal;
}
parameters {
  vector[n_missing] y_missing;
  vector[n_alpha] alpha_raw;
  vector[n_mu] mu;
  vector<lower=0,upper=s_tau>[n_tau] tau;
  vector[n_delta] delta;
  vector[n_beta] beta;
  array[n_study] vector<lower=0,upper=s_sigma>[n_rep] sigma;
  array[n_lambda_current] cholesky_factor_corr[n_rep] lambda_current;
  array[n_lambda_historical] cholesky_factor_corr[n_rep] lambda_historical;
  vector<lower=-1,upper=1>[n_rho_current] rho_current;
  vector<lower=-1,upper=1>[n_rho_historical] rho_historical;
}
transformed parameters {
  vector[n_alpha] alpha;
  matrix[n_rep, n_patient] epsilon;
  array[n_study] matrix[n_rep, n_rep] covariance_cholesky;
  {
    int index;
    int last_visit;
    int first_visit;
    int col_i;
    int row_i;
    int col_n;
    int row_n;
    vector[n_alpha + 1] alpha_latent;
    vector[n_delta + 1] delta_latent;
    vector[n_patient] x_beta_vector;
    vector[n_observe] means;
    matrix[n_rep, n_patient] means_matrix;
    vector[n_observe] y_imputed;
    matrix[n_rep, n_patient] y_matrix;
  
    // Non-centered and latent parameters
    if (model_type == 3) {
      alpha = mu[alpha_rep_index] + (tau[alpha_rep_index] .* alpha_raw);
    } else {
      alpha = alpha_raw;
    }
    alpha_latent[1] = 0;
    alpha_latent[2:(n_alpha + 1)] = alpha;
    delta_latent[1] = 0;
    delta_latent[2:(n_delta + 1)] = delta;
    
    // Multiply beta and x_beta in blocks.
    x_beta_vector = rep_vector(0, n_patient);
    if (n_beta > 0) {
      for (i in 1:n_study_x_beta) {
        col_i = x_beta_col_index[i];
        row_i = x_beta_row_index[i];
        col_n = x_beta_col_n[i];
        row_n = x_beta_row_n[i];
        x_beta_vector[col_i:(col_i + col_n - 1)] = (
          beta[row_i:(row_i + row_n - 1)]' *
          block(x_beta, row_i, col_i, row_n, col_n)
        )';
      }
    }
    
    // Mean response
    if (n_beta > 0) {
      means = alpha_latent[alpha_data_index] +
      delta_latent[delta_data_index] +
      x_beta_vector[index_patient];
    } else {
      means = alpha_latent[alpha_data_index] + delta_latent[delta_data_index];
    }
    
    // Cholesky factor of the covariance matrix for the current study.
    if (covariance_current == covariance_unstructured) {
      covariance_cholesky[n_study] = diag_pre_multiply(
        sigma[n_study],
        lambda_current[1]
      );
    } else if (covariance_current == covariance_ar1) {
      covariance_cholesky[n_study] = diag_pre_multiply(
        sigma[n_study],
        ar1_cholesky(rho_current[1], n_rep)
      );
    } else {
      covariance_cholesky[n_study] = diag_matrix(sigma[n_study]);
    }
    
    // Cholesky factor of the covariance matrix for historical studies.
    for (i in 1:(n_study - 1)) {
      if (covariance_historical == covariance_unstructured) {
        covariance_cholesky[i] = diag_pre_multiply(
          sigma[i],
          lambda_historical[i]
        );
      } else if (covariance_historical == covariance_ar1) {
        covariance_cholesky[i] = diag_pre_multiply(
          sigma[i],
          ar1_cholesky(rho_historical[i], n_rep)
        );
      } else {
        covariance_cholesky[i] = diag_matrix(sigma[i]);
      }
    }
    
    // Impute missing values with parameters.
    for (observation in 1:n_observe) {
      y_imputed[observation] = missing[observation] == 1 ?
      y_missing[count_missing[observation]] :
      y[observation];
    }
    
    // Compute matrices to vectorize the likelihood.
    index = 1;
    for (patient in 1:n_patient) {
      y_matrix[:,patient] = segment(y_imputed, index, n_rep);
      means_matrix[:,patient] = segment(means, index, n_rep);
      index += n_rep;
    }
    epsilon = y_matrix - means_matrix;
  }
}
model {
  // Likelihood
  int start;
  int end;
  int length;
  for (study in 1:n_study) {
    start = index_patient_study[study];
    length = n_patient_study[study];
    end = start + length - 1;
    target += sum(
      log_normal_kernel(
        covariance_cholesky[study],
        epsilon[:, start:end]
      )
    );
  }

  // Priors
  if (model_type == 3) {
    alpha_raw ~ std_normal();
    mu ~ normal(0, s_mu);
    tau ~ uniform(0, s_tau);
  } else {
    alpha_raw ~ normal(0, s_alpha);
  }
  delta ~ normal(0, s_delta);
  beta ~ normal(0, s_beta);
  for (i in 1:n_study) {
    sigma[i] ~ uniform(0, s_sigma);
  }
  if (covariance_current == covariance_unstructured) {
    lambda_current[1] ~ lkj_corr_cholesky(s_lambda);
  } else if (covariance_current == covariance_ar1) {
    rho_current ~ uniform(-1, 1);
  }
  if (covariance_historical == covariance_unstructured) {
    for (i in 1:n_lambda_historical) {
      lambda_historical[i] ~ lkj_corr_cholesky(s_lambda);
    }
  } else if (covariance_historical == covariance_ar1) {
    rho_historical ~ uniform(-1, 1);
  }
}
