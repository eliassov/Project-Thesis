data {
  int<lower=1> No; 
  int<lower=1> Nt; 
  vector[Nt] Y[No]; // Array of vectors for multi_normal
  
  int<lower=1> K;
  matrix[No, K] X;

  int<lower=1> Na; 
  int animal[No]; 
  cov_matrix[Na] A; 
}

transformed data {
  matrix[Na, Na] LA = cholesky_decompose(A);
}

parameters {
  // Fixed effects
  matrix[K, Nt] beta;    

  // --- LATENT VARIABLE (Genetics + PE) ---
  // We keep this structure as requested
  vector[Nt-1] lambda_free;
  real<lower=0> lambda1;
  
  // Proportions of variance for the LV
  // prop_var[1] = Genetic, prop_var[2] = Perm Env
  simplex[2] prop_var; 

  vector[Na] psi_a_std; 
  vector[Na] psi_e_std;

  // --- RESIDUALS (The Fix) ---
  // Instead of an LV, we use a standard covariance matrix
  // This captures the "epsilon" correlation without the identification issues
  vector<lower=0>[Nt] sigma_resid; 
  cholesky_factor_corr[Nt] L_resid; 
}

transformed parameters {
  vector[Nt] lambda;
  vector[Na] psi_ind; 
  
  // Construct Loadings
  lambda[1] = lambda1;
  for(t in 2:Nt) lambda[t] = lambda_free[t-1];

  // Construct Individual LV (Gen + PE)
  // Note: We sum them here to save time in the loop
  psi_ind = (sqrt(prop_var[1]) * (LA * psi_a_std)) + 
            (sqrt(prop_var[2]) * psi_e_std);
}

model {
  // Priors
  to_vector(beta) ~ normal(0, 1);
  lambda1 ~ lognormal(0, 1);
  lambda_free ~ normal(0, 1);
  
  psi_a_std ~ normal(0, 1);
  psi_e_std ~ normal(0, 1);
  
  sigma_resid ~ normal(0, 3);
  L_resid ~ lkj_corr_cholesky(2.0); 

  // Likelihood
  // We pre-calculate the residual covariance Cholesky
  matrix[Nt, Nt] L_Sigma = diag_pre_multiply(sigma_resid, L_resid);

  vector[Nt] mu[No];
  for(i in 1:No){
      // The Mean is Fixed Effects + The Latent Variable Score * Loadings
      mu[i] = (X[i] * beta)' + (psi_ind[animal[i]] * lambda); 
  }
  
  // This effectively models: Y ~ MultiNormal(Mean, Residual_Covariance)
  Y ~ multi_normal_cholesky(mu, L_Sigma);
}

generated quantities {
  // --- OUTPUTS FOR YOUR SUPERVISOR ---
  
  // 1. Variances of the Latent Variable components
  real var_LV_genetic = prop_var[1];
  real var_LV_perm_env = prop_var[2];
  
  // 2. Residual Covariance Matrix (R)
  // This corresponds to the "epsilon" part, just unstructured
  matrix[Nt, Nt] R_mat = multiply_lower_tri_self_transpose(diag_pre_multiply(sigma_resid, L_resid));
  
  // 3. Heritability of traits (Standard formula)
  vector[Nt] h2_traits; 
  for (t in 1:Nt) {
    // Total Variance = Loading^2 * 1.0 + Residual_Variance_t
    real total_var = square(lambda[t]) + R_mat[t,t];
    h2_traits[t] = (square(lambda[t]) * var_LV_genetic) / total_var;
  }
}
