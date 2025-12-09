
data {
  int<lower=1> No;      
  matrix[No, 2] Y;      
  
  int<lower=1> K;       // Number of predictors (2 for intercept + sex)
  matrix[No, K] X;      // Design matrix
  
  array[No] int animal; 
  int<lower=1> Na;      
  cov_matrix[Na] A;     
}


transformed data {
  matrix[Na, Na] LA;  // Cholesky factor of A
  LA = cholesky_decompose(A);
}


parameters {
  // Matrix of fixed effects (K predictors x 2 traits)
  // beta[1,1] = intercept trait 1, beta[2,1] = sex slope trait 1, etc.
  matrix[K, 2] beta;    

  cholesky_factor_corr[2] L_A; 
  vector<lower=0>[2] sd_A; 
  cholesky_factor_corr[2] L_E; 
  vector<lower=0>[2] sd_E; 
  cholesky_factor_corr[2] L_R; 
  vector<lower=0>[2] sd_R; 
  matrix[2, Na] alpha_std; 
  matrix[2, Na] gamma_std; 
}




transformed parameters {
  matrix[2, 2] Sigma_A;  // Genetic covariance matrix
  matrix[2, 2] Sigma_E;  // Permanent covariance matrix
  matrix[2, 2] Sigma_R;  // Residual covariance matrix
  matrix[2, Na] alpha;  // Additive genetic effects (breeding values)
  matrix[2, Na] gamma;  // Raw permanent effects

  Sigma_A = diag_pre_multiply(sd_A, L_A) * diag_pre_multiply(sd_A, L_A)';
  Sigma_E = diag_pre_multiply(sd_E, L_E) * diag_pre_multiply(sd_E, L_E)';
  Sigma_R = diag_pre_multiply(sd_R, L_R) * diag_pre_multiply(sd_R, L_R)';

  alpha = diag_pre_multiply(sd_A, L_A) * alpha_std * LA';
  gamma = diag_pre_multiply(sd_E, L_E) * gamma_std;
}




model {
  // Calculate fixed effects: X * beta -> (No x K) * (K x 2) = (No x 2)
  matrix[No, 2] mu_fixed; // Holds the fixed effect part for every observation
  
  // Priors
  to_vector(beta) ~ normal(0, 100); // Weakly informative prior for all betas
  
  L_A ~ lkj_corr_cholesky(1);
  sd_A ~ normal(0, 10);
  L_E ~ lkj_corr_cholesky(1);
  sd_E ~ normal(0, 10);
  L_R ~ lkj_corr_cholesky(1);
  sd_R ~ normal(0, 10);
  to_vector(alpha_std) ~ normal(0, 1);
  to_vector(gamma_std) ~ normal(0, 1);

  mu_fixed = X * beta;

  for (o in 1:No) {
    vector[2] mu_total;
    
    // Fixed effects + animal random effect + perm env random effect
    mu_total[1] = mu_fixed[o, 1] + alpha[1, animal[o]] + gamma[1, animal[o]];
    mu_total[2] = mu_fixed[o, 2] + alpha[2, animal[o]] + gamma[2, animal[o]];
    
    Y[o] ~ multi_normal(mu_total, Sigma_R);
  }
}



generated quantities {
  corr_matrix[2] cor_A = multiply_lower_tri_self_transpose(L_A);
  corr_matrix[2] cor_E = multiply_lower_tri_self_transpose(L_E);
  corr_matrix[2] cor_R = multiply_lower_tri_self_transpose(L_R);
  real heritability1 = Sigma_A[1,1] / (Sigma_A[1,1] + Sigma_E[1,1] + Sigma_R[1,1]);
  real heritability2 = Sigma_A[2,2] / (Sigma_A[2,2] + Sigma_E[2,2] + Sigma_R[2,2]);
  
}
