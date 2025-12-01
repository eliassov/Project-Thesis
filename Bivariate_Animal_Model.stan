// data {
//   int<lower=1> No;  // Number of observations
//   matrix[No, 2] Y;  // Trait observations (2 columns for bivariate)
//   array[No] int animal;  // Individual identity for each observation
//   int<lower=1> Na;  // Number of individuals
//   cov_matrix[Na] A;  // Additive-genetic relationship matrix
// }


data {
  int<lower=1> No;      
  matrix[No, 2] Y;      
  
  // NEW: Fixed effects structure
  int<lower=1> K;       // Number of predictors (e.g. 2: Intercept + Sex)
  matrix[No, K] X;      // Design matrix
  
  array[No] int animal; 
  int<lower=1> Na;      
  cov_matrix[Na] A;     
}


transformed data {
  vector[2] mean_Y;  // Sample means of Y
  vector[2] sd_Y;    // Sample SDs of Y
  matrix[Na, Na] LA;  // Cholesky factor of A

  for (k in 1:2) {
    mean_Y[k] = mean(Y[, k]);
    sd_Y[k] = sd(Y[, k]);
  }

  LA = cholesky_decompose(A);
}

// parameters {
//   vector[2] mu;  // Overall intercepts
//   cholesky_factor_corr[2] L_A;  // Cholesky for genetic correlations
//   vector<lower=0>[2] sd_A;  // SDs of additive genetic effects
//   cholesky_factor_corr[2] L_E;  // Cholesky for permanent correlations
//   vector<lower=0>[2] sd_E;  // SDs of permanent effects
//   cholesky_factor_corr[2] L_R;  // Cholesky for residual correlations
//   vector<lower=0>[2] sd_R;  // SDs of residuals
//   matrix[2, Na] alpha_std;  // Precursors for additive genetic effects
//   matrix[2, Na] gamma_std;  // Standardised permanent effects
// }

parameters {
  // NEW: Matrix of fixed effects (K predictors x 2 traits)
  // beta[1,1] = Intercept Trait 1, beta[2,1] = Sex Slope Trait 1, etc.
  matrix[K, 2] beta;    

  // ... (Correlation/SD parameters remain unchanged) ...
  cholesky_factor_corr[2] L_A; 
  vector<lower=0>[2] sd_A; 
  cholesky_factor_corr[2] L_E; 
  vector<lower=0>[2] sd_E; 
  cholesky_factor_corr[2] L_R; 
  vector<lower=0>[2] sd_R; 
  matrix[2, Na] alpha_std; 
  matrix[2, Na] gamma_std; 
}




// 
// transformed parameters {
//   matrix[2, 2] Sigma_A;  // Genetic covariance matrix
//   matrix[2, 2] Sigma_E;  // Permanent covariance matrix
//   matrix[2, 2] Sigma_R;  // Residual covariance matrix
//   matrix[2, Na] alpha;  // Raw additive genetic effects
//   matrix[2, Na] gamma;  // Raw permanent effects
// 
//   Sigma_A = diag_pre_multiply(sd_A, L_A) * diag_pre_multiply(sd_A, L_A)';
//   Sigma_E = diag_pre_multiply(sd_E, L_E) * diag_pre_multiply(sd_E, L_E)';
//   Sigma_R = diag_pre_multiply(sd_R, L_R) * diag_pre_multiply(sd_R, L_R)';
// 
//   alpha = diag_pre_multiply(sd_A, L_A) * alpha_std * LA';
//   gamma = diag_pre_multiply(sd_E, L_E) * gamma_std;
// }



transformed parameters {
  matrix[2, 2] Sigma_A;  // Genetic covariance matrix
  matrix[2, 2] Sigma_E;  // Permanent covariance matrix
  matrix[2, 2] Sigma_R;  // Residual covariance matrix
  matrix[2, Na] alpha;  // Raw additive genetic effects
  matrix[2, Na] gamma;  // Raw permanent effects  matrix[2, 2] Sigma_A = diag_pre_multiply(sd_A, L_A) * diag_pre_multiply(sd_A, L_A)';

  Sigma_A = diag_pre_multiply(sd_A, L_A) * diag_pre_multiply(sd_A, L_A)';
  Sigma_E = diag_pre_multiply(sd_E, L_E) * diag_pre_multiply(sd_E, L_E)';
  Sigma_R = diag_pre_multiply(sd_R, L_R) * diag_pre_multiply(sd_R, L_R)';

  alpha = diag_pre_multiply(sd_A, L_A) * alpha_std * LA';
  gamma = diag_pre_multiply(sd_E, L_E) * gamma_std;
}



// 
// model {
//   vector[No] eta1;  // Expected for trait 1
//   vector[No] eta2;  // Expected for trait 2
// 
//   mu[1] ~ normal(76, 38); // wing length prior (based on intuition)
//   mu[2] ~ normal(11, 5.5); // beak length prior (also based on intuition)
// 
//   L_A ~ lkj_corr_cholesky(1);
//   // sd_A ~ exponential(1 ./ sd_Y);
//   sd_A ~ exponential(7);
//   
//   L_E ~ lkj_corr_cholesky(1);
//   // sd_E ~ exponential(1 ./ sd_Y);
//   sd_E ~ exponential(7);
//   
//   L_R ~ lkj_corr_cholesky(1);
//   // sd_R ~ exponential(1 ./ sd_Y);
//   sd_R ~ exponential(7);
// 
//   to_vector(alpha_std) ~ normal(0, 1);
//   to_vector(gamma_std) ~ normal(0, 1);
// 
//   for (o in 1:No) {
//     eta1[o] = mu[1] + alpha[1, animal[o]] + gamma[1, animal[o]];
//     eta2[o] = mu[2] + alpha[2, animal[o]] + gamma[2, animal[o]];
//   }
// 
//   for (o in 1:No) {
//     Y[o,] ~ multi_normal([eta1[o], eta2[o]], Sigma_R);
//   }
// }


model {
  matrix[No, 2] mu_fixed; // Holds the fixed effect part for every observation
  
  // Priors
  to_vector(beta) ~ normal(0, 100); // Weakly informative prior for all betas
  
  // ... (Other priors unchanged) ...
  L_A ~ lkj_corr_cholesky(1);
  sd_A ~ exponential(7);
  L_E ~ lkj_corr_cholesky(1);
  sd_E ~ exponential(7);
  L_R ~ lkj_corr_cholesky(1);
  sd_R ~ exponential(7);
  to_vector(alpha_std) ~ normal(0, 1);
  to_vector(gamma_std) ~ normal(0, 1);

  // Likelihood
  // Calculate fixed effects: X * beta -> (No x K) * (K x 2) = (No x 2)
  mu_fixed = X * beta;

  for (o in 1:No) {
    vector[2] mu_total;
    // Add fixed effects + animal random effect + perm env random effect
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
