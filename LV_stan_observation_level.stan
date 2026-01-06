data {
  int<lower=1> No; // Number of Observations
  int<lower=1> Nt; // Number of Traits (Now 3)
  matrix[No, Nt] Y; 
    
  int<lower=1> K;
  matrix[No, K] X;

  array[No] int animal; 
  int<lower=1> Na; 
  cov_matrix[Na] A; 
}

transformed data {
  matrix[Na, Na] LA = cholesky_decompose(A);
}

parameters {
  // Fixed effects (Matrix is K x 3 now)
  matrix[K, Nt] beta;    

  // Loadings
  // We identify the sign by forcing lambda1 > 0
  vector[Nt-1] lambda_free;
  real<lower=0> lambda1;
    
  // Unique Residual SDs (The "Uniquenesses")
  vector<lower=0>[Nt] sd_R; 

  // Variance Partitioning of the Latent Variable
  simplex[3] prop_var; 

  // Standardized random effects
  vector[Na] psi_a_std; 
  vector[Na] psi_e_std; 
  vector[No] psi_r_std; 
}

transformed parameters {
  vector[Nt] lambda;
  
  // Random Effect Vectors (scaled)
  vector[Na] psi_a;
  vector[Na] psi_e;
  vector[No] psi_r;
  
  // Total Latent Variable Score
  vector[No] psi_total;
  
  matrix[No, Nt] eta;
  matrix[No, Nt] fixed_part;

  // 1. Construct lambda vector
  lambda[1] = lambda1;
  for(t in 2:Nt) lambda[t] = lambda_free[t-1];

  // 2. Construct Latent Variable Components (Sum of vars = 1.0)
  psi_a = sqrt(prop_var[1]) * (LA * psi_a_std);
  psi_e = sqrt(prop_var[2]) * psi_e_std;
  psi_r = sqrt(prop_var[3]) * psi_r_std; 
    
  // 3. Construct Total Latent Score
  psi_total = psi_a[animal] + psi_e[animal] + psi_r;

  // 4. Linear Predictor
  fixed_part = X * beta;
  
  // Broadcast factor score across traits
  eta = fixed_part + psi_total * lambda';
}

model {
  // Priors
  to_vector(beta) ~ normal(0, 1);
  lambda1 ~ lognormal(0, 1);
  lambda_free ~ normal(0, 1);
  sd_R ~ normal(0, 10); 
  
  psi_a_std ~ normal(0, 1);
  psi_e_std ~ normal(0, 1);
  psi_r_std ~ normal(0, 1); 
    
  // Likelihood
  // Note: This assumes independent residuals (diagonal R matrix).
  // If you wanted correlated residuals beyond the factor, you'd need a full cov_matrix for R.
  // With 3 traits, a 1-factor model + diagonal residuals is standard.
  to_vector(Y) ~ normal(to_vector(eta), to_vector(rep_matrix(sd_R', No)));
}

generated quantities {
  // Recover variance components of the Latent Variable itself
  real var_LV_genetic = prop_var[1];
  real var_LV_perm_env = prop_var[2];
  real var_LV_temp_env = prop_var[3];
  
  vector[Nt] h2_traits; 
  vector[Nt] pe2_traits;
  
  for (t in 1:Nt) {
    real total_var = square(lambda[t]) * 1.0 + square(sd_R[t]);
    
    // Trait Heritability (via the Factor)
    h2_traits[t] = (square(lambda[t]) * var_LV_genetic) / total_var;
    
    // Trait PE (via the Factor)
    pe2_traits[t] = (square(lambda[t]) * var_LV_perm_env) / total_var;
  }
}
