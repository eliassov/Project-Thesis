data {
  int<lower=1> No; // Number of observations
  int<lower=1> Nt; // Number of traits 
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
  matrix[K, Nt] beta;    

  // Loadings
  // We identify the sign by forcing lambda1 > 0
  vector[Nt-1] lambda_free;
  real<lower=0> lambda1;
    
  vector<lower=0>[Nt] sd_R; 

  // Variance partitioning of the latent variable
  simplex[3] prop_var; 

  // Standardized random effects
  vector[Na] psi_a_std; 
  vector[Na] psi_e_std; 
  vector[No] psi_r_std; 
}

transformed parameters {
  vector[Nt] lambda;
  
  // Random effect vectors (scaled)
  vector[Na] psi_a;
  vector[Na] psi_e;
  vector[No] psi_r;
  
  vector[No] psi_total;
  
  matrix[No, Nt] eta;
  matrix[No, Nt] fixed_part;

  lambda[1] = lambda1;
  for(t in 2:Nt) lambda[t] = lambda_free[t-1];

  psi_a = sqrt(prop_var[1]) * (LA * psi_a_std);
  psi_e = sqrt(prop_var[2]) * psi_e_std;
  psi_r = sqrt(prop_var[3]) * psi_r_std; 
    
  psi_total = psi_a[animal] + psi_e[animal] + psi_r;

  fixed_part = X * beta;
  
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
    

  to_vector(Y) ~ normal(to_vector(eta), to_vector(rep_matrix(sd_R', No)));
}

generated quantities {
  real var_LV_genetic = prop_var[1];
  real var_LV_perm_env = prop_var[2];
  real var_LV_temp_env = prop_var[3];
  
  vector[Nt] h2_traits; 
  vector[Nt] pe2_traits;
  
  for (t in 1:Nt) {
    real total_var = square(lambda[t]) * 1.0 + square(sd_R[t]);
    
    // Trait heritability (via the factor)
    h2_traits[t] = (square(lambda[t]) * var_LV_genetic) / total_var;
    
    // Trait PE proportion (via the factor)
    pe2_traits[t] = (square(lambda[t]) * var_LV_perm_env) / total_var;
  }
}
