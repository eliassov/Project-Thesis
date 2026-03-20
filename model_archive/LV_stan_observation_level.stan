data {
  int<lower=1> No; // Number of observations
  int<lower=1> Nt; // Number of traits 
  int<lower=1> Nlv; // Number of latent variables
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
  // vector[Nt-1] lambda_free;
  // real<lower=0> lambda1;
  
  vector[Nt-1] lambda_1_free;
  real<lower=0> lambda1_1;
  
  vector[Nt-2] lambda_2_free;
  real<lower=0> lambda2_2;

  vector<lower=0>[Nt] sd_R; 

  // Variance partitioning of the latent variable
  // simplex[3] prop_var; 
  
  simplex[3] prop_var_1; 

  simplex[3] prop_var_2;


  // Standardized random effects
  // vector[Na] psi_a_std; 
  // vector[Na] psi_e_std; 
  // vector[No] psi_r_std; 
  
  vector[Na] psi_1_a_std; 
  vector[Na] psi_1_e_std; 
  vector[No] psi_1_r_std; 
  
  vector[Na] psi_2_a_std;
  vector[Na] psi_2_e_std;
  vector[No] psi_2_r_std;
  
}

transformed parameters {
  // vector[Nt] lambda;
  
  matrix[Nt, Nlv] lambda;
  
  // Random effect vectors (scaled)
  // vector[Na] psi_a;
  // vector[Na] psi_e;
  // vector[No] psi_r;
  
  vector[Na] psi_1_a;
  vector[Na] psi_1_e;
  vector[No] psi_1_r;
  
  // vector[Na] psi_2_a;
  // vector[Na] psi_2_e;
  // vector[No] psi_2_r;
  
  // vector[No] psi_total;
  
  vector[No] psi_1_total;

  // vector[No] psi_2_total;

  
  matrix[No, Nt] eta;
  matrix[No, Nt] fixed_part;

  lambda[1,1] = lambda1_1;
  // lambda[2,2] = lambda2_2;
  
  for(t in 2:Nt) lambda[t,1] = lambda_1_free[t-1];
  
  // lambda[1,2] = lambda_2_free[1];
  // for(t in 3:Nt) lambda[t,2] = lambda_2_free[t-1];

  // psi_a = sqrt(prop_var[1]) * (LA * psi_a_std);
  // psi_e = sqrt(prop_var[2]) * psi_e_std;
  // psi_r = sqrt(prop_var[3]) * psi_r_std; 
  //   
  psi_1_a = sqrt(prop_var_1[1]) * (LA * psi_1_a_std);
  psi_1_e = sqrt(prop_var_1[2]) * psi_1_e_std;
  psi_1_r = sqrt(prop_var_1[3]) * psi_1_r_std; 
  // 
  // psi_2_a = sqrt(prop_var_2[1]) * (LA * psi_2_a_std);
  // psi_2_e = sqrt(prop_var_2[2]) * psi_2_e_std;
  // psi_2_r = sqrt(prop_var_2[3]) * psi_2_r_std; 
  
  
    
  // psi_total = psi_a[animal] + psi_e[animal] + psi_r;
  
  psi_1_total = psi_1_a[animal] + psi_1_e[animal] + psi_1_r;

  // psi_2_total = psi_2_a[animal] + psi_2_e[animal] + psi_2_r;
  
  


  fixed_part = X * beta;
  
  // eta = fixed_part + psi_total * lambda';
  
  eta = fixed_part + psi_1_total * lambda[, 1]' + psi_2_total * lambda[, 2]';
    
    
}

model {
  // Priors
  to_vector(beta) ~ normal(0, 1);
  sd_R ~ normal(0, 10);
  // lambda1 ~ lognormal(0, 1);
  // lambda_free ~ normal(0, 1);
  
  lambda1_1 ~ lognormal(0, 1);
  lambda_1_free ~ normal(0, 1);
  
  lambda2_2 ~ lognormal(0, 1);
  lambda_2_free ~ normal(0, 1);
   
  
  // psi_a_std ~ normal(0, 1);
  // psi_e_std ~ normal(0, 1);
  // psi_r_std ~ normal(0, 1); 
  
  psi_1_a_std ~ normal(0, 1);
  psi_1_e_std ~ normal(0, 1);
  psi_1_r_std ~ normal(0, 1); 
  
  psi_2_a_std ~ normal(0, 1);
  psi_2_e_std ~ normal(0, 1);
  psi_2_r_std ~ normal(0, 1); 
    

  to_vector(Y) ~ normal(to_vector(eta), to_vector(rep_matrix(sd_R', No)));
}

generated quantities {
  // real var_LV_genetic = prop_var[1];
  // real var_LV_perm_env = prop_var[2];
  // real var_LV_temp_env = prop_var[3];
  
  real var_LV_1_genetic = prop_var_1[1];
  real var_LV_1_perm_env = prop_var_1[2];
  real var_LV_1_temp_env = prop_var_1[3];
  
  real var_LV_2_genetic = prop_var_2[1];
  real var_LV_2_perm_env = prop_var_2[2];
  real var_LV_2_temp_env = prop_var_2[3];
  
  vector[Nt] h2_traits; 
  // vector[Nt] pe2_traits;
  
  for (t in 1:Nt) {
    // real total_var = square(lambda[t]) * 1.0 + square(sd_R[t]);
    
    real total_var = square(lambda[t,1]) * 1.0 + square(lambda[t,2]) * 1.0 + square(sd_R[t]);

    // Trait heritability (via the factor)
    h2_traits[t] = (square(lambda[t,1]) * var_LV_1_genetic + square(lambda[t,2]) * var_LV_2_genetic) / total_var;
    
    // Trait PE proportion (via the factor)
    // pe2_traits[t] = (square(lambda[t,1]) * var_LV_1_perm_env + square(lambda[t,2]) * var_LV_2_perm_env) / total_var;
  }
}
