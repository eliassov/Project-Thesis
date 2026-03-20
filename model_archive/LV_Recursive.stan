data {
  int<lower=1> No; // Number of observations
  int<lower=1> Nt; // Number of traits 
  int<lower=1> Nlv; // Number of latent variables
  matrix[No, Nt] Y; 
    
  int<lower=1> K;
  matrix[No, K] X;

  int<lower=1> Na; 
  array[No] int<lower=1, upper=Na> animal; 

  array[Na] int<lower=1, upper=Na+1> dam; // Parent IDs (Na + 1 means missing)
  array[Na] int<lower=1, upper=Na+1> sire;
  vector[Na] dii; // Mendelian sampling variance scaling factors


  // cov_matrix[Na] A;  
  
}

// transformed data {
//   matrix[Na, Na] LA = cholesky_decompose(A);
// }

parameters {
  matrix[K, Nt] beta;    

  // Loadings (Lambda)
  
  // We identify the sign by forcing lambda1 > 0
  // vector[Nt-1] lambda_free;
  // real<lower=0> lambda1;
  
  vector[Nt] lambda_1_free;
  // real<lower=0> lambda1_1;
  
  vector[Nt-1] lambda_2_free;
  // real<lower=0> lambda2_2;
    
  vector<lower=0>[Nt] sd_R; 

  // Variance partitioning of the latent variable
  // simplex[3] prop_var; 
  
  // simplex[3] prop_var_1; 

  // simplex[3] prop_var_2;
  
  real<lower=0> sd_1_a;
  real<lower=0> sd_1_pe;
  real<lower=0> sd_1_te;
  
  real<lower=0, upper=1> rho;

  
  real<lower=0> sd_2_pe;
  real<lower=0> sd_2_te;

  // Standardized random effects
  // vector[Na] psi_a_std; 
  // vector[Na] psi_e_std; 
  // vector[No] psi_r_std; 
  
  // vector[Na] psi_1_a_std; 
  vector[Na] w_1_a;
  vector[Na] w_1_pe; 
  vector[No] w_1_te; 
  
  // vector[Na] psi_2_a_std; 
  vector[Na] w_2_a;
  vector[Na] w_2_pe; 
  vector[No] w_2_te; 
  
}

transformed parameters {
  // vector[Nt] lambda;
  
  matrix[Nt, Nlv] lambda;
  
  // Random effect vectors (scaled)
  // vector[Na] psi_a;
  // vector[Na] psi_e;
  // vector[No] psi_r;
  
  vector[Na+1] psi_1_a_aug; // Vector to handle missing parents (index Na + 1)
  vector[Na] psi_1_a;
  vector[Na] psi_1_pe;
  vector[No] psi_1_te;
  
  vector[Na+1] psi_2_a_aug;
  vector[Na] psi_2_a;
  vector[Na] psi_2_pe;
  vector[No] psi_2_te;
  
  // vector[No] psi_total;
  
  vector[No] psi_1_total;

  vector[No] psi_2_total;

  
  matrix[No, Nt] eta;
  matrix[No, Nt] fixed_part;
  
  real mean_lv_1; 
  real mean_lv_2;
  real sd_lv_1;
  real sd_lv_2;
  
  // vector[No] mean_lv_1_vec;
  // vector[No] mean_lv_2_vec;
  // vector[No] sd_lv_1_vec;
  // vector[No] sd_lv_2_vec;
  
  
  // real<lower=0> var_G_1 = prop_var_1[1];
  real<lower=0> var_G_1 = sd_1_a^2; // ????
  
  real<lower=0> var_G_2 = var_G_1 * rho;


  // lambda[1,1] = lambda1_1;
  for(t in 1:Nt) lambda[t,1] = lambda_1_free[t];
  
  lambda[1,2] = 0; // Prevents rotational freedom // THIS NEEDS TO BE UPDATED AS WELL.
  
  // lambda[2,2] = lambda2_2;
  
  for(t in 1:Nt) lambda[t,2] = lambda_2_free[t-1];
  
  
  // real sd_G_1 = sqrt(prop_var_1[1]);
  // real sd_G_2 = sqrt(prop_var_2[1]);
  real sd_G_1 = sqrt(var_G_1);
  real sd_G_2 = sqrt(var_G_2);
  
  
  psi_1_a_aug[Na+1] = 0;
  psi_2_a_aug[Na+1] = 0;
  
  for (i in 1:Na) {
    real parent_contr_1 = 0.5 * psi_1_a_aug[dam[i]] + 0.5 * psi_1_a_aug[sire[i]];
    real parent_contr_2 = 0.5 * psi_2_a_aug[dam[i]] + 0.5 * psi_2_a_aug[sire[i]];
    
    real mendelian_1 = sd_G_1 * sqrt(dii[i]) * w_1_a[i];
    real mendelian_2 = sd_G_2 * sqrt(dii[i]) * w_2_a[i];
    
    psi_1_a_aug[i] = parent_contr_1 + mendelian_1;
    psi_2_a_aug[i] = parent_contr_2 + mendelian_2;
  }
  
  // Extract actual animals (1 to Na)
  psi_1_a = psi_1_a_aug[1:Na];
  psi_2_a = psi_2_a_aug[1:Na];

  
  // psi_a = sqrt(prop_var[1]) * (LA * psi_a_std);
  // psi_e = sqrt(prop_var[2]) * psi_e_std;
  // psi_r = sqrt(prop_var[3]) * psi_r_std; 
  //   
  // psi_1_a = sqrt(prop_var_1[1]) * (LA * psi_1_a_std);
  // psi_1_pe = sqrt(prop_var_1[2]) * w_1_pe;
  psi_1_pe = sd_1_pe * w_1_pe;
  psi_1_te = sd_1_te * w_1_te;
  // psi_1_te = sqrt(prop_var_1[3]) * w_1_te; 
  
  // psi_2_a = sqrt(prop_var_2[1]) * (LA * psi_2_a_std);
  psi_2_pe = sd_2_pe * w_2_pe;  
  psi_2_te = sd_2_te * w_2_te;  
  
  
    
  // psi_total = psi_a[animal] + psi_e[animal] + psi_r;
  
  psi_1_total = psi_1_a[animal] + psi_1_pe[animal] + psi_1_te;

  psi_2_total = psi_2_a[animal] + psi_2_pe[animal] + psi_2_te;
  
  
  // mean_lv_1 = mean(psi_1_total);
  mean_lv_1 = mean(psi_1_a)
  // mean_lv_1_vec = rep_vector(mean_lv_1, No);
  
  // mean_lv_2 = mean(psi_2_total);
  // mean_lv_2_vec = rep_vector(mean_lv_2, No);
  
  sd_lv_1 = sd(psi_1_total);
  // sd_lv_1_vec = rep_vector(sd_lv_1, No);
  
  // si_lv_2 = sd(psi_2_total);
  // sd_lv_2_vec = rep_vector(sd_lv_2, No);
  
  
  
  psi_1_gen_standard = (psi_1_a - mean_lv_1) / sd_lv_1;
  // psi_2_standard = (psi_2_total - mean_lv_2_vec) / sd_lv_2_vec;



  fixed_part = X * beta;
  
  // eta = fixed_part + psi_total * lambda';
  
  // Outer product: latent variable * loadings
  // eta = fixed_part + psi_1_total * lambda[, 1]' + psi_2_total * lambda[, 2]';
  
  eta = fixed_part + psi_1_standard * lambda[, 1]' + psi_2_standard * lambda[, 2]';

    
    
}

model {
  // Priors
  to_vector(beta) ~ normal(0, 1);
  sd_R ~ normal(0, 1);
  // lambda1 ~ lognormal(0, 1);
  // lambda_free ~ normal(0, 1);
  // 
  lambda1_1 ~ lognormal(0, 1);
  // lambda_1_free ~ normal(0, 1);
  
  lambda1_1 ~ normal(0, 1); 
  lambda2_2 ~ normal(0, 1);
  
  // lambda2_2 ~ lognormal(0, 1);
  lambda_2_free ~ normal(0, 1);
   
  w_1_a ~ normal(0, 1);  
  w_2_a ~ normal(0, 1); 
  
  // NEED TO ADD GAMMA PRIOR FOR THE VARIANCE?
  prop_var_1 ~ dirichlet(rep_vector(1.5, 3)); // To avoid getting stuck in zero, as we expect the variance components to be positive at all levels
  // prop_var_2 ~ dirichlet(rep_vector(1.5, 3));
  
  sd_2_pe ~ normal(0,1);
  sd_2_te ~ normal(0,1);
  
  rho ~ beta(1,1); 
  // psi_a_std ~ normal(0, 1);
  // psi_e_std ~ normal(0, 1);
  // psi_r_std ~ normal(0, 1); 
  
  // psi_1_a_std ~ normal(0, 1);
  w_1_pe ~ normal(0, 1);
  w_1_te ~ normal(0, 1); 
  
  // psi_2_a_std ~ normal(0, 1);
  w_2_pe ~ normal(0, 1);
  w_2_te ~ normal(0, 1); 
    

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
  vector[Nt] pe2_traits;
  
  for (t in 1:Nt) {
    // real total_var = square(lambda[t]) * 1.0 + square(sd_R[t]);
    
    real total_var = square(lambda[t,1]) * 1.0 + square(lambda[t,2]) * 1.0 + square(sd_R[t]);

    // Trait heritability (via the factor)
    h2_traits[t] = (square(lambda[t,1]) * var_LV_1_genetic + square(lambda[t,2]) * var_LV_2_genetic) / total_var;
    
    // Trait PE proportion (via the factor)
    pe2_traits[t] = (square(lambda[t,1]) * var_LV_1_perm_env + square(lambda[t,2]) * var_LV_2_perm_env) / total_var;
  }
}

