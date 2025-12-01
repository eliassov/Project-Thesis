



data {
  int<lower=1> No; 
  int<lower=1> Nt; 
  matrix[No, Nt] Y; 
  
  // NEW
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
  // NEW: Matrix of fixed effects (K predictors x Nt traits)
  // Replaces 'vector[Nt] mu'
  matrix[K, Nt] beta;  

  // ... (Lambda and SD parameters unchanged) ...
  vector[Nt-1] lambda_free;
  real<lower=0> lambda1;
  vector<lower=0>[Nt] sd_R; 
  real<lower=0> sd_psi_a; 
  real<lower=0> sd_psi_e; 
  vector[Na] psi_a_std; 
  vector[Na] psi_e_std; 
}




transformed parameters {
  // Construct the Lambda Vector
  vector[Nt] lambda;
  lambda[1] = lambda1;
  lambda[2:Nt] = lambda_free;

  // CONSTRUCT THE LATENT VARIABLE
  // Instead of "psi_std = psi / sd(psi)", we build it to have variance 1 automatically.
  // Genetic part: scales with sqrt(h2)
  // Environmental part: scales with sqrt(1 - h2)
  
  vector[Na] psi_a = sqrt(h2_psi) * (LA * psi_a_std); 
  vector[Na] psi_e = sqrt(1 - h2_psi) * psi_e_std;
  
  // Total Psi. Variance = h2 + (1-h2) = 1. 
  vector[Na] psi = psi_a + psi_e;

  matrix[No, Nt] eta; 
  
  // Calculate fixed part: (No x K) * (K x Nt) = (No x Nt)
  matrix[No, Nt] fixed_part = X * beta;

  for (o in 1:No) {
    for (t in 1:Nt) {
      // Fixed effect + (Loading * Latent Factor)
      eta[o,t] = fixed_part[o,t] + lambda[t] * psi[animal[o]];
    }
  }
}

model {
  // Priors
  to_vector(beta) ~ normal(0, 1); // Standardized scale, so N(0,1) is fine
  
  // ... (Other priors unchanged) ...
  lambda ~ normal(0, 1);
  sd_R ~ exponential(1);
  sd_psi_a ~ exponential(1);
  sd_psi_e ~ exponential(1);
  psi_a_std ~ normal(0, 1);
  psi_e_std ~ normal(0, 1);
  
  // Likelihood
  for (o in 1:No) {
    Y[o] ~ normal(eta[o], sd_R);
  }
}






generated quantities {
  real var_psi_a = square(sd_psi_a);  // LV genetic variance
  real var_psi_e = square(sd_psi_e);  // LV residual variance
  real h2_psi = var_psi_a / (var_psi_a + var_psi_e);  // LV heritability
  vector[Nt] h2_traits;  // Trait-specific heritabilities
  for (t in 1:Nt) {
    h2_traits[t] = lambda[t]^2 * var_psi_a / (lambda[t]^2 * (var_psi_a + var_psi_e) + sd_R[t]^2);
  }
  vector[Nt] h2_traits_no_residual;
  for (t in 1:Nt) {
  real denom = lambda[t]^2 * (var_psi_a + var_psi_e);
  h2_traits_no_residual[t] = (lambda[t]^2 * var_psi_a) / denom;
}
}



// data {
//   int<lower=1> No;  // Number of observations
//   int<lower=1> Nt;  // Number of traits
//   matrix[No, Nt] Y; // Trait observations
//   array[No] int animal;  // Individual identity
//   int<lower=1> Na;  // Number of individuals
//   cov_matrix[Na] A; // Additive-genetic relationship matrix
// }
// transformed data {
//   matrix[Na, Na] LA = cholesky_decompose(A);
// }
// parameters {
//   vector[Nt] mu;  // Trait means
//   
//   // vector<lower=0>[Nt-1] lambda_free;  // Free loadings (lambda[2:Nt])
//   
//   // vector[Nt] lambda;
//   
//   // Restricting one component to be positive and one is free (to avoid label-switching)
//   vector[Nt-1] lambda_free;
//   real<lower=0> lambda1;
//   
//   
//   vector<lower=0>[Nt] sd_R;  // Residual SDs
//   real<lower=0> sd_psi_a;  // LV genetic SD
//   real<lower=0> sd_psi_e;  // LV environmental SD
//   
//   #real<lower=0> sd_psi_r;  // LV residual SD
// 
//   vector[Na] psi_a_std;  // Standardized LV genetic effects
//   vector[Na] psi_e_std;  // Standardized LV environmental effects
//   
//   #vector[Na] psi_r_std;  // Standardized LV residual effects
// 
// }
// transformed parameters {
//   vector[Nt] lambda;
//   lambda[1] = lambda1;
//   for(t in 2:Nt) lambda[t] = lambda_free[t-1]
//   
//   vector[Na] psi_a = sd_psi_a * (LA * psi_a_std);  // LV genetic effects
//   vector[Na] psi_e = sd_psi_e * psi_e_std;  // LV environmental effects
//   
//   vector[Na] psi = psi_a + psi_e;
//   real psi_mean = mean(psi);
//   real psi_sd   = sd(psi);
//   vector[Na] psi_std = (psi - psi_mean) / psi_sd;
//   
//   
//   #vector[Na] psi_r = sd_psi_r * psi_r_std;  // LV environmental effects
//   
//   matrix[No, Nt] eta;  // Expected values
//   
//   // Construct full lambda vector with first element fixed to 1
//   // lambda[1] = 1.0;
//   // lambda[2:Nt] = lambda_free;
//   
//   // Compute expected values
//   for (o in 1:No) {
//     for (t in 1:Nt) {
//       // eta[o,t] = mu[t] + lambda[t] * (psi_a[animal[o]] + psi_e[animal[o]]);
//       eta[o,t] = mu[t] + lambda[t] * psi_std[animal[o]];
//     }
//   }
// }
// model {
//   // Priors
//   mu ~ normal(0, 1);
//   // lambda_free ~ normal(0, 1);
//   lambda ~ normal(0, 1);
//   sd_R ~ exponential(1);
//   sd_psi_a ~ exponential(1);
//   sd_psi_e ~ exponential(1);
//   psi_a_std ~ normal(0, 1);
//   psi_e_std ~ normal(0, 1);
//   
//   // Likelihood
//   for (o in 1:No) {
//     Y[o] ~ normal(eta[o], sd_R);
//   }
// }





// transformed parameters {
//   vector[Nt] lambda;
//   lambda[1] = lambda1;
//   for(t in 2:Nt) lambda[t] = lambda_free[t-1];
//   
//   vector[Na] psi_a = sd_psi_a * (LA * psi_a_std); 
//   vector[Na] psi_e = sd_psi_e * psi_e_std; 
//   
//   // Note: Standardizing psi inside the model is tricky for inference. 
//   // Usually better to work with raw psi_a + psi_e directly.
//   vector[Na] psi = psi_a + psi_e;

