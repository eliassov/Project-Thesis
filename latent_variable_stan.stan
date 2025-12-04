

data {
  int<lower=1> No; 
  int<lower=1> Nt; 
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
  // Fixed effects
  matrix[K, Nt] beta;  

  // Loadings (All estimated, first constrained positive)
  vector[Nt-1] lambda_free;
  real<lower=0> lambda1;
  
  // Residual SDs for traits
  vector<lower=0>[Nt] sd_R; 

  // We estimate h2 directly. This constrains total var(psi) = 1
  real<lower=0, upper=1> h2_psi; 

  // Standardized random effects
  vector[Na] psi_a_std; 
  vector[Na] psi_e_std; 
}

transformed parameters {
  vector[Nt] lambda;
  vector[Na] psi_a;
  vector[Na] psi_e;
  vector[Na] psi;
  matrix[No, Nt] eta;
  matrix[No, Nt] fixed_part;

  // Construct lambda
  lambda[1] = lambda1;
  for(t in 2:Nt) lambda[t] = lambda_free[t-1];

  // Construct latent variable (Using h2_psi)
  // Genetic part scales with sqrt(h2)
  // Environmental part scales with sqrt(1-h2)
  psi_a = sqrt(h2_psi) * (LA * psi_a_std); 
  psi_e = sqrt(1 - h2_psi) * psi_e_std; 
  
  // Total psi has variance = h2 + (1-h2) = 1
  psi = psi_a + psi_e;

  // Compute expected values
  fixed_part = X * beta;

  for (o in 1:No) {
    for (t in 1:Nt) {
      eta[o,t] = fixed_part[o,t] + lambda[t] * psi[animal[o]];
    }
  }
}

model {
  // Priors
  to_vector(beta) ~ normal(0, 1);
  lambda1 ~ normal(0, 1);
  lambda_free ~ normal(0, 1);
  sd_R ~ exponential(1);
  
  // h2_psi has an implicit uniform(0,1) prior

  psi_a_std ~ normal(0, 1);
  psi_e_std ~ normal(0, 1);
  
  for (o in 1:No) {
    Y[o] ~ normal(eta[o], sd_R);
  }
}

generated quantities {
  real var_psi_a = h2_psi;
  real var_psi_e = 1 - h2_psi;
  
  real sd_psi_a = sqrt(var_psi_a);
  real sd_psi_e = sqrt(var_psi_e);

  vector[Nt] h2_traits; 
  // vector[Nt] h2_traits_no_residual; # Is just the heritability of LV
  
  for (t in 1:Nt) {
    // V_A(trait) = lambda^2 * V_A(psi)
    // V_P(trait) = lambda^2 * V_P(psi) + V_R(trait)
    // Since V_P(psi) = 1:
    real va_trait = square(lambda[t]) * h2_psi;
    real vp_trait = square(lambda[t]) * 1.0 + square(sd_R[t]);
    
    h2_traits[t] = va_trait / vp_trait;
    
    // Heritability ignoring residual (ratio of individual variance that is genetic) 
    // (but this is just the heritability for the LV regardless of the trait)
    // h2_traits_no_residual[t] = va_trait / square(lambda[t]);
  }
}

