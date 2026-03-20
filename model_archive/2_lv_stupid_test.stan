data {
  int<lower=1> No; // Number of observations
  int<lower=1> Nt; // Number of traits (5)
  matrix[No, Nt] Y; 
  
  int<lower=1> K;
  matrix[No, K] X;

  // Pedigree
  int<lower=1> Na; 
  array[No] int<lower=1, upper=Na> animal; 
  array[Na] int<lower=1, upper=Na+1> dam; 
  array[Na] int<lower=1, upper=Na+1> sire;
  vector[Na] dii; 
}

transformed data {
  int N_factors = 2; // Hard-coded for this model
}

parameters {
  matrix[K, Nt] beta;      
  vector<lower=0>[Nt] sd_R; 

  // --- FACTOR PARAMETERS ---
  // Heritability for EACH factor (Vector of length 2)
  vector<lower=0, upper=1>[N_factors] h2_psi;
  // vector<lower=0, upper=1>[N_factors] prop_var; 


  // --- LOADINGS ---
  // Factor 1: Wing (Anchor), others free
  real<lower=0> L1_anchor; 
  vector[Nt-1]  L1_free;

  // Factor 2: Wing (FIXED 0), Beak (Anchor), others free
  real<lower=0> L2_anchor;
  vector[Nt-2]  L2_free;

  // --- RANDOM VECTORS (Matrix form) ---
  // [2 x Na] matrices for vectorized speed
  matrix[N_factors, Na] w_a; 
  matrix[N_factors, Na] w_pe; 
  
  // matrix[N_factors, No] w_te; 

}

transformed parameters {
  // --- 1. RECURSIVE LOOP (Vectorized for 2 factors) ---
  matrix[N_factors, Na] u_a;
  
  {
    matrix[N_factors, Na+1] aug;
    aug[, Na+1] = rep_vector(0, N_factors); // Ghost parent = 0

    for (i in 1:Na) {
      // Vectorized: Updates both factors for Animal 'i' at once
      aug[, i] = 0.5 * aug[, dam[i]] + 0.5 * aug[, sire[i]] 
                 + sqrt(dii[i]) * col(w_a, i);
    }
    u_a = aug[, 1:Na];
  }

  // --- 2. CONSTRUCT FACTORS (Unit Variance) ---
  // We build a [Na x 2] matrix of "True" Latent Values
  
  matrix[Na, N_factors] PSI;
  // matrix[No, N_factors] PSI;

  
  for (k in 1:N_factors) {
     // Row k of u_a is Factor k for all animals
     // Transpose to get [Na] vector
     PSI[, k] = sqrt(h2_psi[k]) * u_a[k, ]' + sqrt(1 - h2_psi[k]) * w_pe[k, ]';
  }

  // --- 3. CONSTRUCT LOADING MATRIX (The Identification) ---
  matrix[Nt, N_factors] Lambda;
  
  // -- Column 1: "General Size" (Anchored on Wing) --
  Lambda[1, 1] = L1_anchor;      // Trait 1 (Wing)
  for(j in 2:Nt) {
    Lambda[j, 1] = L1_free[j-1]; // Traits 2-5
  }

  // -- Column 2: "Shape" (Anchored on Beak) --
  Lambda[1, 2] = 0;              // CONSTRAINT: Factor 2 cannot affect Wing
  Lambda[2, 2] = L2_anchor;      // Trait 2 (Beak)
  for(j in 3:Nt) {
    Lambda[j, 2] = L2_free[j-2]; // Traits 3-5
  }

  // --- 4. PREDICTOR ---
  matrix[No, Nt] mu;
  
  // Matrix Math: PSI is [No x 2], Lambda' is [2 x Nt] -> Result [No x Nt]
  mu = X * beta + PSI[animal, ] * Lambda';
}

model {
  // --- PRIORS ---
  to_vector(beta) ~ normal(0, 5);
  sd_R ~ normal(0, 1);
  
  // Factor Priors
  // Note: Anchors must be positive to fix sign
  L1_anchor ~ normal(1.0, 0.5); 
  L2_anchor ~ normal(1.0, 0.5);
  
  L1_free ~ normal(0, 1);
  L2_free ~ normal(0, 1);

  // Noise Priors
  to_vector(w_a) ~ normal(0, 1);
  to_vector(w_pe) ~ normal(0, 1);

  // --- LIKELIHOOD ---
  to_vector(Y) ~ normal(to_vector(mu), to_vector(rep_matrix(sd_R', No)));
}

generated quantities {
  // Diagnostics
  vector[N_factors] diag_var_genetic;
  vector[N_factors] diag_var_total;
  
  for(k in 1:N_factors) {
    diag_var_genetic[k] = variance(u_a[k, ]');
    diag_var_total[k]   = variance(PSI[, k]);
  }
  
  // Calculate Heritabilities for Traits
  vector[Nt] h2_trait_total;
  
  for(j in 1:Nt) {
     real var_genetic = square(Lambda[j, 1]) * h2_psi[1] + square(Lambda[j, 2]) * h2_psi[2];
     real var_total   = square(Lambda[j, 1]) + square(Lambda[j, 2]) + square(sd_R[j]);
     h2_trait_total[j] = var_genetic / var_total;
  }
}

