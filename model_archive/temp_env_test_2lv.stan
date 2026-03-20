
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
  int N_factors = 2; 
}

parameters {
  matrix[K, Nt] beta;       
  vector<lower=0>[Nt] sd_R; 

  // --- FACTOR VARIANCE DECOMPOSITION ---
  // A simplex ensures these 3 numbers always sum to 1.0
  // Index 1 = Genetic (h2)
  // Index 2 = Permanent Env (pe2)
  // Index 3 = Temporary Env (te2)
  array[N_factors] simplex[3] var_decomp; 

  // --- LOADINGS ---
  real<lower=0> L1_anchor; 
  vector[Nt-1]  L1_free;

  real<lower=0> L2_anchor;
  vector[Nt-2]  L2_free;

  // --- RANDOM VECTORS ---
  // Animal Level (Linked to Pedigree)
  matrix[N_factors, Na] w_a; 
  // Animal Level (Permanent Environment - e.g. Nest effect/Maternal)
  matrix[N_factors, Na] w_pe; 
  // Observation Level (Temporary Environment - Day-to-day fluctuation)
  matrix[N_factors, No] w_te; 
}

transformed parameters {
  // --- 1. RECURSIVE LOOP (Breeding Values) ---
  matrix[N_factors, Na] u_a;
  
  {
    matrix[N_factors, Na+1] aug;
    aug[, Na+1] = rep_vector(0, N_factors); 

    for (i in 1:Na) {
      aug[, i] = 0.5 * aug[, dam[i]] + 0.5 * aug[, sire[i]] 
               + sqrt(dii[i]) * col(w_a, i);
    }
    u_a = aug[, 1:Na];
  }

  // --- 2. CONSTRUCT FACTORS (Unit Variance) ---
  // PSI is [No x 2]. It is the Latent Score for every single observation.
  // It is composed of 3 weighted layers.
  matrix[No, N_factors] PSI;

  for (i in 1:No) {
    for (k in 1:N_factors) {
      // The "Simplex" solution ensures the weights squared sum to 1.0
      PSI[i, k] = sqrt(var_decomp[k, 1]) * u_a[k, animal[i]]   // Genetic (Use Animal ID)
                + sqrt(var_decomp[k, 2]) * w_pe[k, animal[i]]  // PE (Use Animal ID)
                + sqrt(var_decomp[k, 3]) * w_te[k, i];         // TE (Use Observation ID)
    }
  }

  // --- 3. CONSTRUCT LOADING MATRIX ---
  matrix[Nt, N_factors] Lambda;
  
  // Factor 1 (Wing Anchor)
  Lambda[1, 1] = L1_anchor;     
  for(j in 2:Nt) Lambda[j, 1] = L1_free[j-1];

  // Factor 2 (Beak Anchor)
  Lambda[1, 2] = 0;             
  Lambda[2, 2] = L2_anchor;     
  for(j in 3:Nt) Lambda[j, 2] = L2_free[j-2];

  // --- 4. PREDICTOR ---
  matrix[No, Nt] mu;
  // PSI is [No x 2], Lambda' is [2 x 5] -> Result [No x 5]
  mu = X * beta + PSI * Lambda';
}

model {
  // --- PRIORS ---
  to_vector(beta) ~ normal(0, 5);
  sd_R ~ normal(0, 1);
  
  // Factor Priors (Weakly informative)
  // Simplexes define their own uniform prior over the valid space automatically
  
  L1_anchor ~ normal(1.0, 0.5); 
  L2_anchor ~ normal(1.0, 0.5);
  L1_free ~ normal(0, 1);
  L2_free ~ normal(0, 1);

  // Noise Priors (Standard Normal)
  to_vector(w_a) ~ normal(0, 1);
  to_vector(w_pe) ~ normal(0, 1);
  to_vector(w_te) ~ normal(0, 1);
  
  for (k in 1:N_factors) {
      var_decomp[k] ~ dirichlet(rep_vector(1.5, 3));
  }

  // --- LIKELIHOOD ---
  // Vectorized for speed
  to_vector(Y) ~ normal(to_vector(mu), to_vector(rep_matrix(sd_R', No)));
}

generated quantities {
  // --- DIAGNOSTICS & OUTPUTS ---
  
  // 1. Extract the Variance Components for the Latent Variable (0-1 scale)
  vector[N_factors] h2_psi; // Heritability of the Factor
  vector[N_factors] pe2_psi; // PE of the Factor
  vector[N_factors] te2_psi; // TE of the Factor

  for(k in 1:N_factors) {
     h2_psi[k]  = var_decomp[k, 1];
     pe2_psi[k] = var_decomp[k, 2];
     te2_psi[k] = var_decomp[k, 3];
  }

  // 2. Trait-Level Heritabilities
  // Now we account for all 3 layers of the factor + Residual
  vector[Nt] h2_trait_total;
  
  for(j in 1:Nt) {
      // Genetic Variance = Loading^2 * Factor_Heritability
      real var_G = square(Lambda[j, 1]) * h2_psi[1] 
                 + square(Lambda[j, 2]) * h2_psi[2];
      
      // Total Variance = Loading^2 * (1.0) + Residual^2
      // (Because Factor Variance is constrained to 1.0, Total Factor Var = Loading^2)
      real var_Tot = square(Lambda[j, 1]) 
                   + square(Lambda[j, 2]) 
                   + square(sd_R[j]);
      
      h2_trait_total[j] = var_G / var_Tot;
  }
}

