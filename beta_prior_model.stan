
data {
  int<lower=1> No; // Number of observations
  int<lower=1> Nt; // Number of traits 
  int<lower=1> Nlv; // Number of latent variables (2)
  matrix[No, Nt] Y; 
    
  int<lower=1> K;
  matrix[No, K] X;

  int<lower=1> Na; 
  array[No] int<lower=1, upper=Na> animal; 

  array[Na] int<lower=1, upper=Na+1> dam; 
  array[Na] int<lower=1, upper=Na+1> sire;
  vector[Na] dii; // Mendelian sampling variance scaling factors
}

parameters {
  // --- FIXED EFFECTS ---
  matrix[K, Nt] beta;    
  vector<lower=0>[Nt] sd_R; 

  // --- FACTOR 1 PARAMETERS (The "Big" Factor) ---
  // We estimate components directly (No simplex!)
  real<lower=0> sd_G1; 
  real<lower=0> sd_PE1;
  real<lower=0> sd_TE1;

  // --- FACTOR 2 PARAMETERS (The "Small" Factor) ---
  // Bob's Ordering: Genetic variance is a fraction of Factor 1
  real<lower=0, upper=1> rho; 
  real<lower=0> sd_G2;
  real<lower=0> sd_PE2;
  real<lower=0> sd_TE2;

  // --- LOADINGS (Lambdas) ---
  // Factor 1: Trait 1 is fixed to 1.0, others are free
  vector[Nt] lambda_1;
  
  // Factor 2: Trait 1 fixed to 0.0, Trait 2 fixed to 1.0
  // So we need free loadings for Trait 2 (on Fac 1) and Traits 3,4,5... (on Fac 2)
  vector[Nt] lambda_2;

  // --- RAW RANDOM VECTORS (Unit Normal Noise) ---
  // These drive the recursive loop and environmental noise
  vector[Na] w_1_a; 
  vector[Na] w_2_a;
  
  vector[Na] w_1_pe; 
  vector[Na] w_2_pe;
  
  vector[No] w_1_te; 
  vector[No] w_2_te; 
}

transformed parameters {
  // --- 1. RECURSIVE LOOP (Generate Raw Pedigree Shape) ---
  vector[Na] psi_1_a;
  vector[Na] psi_2_a;
  
  {
    // Augment vectors to handle missing parents (Index Na+1 = 0)
    vector[Na+1] aug_1;
    vector[Na+1] aug_2;
    aug_1[Na+1] = 0; 
    aug_2[Na+1] = 0;

    for (i in 1:Na) {
      // Factor 1: Just accumulate structure (Assume scale = 1.0 for now)
      real par1 = 0.5 * aug_1[dam[i]] + 0.5 * aug_1[sire[i]];
      aug_1[i] = par1 + sqrt(dii[i]) * w_1_a[i]; 

      // Factor 2
      real par2 = 0.5 * aug_2[dam[i]] + 0.5 * aug_2[sire[i]];
      aug_2[i] = par2 + sqrt(dii[i]) * w_2_a[i];
    }
    // Extract the actual animals
    psi_1_a = aug_1[1:Na];
    psi_2_a = aug_2[1:Na];
  }

  // --- 2. BOB'S STANDARDIZATION (The Z-Star) ---
  // We force these vectors to be exactly N(0,1), removing A-matrix inflation


  // --- 3. SCALE AND COMBINE ---
  vector[No] psi_1_total;
  vector[No] psi_2_total;
  
  psi_1_total = (sd_G1  * psi_1_a[animal]) 
              + (sd_PE1 * w_1_pe[animal]) 
              + (sd_TE1 * w_1_te);

  psi_2_total = (sd_G2  * psi_2_a[animal]) 
              + (sd_PE2 * w_2_pe[animal]) 
              + (sd_TE2 * w_2_te);
  
  vector[Na] psi_1_std = (psi_1_total - mean(psi_1_total)) / sd(psi_1_total);
  vector[Na] psi_2_std = (psi_2_a_raw - mean(psi_2_total)) / sd(psi_2_total);
  
  // Define Genetic SD for Factor 2 based on Rho (Bob's Prior)
  // Construct the Total Latent Variable for each observation
  // Formula: (Scale_G * Gen) + (Scale_PE * PE) + (Scale_TE * TE)
  // Note: PE and TE are already N(0,1) from parameters, so we just multiply.
  real var_1 = square(sd_G1) + square(sd_PE1) + square(sd_TE1);
  
  real tot_var_2 = rho * var_1;
  
  real var_2 = rho; // As the variance of the standardized first lv is 1??
  
  simplex[3] prop_2; 
  
  sd_G2 = sqrt(prop_2[1] * var_2);
  sd_PE2 = sqrt(prop_2[2] * var_2);
  sd_TE2 = sqrt(prop_2[3] * var_2);
 
 

  // --- 4. CONSTRUCT LAMBDA MATRIX ---
  matrix[Nt, 2] lambda;
  
  // -- FACTOR 1 --
  

  // -- FACTOR 2 --
  

  // --- 5. LINEAR PREDICTOR ---
  matrix[No, Nt] mu;
  matrix[No, Nt] fixed = X * beta;
  
  // Standard Factor Model Equation: Y = Fixed + (Factor * Loadings)
  mu = fixed + psi_1_std * lambda[,1]' + psi_2_std * lambda[,2]';
}

model {
  // --- PRIORS ---
  
  // 1. Variances (Half-Normal on SDs is robust and standard)
  sd_G1 ~ normal(0, 1);
  sd_PE1 ~ normal(0, 1);
  sd_TE1 ~ normal(0, 1);
  
  // Factor 2 (Bob's Beta Ordering)
  rho ~ beta(1, 1); 
  
  sd_G2 ~ normal(0, 1);
  sd_PE2 ~ normal(0, 1);
  sd_TE2 ~ normal(0, 1);

  // 2. Loadings (Standard Normal)
  lambda_1 ~ normal(0, 1);
  lambda_2 ~ normal(0, 1);
  sd_R ~ normal(0, 1);
  to_vector(beta) ~ normal(0, 5);

  // 3. Raw Random Vectors (Must be N(0,1) for the loop/standardization to work)
  w_1_a ~ normal(0, 1);
  w_2_a ~ normal(0, 1);
  
  w_1_pe ~ normal(0, 1);
  w_2_pe ~ normal(0, 1);
  
  w_1_te ~ normal(0, 1);
  w_2_te ~ normal(0, 1);

  // --- LIKELIHOOD ---
  // Vectorized for speed
  to_vector(Y) ~ normal(to_vector(mu), to_vector(rep_matrix(sd_R', No)));
}

generated quantities {
  // Re-calculate Variances for output
  real var_G1  = square(sd_G1);
  real var_PE1 = square(sd_PE1);
  real var_TE1 = square(sd_TE1);
  real var_LV1_total = var_G1 + var_PE1 + var_TE1;

  real var_G2  = square(sd_G2);
  real var_PE2 = square(sd_PE2);
  real var_TE2 = square(sd_TE2);
  real var_LV2_total = var_G2 + var_PE2 + var_TE2;

  vector[Nt] h2_traits; 
  // vector[Nt] pe2_traits;
  
  for (t in 1:Nt) {
    // Total Trait Variance = (Loading^2 * FactorVar) + Residual
    real var_trait_genetic = square(lambda[t,1]) * var_G1  + square(lambda[t,2]) * var_G2;
    real var_trait_total   = square(lambda[t,1]) * var_LV1_total + square(lambda[t,2]) * var_LV2_total + square(sd_R[t]);

    // Heritability calculation
    h2_traits[t] = var_trait_genetic / var_trait_total;
  }
}


