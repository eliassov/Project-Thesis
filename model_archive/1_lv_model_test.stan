data {
  int<lower=1> No; // Number of observations
  int<lower=1> Nt; // Number of traits 
  int<lower=1> Nlv; // Number of latent variables (Passed as 1 now)
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

  // --- FACTOR 1 PARAMETERS (The Only Factor) ---
  // real<lower=0> sd_G1; 
  // real<lower=0> sd_PE1;
  // real<lower=0> sd_TE1;

  // --- LOADINGS (Lambdas) ---
  // Factor 1: Trait 1 is fixed to 1.0, others are free
  vector[Nt-1] lambda_1_free;
  real<lower=0> lambda_1_1;

  
  // --- RAW RANDOM VECTORS (Unit Normal Noise) ---
  vector[Na] w_1_a; 
  vector[Na] w_1_pe; 
  vector[No] w_1_te; 
  
  simplex[3] prop_var_1;
}

transformed parameters {
  // --- 1. RECURSIVE LOOP (Generate Raw Pedigree Shape) ---
  vector[Na] psi_1_a_raw;
  
  {
    // Augment vectors to handle missing parents (Index Na+1 = 0)
    vector[Na+1] aug_1;
    aug_1[Na+1] = 0; 

    for (i in 1:Na) {
      // Just accumulate structure (Assume scale = 1.0 for now)
      real par1 = 0.5 * aug_1[dam[i]] + 0.5 * aug_1[sire[i]];
      aug_1[i] = par1 + sqrt(dii[i]) * w_1_a[i]; 
    }
    // Extract the actual animals
    psi_1_a_raw = aug_1[1:Na];
  }

  // --- 2. BOB'S STANDARDIZATION ---
  // Force vector to be N(0,1)
  // vector[Na] psi_1_a_std = (psi_1_a_raw - mean(psi_1_a_raw)) / sd(psi_1_a_raw);

  // --- 3. SCALE AND COMBINE ---
  vector[No] psi_1_total;
  
  // Formula: (Scale_G * Gen) + (Scale_PE * PE) + (Scale_TE * TE)
  // psi_1_total = (sd_G1  * psi_1_a_std[animal]) 
  //             + (sd_PE1 * w_1_pe[animal]) 
  //             + (sd_TE1 * w_1_te);
              
  psi_1_total = (sqrt(prop_var_1[1])  * psi_1_a_raw[animal]) 
  + (sqrt(prop_var_1[2]) * w_1_pe[animal]) 
  + (sqrt(prop_var_1[3])* w_1_te);

  // --- 4. CONSTRUCT LAMBDA VECTOR ---
  vector[Nt] lambda;
  
  // lambda[1] = 1.0;                  // ANCHOR: Scale Fixed
  lambda[1] = lambda_1_1;
  for(t in 2:Nt) lambda[t] = lambda_1_free[t-1];

  // --- 5. LINEAR PREDICTOR ---
  matrix[No, Nt] mu;
  matrix[No, Nt] fixed = X * beta;
  
  // Standard Factor Model Equation: Y = Fixed + (Factor * Loadings)
  // We use rank-1 outer product here: vector * vector_transpose
  mu = fixed + psi_1_total * lambda'; 
}

model {
  // --- PRIORS ---
  // Variances 
  // sd_G1 ~ normal(0, 1);
  // sd_PE1 ~ normal(0, 1);
  // sd_TE1 ~ normal(0, 1);
  
  prop_var_1 ~ dirichlet(rep_vector(1.5, 3));
  
  // Loadings
  lambda_1_free ~ normal(0, 1);
  lambda_1_1 ~ normal(0, 1);
  sd_R ~ normal(0, 1);
  to_vector(beta) ~ normal(0, 5);

  // Raw Random Vectors 
  w_1_a ~ normal(0, 1);
  w_1_pe ~ normal(0, 1);
  w_1_te ~ normal(0, 1);

  // --- LIKELIHOOD ---
  to_vector(Y) ~ normal(to_vector(mu), to_vector(rep_matrix(sd_R', No)));
}

generated quantities {
  // Re-calculate Variances for output
  // real var_G1  = square(sd_G1);
  // real var_PE1 = square(sd_PE1);
  // real var_TE1 = square(sd_TE1);
  real var_G1  = prop_var_1[1];
  real var_PE1 = prop_var_1[2];
  real var_TE1 = prop_var_1[3];
  real var_LV1_total = var_G1 + var_PE1 + var_TE1;

  vector[Nt] h2_traits; 
  
  for (t in 1:Nt) {
    // Total Trait Variance = (Loading^2 * FactorVar) + Residual
    // Note: Since we only have 1 factor, the math is simpler
    real var_trait_genetic = square(lambda[t]) * var_G1;
    real var_trait_total   = square(lambda[t]) * var_LV1_total + square(sd_R[t]);

    // Heritability calculation
    h2_traits[t] = var_trait_genetic / var_trait_total;
  }
}