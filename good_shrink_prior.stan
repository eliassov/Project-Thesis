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
  vector[Na] dii; 
}

parameters {
  matrix[K, Nt] beta;     
  vector<lower=0>[Nt] sd_R; // Residual (TE)

  // --- FACTOR PARAMETERS ---
  real<lower=0, upper=1> rho; // Shrinkage for Factor 2
  vector<lower=0, upper=1>[Nlv] h2_lv; // Heritability of the FACTORS

  // --- LOADINGS ---
  matrix[Nt, Nlv] lambda_raw;

  // --- LATENT EFFECTS ---
  matrix[Nlv, Na] w_a;  // Genetic noise
  matrix[Nlv, Na] w_pe; // Permanent Env noise
}

transformed parameters {
  // --- 1. PEDIGREE RECURSION (Build G) ---
  matrix[Na, Nlv] G_effect; 
  {
    matrix[Na+1, Nlv] aug_a;
    aug_a[Na+1] = rep_row_vector(0, Nlv);
    for (i in 1:Na) {
       aug_a[i] = 0.5 * aug_a[dam[i]] + 0.5 * aug_a[sire[i]] 
                  + sqrt(dii[i]) * w_a[, i]'; 
    }
    G_effect = aug_a[1:Na];
  }

  // --- 2. CONSTRUCT THE LATENT FACTOR (G + PE) ---
  // We mix G and PE to keep the Factor variance approx 1.0
  // Factor = sqrt(h2)*G + sqrt(1-h2)*PE
  matrix[Na, Nlv] LV;
  for(k in 1:Nlv) {
    LV[, k] = sqrt(h2_lv[k]) * G_effect[, k] 
            + sqrt(1 - h2_lv[k]) * w_pe[k, ]';
  }

  // --- 3. APPLY SHRINKAGE TO LOADINGS ---
  matrix[Nt, Nlv] Lambda;
  Lambda[, 1] = lambda_raw[, 1];           // Factor 1: Full size
  Lambda[, 2] = lambda_raw[, 2] * sqrt(rho); // Factor 2: Shrunk

  // --- 4. LINEAR PREDICTOR ---
  matrix[No, Nt] mu;
  matrix[No, Nt] fixed = X * beta;
  
  for(i in 1:No){
    for(t in 1:Nt){
      // Y = Fixed + Factor1*L1 + Factor2*L2 + Error
      mu[i,t] = fixed[i,t] 
                + LV[animal[i], 1] * Lambda[t, 1] 
                + LV[animal[i], 2] * Lambda[t, 2];
    }
  }
}

model {
  // Priors
  rho ~ beta(1.7, 2.5);
  h2_lv ~ beta(2, 2); // Prior for factor heritability
  
  to_vector(lambda_raw) ~ std_normal();
  sd_R ~ normal(0, 1);
  to_vector(beta) ~ normal(0, 5);

  // Generative Noise (Standard Unit)
  to_vector(w_a) ~ std_normal();
  to_vector(w_pe) ~ std_normal();

  to_vector(Y) ~ normal(to_vector(mu), to_vector(rep_matrix(sd_R', No)));
}




generated quantities {
  // --- 1. VARIANCE PARTITIONING (Per Trait) ---
  vector[Nt] var_G_trait;   // Total Genetic Variance per trait
  vector[Nt] var_PE_trait;  // Total Perm-Env Variance per trait
  vector[Nt] var_R_trait;   // Residual Variance (TE)
  vector[Nt] var_P_total;   // Total Phenotypic Variance
  
  vector[Nt] h2_trait;      // Narrow-sense Heritability (h^2)
  vector[Nt] pe2_trait;     // Permanent Environment Ratio (pe^2)

  for (t in 1:Nt) {
    real g_sum = 0;
    real pe_sum = 0;
    
    // Sum variance contributions across all factors
    for (k in 1:Nlv) {
       // USE THE ALREADY CALCULATED LAMBDA
       // This Lambda already contains the sqrt(rho) shrinkage for factor 2
       
       // Genetic part: Loading^2 * Factor_Heritability
       g_sum += square(Lambda[t, k]) * h2_lv[k]; 
       
       // PE part: Loading^2 * (1 - Factor_Heritability)
       pe_sum += square(Lambda[t, k]) * (1 - h2_lv[k]);
    }
    
    var_G_trait[t] = g_sum;
    var_PE_trait[t] = pe_sum;
    var_R_trait[t] = square(sd_R[t]);
    
    // Total Variance
    var_P_total[t] = var_G_trait[t] + var_PE_trait[t] + var_R_trait[t];
    
    // Ratios
    h2_trait[t] = var_G_trait[t] / var_P_total[t];
    pe2_trait[t] = var_PE_trait[t] / var_P_total[t];
  }
}

