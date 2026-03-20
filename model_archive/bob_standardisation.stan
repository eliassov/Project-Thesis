data {
  int<lower=1> No; 
  int<lower=1> Nt; 
  int<lower=1> Nlv; // 2
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
  vector<lower=0>[Nt] sd_R;

  // --- FACTOR PARAMETERS ---
  real<lower=0, upper=1> rho; // Shrinkage for Factor 2
  
  // Weights for G vs PE in the raw factor construction
  // We need these to determine "how much" of the raw variance is genetic
  vector<lower=0, upper=1>[Nlv] h2_lv_raw; 

  // --- LOADINGS ---
  matrix[Nt, Nlv] lambda_raw;

  // --- RAW NOISE (Must be anchored!) ---
  matrix[Nlv, Na] w_a;  
  matrix[Nlv, Na] w_pe; 
}

transformed parameters {
  // --- 1. PEDIGREE RECURSION (Build G) ---
  matrix[Na, Nlv] G_raw; 
  {
    matrix[Na+1, Nlv] aug_a;
    aug_a[Na+1] = rep_row_vector(0, Nlv);
    for (i in 1:Na) {
       aug_a[i] = 0.5 * aug_a[dam[i]] + 0.5 * aug_a[sire[i]] 
                  + sqrt(dii[i]) * w_a[, i]'; 
    }
    G_raw = aug_a[1:Na];
  }

  // --- 2. CONSTRUCT RAW FACTOR ---
  // We mix G and PE based on h2_lv_raw parameter
  matrix[Na, Nlv] LV_raw;
  for(k in 1:Nlv) {
      LV_raw[, k] = sqrt(h2_lv_raw[k]) * G_raw[, k] 
                  + sqrt(1 - h2_lv_raw[k]) * w_pe[k, ]';
  }

  // --- 3. HARD STANDARDIZATION (Bob's Logic) ---
  // This forces the factor to have exactly mean=0, sd=1 in this sample
  matrix[Na, Nlv] LV;
  for(k in 1:Nlv){
    real m = mean(LV_raw[, k]);
    real s = sd(LV_raw[, k]);
    
    // The Standardization
    LV[, k] = (LV_raw[, k] - m) / s; 
  }

  // --- 4. LOADINGS & PREDICTOR ---
  matrix[Nt, Nlv] Lambda;
  Lambda[, 1] = lambda_raw[, 1]; 
  Lambda[, 2] = lambda_raw[, 2] * sqrt(rho); // Shrinkage

  matrix[No, Nt] mu;
  matrix[No, Nt] fixed = X * beta;
  
  for(i in 1:No){
    for(t in 1:Nt){
      mu[i,t] = fixed[i,t] 
                + LV[animal[i], 1] * Lambda[t, 1] 
                + LV[animal[i], 2] * Lambda[t, 2];
    }
  }
}

model {
  // --- PRIORS ---
  rho ~ beta(1, 1);
  h2_lv_raw ~ beta(2, 2);
  
  to_vector(lambda_raw) ~ std_normal();
  sd_R ~ normal(0, 1);
  to_vector(beta) ~ normal(0, 5);

  // --- THE CRITICAL PART ---
  // Even though we standardize LV later, w_a and w_pe MUST have priors.
  // This anchors the "Raw" scale so the sampler doesn't drift to infinity.
  to_vector(w_a) ~ std_normal();
  to_vector(w_pe) ~ std_normal();

  // --- LIKELIHOOD ---
  // Uses the standardized LV via 'mu'
  to_vector(Y) ~ normal(to_vector(mu), to_vector(rep_matrix(sd_R', No)));
}

