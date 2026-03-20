data {
  int<lower=1> No; 
  int<lower=1> Nt;  
  int<lower=1> Nlv; 
  matrix[No, Nt] Y; 
  int<lower=1> K;
  matrix[No, K] X;
  int<lower=1> Na; 
  array[No] int<lower=1, upper=Na> animal; 
  array[Na] int<lower=1, upper=Na+1> dam; 
  array[Na] int<lower=1, upper=Na+1> sire;
  vector[Na] dii; 
  
  // --- NEW: RANDOM EFFECTS DATA ---
  int<lower=1> N_year;
  array[No] int<lower=1, upper=N_year> year_idx;
  int<lower=1> N_init;
  array[No] int<lower=1, upper=N_init> init_idx;
}

parameters {
  matrix[K, Nt] beta;      
  vector<lower=0>[Nt] sd_R; 

  // --- FACTOR PARAMETERS ---
  real<lower=0, upper=1> rho; 
  vector<lower=0, upper=1>[Nlv] h2_lv; 

  // --- LOADINGS ---
  matrix[Nt, Nlv] lambda_raw;

  // --- LATENT EFFECTS ---
  matrix[Nlv, Na] w_a;  
  matrix[Nlv, Na] w_pe; 
  
  // --- NEW: RANDOM EFFECTS PARAMETERS (Non-centered) ---
  matrix[N_year, Nt] z_year;
  vector<lower=0>[Nt] sd_year;
  
  matrix[N_init, Nt] z_init;
  vector<lower=0>[Nt] sd_init;
}

transformed parameters {
  // --- 1. PEDIGREE RECURSION ---
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

  // --- 2. CONSTRUCT THE LATENT FACTOR ---
  matrix[Na, Nlv] LV;
  for(k in 1:Nlv) {
    LV[, k] = sqrt(h2_lv[k]) * G_effect[, k] 
            + sqrt(1 - h2_lv[k]) * w_pe[k, ]';
  }

  // --- 3. APPLY SHRINKAGE ---
  matrix[Nt, Nlv] Lambda;
  Lambda[, 1] = lambda_raw[, 1];            
  Lambda[, 2] = lambda_raw[, 2] * sqrt(rho); 

  // --- 4. LINEAR PREDICTOR (UPDATED) ---
  matrix[No, Nt] mu;
  matrix[No, Nt] fixed = X * beta;
  
  for(i in 1:No){
    for(t in 1:Nt){
      // Y = Fixed + Factors + Year + Measurer + Error
      mu[i,t] = fixed[i,t] 
                + LV[animal[i], 1] * Lambda[t, 1] 
                + LV[animal[i], 2] * Lambda[t, 2]
                + z_year[year_idx[i], t] * sd_year[t]
                + z_init[init_idx[i], t] * sd_init[t];
    }
  }
}

model {
  // Priors
  rho ~ beta(1.7, 2.5);
  h2_lv ~ beta(2, 2); 
  
  to_vector(lambda_raw) ~ std_normal();
  sd_R ~ normal(0, 1);
  to_vector(beta) ~ normal(0, 5);

  // Generative Noise 
  to_vector(w_a) ~ std_normal();
  to_vector(w_pe) ~ std_normal();
  
  // --- NEW: RANDOM EFFECTS PRIORS ---
  to_vector(z_year) ~ std_normal();
  to_vector(z_init) ~ std_normal();
  sd_year ~ exponential(2); 
  sd_init ~ exponential(2);

  // Likelihood
  to_vector(Y) ~ normal(to_vector(mu), to_vector(rep_matrix(sd_R', No)));
}

generated quantities {
  // --- VARIANCE PARTITIONING (UPDATED) ---
  vector[Nt] var_G_trait;   
  vector[Nt] var_PE_trait;  
  vector[Nt] var_year_trait; // NEW
  vector[Nt] var_init_trait; // NEW
  vector[Nt] var_R_trait;   
  vector[Nt] var_P_total;   
  
  vector[Nt] h2_trait;      
  vector[Nt] pe2_trait;     

  for (t in 1:Nt) {
    real g_sum = 0;
    real pe_sum = 0;
    
    for (k in 1:Nlv) {
       g_sum += square(Lambda[t, k]) * h2_lv[k]; 
       pe_sum += square(Lambda[t, k]) * (1 - h2_lv[k]);
    }
    
    var_G_trait[t] = g_sum;
    var_PE_trait[t] = pe_sum;
    var_year_trait[t] = square(sd_year[t]);
    var_init_trait[t] = square(sd_init[t]);
    var_R_trait[t] = square(sd_R[t]);
    
    // Total Variance now includes Year and Init
    var_P_total[t] = var_G_trait[t] + var_PE_trait[t] + var_year_trait[t] + var_init_trait[t] + var_R_trait[t];
    
    h2_trait[t] = var_G_trait[t] / var_P_total[t];
    pe2_trait[t] = var_PE_trait[t] / var_P_total[t];
  }
}

