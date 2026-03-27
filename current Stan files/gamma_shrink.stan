
data {
  // Dimensions
  int<lower=1> No_morph; 
  int<lower=1> No_repro;
  int<lower=1> No_surv;
  int<lower=1> Nt;  
  int<lower=1> Nlv; 
  int<lower=1> Na;  

  int<lower=0> K_morph;
  int<lower=0> K_repro;
  int<lower=0> K_surv;
  
  matrix[No_morph, K_morph] X_morph;
  matrix[No_repro, K_repro] X_repro;
  matrix[No_surv, K_surv] X_surv;
  
  // Responses
  matrix[No_morph, Nt] Y_morph; 
  array[No_repro] int<lower=0> Y_repro; // Poisson count
  array[No_surv] int<lower=0, upper=1> Y_surv; // Bernoulli 0/1
  
  // Pedigree & Links
  array[No_morph] int<lower=1, upper=Na> animal_morph; 
  array[No_repro] int<lower=1, upper=Na> animal_repro;
  array[No_surv]  int<lower=1, upper=Na> animal_surv;
  
  array[Na] int<lower=1, upper=Na+1> dam; 
  array[Na] int<lower=1, upper=Na+1> sire;
  vector[Na] dii; 
  
  // Random Effects Indices
  int<lower=1> N_year;
  array[No_morph] int<lower=1, upper=N_year> year_morph;
  array[No_repro] int<lower=1, upper=N_year> year_repro;
  array[No_surv]  int<lower=1, upper=N_year> year_surv;
  
  int<lower=1> N_init;
  array[No_morph] int<lower=1, upper=N_init> init_morph;
}


transformed data {
  real<lower=0> nu = 3.0; // Hyperparameter from the paper 
}


parameters {
  // Just Intercepts (since we dropped the X matrix)
  matrix[K_morph, Nt] beta_morph;
  vector[K_repro] beta_repro;
  vector[K_surv] beta_surv;

  // Fitness loadings/gradients 
  vector[Nlv] gamma_repro_raw;
  vector[Nlv] gamma_surv_raw;

  // Variances
  vector<lower=0>[Nt] sd_R; 
  vector<lower=0, upper=1>[Nlv] h2_lv; 
  
  matrix[Nt, Nlv] lambda_raw;

  // Latent Effects
  matrix[Na, Nlv] w_a;  
  matrix[Na, Nlv] w_pe;
  
  // Random Effects (Year and Measurer)
  matrix[N_year, Nt] z_year_morph;
  vector<lower=0>[Nt] sd_year_morph;
  
  vector[N_year] z_year_repro;
  real<lower=0> sd_year_repro;
  
  vector[N_year] z_year_surv;
  real<lower=0> sd_year_surv;
  
  matrix[N_init, Nt] z_init_morph;
  vector<lower=0>[Nt] sd_init_morph;
  
  matrix[Nt, Nlv] phi; // Local shrinkage (per trait and LV)
  
  vector<lower=0>[Nlv] phi_r; // Shrinkage for reproduction
  vector<lower=0>[Nlv] phi_s; // Shrinkage for survival 
  
  vector<lower=0>[Nlv] delta; // Global shrinkage (per LV)
  
  // Shrinkage hyperparameterers
  real<lower=0> a_1; 
  real<lower=0> a_2;
  
  
}

transformed parameters {

  matrix[Na, Nlv] G_effect; 
  matrix[Na, Nlv] LV;
  vector<lower=0>[Nlv] tau;
  matrix[Nt, Nlv] Lambda;
  vector[Nlv] gamma_repro;
  vector[Nlv] gamma_surv;


  // Pedigree Recursion
  {
    matrix[Na+1, Nlv] aug_a;
    aug_a[Na+1] = rep_row_vector(0, Nlv);
    
    for (i in 1:Na) {
       aug_a[i] = 0.5 * aug_a[dam[i]] + 0.5 * aug_a[sire[i]] 
                  + sqrt(dii[i]) * w_a[i]; 
    }
    G_effect = aug_a[1:Na];
  }

  // Latent Factor
  LV = diag_post_multiply(G_effect, sqrt(h2_lv)) 
       + diag_post_multiply(w_pe, sqrt(1 - h2_lv));

  // Global Shrinkage
  tau[1] = delta[1];
  for(i in 2:Nlv) {
    tau[i] = tau[i-1] * delta[i];
  }
  
  // Local shrinkage and fitness gradients
  for(t in 1:Nt) {
    Lambda[t, ] = lambda_raw[t, ] ./ sqrt(phi[t, ] .* tau'); 
  }
  gamma_repro = gamma_repro_raw ./ sqrt(phi_r .* tau);
  gamma_surv  = gamma_surv_raw  ./ sqrt(phi_s .* tau);
}


model {
  matrix[N_year, Nt] year_eff;
  matrix[N_init, Nt] init_eff;
  matrix[No_morph, Nt] mu_morph;
  vector[No_repro] log_lambda;
  vector[No_surv] logit_p;

  // --- PRIORS ---
  h2_lv ~ beta(2.5, 2.5); 
  sd_R ~ normal(0, 1);
  to_vector(lambda_raw) ~ normal(0,1);
  to_vector(beta_morph) ~ normal(0,1);
  beta_repro ~ normal(0,1);
  beta_surv ~ normal(0,1);
  gamma_repro_raw ~ normal(0, 1);
  gamma_surv_raw ~ normal(0, 1);

  to_vector(w_a) ~ std_normal();
  to_vector(w_pe) ~ std_normal();
  to_vector(z_year_morph) ~ std_normal();
  z_year_repro ~ std_normal();
  z_year_surv ~ std_normal();
  to_vector(z_init_morph) ~ std_normal();
  
  sd_year_morph ~ exponential(2); 
  sd_year_repro ~ exponential(2); 
  sd_year_surv ~ exponential(2); 
  sd_init_morph ~ exponential(2);
  
  a_1 ~ gamma(2,1);
  a_2 ~ gamma(2,1);
  
  delta[1] ~ gamma(a_1,1);
  for(h in 2:Nlv) {
    delta[h] ~ gamma(a_2,1);
  }
  
  to_vector(phi) ~ gamma(nu/2, nu/2);
  phi_r ~ gamma(nu/2, nu/2);
  phi_s ~ gamma(nu/2, nu/2);


  // --- LIKELIHOOD CALCULATIONS ---
  
  // Scale the random effects
  year_eff = diag_post_multiply(z_year_morph, sd_year_morph);
  init_eff = diag_post_multiply(z_init_morph, sd_init_morph);
  
  // Morphology
  mu_morph = X_morph * beta_morph 
             + LV[animal_morph] * Lambda' 
             + year_eff[year_morph] 
             + init_eff[init_morph];
             
  for(t in 1:Nt) {
    Y_morph[, t] ~ normal(mu_morph[, t], sd_R[t]);
  }
    
  // Reproduction
  log_lambda = X_repro * beta_repro
               + LV[animal_repro] * gamma_repro
               + (z_year_repro * sd_year_repro)[year_repro];
               
  Y_repro ~ poisson_log(log_lambda);
  
  // Survival
  logit_p = X_surv * beta_surv
            + LV[animal_surv] * gamma_surv
            + (z_year_surv * sd_year_surv)[year_surv];
            
  Y_surv ~ bernoulli_logit(logit_p);
}


generated quantities {
  // --- 1. MORPHOLOGY: Defining the Factors ---
  vector[Nlv] morph_var_explained; 
  vector[Nlv] morph_prop_var_explained; 
  real total_morph_var = 0;

  // --- 2. FITNESS: The Downstream Effects ---
  // (Since the LV has a variance of 1, the variance it explains 
  // on the linear predictor scale is simply the squared coefficient)
  vector[Nlv] repro_var_explained;
  vector[Nlv] surv_var_explained;

  // Calculate absolute variances
  for (h in 1:Nlv) {
    // A. Morphology (Sum of squared loadings for this factor)
    morph_var_explained[h] = 0; 
    for (t in 1:Nt) {
      morph_var_explained[h] += square(Lambda[t, h]);
    }
    total_morph_var += morph_var_explained[h]; 
    
    // B. Fitness (Squared selection gradients) (need to verify if this makes sense)
    repro_var_explained[h] = square(gamma_repro[h]);
    surv_var_explained[h]  = square(gamma_surv[h]);
  }
  
  // Calculate relative proportions (only for morphology)
  for (h in 1:Nlv) {
    morph_prop_var_explained[h] = morph_var_explained[h] / total_morph_var;
  }
}

