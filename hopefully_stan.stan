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

parameters {
  // Intercepts (Since we dropped the X matrix)
  
  matrix[K_morph, Nt] beta_morph;
  vector[K_repro] beta_repro;
  vector[K_surv] beta_surv;

  // Selection Gradients (The link between Size and Fitness)
  vector[Nlv] gamma_repro;
  vector[Nlv] gamma_surv;

  // Variances
  vector<lower=0>[Nt] sd_R; 
  real<lower=0, upper=1> rho;
  vector<lower=0, upper=1>[Nlv] h2_lv; 
  
  matrix[Nt, Nlv] lambda_raw;
  
  real<lower=0> lambda_1_1;
  vector[Nt-1] lambda_free_1;
  real<lower=0> lambda_2_2;
  vector[Nt-1] lambda_free_2;


  // Latent Effects
  matrix[Nlv, Na] w_a;  
  matrix[Nlv, Na] w_pe; 
  
  // Random Effects (Year and Measurer)
  matrix[N_year, Nt] z_year_morph;
  vector<lower=0>[Nt] sd_year_morph;
  
  vector[N_year] z_year_repro;
  real<lower=0> sd_year_repro;
  
  vector[N_year] z_year_surv;
  real<lower=0> sd_year_surv;
  
  matrix[N_init, Nt] z_init_morph;
  vector<lower=0>[Nt] sd_init_morph;
}

transformed parameters {
  // 1. Pedigree Recursion
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

  // 2. Latent Factor
  matrix[Na, Nlv] LV;
  for(k in 1:Nlv) {
    LV[, k] = sqrt(h2_lv[k]) * G_effect[, k] + sqrt(1 - h2_lv[k]) * w_pe[k, ]';
  }


  matrix[Nt, Nlv] Lambda = rep_matrix(0.0, Nt, Nlv);
  
  Lambda[1,1] = lambda_1_1;
  
  for (t in 2:Nt) {
    Lambda[t, 1] = lambda_free_1[t-1];
  }
  
  Lambda[2,2] = lambda_2_2;
  Lambda[1,2] = lambda_free_2[1];
  
  for (t in 3:Nt) {
    Lambda[t, 2] = lambda_free_2[t-2];
  }

}

model {
  // Priors
  // rho ~ beta(1, 2);
  h2_lv ~ beta(2.5, 2.5); 
  lambda_1_1 ~ normal(0, 1.0);
  lambda_free_1 ~ normal(0, 1.0);
  
  lambda_2_2 ~ normal(0, 0.2);
  lambda_free_2 ~ normal(0, 0.2);

  sd_R ~ normal(0, 1);
  
  to_vector(beta_morph) ~ normal(0,1);
  beta_repro ~ normal(0,1);
  beta_surv ~ normal(0,1);
  
  gamma_repro ~ normal(0, 1);
  gamma_surv ~ normal(0, 1);

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

  // --- LIKELIHOODS ---
  
  // 1. Morphology
  for(i in 1:No_morph){
    for(t in 1:Nt){
      real mu =  
                dot_product(X_morph[i], beta_morph[, t])
                + LV[animal_morph[i], 1] * Lambda[t, 1] 
                + LV[animal_morph[i], 2] * Lambda[t, 2]
                + z_year_morph[year_morph[i], t] * sd_year_morph[t]
                + z_init_morph[init_morph[i], t] * sd_init_morph[t];
      Y_morph[i, t] ~ normal(mu, sd_R[t]);
    }
  }
  
// // 2. Reproduction (Poisson)
  for(i in 1:No_repro){
    real log_lambda =  dot_product(X_repro[i], beta_repro)
                      + gamma_repro[1] * LV[animal_repro[i], 1]
                      + gamma_repro[2] * LV[animal_repro[i], 2]
                      + z_year_repro[year_repro[i]] * sd_year_repro;
    Y_repro[i] ~ poisson_log(log_lambda);
  }

  // 3. Survival (Bernoulli)
  for(i in 1:No_surv){
    real logit_p = dot_product(X_surv[i], beta_surv)
                    + gamma_surv[1] * LV[animal_surv[i], 1]
                   + gamma_surv[2] * LV[animal_surv[i], 2]
                   + z_year_surv[year_surv[i]] * sd_year_surv;
    Y_surv[i] ~ bernoulli_logit(logit_p);
  }

}


generated quantities {
  real rho;
  real var_lv1 = 0;
  real var_lv2 = 0; 
  
  for (t in 1:Nt) {
    var_lv1 += square(Lambda[t,1]);
    var_lv2 += square(Lambda[t,2]);
  }
  
  rho = var_lv2 / var_lv1;
}

