data {
    
    int<lower = 1> No;  // Number of observations
    vector[No] Y;       // Trait observations (response variable)
    array[No] int animal;     // Individual identity for each observation
    int<lower = 1> Na;  // Number of individuals
    cov_matrix[Na] A;   // Additive-genetic relationship matrix
    
}


transformed data {
  
  real mean_Y;  // Sample mean of observed Y
  real sd_Y;    // Sample standard deviation of Y
  matrix[Na, Na] LA;  // Lower-triangular cholesky factor of A

  mean_Y = mean(Y);
  sd_Y = sd(Y);

  LA = cholesky_decompose(A);
  
}


parameters {
  
  real mu;  // Overall intercept for centered predictors
  real<lower=0> sd_A;  // Standard deviation of the additive genetic effects
  real<lower=0> sd_E;  // Standard deviation of the permanent individual effects
  real<lower=0> sd_R;  // Standard deviation of the temporary residual effects
  vector[Na] alpha_std;  // Precursors for the additive genetic effects
  vector[Na] gamma_std;  // Standardised permanent individual effects
  
}

transformed parameters {
  
  real<lower = 0> var_A;  // Additive genetic variance
  real<lower = 0> var_E;  // Permanent individual variance
  real<lower = 0> var_R;  // Temporary residual variance
  vector[Na] alpha;  // Raw additive genetic effects (breeding values)
  vector[Na] gamma;  // Raw permanent individual effects
  
  var_A = square(sd_A);
  var_E = square(sd_E);
  var_R = square(sd_R);
  
  alpha = sd_A * (LA * alpha_std);
  gamma = sd_E * gamma_std;
  
}

model {
  
  vector[No] eta; // Expected phenotype for each observation
  
  mu ~ normal(76, 38); // WHAT TO PUT HERE??
  
  sd_A ~ exponential(7);
  // WHY DOES THE STANDARD DEVIATION FOLLOW AN EXPONENTIAL DISTRIBUTION? Look into the distribution stuff

  sd_E ~ exponential(7);

  sd_A ~ exponential(7);

  // Alternative priors for sd include e.g. normal, student_t (and gamma?)
  
  alpha_std ~ normal(0, 1);  // Implies alpha ~ MVN(0, var_A * A)
  gamma_std ~ normal(0, 1);  // Implies alpha ~ normal(0, sd_I)
  
  for (o in 1:No)
    eta[o] = mu + alpha[animal[o]] + gamma[animal[o]];
    
  Y ~ normal(eta, sd_R); // Or target += normal_lpdf(Y | eta, sd_R);
  
}

generated quantities {
  
  real<lower=0> var_P;  // Phenotypic variance (i.e. total variance)
  real<lower=0> evolvability;  // Mean standardized additive genetic variance
  real<lower=0> heritability;  // Proportion of total variance due to additive genetic variation
  
  var_P = var_A + var_E + var_R;
  
  evolvability = var_A / square(mu);
  heritability = var_A / var_P;
  
}
