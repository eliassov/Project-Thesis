
data {
    int<lower = 1> No;        // Number of observations
    vector[No] Y;             // Trait observations
    
    int<lower = 1> K;         // Number of fixed effect predictors (intercept + sex)
    matrix[No, K] X;          // Design matrix 

    array[No] int animal;     // Individual ID map
    int<lower = 1> Na;        // Number of individuals
    cov_matrix[Na] A;         // Relatedness matrix
}

transformed data {
    matrix[Na, Na] LA;
    LA = cholesky_decompose(A);
}

parameters {

    vector[K] beta;           

    real<lower=0> sd_A;       
    real<lower=0> sd_E;       
    real<lower=0> sd_R;       
    vector[Na] alpha_std;     
    vector[Na] gamma_std;     
}

transformed parameters {
    real<lower = 0> var_A;
    real<lower = 0> var_E;
    real<lower = 0> var_R;
    vector[Na] alpha;
    vector[Na] gamma;

    var_A = square(sd_A);
    var_E = square(sd_E);
    var_R = square(sd_R);

    alpha = sd_A * (LA * alpha_std);
    gamma = sd_E * gamma_std;
}

model {
    vector[No] eta;

    // Priors
    // Using a generic normal prior for both coefficients
    beta ~ normal(0, 100); 

    sd_A ~ exponential(7); 
    sd_E ~ exponential(7);
    sd_R ~ exponential(7);

    alpha_std ~ normal(0, 1);
    gamma_std ~ normal(0, 1);

    vector[No] fixed_part = X * beta;
    for (o in 1:No) {
        eta[o] = fixed_part[o] + alpha[animal[o]] + gamma[animal[o]];
    }

    Y ~ normal(eta, sd_R);
}

generated quantities {
    real<lower=0> var_P;
    real<lower=0> heritability;

    var_P = var_A + var_E + var_R;
    heritability = var_A / var_P;
}



