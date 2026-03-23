# ==============================================================================
# SCRIPT 2: HPC STAN EXECUTION
# ==============================================================================

library(rstan)

# 1. Load the pre-processed data list
dataset_joint <- readRDS("data/clean_stan_data.rds")

# 2. Define the Initialization Function
init_fn_joint <- function() {
  L_init <- matrix(0, nrow = dataset_joint$Nt, ncol = dataset_joint$Nlv)
  L_init[, 1] <- rnorm(dataset_joint$Nt, 0.5, 0.1) # Size (mostly positive)
  L_init[, 2] <- rnorm(dataset_joint$Nt, 0, 0.05)   # Shape (small noise)
  
  list(
    Lambda = L_init,
    h2_lv = as.array(runif(dataset_joint$Nlv, 0.4, 0.6)),
    sd_R = as.array(runif(dataset_joint$Nt, 0.5, 1.0)),
    beta_morph = matrix(rnorm(dataset_joint$K_morph * dataset_joint$Nt, 0, 0.1), 
                        nrow = dataset_joint$K_morph, ncol = dataset_joint$Nt),
    beta_repro = as.array(rnorm(dataset_joint$K_repro, 0, 0.1)),
    beta_surv  = as.array(rnorm(dataset_joint$K_surv, 0, 0.1)),
    gamma_repro = as.array(rnorm(dataset_joint$Nlv, 0, 0.1)),
    gamma_surv  = as.array(rnorm(dataset_joint$Nlv, 0, 0.1)),
    w_a = matrix(rnorm(dataset_joint$Nlv * dataset_joint$Na, 0, 0.1), 
                 nrow = dataset_joint$Nlv, ncol = dataset_joint$Na),
    w_pe = matrix(rnorm(dataset_joint$Nlv * dataset_joint$Na, 0, 0.1), 
                  nrow = dataset_joint$Nlv, ncol = dataset_joint$Na),
    sd_year_morph = as.array(runif(dataset_joint$Nt, 0.05, 0.2)),
    sd_init_morph = as.array(runif(dataset_joint$Nt, 0.05, 0.2)),
    sd_year_repro = runif(1, 0.05, 0.2),
    sd_year_surv  = runif(1, 0.05, 0.2),
    z_year_morph = matrix(rnorm(dataset_joint$N_year * dataset_joint$Nt, 0, 0.1), 
                          nrow = dataset_joint$N_year, ncol = dataset_joint$Nt),
    z_year_repro = as.array(rnorm(dataset_joint$N_year, 0, 0.1)),
    z_year_surv  = as.array(rnorm(dataset_joint$N_year, 0, 0.1)),
    z_init_morph = matrix(rnorm(dataset_joint$N_init * dataset_joint$Nt, 0, 0.1), 
                          nrow = dataset_joint$N_init, ncol = dataset_joint$Nt)
  )
}

# 3. HPC Stan Configuration
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# 4. Run the Model
out_joint <- stan(
  file = "models/baseline_2lv.stan",
  data = dataset_joint,
  init = init_fn_joint,        
  chains = 4, 
  pars = c("beta_morph", "beta_repro", "beta_surv",
           "gamma_repro", "gamma_surv", "Lambda", "sd_R", 
           "sd_year_morph", "sd_year_repro", "sd_year_surv", "sd_init_morph", 
           "h2_lv", "rho"),
  control = list(adapt_delta = 0.95), 
  iter = 20000,          
  warmup = 10000
)

# 5. Save output directly to the outputs folder
saveRDS(out_joint, 'outputs/Output_hopefully_no_rho.rds')