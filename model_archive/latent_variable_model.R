# SETUP AND LIBRARIES

personal_lib <- "~/R/library"
if (!dir.exists(personal_lib)) dir.create(personal_lib, recursive = TRUE)
.libPaths(c(personal_lib, .libPaths()))

required_packages <- c("rstan", "tidyverse", "nadiv") 
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, lib = personal_lib, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# Set working directory
script_path <- NULL
cmd_args <- commandArgs(trailingOnly = FALSE)
if (any(grepl("^--file=", cmd_args))) {
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  script_path <- dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  script_path <- getwd()
}
setwd(script_path)
cat("Working directory set to:", getwd(), "\n")


# LOAD AND PREPROCESS PEDIGREE

pedData <- read.csv("pedigree_data.txt", sep = " ", header = TRUE, na.strings = "NA")

pedigree_clean <- pedData %>%
  dplyr::select(ringnr, dam, sire) %>%
  mutate(
    across(c(ringnr, dam, sire), ~ as.character(trimws(gsub("_.*$", "", .)))),
    dam = if_else(is.na(dam) | dam == "" | dam == "NA", NA_character_, dam),
    sire = if_else(is.na(sire) | sire == "" | sire == "NA", NA_character_, sire)
  )

# Add external parents
all_parents <- unique(c(pedigree_clean$dam, pedigree_clean$sire))
external_parents <- setdiff(all_parents, c(pedigree_clean$ringnr, NA_character_))

if (length(external_parents) > 0) {
  founder_df <- data.frame(ringnr = external_parents, dam = NA_character_, sire = NA_character_)
  pedigree_full <- bind_rows(founder_df, pedigree_clean) %>% distinct(ringnr, .keep_all = TRUE)
} else {
  pedigree_full <- pedigree_clean
}

pedigree_prepped <- prepPed(pedigree_full, gender = NULL, check = TRUE)


# PREPARE TRAIT DATA (LATENT VARIABLE)

traitData <- read.csv("morphology.txt", sep = ";", header = TRUE, 
                      na.strings = c("NA", ""), fileEncoding = "Windows-1252", 
                      stringsAsFactors = FALSE, colClasses = c(ringnr = "character"))

# Filter and clean data: Filter for Hestmannøy + both traits + valid sex
trait_subset_temp <- traitData %>%
  filter(sted_r == "hestmannøy", 
         !is.na(ving_h), 
         !is.na(nebb_l)) %>%
  
  mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr)))) %>%
  
  # Sex logic
  mutate(
    raw_sex = if_else(!is.na(scriptsex) & scriptsex != "u", scriptsex, fieldsex) # prioritize scriptsex over fieldsex
  ) %>%
  mutate(
    sex = case_when(
      raw_sex %in% c("m", "pm") ~ "m",
      raw_sex %in% c("f", "pf") ~ "f",
      TRUE ~ NA_character_ 
    )
  ) %>%
  filter(!is.na(sex))






# ALIGN DATA WITH PEDIGREE

valid_ids <- intersect(trait_subset_temp$ringnr, pedigree_prepped$ringnr)

# Create final data subset
trait_subset_lv <- trait_subset_temp %>%
  filter(ringnr %in% valid_ids)

# Standardize traits (Store mean/sd for back-transformation)
mean_wing <- mean(trait_subset_lv$ving_h, na.rm = TRUE)
sd_wing   <- sd(trait_subset_lv$ving_h, na.rm = TRUE)
mean_beak <- mean(trait_subset_lv$nebb_l, na.rm = TRUE)
sd_beak   <- sd(trait_subset_lv$nebb_l, na.rm = TRUE)

trait_subset_lv <- trait_subset_lv %>%
  mutate(
    ving_h_std = (ving_h - mean_wing) / sd_wing,
    nebb_l_std = (nebb_l - mean_beak) / sd_beak
  )

# BUILD A MATRIX (Subsetted)
get_ancestors <- function(ids, pedigree) {
  ancestors <- ids
  to_process <- ids
  while (length(to_process) > 0) {
    parents <- unique(c(pedigree$dam[pedigree$ringnr %in% to_process], 
                        pedigree$sire[pedigree$ringnr %in% to_process]))
    parents <- parents[!parents %in% ancestors & !is.na(parents)]
    ancestors <- c(ancestors, parents)
    to_process <- parents
  }
  return(ancestors)
}

all_relevant_ids <- get_ancestors(valid_ids, pedigree_prepped)

smaller_pedigree <- pedigree_prepped %>%
  filter(ringnr %in% all_relevant_ids) %>%
  arrange(match(ringnr, all_relevant_ids))

smaller_ped_prepped <- prepPed(smaller_pedigree, gender = NULL, check = TRUE)
A_small_nadiv <- makeA(smaller_ped_prepped)

# Subset A to valid phenotyped individuals
A_hestmannoy <- A_small_nadiv[valid_ids, valid_ids]

# Safety check
if (!all(unique(trait_subset_lv$ringnr) %in% rownames(A_hestmannoy))) {
  stop("Error: Data contains individuals not in the A-matrix!")
}
cat("Data and Matrix aligned. No =", nrow(trait_subset_lv), "Na =", nrow(A_hestmannoy), "\n")


# PREPARE STAN DATA

# Map animals to integers 1:Na based on the sorted matrix
matrix_ids <- rownames(A_hestmannoy)
animal_map <- setNames(seq_len(length(matrix_ids)), matrix_ids)
animal_idx <- animal_map[trait_subset_lv$ringnr]

Y_mat <- as.matrix(trait_subset_lv[, c("ving_h_std", "nebb_l_std")])
X_mat <- model.matrix(~ sex, data = trait_subset_lv)
K <- ncol(X_mat)

dataset_lv <- list(
  No = nrow(trait_subset_lv),
  Nt = 2,
  Y = Y_mat,
  X = X_mat,
  K = K,
  animal = animal_idx,
  Na = length(matrix_ids),
  A = as.matrix(A_hestmannoy)
)




# DIAGNOSTICS LOGGING

sink("Output_LV_PreRun_Diagnostics.txt") # Redirect output to file

cat("=== DATA & MATRIX DIAGNOSTICS ===\n")
cat("TimeStamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("1. Data Dimensions:\n")
cat("   Observations (No):", nrow(trait_subset_lv), "\n")
cat("   Individuals (Na): ", nrow(A_hestmannoy), "\n")
cat("   Traits (Nt):      ", 2, "\n\n")

cat("2. Fixed Effects (X Matrix):\n")
cat("   Predictors (K):   ", ncol(X_mat), "\n")
cat("   Colnames:         ", paste(colnames(X_mat), collapse=", "), "\n\n")

cat("3. ID Alignment Check:\n")
# Check if all data IDs are in matrix
ids_in_matrix <- unique(trait_subset_lv$ringnr) %in% rownames(A_hestmannoy)
cat("   All Data IDs in Matrix? ", all(ids_in_matrix), "\n")
if(!all(ids_in_matrix)) {
  cat("   MISSING IDs: ", head(unique(trait_subset_lv$ringnr)[!ids_in_matrix]), "...\n")
}

# Check matrix properties
cat("\n4. Matrix Properties:\n")
cat("   Symmetric?        ", isSymmetric(A_hestmannoy), "\n")
cat("   Positive Definite?", all(eigen(A_hestmannoy, only.values=TRUE)$values > -1e-6), "\n")
cat("   Mean Diagonal:    ", mean(diag(A_hestmannoy)), "(Expected ~1.0)\n")

sink() # Turn off redirection (restore console output)




# RUN STAN
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

nc <- 4
nw <- 5000
ni <- 10000
nt <- 1

tomonitor_lv <- c(
  "beta",             # Fixed effects (intercepts and sex)
  "lambda",           # Factor loadings for traits
  "var_LV_genetic",   # Proportion of LV variance that is genetic (h2_psi)
  "var_LV_perm_env",  # Proportion of LV variance that is permanent environmental
  "R_mat",            # The full 2x2 residual covariance matrix
  "h2_traits",        # Reconstructed heritabilities for traits
  "psi_ind"           # Individual latent scores (for caterpillar plots)
)

out_lv <- stan(file = 'lv_corr_resid.stan',
               data = dataset_lv,
               pars = tomonitor_lv,  # Optional, but it doesn't take too much longer to monitor all parameters, and then you get breeding values and environmental random effects
               chains = nc, iter = ni, warmup = nw, thin = nt,
               open_progress = FALSE,
               seed = 123)

saveRDS(out_lv, 'Output_LV_corr_resid.rds')
summ_lv <- summary(out_lv)$summary
write.csv(as.data.frame(summ_lv), file = "Output_LV_Summary_corr_resid.csv", row.names = TRUE)




# BACK TRANSFORMATION 
samples <- rstan::extract(out_lv)

# --- 1. Fixed Effects (Means) ---
mu_female_wing <- samples$beta[,,1][,1] * sd_wing + mean_wing
eff_male_wing  <- samples$beta[,,1][,2] * sd_wing 
mu_female_beak <- samples$beta[,,2][,1] * sd_beak + mean_beak
eff_male_beak  <- samples$beta[,,2][,2] * sd_beak

# --- 2. Loadings ---
loading_wing <- samples$lambda[,1]
loading_beak <- samples$lambda[,2]

# --- 3. Variance Components (Scaled back to mm^2) ---

# A. Additive Genetic Variance (via Latent Variable)
VA_wing <- samples$lambda[,1]^2 * samples$var_LV_genetic * sd_wing^2
VA_beak <- samples$lambda[,2]^2 * samples$var_LV_genetic * sd_beak^2

# B. Permanent Environmental Variance (via Latent Variable)
V_PE_wing <- samples$lambda[,1]^2 * samples$var_LV_perm_env * sd_wing^2
V_PE_beak <- samples$lambda[,2]^2 * samples$var_LV_perm_env * sd_beak^2

# C. Residual Variance (from R_mat)
# R_mat contains the total temporary environmental variance (shared + unique)
V_R_wing <- samples$R_mat[, 1, 1] * sd_wing^2
V_R_beak <- samples$R_mat[, 2, 2] * sd_beak^2

# D. Total Phenotypic Variance
VP_wing <- VA_wing + V_PE_wing + V_R_wing
VP_beak <- VA_beak + V_PE_beak + V_R_beak

# --- 4. Ratios and Correlations ---
# Trait-specific ratios
h2_wing_calc  <- VA_wing / VP_wing
h2_beak_calc  <- VA_beak / VP_beak
pe2_wing_calc <- V_PE_wing / VP_wing
pe2_beak_calc <- V_PE_beak / VP_beak

# Residual Correlation (Environmental correlation not captured by the LV)
res_corr <- samples$R_mat[, 1, 2] / sqrt(samples$R_mat[, 1, 1] * samples$R_mat[, 2, 2])

# --- 5. Compile Results ---
get_summary <- function(param_vector) {
  c(mean = mean(param_vector, na.rm = TRUE), 
    sd = sd(param_vector, na.rm = TRUE), 
    L95 = quantile(param_vector, 0.025, na.rm = TRUE), 
    U95 = quantile(param_vector, 0.975, na.rm = TRUE))
}

params_list <- list(
  Female_Mean_Wing = mu_female_wing,
  Male_Effect_Wing = eff_male_wing,
  Female_Mean_Beak = mu_female_beak,
  Male_Effect_Beak = eff_male_beak,
  
  Loading_Wing = loading_wing,
  Loading_Beak = loading_beak,
  
  LV_Prop_Genetic = samples$var_LV_genetic,
  LV_Prop_PermEnv = samples$var_LV_perm_env,
  Residual_Correlation = res_corr,
  
  VA_Wing   = VA_wing,
  V_PE_Wing = V_PE_wing,
  V_R_Wing  = V_R_wing,
  VP_Wing   = VP_wing,
  
  VA_Beak   = VA_beak,
  V_PE_Beak = V_PE_beak,
  V_R_Beak  = V_R_beak,
  VP_Beak   = VP_beak,
  
  h2_Wing = h2_wing_calc,
  h2_Beak = h2_beak_calc,
  pe2_Wing = pe2_wing_calc,
  pe2_Beak = pe2_beak_calc
)

posterior_summary <- do.call(rbind, lapply(params_list, get_summary))
posterior_summary <- as.data.frame(posterior_summary)
write.csv(posterior_summary, "Output_LV_Final_Results_corr_resid.csv", row.names = TRUE)
print(posterior_summary)




# # BACK TRANSFORMATION 
# samples <- rstan::extract(out_lv)
# 
# # Fixed effects (means):
# # beta[sample, predictor, trait]
# # predictor 1 = Intercept (female mean)
# # predictor 2 = Sexm (male effect)
# 
# # Trait 1: wing length
# mu_female_wing <- samples$beta[,,1][,1] * sd_wing + mean_wing
# eff_male_wing  <- samples$beta[,,1][,2] * sd_wing # Differences scale with SD only
# 
# # Trait 2: beak length
# mu_female_beak <- samples$beta[,,2][,1] * sd_beak + mean_beak
# eff_male_beak  <- samples$beta[,,2][,2] * sd_beak
# 
# # Loadings
# loading_wing <- samples$lambda[,1]
# loading_beak <- samples$lambda[,2]
# 
# # Variance components (unscaled)
# # h2_psi approach puts variance into lambda scaling
# VA_wing   <- samples$lambda[,1]^2 * samples$var_LV_genetic * sd_wing^2
# VA_beak   <- samples$lambda[,2]^2 * samples$var_LV_genetic * sd_beak^2
# 
# V_PE_wing   <- samples$lambda[,1]^2 * samples$var_LV_perm_env * sd_wing^2
# V_PE_beak   <- samples$lambda[,2]^2 * samples$var_LV_perm_env * sd_beak^2
# 
# V_TE_wing <- samples$lambda[,1]^2 * samples$var_LV_temp_env * sd_wing^2
# V_TE_beak <- samples$lambda[,2]^2 * samples$var_LV_temp_env * sd_beak^2
# 
# V_R_wing   <- samples$var_R[,1] * sd_wing^2
# V_R_beak   <- samples$var_R[,2] * sd_beak^2
# 
# # E. Total Phenotypic Variance (Reconstructed)
# VP_wing <- VA_wing + V_PE_wing + V_TE_wing + V_R_wing
# VP_beak <- VA_beak + V_PE_beak + V_TE_beak + V_R_beak
# 
# # --- 4. Ratios (Heritability & PE Ratio) ---
# # We calculate these using the reconstructed totals to be safe
# h2_wing_calc <- VA_wing / VP_wing
# h2_beak_calc <- VA_beak / VP_beak
# 
# pe2_wing_calc <- V_PE_wing / VP_wing
# pe2_beak_calc <- V_PE_beak / VP_beak
# 
# # --- 5. Compile Results ---
# # Define a helper function for summary stats to keep code clean
# get_summary <- function(param_vector) {
#   c(mean = mean(param_vector), 
#     sd = sd(param_vector), 
#     L95 = quantile(param_vector, 0.025), 
#     U95 = quantile(param_vector, 0.975))
# }
# 
# # Create list of all parameters we want to summarize
# params_list <- list(
#   # Fixed Effects
#   Female_Mean_Wing = mu_female_wing,
#   Male_Effect_Wing = eff_male_wing,
#   Female_Mean_Beak = mu_female_beak,
#   Male_Effect_Beak = eff_male_beak,
#   
#   # Latent Variable Structure
#   Loading_Wing = loading_wing,
#   Loading_Beak = loading_beak,
#   LV_Prop_Genetic = samples$var_LV_genetic,
#   LV_Prop_PermEnv = samples$var_LV_perm_env,
#   LV_Prop_TempEnv = samples$var_LV_temp_env,
#   
#   # Variance Components (Wing)
#   VA_Wing   = VA_wing,
#   V_PE_Wing = V_PE_wing,
#   V_TE_Wing = V_TE_wing, # Shared residual
#   V_R_Wing  = V_R_wing,  # Unique residual
#   VP_Wing   = VP_wing,   # Total
#   
#   # Variance Components (Beak)
#   VA_Beak   = VA_beak,
#   V_PE_Beak = V_PE_beak,
#   V_TE_Beak = V_TE_beak,
#   V_R_Beak  = V_R_beak,
#   VP_Beak   = VP_beak,
#   
#   # Ratios
#   h2_Wing = h2_wing_calc,
#   h2_Beak = h2_beak_calc,
#   pe2_Wing = pe2_wing_calc,
#   pe2_Beak = pe2_beak_calc
# )
# 
# # Apply summary function and bind into dataframe
# posterior_summary <- do.call(rbind, lapply(params_list, get_summary))
# posterior_summary <- as.data.frame(posterior_summary)
# write.csv(posterior_summary, "Output_LV_Final_Results_corr_resid_.csv", row.names = TRUE)
# print(posterior_summary)
