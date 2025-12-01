# ==============================================================================
# 1. SETUP & LIBRARIES
# ==============================================================================
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

# ==============================================================================
# 2. LOAD & PREPROCESS PEDIGREE
# ==============================================================================
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

# ==============================================================================
# 3. PREPARE TRAIT DATA (LATENT VARIABLE)
# ==============================================================================
traitData <- read.csv("morphology.txt", sep = ";", header = TRUE, 
                      na.strings = c("NA", ""), fileEncoding = "Windows-1252", 
                      stringsAsFactors = FALSE, colClasses = c(ringnr = "character"))

# Step A: Filter and Clean Data FIRST
# Filter for Hestmannøy + BOTH Traits + Valid Sex
trait_subset_temp <- traitData %>%
  filter(sted_r == "hestmannøy", 
         !is.na(ving_h), 
         !is.na(nebb_l)) %>%
  
  mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr)))) %>%
  
  # Sex Logic
  mutate(
    raw_sex = if_else(!is.na(scriptsex) & scriptsex != "u", scriptsex, fieldsex)
  ) %>%
  mutate(
    sex = case_when(
      raw_sex %in% c("m", "pm") ~ "m",
      raw_sex %in% c("f", "pf") ~ "f",
      TRUE ~ NA_character_ 
    )
  ) %>%
  filter(!is.na(sex))

# ==============================================================================
# 4. ALIGN DATA WITH PEDIGREE
# ==============================================================================
valid_ids <- intersect(trait_subset_temp$ringnr, pedigree_prepped$ringnr)

# Create FINAL data subset sorted by ID
trait_subset_lv <- trait_subset_temp %>%
  filter(ringnr %in% valid_ids) %>%
  arrange(ringnr)

# Standardize Traits (Store mean/sd for back-transformation)
mean_wing <- mean(trait_subset_lv$ving_h)
sd_wing   <- sd(trait_subset_lv$ving_h)
mean_beak <- mean(trait_subset_lv$nebb_l)
sd_beak   <- sd(trait_subset_lv$nebb_l)

trait_subset_lv <- trait_subset_lv %>%
  mutate(
    ving_h_std = (ving_h - mean_wing) / sd_wing,
    nebb_l_std = (nebb_l - mean_beak) / sd_beak
  )

# ==============================================================================
# 5. BUILD A-MATRIX (Subsetted)
# ==============================================================================
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

A_lv <- A_small_nadiv[valid_ids, valid_ids]
A_lv <- A_lv[order(rownames(A_lv)), order(colnames(A_lv))]

stopifnot(all(rownames(A_lv) == trait_subset_lv$ringnr))

# ==============================================================================
# 6. PREPARE STAN DATA
# ==============================================================================

animal_map <- setNames(seq_len(length(valid_ids)), valid_ids)
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
  Na = length(valid_ids),
  A = as.matrix(A_lv)
)

# ==============================================================================
# 7. RUN STAN
# ==============================================================================
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

nc <- 4
nw <- 4000
ni <- 8000
nt <- 1

tomonitor_lv <- c(
  "beta", "lambda", "sd_psi_a", "sd_psi_e", "sd_R", 
  "var_psi_a", "var_psi_e", "h2_psi", "h2_traits", "h2_traits_no_residual"
)

out_lv <- stan(file = 'latent_variable_stan.stan',
               data = dataset_lv,
               pars = tomonitor_lv,
               chains = nc, iter = ni, warmup = nw, thin = nt,
               open_progress = FALSE,
               seed = 123)

saveRDS(out_lv, 'Output_LV_Animal_Model.rds')
summ_lv <- summary(out_lv)$summary
write.csv(as.data.frame(summ_lv), file = "Output_LV_Summary.csv", row.names = TRUE)

# ==============================================================================
# 8. BACK-TRANSFORMATION (Updated for Fixed Effects)
# ==============================================================================
samples <- rstan::extract(out_lv)

# Flip sign if lambda[1] negative (though we constrain it, good practice)
# Since we constrain lambda1 > 0, this step is technically redundant but safe
# flip <- samples$lambda[,1] < 0
# samples$lambda[flip, ] <- -samples$lambda[flip, ]

# --- A. Fixed Effects (Means) ---
# beta[sample, predictor, trait]
# predictor 1 = Intercept (Female Mean)
# predictor 2 = Sexm (Male Effect)

# Trait 1: Wing Length
mu_female_wing <- samples$beta[,,1][,1] * sd_wing + mean_wing
eff_male_wing  <- samples$beta[,,1][,2] * sd_wing # Differences scale with SD only

# Trait 2: Beak Length
mu_female_beak <- samples$beta[,,2][,1] * sd_beak + mean_beak
eff_male_beak  <- samples$beta[,,2][,2] * sd_beak

# --- B. Loadings ---
loading_wing <- samples$lambda[,1]
loading_beak <- samples$lambda[,2]

# --- C. Variance Components (Unscaled) ---
# Note: h2_psi approach puts variance into lambda scaling
VA_wing   <- samples$lambda[,1]^2 * samples$var_psi_a * sd_wing^2
VA_beak   <- samples$lambda[,2]^2 * samples$var_psi_a * sd_beak^2
VR_wing   <- samples$sd_R[,1]^2 * sd_wing^2
VR_beak   <- samples$sd_R[,2]^2 * sd_beak^2

# --- D. Compile Results ---
posterior_summary <- data.frame(
  Parameter = c(
    "Female_Mean_Wing", "Male_Effect_Wing",
    "Female_Mean_Beak", "Male_Effect_Beak",
    "Loading_Wing", "Loading_Beak",
    "VA_Wing", "VA_Beak",
    "VR_Wing", "VR_Beak",
    "h2_Wing", "h2_Beak",
    "h2_Latent_Size", 
    "h2_Wing_No_Resid, h2_Beak_No_Resid"
  ),
  Mean = c(
    mean(mu_female_wing), mean(eff_male_wing),
    mean(mu_female_beak), mean(eff_male_beak),
    mean(loading_wing), mean(loading_beak),
    mean(VA_wing), mean(VA_beak),
    mean(VR_wing), mean(VR_beak),
    mean(samples$h2_traits[,1]), mean(samples$h2_traits[,2]),
    mean(samples$h2_psi),
    mean(samples$h2_traits_no_residual[,1]), mean(samples$h2_traits_no_residual[,2])
  ),
  SD = c(
    sd(mu_female_wing), sd(eff_male_wing),
    sd(mu_female_beak), sd(eff_male_beak),
    sd(loading_wing), sd(loading_beak),
    sd(VA_wing), sd(VA_beak),
    sd(VR_wing), sd(VR_beak),
    sd(samples$h2_traits[,1]), sd(samples$h2_traits[,2]),
    sd(samples$h2_psi),
    sd(samples$h2_traits_no_residual[,1]), sd(samples$h2_traits_no_residual[,2])
  )
)

# Add Credible Intervals
posterior_summary$Lower_95 <- c(
  quantile(mu_female_wing, 0.025), quantile(eff_male_wing, 0.025),
  quantile(mu_female_beak, 0.025), quantile(eff_male_beak, 0.025),
  quantile(loading_wing, 0.025), quantile(loading_beak, 0.025),
  quantile(VA_wing, 0.025), quantile(VA_beak, 0.025),
  quantile(VR_wing, 0.025), quantile(VR_beak, 0.025),
  quantile(samples$h2_traits[,1], 0.025), quantile(samples$h2_traits[,2], 0.025),
  quantile(samples$h2_psi, 0.025),
  quantile(samples$h2_traits_no_residual[,1], 0.025), quantile(samples$h2_traits_no_residual[,2], 0.025)
)
posterior_summary$Upper_95 <- c(
  quantile(mu_female_wing, 0.975), quantile(eff_male_wing, 0.975),
  quantile(mu_female_beak, 0.975), quantile(eff_male_beak, 0.975),
  quantile(loading_wing, 0.975), quantile(loading_beak, 0.975),
  quantile(VA_wing, 0.975), quantile(VA_beak, 0.975),
  quantile(VR_wing, 0.975), quantile(VR_beak, 0.975),
  quantile(samples$h2_traits[,1], 0.975), quantile(samples$h2_traits[,2], 0.975),
  quantile(samples$h2_psi, 0.975),
  quantile(samples$h2_traits_no_residual[,1], 0.975), quantile(samples$h2_traits_no_residual[,2], 0.975)
)

write.csv(posterior_summary, "Output_LV_Final_Results.csv", row.names = TRUE)
print(posterior_summary)






# # Set up personal library
# personal_lib <- "~/R/library"
# if (!dir.exists(personal_lib)) {
#   dir.create(personal_lib, recursive = TRUE)
# }
# .libPaths(c(personal_lib, .libPaths()))
# 
# # Install missing packages
# required_packages <- c("rstan", "tidyverse", "nadiv")
# for (pkg in required_packages) {
#   if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
#     install.packages(pkg, lib = personal_lib, repos = "https://cloud.r-project.org")
#     library(pkg, character.only = TRUE)
#   }
# }
# 
# 
# 
# script_path <- NULL
# cmd_args <- commandArgs(trailingOnly = FALSE)
# if (any(grepl("^--file=", cmd_args))) {
#   file_arg <- grep("^--file=", cmd_args, value = TRUE)
#   script_path <- sub("^--file=", "", file_arg)
#   script_path <- normalizePath(script_path)
# }
# if (is.null(script_path)) {
#   src_file <- tryCatch(getSrcFilename(function(x) x, full.names = TRUE),
#                        error = function(e) NULL)
#   if (!is.null(src_file) && file.exists(src_file)) script_path <- src_file
# }
# if (is.null(script_path) || !file.exists(script_path)) {
#   script_path <- getwd()
#   message("Warning: Using getwd(): ", script_path)
# } else {
#   script_path <- dirname(script_path)
#   message("Success: Script directory: ", script_path)
# }
# setwd(script_path)
# cat("Working directory set to:", getwd(), "\n")
# 
# 
# 
# # Load data
# traitData <- read.csv("morphology.txt",
#                       sep = ";",
#                       header = TRUE,
#                       na.strings = c("NA", ""),
#                       fileEncoding = "Windows-1252",
#                       stringsAsFactors = FALSE,
#                       colClasses = c(ringnr = "character"))
# 
# pedData <- read.csv("pedigree_data.txt",
#                     sep = " ",
#                     header = TRUE,
#                     na.strings = "NA")
# 
# # Preprocess pedigree
# pedigree_clean <- pedData %>%
#   dplyr::select(ringnr, dam, sire) %>%
#   mutate(
#     across(c(ringnr, dam, sire), ~ as.character(trimws(gsub("_.*$", "", .)))),
#     dam = if_else(is.na(dam) | dam == "" | dam == "NA", NA_character_, dam),
#     sire = if_else(is.na(sire) | sire == "" | sire == "NA", NA_character_, sire)
#   )
# 
# # Add external parents
# all_parents <- unique(c(pedigree_clean$dam, pedigree_clean$sire))
# external_parents <- setdiff(all_parents, c(pedigree_clean$ringnr, NA_character_))
# 
# if (length(external_parents) > 0) {
#   founder_df <- data.frame(
#     ringnr = external_parents,
#     dam = NA_character_,
#     sire = NA_character_
#   )
#   pedigree_full <- bind_rows(founder_df, pedigree_clean) %>%
#     distinct(ringnr, .keep_all = TRUE)
# } else {
#   pedigree_full <- pedigree_clean
# }
# 
# # Prepare pedigree and compute full A matrix
# pedigree_prepped <- prepPed(pedigree_full, gender = NULL, check = TRUE)
# A_nadiv <- makeA(pedigree_prepped)
# 
# 
# 
# 
# 
# 
# 
# 
# # Step A: Filter and Clean Data FIRST
# # We filter for Island + BOTH Traits + Valid Sex
# trait_subset_temp <- traitData %>%
#   filter(sted_r == "hestmannøy", 
#          !is.na(ving_h), 
#          !is.na(nebb_l)) %>%  # Bivariate: Must have both traits
#   
#   mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr)))) %>%
#   
#   # Logic to clean Sex Column (Priority: scriptsex > fieldsex)
#   mutate(
#     raw_sex = if_else(!is.na(scriptsex) & scriptsex != "u", scriptsex, fieldsex)
#   ) %>%
#   mutate(
#     sex = case_when(
#       raw_sex %in% c("m", "pm") ~ "m",
#       raw_sex %in% c("f", "pf") ~ "f",
#       TRUE ~ NA_character_ 
#     )
#   ) %>%
#   filter(!is.na(sex)) # Drop individuals with unknown sex
# 
# # ==============================================================================
# # 4. ALIGN DATA WITH PEDIGREE
# # ==============================================================================
# 
# # Identify valid IDs (intersection of clean data and pedigree)
# valid_ids <- intersect(trait_subset_temp$ringnr, pedigree_prepped$ringnr)
# 
# # Create FINAL data subset sorted by ID
# trait_subset <- trait_subset_temp %>%
#   filter(ringnr %in% valid_ids) %>%
#   arrange(ringnr)
# 
# # ==============================================================================
# # 5. BUILD A-MATRIX (Subsetted)
# # ==============================================================================
# 
# get_ancestors <- function(ids, pedigree) {
#   ancestors <- ids
#   to_process <- ids
#   while (length(to_process) > 0) {
#     parents <- unique(c(pedigree$dam[pedigree$ringnr %in% to_process], 
#                         pedigree$sire[pedigree$ringnr %in% to_process]))
#     parents <- parents[!parents %in% ancestors & !is.na(parents)]
#     ancestors <- c(ancestors, parents)
#     to_process <- parents
#   }
#   return(ancestors)
# }
# 
# # Build small pedigree containing only relevant ancestors
# all_relevant_ids <- get_ancestors(valid_ids, pedigree_prepped)
# 
# smaller_pedigree <- pedigree_prepped %>%
#   filter(ringnr %in% all_relevant_ids) %>%
#   arrange(match(ringnr, all_relevant_ids))
# 
# # Calculate A-matrix
# smaller_ped_prepped <- prepPed(smaller_pedigree, gender = NULL, check = TRUE)
# A_small_nadiv <- makeA(smaller_ped_prepped)
# 
# # Subset A to valid phenotyped individuals AND sort to match data
# A_bivariate <- A_small_nadiv[valid_ids, valid_ids]
# A_bivariate <- A_bivariate[order(rownames(A_bivariate)), order(colnames(A_bivariate))]
# 
# # Safety Check
# stopifnot(all(rownames(A_bivariate) == trait_subset$ringnr))
# cat("Data and Matrix aligned. No =", nrow(trait_subset), "Na =", length(valid_ids), "\n")
# 
# # =================
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Filter for Hestmannøy birds with traits
# phenotyped_ids <- traitData %>%
#   filter(sted_r == "hestmannøy" & !is.na(ving_h) & !is.na(nebb_l)) %>%
#   mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr)))) %>%
#   pull(ringnr) %>%
#   unique()
# 
# valid_phenotyped_ids <- phenotyped_ids[phenotyped_ids %in% pedigree_prepped$ringnr]
# 
# # Get ancestors
# get_ancestors <- function(ids, pedigree) {
#   ancestors <- ids
#   to_process <- ids
#   while (length(to_process) > 0) {
#     parents <- unique(c(pedigree$dam[pedigree$ringnr %in% to_process],
#                         pedigree$sire[pedigree$ringnr %in% to_process]))
#     parents <- parents[!parents %in% ancestors & !is.na(parents)]
#     ancestors <- c(ancestors, parents)
#     to_process <- parents
#   }
#   return(ancestors)
# }
# 
# all_relevant_ids <- get_ancestors(valid_phenotyped_ids, pedigree_prepped)
# 
# # Create smaller pedigree
# smaller_pedigree <- pedigree_prepped %>%
#   filter(ringnr %in% all_relevant_ids) %>%
#   arrange(match(ringnr, all_relevant_ids))
# 
# # Compute smaller A matrix
# smaller_ped_prepped <- prepPed(smaller_pedigree, gender = NULL, check = TRUE)
# A_small_nadiv <- makeA(smaller_ped_prepped)
# 
# # Subset A to phenotyped individuals
# A_hestmannoy <- A_small_nadiv[valid_phenotyped_ids, valid_phenotyped_ids]
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # ... (Data cleaning and A-matrix steps) ...
# 
# # 1. Standardization Note
# # Since you are adding fixed effects, standardizing Y by simple mean/sd is slightly risky 
# # because the mean differs by sex. Ideally, standardise by the residual SD of a linear model.
# # BUT, for this project, stick to simple standardization to keep it moving.
# trait_subset_lv <- trait_subset_lv %>%
#   mutate(
#     ving_h_std = (ving_h - mean(ving_h, na.rm=TRUE)) / sd(ving_h, na.rm=TRUE),
#     nebb_l_std = (nebb_l - mean(nebb_l, na.rm=TRUE)) / sd(nebb_l, na.rm=TRUE)
#   )
# 
# # 2. Prepare Matrices
# Y_lv <- as.matrix(trait_subset_lv[, c("ving_h_std", "nebb_l_std")])
# X_mat_lv <- model.matrix(~ sex, data = trait_subset_lv)
# K_lv <- ncol(X_mat_lv)
# 
# # 3. Bundle Data
# dataset_lv <- list(
#   No = nrow(trait_subset_lv),
#   Nt = 2,
#   Y = Y_lv,
#   X = X_mat_lv,       # NEW
#   K = K_lv,           # NEW
#   animal = animal_lv,
#   Na = length(valid_phenotyped_ids),
#   A = as.matrix(A_hestmannoy)
# )
# 
# # Monitor 'beta' instead of 'mu'
# tomonitor_lv <- c("beta", "lambda", "sd_psi_a", "sd_psi_e", "sd_R", "var_psi_a", "var_psi_e", "h2_psi")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Prepare data for Stan (two traits: ving_h and nebb_l)
# trait_subset_lv <- traitData %>%
#   filter(ringnr %in% valid_phenotyped_ids, sted_r == "hestmannøy", !is.na(ving_h), !is.na(nebb_l)) %>%
#   mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr)))) %>%
#   mutate(
#     ving_h_std = (ving_h - mean(ving_h)) / sd(ving_h),
#     nebb_l_std = (nebb_l - mean(nebb_l)) / sd(nebb_l)
#   )
# 
# No_lv <- nrow(trait_subset_lv)
# Nt_lv <- 2
# Na_lv <- length(valid_phenotyped_ids)
# animal_map_lv <- setNames(seq_len(Na_lv), valid_phenotyped_ids)
# animal_lv <- animal_map_lv[trait_subset_lv$ringnr]
# 
# Y_lv <- as.matrix(trait_subset_lv[, c("ving_h_std", "nebb_l_std")])
# 
# # Bundle data
# dataset_lv <- list(
#   No = No_lv,
#   Nt = Nt_lv,
#   Y = Y_lv,
#   animal = animal_lv,
#   Na = Na_lv,
#   A = as.matrix(A_hestmannoy)
# )
# 
# # Set Stan options
# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())
# 
# # MCMC settings
# nc_lv <- 4
# nw_lv <- 3000
# ni_lv <- 6000
# nt_lv <- 1
# 
# # Parameters to monitor
# tomonitor_lv <- c(
#   "mu",                  # trait means (standardised scale)
#   "lambda",              # loadings 
#   "sd_psi_a",            # LV genetic SD
#   "sd_psi_e",            # LV non-genetic individual SD
#   "sd_R",                # trait-specific residual SDs
#   "var_psi_a",           
#   "var_psi_e",
#   "h2_psi",
#   "h2_traits",           # observed trait heritabilities
#   "h2_traits_no_residual"   # individual-level heritabilities (without residual error)
# )
# 
# # Call Stan
# out_lv <- stan(file = 'latent_variable_stan.stan',
#                data = dataset_lv,
#                pars = tomonitor_lv, # Can be removed to keep all parameters
#                chains = nc_lv, iter = ni_lv, warmup = nw_lv, thin = nt_lv,
#                open_progress = FALSE,
#                refresh = 10,
#                seed = 123)
# 
# 
# # Save results
# saveRDS(out_lv, 'Output_LV_Animal_Model.rds')
# 
# 
# # Write summary to CSV
# summ_lv <- summary(out_lv)$summary
# write.csv(as.data.frame(summ_lv), file = "Output_LV_Summary.csv", row.names = TRUE)
# cat("Summary written to: Output_LV_Summary.csv\n")
# 
# 
# 
# 
# # =============================================================================
# # BACK-TRANSFORMATION 
# 
# fit <- out_lv
# samples <- rstan::extract(fit)
# 
# flip <- samples$lambda[,1] < 0
# samples$lambda[flip, ] <- -samples$lambda[flip, ]
# 
# # Original trait means and SDs (used when standardizing)
# mean_wing <- mean(trait_subset_lv$ving_h, na.rm = TRUE)
# sd_wing   <- sd(trait_subset_lv$ving_h, na.rm = TRUE)
# mean_beak <- mean(trait_subset_lv$nebb_l, na.rm = TRUE)
# sd_beak   <- sd(trait_subset_lv$nebb_l, na.rm = TRUE)
# 
# # 1. Trait means
# mu_wing  <- samples$mu[,1] * sd_wing + mean_wing
# mu_beak  <- samples$mu[,2] * sd_beak + mean_beak
# 
# # 2. Loadings (scale-free)
# loading_wing <- samples$lambda[,1]
# loading_beak <- samples$lambda[,2]
# 
# # 3. Variance components on original scale
# VA_wing   <- samples$lambda[,1]^2 * samples$var_psi_a * sd_wing^2   # additive-genetic
# VA_beak   <- samples$lambda[,2]^2 * samples$var_psi_a * sd_beak^2
# VInd_wing <- samples$lambda[,1]^2 * samples$var_psi_e * sd_wing^2   # individual-level non-genetic
# VInd_beak <- samples$lambda[,2]^2 * samples$var_psi_e * sd_beak^2
# VR_wing   <- samples$sd_R[,1]^2 * sd_wing^2                        # residual
# VR_beak   <- samples$sd_R[,2]^2 * sd_beak^2
# 
# # 4. Final summary table
# posterior_summary <- data.frame(
#   Parameter = c(
#     "mu_wing", "mu_beak",
#     "loading_wing", "loading_beak",
#     "VA_wing", "VA_beak",
#     "VInd_wing", "VInd_beak",
#     "VR_wing", "VR_beak",
#     "h2_wing", "h2_beak",
#     "h2_wing_indiv", "h2_beak_indiv",
#     "h2_psi"
#   ),
#   Mean = c(
#     mean(mu_wing), mean(mu_beak),
#     mean(loading_wing), mean(loading_beak),
#     mean(VA_wing), mean(VA_beak),
#     mean(VInd_wing), mean(VInd_beak),
#     mean(VR_wing), mean(VR_beak),
#     mean(samples$h2_traits[,1]), mean(samples$h2_traits[,2]),
#     mean(samples$h2_traits_no_residual[,1]), mean(samples$h2_traits_no_residual[,2]),
#     mean(samples$h2_psi)
#   ),
#   SD = c(
#     sd(mu_wing), sd(mu_beak),
#     sd(loading_wing), sd(loading_beak),
#     sd(VA_wing), sd(VA_beak),
#     sd(VInd_wing), sd(VInd_beak),
#     sd(VR_wing), sd(VR_beak),
#     sd(samples$h2_traits[,1]), sd(samples$h2_traits[,2]),
#     sd(samples$h2_traits_no_residual[,1]), sd(samples$h2_traits_no_residual[,2]),
#     sd(samples$h2_psi)
#   ),
#   `X2.5%` = c(
#     quantile(mu_wing, 0.025), quantile(mu_beak, 0.025),
#     quantile(loading_wing, 0.025), quantile(loading_beak, 0.025),
#     quantile(VA_wing, 0.025), quantile(VA_beak, 0.025),
#     quantile(VInd_wing, 0.025), quantile(VInd_beak, 0.025),
#     quantile(VR_wing, 0.025), quantile(VR_beak, 0.025),
#     quantile(samples$h2_traits[,1], 0.025), quantile(samples$h2_traits[,2], 0.025),
#     quantile(samples$h2_traits_no_residual[,1], 0.025), quantile(samples$h2_traits_no_residual[,2], 0.025),
#     quantile(samples$h2_psi, 0.025)
#   ),
#   `X97.5%` = c(
#     quantile(mu_wing, 0.975), quantile(mu_beak, 0.975),
#     quantile(loading_wing, 0.975), quantile(loading_beak, 0.975),
#     quantile(VA_wing, 0.975), quantile(VA_beak, 0.975),
#     quantile(VInd_wing, 0.975), quantile(VInd_beak, 0.975),
#     quantile(VR_wing, 0.975), quantile(VR_beak, 0.975),
#     quantile(samples$h2_traits[,1], 0.975), quantile(samples$h2_traits[,2], 0.975),
#     quantile(samples$h2_traits_no_residual[,1], 0.975), quantile(samples$h2_traits_no_residual[,2], 0.975),
#     quantile(samples$h2_psi, 0.975)
#   )
# )
# 
# # Round and save
# posterior_summary[,2:5] <- lapply(posterior_summary[,2:5], round, digits = 4)
# write.csv(posterior_summary, "Output_LV_Final_Results.csv", row.names = FALSE)
# 
# cat("\n=== LATENT VARIABLE MODEL RESULTS (original scale) ===\n")
# print(posterior_summary, row.names = FALSE)
# cat("\nResults saved to: Output_LV_Final_Results.csv\n")