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
A_lv <- A_small_nadiv[valid_ids, valid_ids]

# Safety check
if (!all(unique(trait_subset_lv$ringnr) %in% rownames(A_lv))) {
  stop("Error: Data contains individuals not in the A-matrix!")
}
cat("Data and Matrix aligned. No =", nrow(trait_subset_lv), "Na =", nrow(A_lv), "\n")


# PREPARE STAN DATA

# Map animals to integers 1:Na based on the sorted matrix
matrix_ids <- rownames(A_lv)
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
  A = as.matrix(A_lv)
)




# DIAGNOSTICS LOGGING

sink("Output_LV_PreRun_Diagnostics.txt") # Redirect output to file

cat("=== DATA & MATRIX DIAGNOSTICS ===\n")
cat("TimeStamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("1. Data Dimensions:\n")
cat("   Observations (No):", nrow(trait_subset_lv), "\n")
cat("   Individuals (Na): ", nrow(A_lv), "\n")
cat("   Traits (Nt):      ", 2, "\n\n")

cat("2. Fixed Effects (X Matrix):\n")
cat("   Predictors (K):   ", ncol(X_mat), "\n")
cat("   Colnames:         ", paste(colnames(X_mat), collapse=", "), "\n\n")

cat("3. ID Alignment Check:\n")
# Check if all data IDs are in Matrix
ids_in_matrix <- unique(trait_subset_lv$ringnr) %in% rownames(A_lv)
cat("   All Data IDs in Matrix? ", all(ids_in_matrix), "\n")
if(!all(ids_in_matrix)) {
  cat("   MISSING IDs: ", head(unique(trait_subset_lv$ringnr)[!ids_in_matrix]), "...\n")
}

# Check Matrix Properties
cat("\n4. Matrix Properties:\n")
cat("   Symmetric?        ", isSymmetric(A_lv), "\n")
cat("   Positive Definite?", all(eigen(A_lv, only.values=TRUE)$values > -1e-6), "\n")
cat("   Mean Diagonal:    ", mean(diag(A_lv)), "(Expected ~1.0)\n")

sink() # Turn off redirection (restore console output)




# RUN STAN
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
               #pars = tomonitor_lv,
               chains = nc, iter = ni, warmup = nw, thin = nt,
               open_progress = FALSE,
               seed = 123)

saveRDS(out_lv, 'Output_LV_Animal_Model.rds')
summ_lv <- summary(out_lv)$summary
write.csv(as.data.frame(summ_lv), file = "Output_LV_Summary.csv", row.names = TRUE)


# 8. BACK-TRANSFORMATION 
samples <- rstan::extract(out_lv)

# --- Fixed Effects (Means) ---
# beta[sample, predictor, trait]
# predictor 1 = Intercept (female mean)
# predictor 2 = Sexm (male effect)

# Trait 1: wing length
mu_female_wing <- samples$beta[,,1][,1] * sd_wing + mean_wing
eff_male_wing  <- samples$beta[,,1][,2] * sd_wing # Differences scale with SD only

# Trait 2: beak length
mu_female_beak <- samples$beta[,,2][,1] * sd_beak + mean_beak
eff_male_beak  <- samples$beta[,,2][,2] * sd_beak

# --- Loadings ---
loading_wing <- samples$lambda[,1]
loading_beak <- samples$lambda[,2]

# --- Variance Components (Unscaled) ---
# Note: h2_psi approach puts variance into lambda scaling
VA_wing   <- samples$lambda[,1]^2 * samples$var_psi_a * sd_wing^2
VA_beak   <- samples$lambda[,2]^2 * samples$var_psi_a * sd_beak^2
VR_wing   <- samples$sd_R[,1]^2 * sd_wing^2
VR_beak   <- samples$sd_R[,2]^2 * sd_beak^2

# Compile
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

# Add credible intervals
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
