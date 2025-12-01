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

# CORRECT SAFETY CHECK (Sets, not lengths)
if (!all(unique(trait_subset_lv$ringnr) %in% rownames(A_lv))) {
  stop("Error: Data contains individuals not in the A-matrix!")
}
cat("Data and Matrix aligned. No =", nrow(trait_subset_lv), "Na =", nrow(A_lv), "\n")






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
# 4. ALIGN DATA WITH PEDIGREE
# ==============================================================================
valid_ids <- intersect(trait_subset_temp$ringnr, pedigree_prepped$ringnr)

# Create FINAL data subset
# We do NOT sort the data by ID here (it's not strictly necessary for Stan),
# but we MUST ensure the A-matrix is sorted and the mapping matches.
trait_subset_lv <- trait_subset_temp %>%
  filter(ringnr %in% valid_ids)

# Standardize Traits (Store mean/sd for back-transformation)
mean_wing <- mean(trait_subset_lv$ving_h, na.rm = TRUE)
sd_wing   <- sd(trait_subset_lv$ving_h, na.rm = TRUE)
mean_beak <- mean(trait_subset_lv$nebb_l, na.rm = TRUE)
sd_beak   <- sd(trait_subset_lv$nebb_l, na.rm = TRUE)

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

# Subset A to valid phenotyped individuals
A_lv <- A_small_nadiv[valid_ids, valid_ids]

# CRITICAL: Sort A-matrix alphabetically so indices 1:Na are stable
A_lv <- A_lv[order(rownames(A_lv)), order(colnames(A_lv))]

# CORRECT SAFETY CHECK (Sets, not lengths)
if (!all(unique(trait_subset_lv$ringnr) %in% rownames(A_lv))) {
  stop("Error: Data contains individuals not in the A-matrix!")
}
cat("Data and Matrix aligned. No =", nrow(trait_subset_lv), "Na =", nrow(A_lv), "\n")

# ==============================================================================
# 6. PREPARE STAN DATA
# ==============================================================================

# Map animals to integers 1:Na based on the SORTED Matrix
# This ensures Index 1 in Stan = Row 1 in A_lv
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



# ==============================================================================
# DIAGNOSTICS LOGGING
# ==============================================================================
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















## Generate final table if it didn't work 


library(rstan)
library(tidyverse)

# 1. Load Model
fit <- readRDS("Output_LV_Animal_Model.rds")
samples <- rstan::extract(fit)

# 2. Load Data (ROBUST METHOD)
# We force columns to be numeric to handle potential read errors
# ==============================================================================
# 2. Load Data (CORRECTED)
# ==============================================================================

# We use the exact settings that worked for you before
traitData <- read.csv("morphology.txt", 
                      sep = ";", 
                      header = TRUE,
                      na.strings = c("NA", ""), 
                      fileEncoding = "Windows-1252",  # <--- Critical for "ø"
                      stringsAsFactors = FALSE, 
                      colClasses = c(ringnr = "character"))

# Now apply the cleaning
trait_subset <- traitData %>% 
  # Filter for the island (now it should match correctly)
  filter(sted_r == "hestmannøy") %>%
  
  # Handle decimal commas if they exist (safety step)
  mutate(
    ving_h = as.numeric(gsub(",", ".", ving_h)), 
    nebb_l = as.numeric(gsub(",", ".", nebb_l))
  ) %>%
  
  # Filter for existing trait data
  filter(!is.na(ving_h), !is.na(nebb_l)) %>%
  
  # Apply sex logic
  mutate(
    raw_sex = if_else(!is.na(scriptsex) & scriptsex != "u", scriptsex, fieldsex),
    sex = case_when(
      raw_sex %in% c("m", "pm") ~ "m",
      raw_sex %in% c("f", "pf") ~ "f",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(sex))

# Check again
if(nrow(trait_subset) == 0) {
  print("DEBUG: Unique islands found in data:")
  print(unique(traitData$sted_r))
  stop("Error: Data subset is still empty! Check the island spelling above.")
} else {
  cat("Success! Loaded", nrow(trait_subset), "rows for Hestmannøy.\n")
}
# 3. Calculate Scaling Factors
# Check if we actually have data
if(nrow(trait_subset) == 0) stop("Error: Data subset is empty! Check filtering.")

mean_wing <- mean(trait_subset$ving_h, na.rm = TRUE)
sd_wing   <- sd(trait_subset$ving_h, na.rm = TRUE)
mean_beak <- mean(trait_subset$nebb_l, na.rm = TRUE)
sd_beak   <- sd(trait_subset$nebb_l, na.rm = TRUE)

cat("Scaling Factors Found:\n")
cat("Wing: Mean =", mean_wing, " SD =", sd_wing, "\n")
cat("Beak: Mean =", mean_beak, " SD =", sd_beak, "\n")

# 4. Calculate Posteriors (Back-transformed)

# --- A. Fixed Effects ---
# beta[iteration, predictor, trait]
# Predictor 1 = Intercept (Female), Predictor 2 = Sex (Male Effect)
mu_female_wing <- samples$beta[,,1][,1] * sd_wing + mean_wing
eff_male_wing  <- samples$beta[,,1][,2] * sd_wing

mu_female_beak <- samples$beta[,,2][,1] * sd_beak + mean_beak
eff_male_beak  <- samples$beta[,,2][,2] * sd_beak

# --- B. Loadings ---
loading_wing <- samples$lambda[,1]
loading_beak <- samples$lambda[,2]

# --- C. Variance Components (Original Scale) ---
VA_wing   <- samples$lambda[,1]^2 * samples$h2_psi * sd_wing^2
VA_beak   <- samples$lambda[,2]^2 * samples$h2_psi * sd_beak^2
VR_wing   <- samples$sd_R[,1]^2 * sd_wing^2
VR_beak   <- samples$sd_R[,2]^2 * sd_beak^2

# --- D. Create Data Frame ---
values_list <- list(
  mu_female_wing, eff_male_wing,
  mu_female_beak, eff_male_beak,
  loading_wing, loading_beak,
  VA_wing, VA_beak,
  VR_wing, VR_beak,
  samples$h2_traits[,1], samples$h2_traits[,2],
  samples$h2_traits_no_residual[,1], samples$h2_traits_no_residual[,2],
  samples$h2_psi
)

posterior_summary <- data.frame(
  Parameter = c(
    "Female_Mean_Wing", "Male_Effect_Wing",
    "Female_Mean_Beak", "Male_Effect_Beak",
    "Loading_Wing", "Loading_Beak",
    "VA_Wing", "VA_Beak",
    "VR_Wing", "VR_Beak",
    "h2_Wing", "h2_Beak",
    "h2_Wing_Indiv_Only", "h2_Beak_Indiv_Only",
    "h2_Latent_Size"
  ),
  Mean = sapply(values_list, mean, na.rm=TRUE),
  SD   = sapply(values_list, sd, na.rm=TRUE),
  Lower_95 = sapply(values_list, quantile, probs = 0.025, na.rm=TRUE),
  Upper_95 = sapply(values_list, quantile, probs = 0.975, na.rm=TRUE)
)

# Save and Print
write.csv(posterior_summary, "Output_LV_Final_Results_Fixed.csv", row.names = FALSE)
print(posterior_summary)















# --- Prepare Data for Stacked Plot ---

# 1. Calculate Total Variance for each trait
VP_wing <- VA_wing + VInd_wing + VR_wing # Using vectors from previous step
VP_beak <- VA_beak + VInd_beak + VR_beak

# 2. Calculate Means of Proportions (Variance Components)
# We categorize variance into 3 types: Genetic (Latent), Perm Env (Latent), Residual (Specific)
df_stack <- data.frame(
  Trait = rep(c("Wing Length", "Beak Length"), each = 3),
  Component = rep(c("Genetic (Latent Size)", "Perm. Env (Latent Size)", "Residual (Specific)"), 2),
  Value = c(
    mean(VA_wing / VP_wing), mean(VInd_wing / VP_wing), mean(VR_wing / VP_wing),
    mean(VA_beak / VP_beak), mean(VInd_beak / VP_beak), mean(VR_beak / VP_beak)
  )
)

# Force the order of components (Residual at top, Genetics at bottom)
df_stack$Component <- factor(df_stack$Component, 
                             levels = c("Residual (Specific)", "Perm. Env (Latent Size)", "Genetic (Latent Size)"))

# 3. Plot: Stacked Bar Chart
p_stack <- ggplot(df_stack, aes(x = Trait, y = Value, fill = Component)) +
  geom_col(width = 0.5, color = "black") +
  scale_fill_manual(values = c("Genetic (Latent Size)" = "#D55E00", 
                               "Perm. Env (Latent Size)" = "#E69F00", 
                               "Residual (Specific)" = "#999999")) +
  labs(title = "Proportion of Phenotypic Variance Explained",
       subtitle = "How the Latent 'Body Size' factor drives trait variation",
       y = "Proportion of Variance", x = NULL) +
  theme_minimal() +
  theme(legend.position = "right")

# Save and View
ggsave("Plot_LV_Variance_Partition.png", p_stack, width = 7, height = 5)
print(p_stack)













library(ggplot2)
library(bayesplot)
library(gridExtra)

# Ensure you have the scaling factors from the previous step
# mean_wing, sd_wing, mean_beak, sd_beak

# --- Prepare Data for Plotting ---
# Extract the raw posterior samples for the Male Effect
male_eff_wing_post <- samples$beta[,,1][,2] * sd_wing
male_eff_beak_post <- samples$beta[,,2][,2] * sd_beak

# Extract Heritabilities
h2_wing_post <- samples$h2_traits[,1]
h2_beak_post <- samples$h2_traits[,2]
h2_size_post <- samples$h2_psi

# Combine into a long data frame for ggplot
df_sex <- data.frame(
  Value = c(male_eff_wing_post, male_eff_beak_post),
  Trait = rep(c("Wing Length", "Beak Length"), each = length(male_eff_wing_post))
)

df_h2 <- data.frame(
  Value = c(h2_wing_post, h2_beak_post, h2_size_post),
  Parameter = factor(rep(c("h² Wing", "h² Beak", "h² Latent Size"), each = length(h2_wing_post)),
                     levels = c("h² Wing", "h² Beak", "h² Latent Size"))
)

# --- PLOT 1: Sexual Dimorphism (The Contrast) ---
# We use 'scales="free"' because Wing effect is huge (2.4) and Beak is tiny (-0.07)
p1 <- ggplot(df_sex, aes(x = Value, fill = Trait)) +
  geom_density(alpha = 0.7, color = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  facet_wrap(~Trait, scales = "free") + # Crucial for seeing the small beak effect
  scale_fill_manual(values = c("Wing Length" = "#E69F00", "Beak Length" = "#56B4E9")) +
  labs(title = "Posterior Distribution of Sexual Dimorphism (Male - Female)",
       subtitle = "Males have significantly longer wings, but slightly shorter/equal beaks",
       x = "Effect Size (mm)", y = "Posterior Density") +
  theme_minimal() +
  theme(legend.position = "none")

# --- PLOT 2: Heritabilities (The Latent Power) ---
p2 <- ggplot(df_h2, aes(x = Value, fill = Parameter)) +
  geom_density(alpha = 0.6, color = "white") +
  scale_fill_manual(values = c("h² Wing" = "#999999", "h² Beak" = "#999999", "h² Latent Size" = "#D55E00")) +
  labs(title = "Heritability of Traits vs. Latent 'Body Size'",
       subtitle = "The latent size factor captures the shared genetic signal (h² = 0.42)",
       x = "Heritability (h²)", y = "Density") +
  theme_minimal() +
  theme(legend.position = "bottom")

# --- Save and Show ---
ggsave("Plot_Sex_Dimorphism.png", p1, width = 8, height = 4)
ggsave("Plot_Latent_Heritability.png", p2, width = 6, height = 4)

print(p1)
print(p2)





# ==============================================================================
# PLOT: Factor Loadings (Zoomed In)
# ==============================================================================

df_loadings <- data.frame(
  Trait = c("Wing Length", "Beak Length"),
  Mean = c(mean(loading_wing), mean(loading_beak)),
  Lower = c(quantile(loading_wing, 0.025), quantile(loading_beak, 0.025)),
  Upper = c(quantile(loading_wing, 0.975), quantile(loading_beak, 0.975))
)

p_loadings <- ggplot(df_loadings, aes(x = Trait, y = Mean)) +
  geom_point(size = 4, color = "darkblue") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.15, size = 1, color = "darkblue") +
  labs(title = "Factor Loadings on Latent 'Body Size'",
       subtitle = "Standardized loadings showing correlation with the latent factor",
       y = "Loading Coefficient", x = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(r = 10)))

ggsave("Plot_LV_Loadings_Zoomed.png", p_loadings, width = 5, height = 4)
print(p_loadings)








# ==============================================================================
# RE-CALCULATE VARIANCE VECTORS (To fix the "Object not found" error)
# ==============================================================================
# We need the full posterior vectors to calculate means for the plot

# 1. Genetic Variance (Latent)
VA_wing   <- samples$lambda[,1]^2 * samples$h2_psi * sd_wing^2
VA_beak   <- samples$lambda[,2]^2 * samples$h2_psi * sd_beak^2

# 2. Permanent Environment Variance (Latent) -> THIS WAS MISSING
# The latent variable has variance 1. The part that isn't genetic (h2) is Env (1-h2).
VInd_wing <- samples$lambda[,1]^2 * (1 - samples$h2_psi) * sd_wing^2
VInd_beak <- samples$lambda[,2]^2 * (1 - samples$h2_psi) * sd_beak^2

# 3. Residual Variance (Specific)
VR_wing   <- samples$sd_R[,1]^2 * sd_wing^2
VR_beak   <- samples$sd_R[,2]^2 * sd_beak^2

# ==============================================================================
# PLOT: Stacked Variance Partitioning
# ==============================================================================

# Calculate Totals
VP_wing <- VA_wing + VInd_wing + VR_wing
VP_beak <- VA_beak + VInd_beak + VR_beak

# Create Data Frame for Plotting
df_stack <- data.frame(
  Trait = rep(c("Wing Length", "Beak Length"), each = 3),
  Component = rep(c("Genetic (Latent Size)", "Perm. Env (Latent Size)", "Residual (Specific)"), 2),
  Value = c(
    mean(VA_wing / VP_wing), mean(VInd_wing / VP_wing), mean(VR_wing / VP_wing),
    mean(VA_beak / VP_beak), mean(VInd_beak / VP_beak), mean(VR_beak / VP_beak)
  )
)

# Order the components logically (Residual on top)
df_stack$Component <- factor(df_stack$Component, 
                             levels = c("Residual (Specific)", "Perm. Env (Latent Size)", "Genetic (Latent Size)"))

# Plot
p_stack <- ggplot(df_stack, aes(x = Trait, y = Value, fill = Component)) +
  geom_col(width = 0.5, color = "black", alpha = 0.8) +
  scale_fill_manual(values = c("Genetic (Latent Size)" = "#D55E00", 
                               "Perm. Env (Latent Size)" = "#E69F00", 
                               "Residual (Specific)" = "#999999")) +
  labs(title = "Variance Partitioning (Latent Variable Model)",
       subtitle = "Genetic and Environmental effects of the shared 'Size' factor",
       y = "Proportion of Phenotypic Variance", x = NULL) +
  theme_minimal() +
  theme(legend.position = "right", legend.title = element_blank())

ggsave("Plot_LV_Variance_Partition.png", p_stack, width = 7, height = 5)
print(p_stack)























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