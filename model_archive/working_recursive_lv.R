# ==============================================================================
# 1. SETUP & LIBRARIES
# ==============================================================================
personal_lib <- "~/R/library"
if (!dir.exists(personal_lib)) dir.create(personal_lib, recursive = TRUE)
.libPaths(c(personal_lib, .libPaths()))

required_packages <- c("rstan", "tidyverse", "nadiv", "MCMCglmm") 
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, lib = personal_lib, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# ==============================================================================
# 2. LOAD & CLEAN DATA
# ==============================================================================
# Pedigree
pedData <- read.csv("pedigree_data.txt", sep = " ", header = TRUE, na.strings = "NA")
pedigree_clean <- pedData %>%
  dplyr::select(ringnr, dam, sire) %>%
  mutate(across(c(ringnr, dam, sire), ~ as.character(trimws(gsub("_.*$", "", .)))),
         dam = if_else(is.na(dam) | dam == "" | dam == "NA", NA_character_, dam),
         sire = if_else(is.na(sire) | sire == "" | sire == "NA", NA_character_, sire))

# External Parents
all_parents <- unique(c(pedigree_clean$dam, pedigree_clean$sire))
external_parents <- setdiff(all_parents, c(pedigree_clean$ringnr, NA_character_))
if (length(external_parents) > 0) {
  founder_df <- data.frame(ringnr = external_parents, dam = NA_character_, sire = NA_character_)
  pedigree_full <- bind_rows(founder_df, pedigree_clean) %>% distinct(ringnr, .keep_all = TRUE)
} else {
  pedigree_full <- pedigree_clean
}

pedigree_prepped <- prepPed(pedigree_full, gender = NULL, check = TRUE)




# Trait Data
traitData <- read.csv("morphology.txt", sep = ";", header = TRUE, 
                      na.strings = c("NA", ""), fileEncoding = "Windows-1252", 
                      stringsAsFactors = FALSE, colClasses = c(ringnr = "character"))

# # Filter & Standardize
# trait_subset_temp <- traitData %>%
#   filter(sted_r == "hestmannøy", 
#          !is.na(ving_h), !is.na(nebb_l), !is.na(tars_h), !is.na(nebb_h), !is.na(vekt)) %>%
#   mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr)))) %>%
#   mutate(raw_sex = if_else(!is.na(scriptsex) & scriptsex != "u", scriptsex, fieldsex)) %>%
#   mutate(sex = case_when(raw_sex %in% c("m", "pm") ~ "m", raw_sex %in% c("f", "pf") ~ "f", TRUE ~ NA_character_)) %>%
#   filter(!is.na(sex))


# ==============================================================================
# 2. LOAD & CLEAN DATA (UPDATED)
# ==============================================================================

# ==============================================================================
# 2. LOAD & CLEAN DATA (OPTIMIZED ORDER)
# ==============================================================================

# A. LOAD & LIGHT FORMATTING (Crucial for matching)
trait_raw <- traitData %>%
  filter(sted_r == "hestmannøy", 
         !is.na(ving_h), !is.na(nebb_l), !is.na(tars_h), !is.na(nebb_h), !is.na(vekt)) %>%
  # Fix Ring Numbers NOW so they match the pedigree
  mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr))))

# B. ALIGN WITH PEDIGREE (The "Gatekeeper" Step)
# We do this EARLY so we only clean/analyze relevant birds
valid_ids <- intersect(trait_raw$ringnr, pedigree_prepped$ringnr)
trait_subset_temp <- trait_raw %>% filter(ringnr %in% valid_ids)

cat("Birds in Raw File:    ", nrow(trait_raw), "\n")
cat("Birds with Pedigree:  ", nrow(trait_subset_temp), "\n")

# C. DEEP CLEANING (Sex, Age, Bad IDs) on the SUBSET
# 1. Define Known Bad IDs (Impossible values)
bad_ids <- c("8N27558", "8N42527", "8N13933", "8M71874", "8N87712", "8L89527")

trait_subset_temp <- trait_subset_temp %>%
  # Remove Known Bad IDs
  filter(!ringnr %in% bad_ids) %>%
  
  # SEX LOGIC
  mutate(raw_sex = if_else(!is.na(scriptsex) & scriptsex != "u", scriptsex, fieldsex)) %>%
  mutate(sex = case_when(
    raw_sex %in% c("m", "pm") ~ "m", 
    raw_sex %in% c("f", "pf") ~ "f", 
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(sex)) %>%
  
  # AGE LOGIC
  mutate(
    s_age = trimws(as.character(scriptage)),
    f_age = trimws(as.character(fieldage)),
    
    age_class = case_when(
      grepl("^1K", s_age) ~ "juvenile",
      grepl("^[2-9]K", s_age) | grepl("^1[0-9]K", s_age) ~ "adult",
      f_age %in% c("j", "pj") ~ "juvenile",
      f_age %in% c("a", "pa") ~ "adult",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(age_class))

# D. STATISTICAL OUTLIER REMOVAL (Residual Check)
# Now the residuals are calculated strictly on the birds entering the model
outlier_ids_to_kill <- c()
traits_to_check <- c("ving_h", "tars_h", "nebb_l", "nebb_h", "vekt")

for(t in traits_to_check) {
  # Fit Model
  f <- as.formula(paste(t, "~ sex + age_class"))
  fit <- lm(f, data = trait_subset_temp, na.action = na.exclude)
  
  # Check Residuals
  resids <- rstudent(fit)
  cutoff <- if(t == "vekt") 6.0 else 4.0 # Loose for weight, strict for bones
  
  bad_rows <- which(abs(resids) > cutoff)
  
  if(length(bad_rows) > 0) {
    new_bad_ids <- trait_subset_temp$ringnr[bad_rows]
    outlier_ids_to_kill <- c(outlier_ids_to_kill, new_bad_ids)
    cat("Trait:", t, "- Removing", length(new_bad_ids), "IDs\n")
  }
}

# E. FINAL DATASET
trait_subset_lv <- trait_subset_temp %>%
  filter(!ringnr %in% unique(outlier_ids_to_kill))

cat("Final Clean N: ", nrow(trait_subset_lv), "\n")


# Load necessary library if not already loaded


# ==============================================================================
# OUTLIER DETECTION: RESIDUAL CHECK (Post-Age/Sex Correction)
# ==============================================================================


# # 2. Setup storage for outliers
# traits <- c("ving_h", "tars_h", "nebb_l", "nebb_h", "vekt")
# outlier_list <- data.frame()
# 
# 
# 
# 
# ids_to_remove <- all_outliers %>%
#   filter(
#     (Trait %in% c("ving_h", "tars_h", "nebb_l", "nebb_h") & abs(Residual_Score) > 4) |
#       (Trait == "vekt" & abs(Residual_Score) > 6)
#   ) %>%
#   pull(ringnr) %>%
#   unique()
# 
# cat("Removing", length(ids_to_remove), "birds based on residual checks.\n")
# 
# # Apply the filter
# trait_subset_lv <- trait_subset_clean %>%
#   filter(!ringnr %in% ids_to_remove)
# 
# 
# 
# # 4. Show the "Real" Outliers
# if(nrow(outlier_list) > 0) {
#   cat("\nPOTENTIAL OUTLIERS DETECTED (Residual > 4):\n")
#   print(outlier_list %>% arrange(desc(abs(Residual_Score))))
# } else {
#   cat("\nNo extreme outliers found (all residuals < 4).\n")
# }
# 
# # 5. Visual Check (Optional but Recommended)
# # Plot the residuals to see if they look "Normal" (Bell curve)
# par(mfrow=c(2,3))
# for(t in traits) {
#   f <- as.formula(paste(t, "~ sex + age_class"))
#   fit <- lm(f, data = trait_subset_clean, na.action = na.exclude)
#   hist(rstudent(fit), main=paste("Resids:", t), xlab="SD from Mean", breaks=30, col="lightblue")
#   abline(v=c(-4, 4), col="red", lty=2)
# }





# Z-score standardization
trait_subset_lv <- trait_subset_lv %>%
  mutate(
    ving_h_std = scale(ving_h),
    nebb_l_std = scale(nebb_l),
    tars_h_std = scale(tars_h),
    nebb_h_std = scale(nebb_h),
    vekt_std   = scale(vekt)
  )



# 1. Define the "Kill List" of specific IDs identified as errors
# (Beak errors, Tarsus error, Molting/Tiny wings)
# bad_ids <- c(
#   "8N27558", "8N42527", # The 23mm beaks (Z > 11)
#   "8N13933",            # The 12mm tarsus (Z < -8)
#   "8M71874", "8N87712", # The 58mm wings (Z < -7)
#   "8L89527"             # Very small overall (Z < -5), likely juv/error
# )

# # 2. Filter them out
# trait_subset_clean <- trait_subset_lv %>%
#   filter(!ringnr %in% bad_ids)
# 

# trait_subset_lv <- trait_subset_clean


# z_threshold <- 3.5 
# 
# outliers_df <- trait_subset_lv %>%
#   # Calculate Z-scores for all 5 traits
#   mutate(
#     z_ving = (ving_h - mean(ving_h, na.rm=TRUE)) / sd(ving_h, na.rm=TRUE),
#     z_nebb_l = (nebb_l - mean(nebb_l, na.rm=TRUE)) / sd(nebb_l, na.rm=TRUE),
#     z_tars = (tars_h - mean(tars_h, na.rm=TRUE)) / sd(tars_h, na.rm=TRUE),
#     z_nebb_h = (nebb_h - mean(nebb_h, na.rm=TRUE)) / sd(nebb_h, na.rm=TRUE),
#     z_vekt = (vekt - mean(vekt, na.rm=TRUE)) / sd(vekt, na.rm=TRUE)
#   ) %>%
#   # Filter rows where ANY trait exceeds the threshold
#   filter(
#     abs(z_ving) > z_threshold |
#       abs(z_nebb_l) > z_threshold |
#       abs(z_tars) > z_threshold |
#       abs(z_nebb_h) > z_threshold |
#       abs(z_vekt) > z_threshold
#   ) %>%
#   # Select just the ID and the raw values to inspect
#   select(ringnr, ving_h, nebb_l, tars_h, nebb_h, vekt, z_ving, z_nebb_l, z_tars, z_nebb_h, z_vekt)
# 
# # Print the list
# print(outliers_df)
# 
# 
# library(tidyverse)
# 
# # Define your trait columns (using the raw names from your code)
# traits <- c("ving_h_std", "tars_h_std", "nebb_l_std", "nebb_h_std", "vekt_std")
# 
# # Reshape data to long format for easy plotting
# long_data <- trait_subset_lv %>%
#   pivot_longer(cols = all_of(traits), names_to = "Trait", values_to = "Value")
# 
# # Create the plot
# ggplot(long_data, aes(x = Trait, y = Value)) +
#   geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 3) +
#   facet_wrap(~Trait, scales = "free") +  # "free" scales are crucial since units differ
#   theme_bw() +
#   labs(title = "Outlier Detection: Raw Trait Values", y = "Measurement")


cat("Data Loaded. N Rows:", nrow(trait_subset_lv), "\n")


# ==============================================================================
# 3. PEDIGREE PROCESSING (Robust String-to-Integer Mapping)
# ==============================================================================
pedigree_working <- pedigree_prepped
colnames(pedigree_working) <- c("id", "dam", "sire") 
pedigree_working$id   <- as.character(pedigree_working$id)
pedigree_working$dam  <- as.character(pedigree_working$dam)
pedigree_working$sire <- as.character(pedigree_working$sire)

# Prune
ped_pruned <- MCMCglmm::prunePed(pedigree_working, keep = trait_subset_lv$ringnr)

# Inverse A 
A_inv_obj <- MCMCglmm::inverseA(ped_pruned, nodes = "ALL", scale = FALSE)

# Recover Ordered IDs
ordered_ids <- rownames(A_inv_obj$Ainv)

# Create Ordered Pedigree DF
# We assume columns 2 and 3 might be STRINGS now
ped_ordered <- data.frame(
  id = ordered_ids,
  dam = as.character(A_inv_obj$pedigree[, 2]), 
  sire = as.character(A_inv_obj$pedigree[, 3]),
  stringsAsFactors = FALSE
)

dii_vector <- A_inv_obj$dii
Na <- nrow(ped_ordered)

# *** CRITICAL FIX: MAP PARENT STRINGS TO ROW INTEGERS ***
# We look up the dam's name in the 'id' column to find her Row Number.
dam_idx_mapped  <- match(ped_ordered$dam, ped_ordered$id)
sire_idx_mapped <- match(ped_ordered$sire, ped_ordered$id)

# Handle Missing Parents
# match() returns NA if the parent is NA or not found.
# We map these NAs to the dummy index (Na + 1).
dam_idx_mapped[is.na(dam_idx_mapped)]   <- Na + 1
sire_idx_mapped[is.na(sire_idx_mapped)] <- Na + 1

# Create Map
id_map <- setNames(1:Na, ped_ordered$id)

# ==============================================================================
# 4. DATA MAPPING & SAFETY CHECK
# ==============================================================================
animal_idx <- id_map[trait_subset_lv$ringnr]
missing_mask <- is.na(animal_idx)

if (sum(missing_mask) > 0) {
  cat("WARNING: Removing", sum(missing_mask), "rows. ID mismatch persists.\n")
  trait_subset_lv <- trait_subset_lv[!missing_mask, ]
  animal_idx <- id_map[trait_subset_lv$ringnr]
}

# Final Matrices
Y_mat <- as.matrix(trait_subset_lv[, c(
                                        "tars_h_std",
                                        "ving_h_std",
                                        "nebb_l_std", 
                                        "nebb_h_std",
                                        "vekt_std"
                                       )])
# X_mat <- model.matrix(~ sex, data = trait_subset_lv)
X_mat <- model.matrix(~ sex + age_class, data = trait_subset_lv)


cat("==============================================\n")
cat("       FINAL DATA VERIFICATION CHECK          \n")
cat("==============================================\n")
cat("1. Total Data Points (N_obs):     ", nrow(trait_subset_lv), "\n")
cat("   (This is the number of rows Stan analyzes)\n\n")

cat("2. Unique Phenotyped Birds:       ", length(unique(trait_subset_lv$ringnr)), "\n")
cat("   (This is the actual number of birds you caught)\n\n")

cat("3. Total Pedigree Size (Na):      ", Na, "\n")
cat("   (This is Phenotyped Birds + Ancestors)\n")
cat("==============================================\n")

# ==============================================================================
# 5. STAN EXECUTION
# ==============================================================================
# dataset_lv_final <- list(
#   No = nrow(trait_subset_lv),
#   Nt = 5, 
#   Nlv = 2,
#   Y = Y_mat,
#   X = X_mat,
#   K = ncol(X_mat),
#   animal = as.array(as.integer(animal_idx)),
#   Na = Na,
#   
#   # USE THE MAPPED INDICES
#   dam = as.array(as.integer(dam_idx_mapped)), 
#   sire = as.array(as.integer(sire_idx_mapped)),
#   
#   dii = as.array(dii_vector)
# )


dataset_lv_final <- list(
  No = nrow(trait_subset_lv),
  Nt = 5, 
  Nlv = 2,
  Y = Y_mat,
  X = X_mat,
  K = ncol(X_mat),
  animal = as.array(as.integer(animal_idx)),
  Na = Na,
  dam = as.array(as.integer(dam_idx_mapped)), 
  sire = as.array(as.integer(sire_idx_mapped)),
  dii = as.array(dii_vector)
)


# Check raw correlations among your 5 traits
# cor(dataset_lv_final$Y)



rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# 
# tomonitor_lv <- c(
#   "beta", "lambda", 
#   "var_G1", "var_PE1", "var_TE1", "var_LV1_total",
#   "var_G2", "var_PE2", "var_TE2", "var_LV2_total",
#   "sd_R", "h2_traits", "lp__"
# )


# tomonitor_lv <- c(
#   "beta", "lambda", 
#   "var_G1", "var_PE1", "var_TE1", "var_LV1_total",
#   "sd_R", "h2_traits", "lp__"
# )

tomonitor_lv <- c(
  # --- The Basics ---
  "beta",         # Fixed effects
  "Lambda",       # Factor loadings (reconstructed)
  "sd_R",         # Residual standard deviations
  
  # --- The New Key Parameter ---
  # "h2_psi",       # Heritability of the Latent Factor (0 to 1)
  "h2_lv",
  "rho",
  "LV",
  "h2_trait",
  "pe2_trait",
  "var_G_trait",
  "var_PE_trait",
  "var_R_trait",
  "var_P_total",
  # "te2_psi",
  
  # --- Outputs & Diagnostics ---
  # "h2_trait_total",            # Heritability of the observed traits
  # "diag_var_genetic", # Check this is ~1.0 (Proves recursion worked)
  # "diag_var_total",   # Check this is ~1.0 (Proves factor scaling worked)
  
  "lp__"          # Log-probability (for checking convergence)
)


# 1. Verify Sort Order
# If this errors, use MCMCglmm::orderPed() to fix pedigree_clean BEFORE processing
if (any(dam_idx_mapped[dam_idx_mapped <= Na] >= (1:Na)[dam_idx_mapped <= Na])) {
  stop("CRITICAL ERROR: A Dam appears in the list AFTER her child. Recursion will fail.")
}
if (any(sire_idx_mapped[sire_idx_mapped <= Na] >= (1:Na)[sire_idx_mapped <= Na])) {
  stop("CRITICAL ERROR: A Sire appears in the list AFTER his child.")
}

# 2. Initialization Function
# Use this to prevent the "Zero Variance" trap
# init_fn <- function() {
#   list(
#     h2_psi = runif(1, 0.2, 0.8),         # Start with a healthy heritability
#     lambda_anchor = runif(1, 0.5, 1.0),  # Start positive
#     lambda_free = rnorm(dataset_lv_final$Nt - 1, 0, 0.5),
#     sd_R = runif(dataset_lv_final$Nt, 0.5, 1.0),
#     w_a = rnorm(dataset_lv_final$Na, 0, 0.1),
#     w_pe = rnorm(dataset_lv_final$Na, 0, 0.1)
#   )
# }

# 
# init_fn_2factor <- function() {
#   list(
#     h2_psi = runif(2, 0.3, 0.7),
#     
#     L1_anchor = runif(1, 0.5, 1.0),
#     L1_free = rnorm(4, 0, 0.2), # 4 traits left
#     
#     L2_anchor = runif(1, 0.5, 1.0),
#     L2_free = rnorm(3, 0, 0.2), # 3 traits left (5 - 1 fixed - 1 anchor)
#     
#     w_a = matrix(rnorm(2 * dataset_lv_final$Na, 0, 0.1), nrow=2),
#     w_pe = matrix(rnorm(2 * dataset_lv_final$Na, 0, 0.1), nrow=2),
#     
#     sd_R = runif(5, 0.5, 1.0)
#   )
# }


init_fn_shrinkage <- function() {
  # Nlv = 2, Nt = 5
  
  # 1. LOADINGS
  # We initialize Trait 1 (Wing?) on Factor 1 to be POSITIVE (0.5 to 1.5).
  # This encourages the sampler to orient Factor 1 as "Size" (Positive correlation).
  L_init <- matrix(rnorm(dataset_lv_final$Nt * dataset_lv_final$Nlv, 0, 0.1), 
                   nrow = dataset_lv_final$Nt, 
                   ncol = dataset_lv_final$Nlv)
  
  L_init[1, 1] <- runif(1, 0.8, 1.1) # Nudge Factor 1 to be positive for Trait 1
  L_init[1, 2] <- runif(1, -0.1, 0.1) # Start Factor 2 small for Trait 1
  
  list(
    # --- Anchors ---
    lambda_raw = L_init, 
    
    # --- Variances ---
    rho = runif(1, 0.4, 0.6),        # Start in the middle (don't start at 0!)
    h2_lv = runif(2, 0.4, 0.6),      # Start split evenly between G and PE
    sd_R = runif(dataset_lv_final$Nt, 0.8, 1.2), # Data is Z-scored, so Resids ~ 1.0
    
    # --- Fixed Effects ---
    beta = matrix(rnorm(dataset_lv_final$K * dataset_lv_final$Nt, 0, 1),
                  nrow = dataset_lv_final$K),
    
    # --- Latent Vectors ---
    # Start these near zero so the first few steps aren't crazy
    w_a = matrix(rnorm(dataset_lv_final$Nlv * dataset_lv_final$Na, 0, 0.1),
                 nrow = dataset_lv_final$Nlv),
    w_pe = matrix(rnorm(dataset_lv_final$Nlv * dataset_lv_final$Na, 0, 0.1),
                  nrow = dataset_lv_final$Nlv)
  )
}

# 3. Run
out_lv <- stan(
  file = 'good_shrink_prior.stan', 
  data = dataset_lv_final,
  init = init_fn_shrinkage,       # <--- Don't forget this!
  chains = 4, 
  pars = tomonitor_lv,
  control = list(adapt_delta = 0.9),
  iter = 8000,          
  warmup = 4000
)

# out_lv <- stan(
#   file = '1_lv_model_test.stan', 
#   data = dataset_lv_final,
#   pars = tomonitor_lv,
#   chains = 4, iter = 10000, warmup = 5000, thin = 1,
#   seed = 123
#   # ,
#   # control = list(
#   #   adapt_delta = 0.99,   # Forces smaller steps (Fixes divergences)
#   #   max_treedepth = 15    # Allows more steps (Fixes efficiency)
#   # )
# )
# 
saveRDS(out_lv, 'Output_shrink_prior.rds')




# 
# # Load the output if not already in memory
out_lv <- readRDS('Output_shrink_prior_long.rds')
# 
# # Print Summary of Key Parameters
print(out_lv,
      pars = c("Lambda",       # Loadings on the factor
               # --- The Basics ---
               "beta",         # Fixed effects
               "sd_R",         # Residual standard deviations
               
               # --- The New Key Parameter ---
               # "h2_psi",       # Heritability of the Latent Factor (0 to 1)
               "h2_lv",
               "rho",
               "h2_trait",
               "pe2_trait",
               "var_G_trait",
               "var_PE_trait",
               "var_R_trait",
               "var_P_total",
               # "te2_psi",
               
               # --- Outputs & Diagnostics ---
               # "h2_trait_total",            # Heritability of the observed traits
               # "diag_var_genetic", # Check this is ~1.0 (Proves recursion worked)
               # "diag_var_total",   # Check this is ~1.0 (Proves factor scaling worked)
               
               "lp__"    ),      # Log-probability (for checking convergence)
      probs = c(0.025, 0.975),
      digits_summary = 3)


library(rstan)
# 
rstan::traceplot(out_lv, pars = c("Lambda[1,1]", "Lambda[1,2]", "h2_trait[1]", "h2_trait[2]"), inc_warmup = FALSE)
# rstan::traceplot(out_lv, pars = c("te2_psi[1]", "te2_psi[2]", "Lambda[1,2]", "Lambda[2,2]", "Lambda[3,2]", "Lambda[4,2]", "Lambda[5,2]"), inc_warmup = FALSE)
# rstan::traceplot(out_lv, pars = c("sd_R[1]", "sd_R[2]", "sd_R[3]", "sd_R[4]", "sd_R[5]", "h2_trait_total"), inc_warmup = FALSE)



library(rstan)

# 1. Extract the 3D array of draws [iterations, chains, parameters]
draws <- as.array(out_lv)

# 2. Loop through the 4 chains
for (chain in 1:4) {
  
  # Check the mean of Beak Height on Factor 2 (Lambda[4,2]) for this specific chain
  chain_mean <- mean(draws[, chain, "Lambda[4,2]"])
  
  # If the chain found the "negative" mode, we flip it!
  if (chain_mean < 0) {
    cat("Flipping signs for Factor 2 in Chain", chain, "\n")
    
    # Multiply all Factor 2 loadings by -1 for this chain
    draws[, chain, "Lambda[1,2]"] <- draws[, chain, "Lambda[1,2]"] * -1
    draws[, chain, "Lambda[2,2]"] <- draws[, chain, "Lambda[2,2]"] * -1
    draws[, chain, "Lambda[3,2]"] <- draws[, chain, "Lambda[3,2]"] * -1
    draws[, chain, "Lambda[4,2]"] <- draws[, chain, "Lambda[4,2]"] * -1
    draws[, chain, "Lambda[5,2]"] <- draws[, chain, "Lambda[5,2]"] * -1
    
    # CRITICAL NOTE: If you saved "LV" (Latent Variables) in this run, 
    # you technically need to multiply LV[, 2] by -1 for this chain as well, 
    # because Predicted = LV * Lambda. If you flip one, you must flip the other!
  }
}





# 3. Create a new, corrected summary table for just the Lambdas
corrected_summary <- rstan::monitor(draws, warmup = 0, print = FALSE)

summary_df <- as.data.frame(corrected_summary)

# 5. View the fixed results!
lambdas_only <- summary_df[grep("Lambda", rownames(summary_df)), ]

# Print the clean, converged table
print(lambdas_only[, c("mean", "sd", "2.5%", "97.5%", "n_eff", "Rhat")], digits = 3)





# 1. Get the summary table from your CORRECTED draws array

# 2. Extract the posterior means for the key matrices
beta_means   <- summary_df[grep("^beta\\[", rownames(summary_df)), "mean"]
Lambda_means <- summary_df[grep("^Lambda\\[", rownames(summary_df)), "mean"]
LV_means     <- summary_df[grep("^LV\\[", rownames(summary_df)), "mean"]

# 3. Reshape them back into matrices
# Ensure these match your actual dimensions! 
# dataset_lv_final$K = number of fixed effects
B_hat  <- matrix(beta_means,   nrow = dataset_lv_final$K,  ncol = dataset_lv_final$Nt)
L_hat  <- matrix(Lambda_means, nrow = dataset_lv_final$Nt, ncol = dataset_lv_final$Nlv)
LV_hat <- matrix(LV_means,     nrow = dataset_lv_final$Na, ncol = dataset_lv_final$Nlv)

# 4. Map the Latent Variables (Na) to the Observations (No)
LV_obs <- LV_hat[dataset_lv_final$animal, ]

# 5. Calculate Predicted Values (Mu)
# Mu = Fixed Effects + (Factor Scores * Loadings)
Mu_hat <- (dataset_lv_final$X %*% B_hat) + (LV_obs %*% t(L_hat))

# 6. Calculate Residuals
# Residuals = Observed - Predicted
Residuals <- dataset_lv_final$Y - Mu_hat




cat("Residual Correlation Matrix:\n")
resid_cor_matrix <- cor(Residuals)
print(round(resid_cor_matrix, 2))




trait_names <- c("Tarsus", "Wing", "Beak L", "Beak H", "Mass")
par(mfrow = c(2, 3)) # Sets up a grid for the plots

for(i in 1:5) {
  qqnorm(Residuals[, i], main = paste("QQ Plot:", trait_names[i]))
  qqline(Residuals[, i], col = "red", lwd = 2)
}




par(mfrow = c(2, 3))

for(i in 1:5) {
  plot(Mu_hat[, i], Residuals[, i], 
       main = paste("Resid vs Fit:", trait_names[i]),
       xlab = "Predicted Value", 
       ylab = "Residual",
       pch = 16, col = rgb(0,0,0,0.2)) # transparent points for dense data
  abline(h = 0, col = "red", lwd = 2)
}

# Residual analysis 
# 
# # 1. Extract the Posterior Means (Point Estimates)
# # We use the means to get a "consensus" prediction
# # (You could also do this for every sample, but that takes forever)
# post_beta   <- summary(out_shrink, pars = "beta")$summary[, "mean"]
# post_Lambda <- summary(out_shrink, pars = "Lambda")$summary[, "mean"]
# post_LV     <- summary(out_shrink, pars = "LV")$summary[, "mean"]
# 
# # Reshape them back into matrices (Stan flattens them)
# # Check dimensions: beta is K x Nt, Lambda is Nt x 2, LV is Na x 2
# B_hat <- matrix(post_beta, nrow = dataset_lv_final$K, ncol = dataset_lv_final$Nt)
# L_hat <- matrix(post_Lambda, nrow = dataset_lv_final$Nt, ncol = dataset_lv_final$Nlv)
# LV_hat <- matrix(post_LV, nrow = dataset_lv_final$Na, ncol = dataset_lv_final$Nlv)
# 
# # 2. Reconstruct the Linear Predictor (Mu)
# # Formula: Mu = X * Beta + LV * Lambda'
# # Note: We need to map LV from 'Animal' ID to 'Observation' ID
# LV_obs <- LV_hat[dataset_lv_final$animal, ] 
# 
# Mu_hat <- (dataset_lv_final$X %*% B_hat) + (LV_obs %*% t(L_hat))
# 
# # 3. Calculate Residuals
# # Residual = Observed - Predicted
# Residuals <- dataset_lv_final$Y - Mu_hat
# 
# # =========================================================
# # ANALYSIS 1: Correlation of Residuals (The TE Check)
# # =========================================================
# # If these correlations are high, it means your model is missing 
# # a "Temporary Environment" correlation layer.
# resid_cor_matrix <- cor(Residuals)
# print(round(resid_cor_matrix, 2))
# 
# # Interpretation: 
# # If off-diagonals are close to 0 (e.g., < 0.2), your current model is fine.
# # If you see 0.5+ correlations, you need that TE layer.
# 
# # =========================================================
# # ANALYSIS 2: Error Analysis (Normality Check)
# # =========================================================
# # Check if residuals are normally distributed (as assumed)
# par(mfrow = c(2, 3)) # Grid layout
# trait_names <- c("Tarsus", "Wing", "Beak L", "Beak H", "Mass")
# 
# for(i in 1:5) {
#   qqnorm(Residuals[, i], main = paste("QQ Plot:", trait_names[i]))
#   qqline(Residuals[, i], col = "red")
# }
# 
# # =========================================================
# # ANALYSIS 3: Residual vs. Fitted (Homoscedasticity)
# # =========================================================
# # Check for "The Cone" shape (variance increasing with size)
# par(mfrow = c(2, 3))
# 
# for(i in 1:5) {
#   plot(Mu_hat[, i], Residuals[, i], 
#        main = paste("Resid vs Fit:", trait_names[i]),
#        xlab = "Predicted Value", ylab = "Residual")
#   abline(h = 0, col = "red")
# }
# 
# 
# 
# 
# 









# # 1. Extract draws as a 3D array: [Iterations, Chains, Parameters]
# fit_array <- as.array(out_lv)
# 
# # 2. Keep only Chains 1, 2, 4 (Index 1, 2, 4 in the second dimension)
# clean_array <- fit_array[, c(1, 2, 4), ]
# 
# # 3. Use rstan's monitor function to get the summary table
# #    warmup = 0 because we already excluded warmup in extraction
# clean_summary <- monitor(clean_array, warmup = 0, print = FALSE)
# 
# # 4. Filter for the parameters you want and print
# pars_of_interest <- c("Lambda", "h2_psi", "pe2_psi", "te2_psi", "sd_R", "h2_trait_total")
# # Grep looks for these patterns in the row names
# rows_to_show <- grep(paste(pars_of_interest, collapse="|"), rownames(clean_summary))
# 
# print(clean_summary[rows_to_show, c("mean", "sd", "2.5%", "50%", "97.5%", "n_eff", "Rhat")], digits = 3)






# 
# print(out_lv, 
#       pars = c("diag_var_genetic", "diag_var_total"),
#       probs = c(0.025, 0.5, 0.975),
#       digits_summary = 3)



# 
# 
# # ==============================================================================
# # FAST RESIDUAL CHECK (No Re-run Needed)
# # ==============================================================================
# 
# # 1. Get the posterior mean of Beta (Fixed Effects)
# beta_summary <- summary(out_lv, pars = "beta")$summary
# beta_mean <- matrix(beta_summary[, "mean"], 
#                     nrow = ncol(dataset_lv_final$X), 
#                     ncol = dataset_lv_final$Nt)
# 
# # 2. Calculate "Marginal Predictions" (Just Sex/Fixed Effects)
# #    Formula: Pred = X * Beta
# pred_marginal <- dataset_lv_final$X %*% beta_mean
# 
# # 3. Calculate "Total Residuals" (Genetics + Environment)
# resid_total <- dataset_lv_final$Y - pred_marginal
# 
# # 4. Plot to check for "The Cone" (Log Transform Needed?)
# trait_names <- c("Wing", "Tarsus", "Beak_L", "Beak_D", "Mass")
# par(mfrow = c(2, 3)) 
# 
# for(i in 1:5) {
#   # Plot Predicted (Marginal) vs Total Residual
#   plot(pred_marginal[, i], resid_total[, i],
#        main = paste(trait_names[i], "- Total Resid"),
#        xlab = "Predicted Value (Sex Effect)", 
#        ylab = "Total Residual (Gen + Env)",
#        pch = 20, col = rgb(0,0,0,0.2))
#   abline(h = 0, col = "red")
# }
# 
# cor(resid_total)
# 
# # ==============================================================================
# # 1. EXTRACT PREDICTED VALUES (Genetics + Fixed Effects)
# # ==============================================================================
# # Get the mean 'mu' for every observation
# mu_summary <- summary(out_lv, pars = "mu")$summary # Use your fit object name
# pred_matrix <- matrix(mu_summary[, "mean"], 
#                       nrow = dataset_lv_final$No, 
#                       ncol = dataset_lv_final$Nt, 
#                       byrow = FALSE)
# 
# # ==============================================================================
# # 2. CALCULATE RESIDUALS
# # ==============================================================================
# # Residual = Observed Data - Predicted Mean
# resid_matrix <- dataset_lv_final$Y - pred_matrix
# 
# # ==============================================================================
# # 3. CHECK 1: "The Cone" (Heteroscedasticity)
# # ==============================================================================
# trait_names <- c("Wing", "Tarsus", "Beak_L", "Beak_D", "Mass")
# par(mfrow = c(2, 3)) 
# 
# for(i in 1:5) {
#   plot(pred_matrix[, i], resid_matrix[, i],
#        main = trait_names[i],
#        xlab = "Predicted Genetic Value", ylab = "Environmental Residual",
#        pch = 20, col = rgb(0,0,0,0.2))
#   abline(h = 0, col = "red")
# }
# # LOOK FOR: A fan shape (getting wider/narrower). 
# # IF FOUND: Log-transform that trait.
# 
# # ==============================================================================
# # 4. CHECK 2: "The Missing Correlation" (Do we need a TE Factor?)
# # ==============================================================================
# # Calculate the correlation of the *residuals*
# res_cor <- cor(resid_matrix)
# print(round(res_cor, 2))
# 
# 
# out_2 <- readRDS('Output_stupid_test_all_traits.rds')
# 
# library(rstan)
# 
# print(out_2, 
#       probs = c(0.025, 0.5, 0.975),
#       digits_summary = 3)
# 
# rstan::traceplot(out_2, pars = c("h2_psi", "lambda[1]", "lambda[2]", "lambda[3]", "lambda[4]", "lambda[5]"), inc_warmup = FALSE)

# 
# # Print summary of the main parameters
# print(out_lv, pars = c(  "beta", "lambda",
#                          "var_G1", "var_PE1", "var_TE1", "var_LV1_total",
#                          "var_G2", "var_PE2", "var_TE2", "var_LV2_total",
#                          "sd_R", "h2_traits", "lp__"),
#       probs = c(0.025, 0.975))


# Print summary of the main parameters
# print(out_lv, pars = c(  "beta", "lambda", 
#                          "var_G1", "var_PE1", "var_TE1", "var_LV1_total",
#                          "sd_R", "h2_traits", "lp__"),
#       probs = c(0.025, 0.975))


# 
# 
# library(dplyr)
# fit_object <- out_lv
# 
# # 1. Get the array of draws: [iterations, chains, parameters]
# draws_array <- as.array(fit_object)
# 
# # 2. Subset to keep only Chains 1, 2, and 3 (Dropping the slow Chain 4)
# # Dimensions: [all_rows, chains 1:3, all_columns]
# draws_subset <- draws_array[, 1:3, ]
# 
# # 3. Calculate diagnostics on just these 3 chains
# # We use the 'posterior' package for this (standard with cmdstanr, or available on CRAN)
# if (!require("posterior")) install.packages("posterior")
# library(posterior)
# 
# # Create a summary table
# summary_subset <- summarise_draws(draws_subset, 
#                                   default_summary_measures(),
#                                   default_convergence_measures())
# 
# # 4. Check the new Rhats for your trouble parameters
# summary_subset %>% 
#   filter(variable %in% c("lambda[2]", "lambda[3]", "var_G1", "var_PE1")) %>%
#   select(variable, mean, sd, rhat, ess_bulk)
