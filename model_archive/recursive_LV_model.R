# SETUP AND LIBRARIES
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

# LOAD DATA
# ------------------------------------------------------------------------------
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


# PREPARE TRAIT DATA (NOW 3 TRAITS)
# ------------------------------------------------------------------------------
traitData <- read.csv("morphology.txt", sep = ";", header = TRUE, 
                      na.strings = c("NA", ""), fileEncoding = "Windows-1252", 
                      stringsAsFactors = FALSE, colClasses = c(ringnr = "character"))

# NOTE: Adjust 'tars_h' below if your column is named 'tars_h_l', 'fot', etc.
trait_subset_temp <- traitData %>%
  filter(sted_r == "hestmannøy", 
         !is.na(ving_h), 
         !is.na(nebb_l),
         !is.na(tars_h),
         !is.na(nebb_h),
         !is.na(vekt)) %>% 
  
  mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr)))) %>%
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

# ALIGN DATA WITH PEDIGREE
# ------------------------------------------------------------------------------
valid_ids <- intersect(trait_subset_temp$ringnr, pedigree_prepped$ringnr)

trait_subset_lv <- trait_subset_temp %>%
  filter(ringnr %in% valid_ids)

n_ind_final <- nrow(trait_subset_lv)
cat("Total number of measurements in analysis:", n_ind_final, "\n")

# Standardize traits (Store mean/sd for all 3 traits)
mean_wing <- mean(trait_subset_lv$ving_h, na.rm = TRUE)
sd_wing   <- sd(trait_subset_lv$ving_h, na.rm = TRUE)

mean_beak <- mean(trait_subset_lv$nebb_l, na.rm = TRUE)
sd_beak   <- sd(trait_subset_lv$nebb_l, na.rm = TRUE)

mean_tars_h <- mean(trait_subset_lv$tars_h, na.rm = TRUE)
sd_tars_h   <- sd(trait_subset_lv$tars_h, na.rm = TRUE)

mean_nebb_h <- mean(trait_subset_lv$nebb_h, na.rm = TRUE)
sd_nebb_h   <- sd(trait_subset_lv$nebb_h, na.rm = TRUE)

mean_vekt <- mean(trait_subset_lv$vekt, na.rm = TRUE)
sd_vekt   <- sd(trait_subset_lv$vekt, na.rm = TRUE)



trait_subset_lv <- trait_subset_lv %>%
  mutate(
    ving_h_std = (ving_h - mean_wing) / sd_wing,
    nebb_l_std = (nebb_l - mean_beak) / sd_beak,
    tars_h_std = (tars_h - mean_tars_h) / sd_tars_h,
    nebb_h_std = (nebb_h - mean_nebb_h) / sd_nebb_h,
    vekt_std   = (vekt - mean_vekt) / sd_vekt
    
  )







# 1. Setup the Working Pedigree
pedigree_working <- pedigree_prepped

# *** CRITICAL FIX: Rename 'ringnr' to 'id' for MCMCglmm ***
# MCMCglmm strictly requires "id", "dam", "sire" (or "animal")
colnames(pedigree_working) <- c("id", "dam", "sire") 

# 2. Force IDs to Character (Safety)
pedigree_working$id   <- as.character(pedigree_working$id)
pedigree_working$dam  <- as.character(pedigree_working$dam)
pedigree_working$sire <- as.character(pedigree_working$sire)


# Use recursive method (new) instead of A matrix

# 3. Prune Pedigree
# Now prune using the EXACT IDs from your data
# Note: We use trait_subset_temp here (the pre-filtered data) to ensure we don't start with 0 rows
ped_pruned <- MCMCglmm::prunePed(pedigree_working, keep = trait_subset_lv$ringnr)

# # Calculate a A-inverse and reordering. Also we don't only want the offspring (therefore "ALL"), and we don't want to assume a variance, we want to estimate it with Stan (therefore "FALSE")
# A_inv_obj <- MCMCglmm::inverseA(ped_pruned, nodes = "ALL", scale = FALSE)













# Calculate A-inverse
A_inv_obj <- MCMCglmm::inverseA(ped_pruned, nodes = "ALL", scale = FALSE)

# *** FIX: Get IDs from the Ainv matrix rownames, NOT the pedigree rownames ***
# MCMCglmm stores the character IDs in the sparse matrix names
ordered_ids <- rownames(A_inv_obj$Ainv) 

# Safety check to ensure we have IDs
if(is.null(ordered_ids)) stop("Error: IDs were lost in inverseA calculation.")

# Create the ped_ordered dataframe using the correct IDs
# Note: A_inv_obj$pedigree[, 2] and [, 3] are already integer indices pointing 
# to the rows of this new ordered list, which is exactly what Stan needs.
ped_ordered <- data.frame(
  id = ordered_ids,
  dam = A_inv_obj$pedigree[, 2],
  sire = A_inv_obj$pedigree[, 3]
)

dii_vector <- A_inv_obj$dii
Na <- nrow(ped_ordered)

# Handle missing parents: set to Na + 1 (dummy)
dam_idx <- ped_ordered$dam
sire_idx <- ped_ordered$sire

# MCMCglmm returns 0 for missing parents in the integer index; map these to Na + 1
dam_idx[is.na(dam_idx) | dam_idx == 0] <- Na + 1
sire_idx[is.na(sire_idx) | sire_idx == 0] <- Na + 1

# Create the map
id_map <- setNames(1:Na, ped_ordered$id)







# 
# 
# 
# # Ordering and extracting the dii vector (relatedness scaling factor)
# ped_ordered <- data.frame(
#   id = rownames(A_inv_obj$pedigree),
#   dam = A_inv_obj$pedigree[, 2],
#   sire = A_inv_obj$pedigree[, 3]
# )
# dii_vector <- A_inv_obj$dii
# 
# Na <- nrow(ped_ordered)
# 
# # Handle missing parents: set to Na + 1 (dummy)
# dam_idx <- ped_ordered$dam
# sire_idx <- ped_ordered$sire
# dam_idx[is.na(dam_idx) | dam_idx == 0] <- Na + 1
# sire_idx[is.na(sire_idx) | sire_idx == 0] <- Na + 1
# 
# id_map <- setNames(1:Na, ped_ordered$id)



# 8. Map Data and Verify
temp_animal_idx <- id_map[trait_subset_lv$ringnr]
missing_mask <- is.na(temp_animal_idx)

# 9. Final Safety Check
if (sum(missing_mask) > 0) {
  cat("WARNING: Removing", sum(missing_mask), "rows. ID mismatch persists.\n")
  trait_subset_lv <- trait_subset_lv[!missing_mask, ]
  animal_idx <- id_map[trait_subset_lv$ringnr]
} else {
  cat("SUCCESS: All", nrow(trait_subset_lv), "data rows mapped to pedigree successfully.\n")
  animal_idx <- temp_animal_idx
}





# animal_idx <- id_map[trait_subset_lv$ringnr]



# # BUILD A MATRIX (Subsetted)
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
# all_relevant_ids <- get_ancestors(valid_ids, pedigree_prepped)
# smaller_pedigree <- pedigree_prepped %>%
#   filter(ringnr %in% all_relevant_ids) %>%
#   arrange(match(ringnr, all_relevant_ids))
# 
# smaller_ped_prepped <- prepPed(smaller_pedigree, gender = NULL, check = TRUE)
# A_small_nadiv <- makeA(smaller_ped_prepped)
# A_hestmannoy <- A_small_nadiv[valid_ids, valid_ids]
# 
# # PREPARE STAN DATA
# # ------------------------------------------------------------------------------
# matrix_ids <- rownames(A_hestmannoy)
# animal_map <- setNames(seq_len(length(matrix_ids)), matrix_ids)
# animal_idx <- animal_map[trait_subset_lv$ringnr]

# Na = length(matrix_ids)

# Y Matrix now has 3 columns
Y_mat <- as.matrix(trait_subset_lv[, c("ving_h_std", "nebb_l_std", "tars_h_std", "nebb_h_std", "vekt_std")])
X_mat <- model.matrix(~ sex, data = trait_subset_lv)
K <- ncol(X_mat)





cat("\n=== DATA COUNTS ===\n")
cat(sprintf("Observations (Rows in Y): %d\n", nrow(trait_subset_lv)))
cat(sprintf("Unique Individuals:       %d\n", length(unique(trait_subset_lv$ringnr))))
cat(sprintf("Total Animals in Ped (Na):%d\n", Na))
cat("===================\n\n")

cat("=== ID MATCHING CHECK (First 5 Rows) ===\n")














# A. Create a precise lookup map: Bird ID (Character) -> Row Number (Integer)
# This ensures we map "8573160" to something like "505"
all_ped_ids <- as.character(ped_ordered$id)
id_map_lookup <- setNames(seq_along(all_ped_ids), all_ped_ids)

# B. Extract parents as characters (to match the lookup keys)
dam_chars  <- as.character(ped_ordered$dam)
sire_chars <- as.character(ped_ordered$sire)

# C. Perform the translation
clean_dam  <- id_map_lookup[dam_chars]
clean_sire <- id_map_lookup[sire_chars]

# D. Handle Missing Parents
# Any parent not found in the pruned pedigree (NA) gets the dummy code (Na + 1)
missing_code <- Na + 1
clean_dam[is.na(clean_dam)]   <- missing_code
clean_sire[is.na(clean_sire)] <- missing_code

# E. Final Integer Enforcement
clean_dam  <- as.integer(clean_dam)
clean_sire <- as.integer(clean_sire)

# F. Diagnostic Check
if(max(clean_dam) > (Na + 1) || max(clean_sire) > (Na + 1)) {
  stop("CRITICAL ERROR: Generated indices are still too large!")
} else {
  cat(sprintf("Success: Parent indices mapped. Max Index: %d (Limit: %d)\n", max(clean_dam), Na + 1))
}

# ==============================================================================
# 2. FINAL DATA LIST & RUN
# ==============================================================================

dataset_lv_final <- list(
  No = nrow(trait_subset_lv),
  Nt = 5, 
  Nlv = 2,
  Y = as.matrix(trait_subset_lv[, c("ving_h_std", "nebb_l_std", "tars_h_std", "nebb_h_std", "vekt_std")]),
  X = model.matrix(~ sex, data = trait_subset_lv),
  K = ncol(model.matrix(~ sex, data = trait_subset_lv)),
  
  animal = as.array(as.integer(animal_idx)),
  Na = Na,
  dam = as.array(clean_dam),   # The corrected integer array
  sire = as.array(clean_sire), # The corrected integer array
  dii = as.array(dii_vector)
)

cat("Data list verified. Launching Stan...\n")

tomonitor_lv <- c(
  "beta", "lambda", 
  "var_LV_1_genetic", "var_LV_1_perm_env", "var_LV_1_temp_env",
  "var_LV_2_genetic", "var_LV_2_perm_env", "var_LV_2_temp_env",
  "sd_R", "h2_traits", "pe2_traits", "lp__"
)

out_lv <- stan(
  file = 'LV_Recursive.stan', 
  data = dataset_lv_final,
  pars = tomonitor_lv,
  chains = 4, iter = 10000, warmup = 5000, thin = 1,
  seed = 123
)

saveRDS(out_lv, 'Output_LV_Recursive.rds')
cat("Done. Sampling started.\n")














# # 1. Grab the first 5 IDs from your data
# orig_ids <- trait_subset_lv$ringnr[1:5]
# 
# # 2. Grab the first 5 integer pointers sending to Stan
# mapped_indices <- animal_idx[1:5]
# 
# # 3. Use those integers to look up the name in the pedigree key
# recovered_ids <- ped_ordered$id[mapped_indices]
# 
# # 4. Show side-by-side
# print(data.frame(
#   Row = 1:5,
#   Original_ID_Data = orig_ids,
#   Stan_Integer_Map = mapped_indices,
#   Recovered_ID_Ped = recovered_ids,
#   MATCH = (orig_ids == recovered_ids)
# ))











# ==============================================================================
# ROBUST INDEX SANITIZATION
# ==============================================================================
# We force the parent columns to be numeric integers. 
# If R sees them as Factors or Characters, this converts them safely.

# # 1. Force to numeric (using as.character first handles Factor issues)
# clean_dam  <- as.numeric(as.character(ped_ordered$dam))
# clean_sire <- as.numeric(as.character(ped_ordered$sire))
# 
# # 2. The "Na" variable is the number of animals. 
# #    We map all Missing (NA) or 0 values to (Na + 1)
# missing_code <- Na + 1
# 
# clean_dam[is.na(clean_dam) | clean_dam == 0]   <- missing_code
# clean_sire[is.na(clean_sire) | clean_sire == 0] <- missing_code
# 
# # 3. Final verification checks
# if (any(is.na(clean_dam)) || any(is.na(clean_sire))) {
#   stop("CRITICAL: Index sanitization failed. Still found NAs in dam/sire vectors.")
# }
# 
# cat("Success: Parent indices sanitized. No NAs found.\n")
# 
# # ==============================================================================
# # FINAL DATA LIST & RUN
# # ==============================================================================
# 
# dataset_lv_final <- list(
#   No = nrow(trait_subset_lv),
#   Nt = 5, 
#   Nlv = 2,
#   # Ensure Y is a matrix
#   Y = as.matrix(trait_subset_lv[, c("ving_h_std", "nebb_l_std", "tars_h_std", "nebb_h_std", "vekt_std")]),
#   X = model.matrix(~ sex, data = trait_subset_lv),
#   K = ncol(model.matrix(~ sex, data = trait_subset_lv)),
#   
#   # Ensure all integer arrays are clean
#   animal = as.array(as.integer(animal_idx)),
#   Na = Na,
#   dam = as.array(as.integer(clean_dam)),   # Use the CLEAN variables
#   sire = as.array(as.integer(clean_sire)), # Use the CLEAN variables
#   dii = as.array(dii_vector)
# )
# 
# # Run Stan
# cat("Launching Stan...\n")
# 
# tomonitor_lv <- c(
#   "beta", "lambda", 
#   "var_LV_1_genetic", "var_LV_1_perm_env", "var_LV_1_temp_env",
#   "var_LV_2_genetic", "var_LV_2_perm_env", "var_LV_2_temp_env",
#   "sd_R", "h2_traits", "pe2_traits", "lp__"
# )
# 
# out_lv <- stan(
#   file = 'LV_Recursive.stan', 
#   data = dataset_lv_final,
#   pars = tomonitor_lv,
#   chains = 4, iter = 10000, warmup = 5000, thin = 1,
#   seed = 123
# )
# 
# saveRDS(out_lv, 'Output_LV_Recursive.rds')
# cat("Done.\n")
# cat("Sampling started.\n")


















# 
# 
# dataset_lv <- list(
#   No = nrow(trait_subset_lv),
#   Nt = 5, 
#   Nlv = 2,
#   Y = Y_mat,
#   X = X_mat,
#   K = K,
#   animal = as.array(animal_idx),
#   Na = Na,
#   dam = as.array(dam_idx),
#   sire = as.array(sire_idx),
#   dii = as.array(dii_vector)
#   # ,
#   # A = as.matrix(A_hestmannoy)
# )
# 
# # RUN STAN
# # ------------------------------------------------------------------------------
# # library(rstan)
# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())
# 
# 
# tomonitor_lv <- c(
#   "beta", 
#   "lambda", 
#   "var_LV_1_genetic", 
#   "var_LV_1_perm_env", 
#   "var_LV_1_temp_env",
#   "var_LV_2_genetic", 
#   "var_LV_2_perm_env", 
#   "var_LV_2_temp_env",
#   "sd_R", 
#   "h2_traits", 
#   "pe2_traits",
#   "lp__"
# )
# 
# 
# out_lv <- stan(file = 'LV_Recursive.stan',
#                data = dataset_lv,
#                pars = tomonitor_lv,
#                chains = 4, iter = 10000, warmup = 5000, thin = 1,
#                seed = 123)
# 
# saveRDS(out_lv, 'Output_LV_Recursive.rds')
# 
# 
# 


# Load the output if not already in memory
out_lv <- readRDS('Output_LV_Recursive.rds')

# Print summary of the main parameters
print(out_lv, pars = c( "beta", 
                        "lambda", 
                        "var_LV_1_genetic", 
                        "var_LV_1_perm_env", 
                        "var_LV_1_temp_env",
                        "var_LV_2_genetic", 
                        "var_LV_2_perm_env", 
                        "var_LV_2_temp_env",
                        "sd_R", 
                        "h2_traits", 
                        "pe2_traits",
                        "lp__"), 
      probs = c(0.025, 0.5, 0.975))

