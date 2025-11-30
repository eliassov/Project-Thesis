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

# Add external parents (founders)
all_parents <- unique(c(pedigree_clean$dam, pedigree_clean$sire))
external_parents <- setdiff(all_parents, c(pedigree_clean$ringnr, NA_character_))

if (length(external_parents) > 0) {
  founder_df <- data.frame(ringnr = external_parents, dam = NA_character_, sire = NA_character_)
  pedigree_full <- bind_rows(founder_df, pedigree_clean) %>% distinct(ringnr, .keep_all = TRUE)
} else {
  pedigree_full <- pedigree_clean
}

# Prepare pedigree for A-matrix calculation
pedigree_prepped <- prepPed(pedigree_full, gender = NULL, check = TRUE)

# ==============================================================================
# 3. PREPARE TRAIT DATA (BIVARIATE)
# ==============================================================================
traitData <- read.csv("morphology.txt", sep = ";", header = TRUE, 
                      na.strings = c("NA", ""), fileEncoding = "Windows-1252", 
                      stringsAsFactors = FALSE, colClasses = c(ringnr = "character"))

# Step A: Filter and Clean Data FIRST
# We filter for Island + BOTH Traits + Valid Sex
trait_subset_temp <- traitData %>%
  filter(sted_r == "hestmannøy", 
         !is.na(ving_h), 
         !is.na(nebb_l)) %>%  # Bivariate: Must have both traits
  
  mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr)))) %>%
  
  # Logic to clean Sex Column (Priority: scriptsex > fieldsex)
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
  filter(!is.na(sex)) # Drop individuals with unknown sex

# ==============================================================================
# 4. ALIGN DATA WITH PEDIGREE
# ==============================================================================

# Identify valid IDs (intersection of clean data and pedigree)
valid_ids <- intersect(trait_subset_temp$ringnr, pedigree_prepped$ringnr)

# Create FINAL data subset sorted by ID
trait_subset <- trait_subset_temp %>%
  filter(ringnr %in% valid_ids) %>%
  arrange(ringnr)

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

# Build small pedigree containing only relevant ancestors
all_relevant_ids <- get_ancestors(valid_ids, pedigree_prepped)

smaller_pedigree <- pedigree_prepped %>%
  filter(ringnr %in% all_relevant_ids) %>%
  arrange(match(ringnr, all_relevant_ids))

# Calculate A-matrix
smaller_ped_prepped <- prepPed(smaller_pedigree, gender = NULL, check = TRUE)
A_small_nadiv <- makeA(smaller_ped_prepped)

# Subset A to valid phenotyped individuals AND sort to match data
A_bivariate <- A_small_nadiv[valid_ids, valid_ids]
A_bivariate <- A_bivariate[order(rownames(A_bivariate)), order(colnames(A_bivariate))]

# Safety Check
stopifnot(all(rownames(A_bivariate) == trait_subset$ringnr))
cat("Data and Matrix aligned. No =", nrow(trait_subset), "Na =", length(valid_ids), "\n")

# ==============================================================================
# 6. PREPARE STAN DATA
# ==============================================================================

# Map animals to integers
animal_map <- setNames(seq_len(length(valid_ids)), valid_ids)
animal_idx <- animal_map[trait_subset$ringnr]

# Prepare Y Matrix (No x 2)
Y_mat <- as.matrix(trait_subset[, c("ving_h", "nebb_l")])

# Prepare Fixed Effects Matrix X (Intercept + Sex)
X_mat <- model.matrix(~ sex, data = trait_subset)
K <- ncol(X_mat)

dataset_bivariate <- list(
  No = nrow(trait_subset),
  Y = Y_mat,
  X = X_mat,       
  K = K,          
  animal = animal_idx,
  Na = length(valid_ids),
  A = as.matrix(A_bivariate)
)

# ==============================================================================
# 7. RUN STAN
# ==============================================================================
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# MCMC settings
nc <- 4
nw <- 3000
ni <- 6000
nt <- 1

# Parameters to monitor (Note: 'beta' instead of 'mu')
tomonitor <- c('beta', 'Sigma_A', 'Sigma_E', 'Sigma_R', 
               'heritability1', 'heritability2', 
               'cor_A', 'cor_E', 'cor_R')

out_bivariate <- stan(file = 'Bivariate_Animal_Model.stan',
                      data = dataset_bivariate, 
                      pars = tomonitor,
                      chains = nc, iter = ni, warmup = nw, thin = nt,
                      open_progress = FALSE,
                      seed = 123)

# Save outputs
saveRDS(out_bivariate, 'Output_Bivariate_Animal_Model.rds')

summ_biv <- summary(out_bivariate)$summary
write.csv(as.data.frame(summ_biv), file = "Output_Bivariate_Summary.csv", row.names = TRUE)

cat("Bivariate analysis complete. Summary saved.\n")







# 
# # Set up personal library
# personal_lib <- "~/R/library"
# if (!dir.exists(personal_lib)) {
#   dir.create(personal_lib, recursive = TRUE)
# }
# .libPaths(c(personal_lib, .libPaths()))
# 
# # Install missing packages to personal library
# required_packages <- c("rstan", "tidyverse", "nadiv") 
# for (pkg in required_packages) {
#   if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
#     install.packages(pkg, lib = personal_lib, repos = "https://cloud.r-project.org")
#     library(pkg, character.only = TRUE)
#   }
# }
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
# 
# traitData <- read.csv("morphology.txt",
#                       sep = ";",
#                       header = TRUE,
#                       na.strings = c("NA", ""),
#                       fileEncoding = "Windows-1252",
#                       stringsAsFactors = FALSE,
#                       colClasses = c(ringnr = "character"))
# 
# 
# pedData <- read.csv("pedigree_data.txt",
#                     sep = " ",
#                     header = TRUE,
#                     na.strings = "NA")
# 
# # Re-preprocess from pedigree_clean
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
# 
# # Run prepPed with NA as missing
# pedigree_prepped <- prepPed(pedigree_full, gender = NULL, check = TRUE)
# 
# 
# 
# # Re-run makeA
# tryCatch({
#   A_nadiv <- makeA(pedigree_prepped)
#   A_nadiv_dense <- as.matrix(A_nadiv)
#   cat("New A matrix computed. Dimensions:", dim(A_nadiv_dense), "\n")
#   cat("Diagonal summary:", summary(diag(A_nadiv_dense)), "\n")
# }, error = function(e) {
#   cat("nadiv Error:", e$message, "\n")
# })
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
# trait_subset_temp <- traitData %>%
#   # A. Initial Filter: Island and Traits
#   filter(sted_r == "hestmannøy", 
#          !is.na(ving_h)) %>%
#   
#   # B. Clean Ring Number (Standardize ID format)
#   mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr)))) %>%
#   
#   # C. Sex Logic: Create raw column, then clean/map it
#   mutate(
#     raw_sex = if_else(!is.na(scriptsex) & scriptsex != "u", scriptsex, fieldsex)
#   ) %>%
#   mutate(
#     sex = case_when(
#       raw_sex %in% c("m", "pm") ~ "m",
#       raw_sex %in% c("f", "pf") ~ "f",
#       TRUE ~ NA_character_ # Maps 'u', NA, and typos to NA
#     )
#   ) %>%
#   
#   # D. FINAL FILTER: Drop rows where sex is missing
#   # This is the critical step that was missing from your ID selection
#   filter(!is.na(sex)) 
# 
# # ==============================================================================
# # 2. DEFINE VALID IDs BASED ON CLEAN DATA
# # ==============================================================================
# 
# # Now we know exactly which birds are left. 
# # We only keep those that are ALSO in the pedigree.
# valid_phenotyped_ids <- intersect(trait_subset_temp$ringnr, pedigree_prepped$ringnr)
# 
# # Now create the final data object using ONLY these valid IDs
# trait_subset <- trait_subset_temp %>%
#   filter(ringnr %in% valid_phenotyped_ids) %>%
#   arrange(ringnr) # Sort to ensure alignment with matrix
# 
# # ==============================================================================
# # 3. BUILD A-MATRIX (Using only the Valid IDs)
# # ==============================================================================
# 
# # Get ancestry for ONLY the valid individuals
# # (Reuse your existing get_ancestors function)
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
# # Subset and Sort A-matrix to match trait_subset
# A_hestmannoy <- A_small_nadiv[valid_phenotyped_ids, valid_phenotyped_ids]
# A_hestmannoy <- A_hestmannoy[order(rownames(A_hestmannoy)), order(colnames(A_hestmannoy))]
# 
# # ==============================================================================
# # 4. PREPARE STAN DATA
# # ==============================================================================
# 
# No <- nrow(trait_subset)
# Na <- length(valid_phenotyped_ids)
# 
# 
# # Check 1: Every animal in the data must exist in the matrix
# if (!all(unique(trait_subset$ringnr) %in% rownames(A_hestmannoy))) {
#   stop("Error: Data contains individuals not in the A-matrix!")
# }
# 
# # Check 2: Every animal in the matrix must appear in the data (optional but good practice)
# if (!all(rownames(A_hestmannoy) %in% trait_subset$ringnr)) {
#   warning("Warning: A-matrix contains individuals with no data observations.")
# }
# 
# # Map animals to integers 1:Na
# animal_map <- setNames(seq_len(Na), valid_phenotyped_ids)
# animal <- animal_map[trait_subset$ringnr]
# 
# # Prepare Y (Response Variable)
# # Note: You usually only need the trait column, not the sex column in Y
# Y <- as.numeric(trait_subset$ving_h) 
# 
# 
# # 1. Prepare Fixed Effects Matrix (X)
# # This creates a matrix with an Intercept column (all 1s) and a Sex column (0/1)
# X_mat <- model.matrix(~ sex, data = trait_subset)
# 
# # Get number of fixed effect predictors (K)
# K <- ncol(X_mat) 
# 
# # 2. Update the Data Bundle
# dataset <- list(
#   No = nrow(trait_subset),       # Number of observations
#   Y = Y,
#   X = X_mat,                     # NEW: Fixed effects design matrix
#   K = K,                         # NEW: Number of fixed effects
#   animal = animal_map[trait_subset$ringnr], # Map IDs to integers
#   Na = length(valid_phenotyped_ids),        # Number of unique individuals
#   A = as.matrix(A_hestmannoy)    # Relationship matrix
# )
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
# # GET ONLY THE BIRDS FROM HESTMANNØY (that are phenotyped)
# 
# # Identify phenotyped individuals meeting criteria
# phenotyped_ids <- traitData %>%
#   filter(sted_r == "hestmannøy" & !is.na(ving_h)) %>%
#   mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr)))) %>%  # Match pedigree preprocessing
#   pull(ringnr) %>%
#   unique()
# 
# # Ensure phenotyped_ids are in pedigree_prepped
# valid_phenotyped_ids <- phenotyped_ids[phenotyped_ids %in% pedigree_prepped$ringnr]
# 
# 
# # Get all ancestors for valid phenotyped IDs
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
# # Subset to valid phenotyped individuals
# A_hestmannoy <- A_small_nadiv[valid_phenotyped_ids, valid_phenotyped_ids]
# 
# 
# # Diagnostics
# # cat("1. Phenotyped individuals (after filtering):", length(valid_phenotyped_ids), "\n")
# # cat("2. Smaller pedigree rows:", nrow(smaller_pedigree), "\n")
# # cat("3. A_small dimensions:", dim(A_hestmannoy), "\n")
# # cat("4. A_small symmetric:", isSymmetric(A_hestmannoy), "\n")
# # cat("5. Diagonal summary (inbreeding):", summary(diag(as.matrix(A_hestmannoy))), "\n")
# # cat("6. All phenotyped IDs in A_small:", all(valid_phenotyped_ids %in% rownames(A_hestmannoy)), "\n")
# # cat("7. Positive semi-definite:", all(eigen(as.matrix(A_hestmannoy), only.values = TRUE)$values >= 0), "\n")
# 
# 
# 
# 
# 
# 
# 
# # ====================
# 
# ## Bivariate animal model
# 
# 
# # Step 1: Identify phenotyped individuals meeting criteria
# phenotyped_ids_bivariate <- traitData %>%
#   filter(sted_r == "hestmannøy",
#          !is.na(ving_h),
#          !is.na(nebb_l)
#   )  %>%
#   mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr)))) %>%  # Match pedigree preprocessing
#   pull(ringnr) %>%
#   unique()
# 
# 
# # Step 2: Ensure phenotyped_ids are in pedigree_prepped
# valid_phenotyped_ids_bivariate <- phenotyped_ids_bivariate[phenotyped_ids_bivariate %in% pedigree_prepped$ringnr]
# 
# 
# 
# A_hestmannoy_bivariate <- A_small_nadiv[valid_phenotyped_ids_bivariate, valid_phenotyped_ids_bivariate]
# 
# 
# 
# 
# 
# # Extract real data for Hestmannøy (two traits: ving_h and nebb_l)
# trait_subset_bivariate <- traitData %>%
#   filter(ringnr %in% valid_phenotyped_ids_bivariate, sted_r == "hestmannøy", !is.na(ving_h), !is.na(nebb_l)) %>%  # Filter for non-missing both traits
#   mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr))))
# 
# No_2 <- nrow(trait_subset_bivariate)
# Na_2 <- length(valid_phenotyped_ids_bivariate)
# animal_map_2 <- setNames(seq_len(Na_2), valid_phenotyped_ids_bivariate)
# animal_2 <- animal_map_2[trait_subset_bivariate$ringnr]
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
# # ... (Data cleaning and A-matrix steps remain the same) ...
# 
# # 1. Prepare Response Matrix Y (No x 2)
# Y_2 <- as.matrix(trait_subset_bivariate[, c("ving_h", "nebb_l")])
# 
# # 2. Prepare Design Matrix X
# # We need an intercept and sex effect. 
# # Standard 'model.matrix' creates: Column 1 = Intercept (1s), Column 2 = Sex (0/1)
# X_mat <- model.matrix(~ sex, data = trait_subset_bivariate)
# K <- ncol(X_mat) # K=2 (Intercept, Sex)
# 
# # 3. Bundle Data for Stan
# dataset_bivariate <- list(
#   No = nrow(trait_subset_bivariate),
#   Y = Y_2,
#   X = X_mat,                     # NEW: Design matrix
#   K = K,                         # NEW: Number of fixed effects
#   animal = animal_2,
#   Na = length(valid_phenotyped_ids_bivariate),
#   A = as.matrix(A_hestmannoy_bivariate)
# )
# 
# # Update monitored parameters to include 'beta'
# tomonitor_bivariate <- c('beta', 'Sigma_A', 'Sigma_E', 'Sigma_R', 'heritability1', 'heritability2', 'cor_A', 'cor_E', 'cor_R')
# 
# # ... (Call Stan as before) ...
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
# 
# 
# 
# 
# 
# 
# 
# # Response matrix (No x 3 for bivariate with sex as fixed effect)
# Y_2 <- matrix(0, No_2, 3)
# Y_2[,1] <- trait_subset_bivariate$ving_h
# Y_2[,2] <- trait_subset_bivariate$nebb_l
# Y_2[,3] <- trait_subset_bivariate$scriptsex
# 
# # Bundle the data
# dataset_bivariate <- list(
#   No = No_2,
#   Y = Y_2,
#   animal = animal_2,
#   Na = Na_2,
#   A = as.matrix(A_hestmannoy_bivariate)  # Dense for Stan
# )
# 
# ## Analyse the data with the bivariate animal model
# library(rstan)
# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())
# 
# # MCMC settings (increased for bivariate convergence)
# nc_2 <- 4
# nw_2 <- 3000
# ni_2 <- 6000
# nt_2 <- 1
# 
# # Parameters to monitor
# tomonitor_bivariate <- c('mu', 'Sigma_A', 'Sigma_E', 'Sigma_R', 'heritability1', 'heritability2', 'cor_A', 'cor_E', 'cor_R')
# 
# # Call Stan from R
# out_bivariate <- stan(file = 'Bivariate_Animal_Model.stan',
#                       data = dataset_bivariate, 
#                       pars = tomonitor_bivariate,
#                       chains = nc_2, iter = ni_2, warmup = nw_2, thin = nt_2,
#                       open_progress = FALSE,
#                       refresh = 10,
#                       seed = 123)
# 
# # Save the results
# saveRDS(out_bivariate, 'Output_Bivariate_Animal_Model.rds')
# 
# 
# summ_biv <- summary(out_bivariate)$summary
# 
# write.csv(
#   as.data.frame(summ_biv),
#   file = "Output_Bivariate_Summary.csv",
#   row.names = TRUE
# )
# 
# cat("Summary written to: Output_Bivariate_Summary.csv\n")
# 
# 
# 
# 
# 
# 
