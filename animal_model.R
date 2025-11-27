

# Set up personal library
personal_lib <- "~/R/library"
if (!dir.exists(personal_lib)) {
  dir.create(personal_lib, recursive = TRUE)
}
.libPaths(c(personal_lib, .libPaths()))

# Install missing packages to personal library
required_packages <- c("rstan", "tidyverse", "nadiv", "shinystan")  # add any other packages you need
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, lib = personal_lib, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# Rest of your script continues here...


#!/usr/bin/env Rscript
## --------------------------------------------------------------
## 1. Auto-detect script folder (Rscript, source, interactive)
## --------------------------------------------------------------
script_path <- NULL
cmd_args <- commandArgs(trailingOnly = FALSE)
if (any(grepl("^--file=", cmd_args))) {
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  script_path <- sub("^--file=", "", file_arg)
  script_path <- normalizePath(script_path)
}
if (is.null(script_path)) {
  src_file <- tryCatch(getSrcFilename(function(x) x, full.names = TRUE),
                       error = function(e) NULL)
  if (!is.null(src_file) && file.exists(src_file)) script_path <- src_file
}
if (is.null(script_path) || !file.exists(script_path)) {
  script_path <- getwd()
  message("Warning: Using getwd(): ", script_path)
} else {
  script_path <- dirname(script_path)
  message("Success: Script directory: ", script_path)
}
setwd(script_path)
cat("Working directory set to:", getwd(), "\n")



library(tidyverse)
library(rstan)
library(nadiv)



traitData <- read.csv("morphology.txt",
                      sep = ";",
                      header = TRUE,
                      na.strings = c("NA", ""),
                      fileEncoding = "Windows-1252",
                      stringsAsFactors = FALSE,
                      colClasses = c(ringnr = "character"))


pedData <- read.csv("pedigree_data.txt",
                    sep = " ",
                    header = TRUE,
                    na.strings = "NA")

# Re-preprocess from pedigree_clean
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
  founder_df <- data.frame(
    ringnr = external_parents,
    dam = NA_character_,
    sire = NA_character_
  )
  pedigree_full <- bind_rows(founder_df, pedigree_clean) %>%
    distinct(ringnr, .keep_all = TRUE)
} else {
  pedigree_full <- pedigree_clean
}

library(nadiv)

# Run prepPed with NA as missing
pedigree_prepped <- prepPed(pedigree_full, gender = NULL, check = TRUE)

# # Check for UNKNOWN in ringnr
# cat("Any ringnr = UNKNOWN:", any(pedigree_prepped$ringnr == "UNKNOWN"), "\n")

# Re-run makeA
tryCatch({
  A_nadiv <- makeA(pedigree_prepped)
  A_nadiv_dense <- as.matrix(A_nadiv)
  cat("New A matrix computed. Dimensions:", dim(A_nadiv_dense), "\n")
  cat("Diagonal summary:", summary(diag(A_nadiv_dense)), "\n")
  cat("Parent-offspring check (8L19506, 8L19505):", A_nadiv_dense["8L19506", "8L19505"], "\n")
}, error = function(e) {
  cat("nadiv Error:", e$message, "\n")
})



# GET ONLY THE BIRDS FROM HESTMANNØY (that are phenotyped)

library(dplyr)
library(nadiv)
library(Matrix)

# Step 1: Identify phenotyped individuals meeting criteria
phenotyped_ids <- traitData %>%
  filter(sted_r == "hestmannøy" & !is.na(ving_h)) %>%
  mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr)))) %>%  # Match pedigree preprocessing
  pull(ringnr) %>%
  unique()

# Step 2: Ensure phenotyped_ids are in pedigree_prepped
valid_phenotyped_ids <- phenotyped_ids[phenotyped_ids %in% pedigree_prepped$ringnr]
if (length(valid_phenotyped_ids) < length(phenotyped_ids)) {
  cat("Warning: Some phenotyped IDs not in pedigree:", setdiff(phenotyped_ids, pedigree_prepped$ringnr), "\n")
}

# Step 3: Get all ancestors for valid phenotyped IDs
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

all_relevant_ids <- get_ancestors(valid_phenotyped_ids, pedigree_prepped)

# Step 4: Create smaller pedigree
smaller_pedigree <- pedigree_prepped %>%
  filter(ringnr %in% all_relevant_ids) %>%
  arrange(match(ringnr, all_relevant_ids))

# Step 5: Compute smaller A matrix
smaller_ped_prepped <- prepPed(smaller_pedigree, gender = NULL, check = TRUE)
A_small_nadiv <- makeA(smaller_ped_prepped)

# Subset to valid phenotyped individuals
A_hestmannoy <- A_small_nadiv[valid_phenotyped_ids, valid_phenotyped_ids]


# Diagnostics
cat("1. Phenotyped individuals (after filtering):", length(valid_phenotyped_ids), "\n")
cat("2. Smaller pedigree rows:", nrow(smaller_pedigree), "\n")
cat("3. A_small dimensions:", dim(A_hestmannoy), "\n")
cat("4. A_small symmetric:", isSymmetric(A_hestmannoy), "\n")
cat("5. Diagonal summary (inbreeding):", summary(diag(as.matrix(A_hestmannoy))), "\n")
cat("6. All phenotyped IDs in A_small:", all(valid_phenotyped_ids %in% rownames(A_hestmannoy)), "\n")
cat("7. Positive semi-definite:", all(eigen(as.matrix(A_hestmannoy), only.values = TRUE)$values >= 0), "\n")
cat("8. Sample parent-offspring (if available, e.g., 8L19506, 8L19505):",
    if ("8L19506" %in% valid_phenotyped_ids & "8L19505" %in% valid_phenotyped_ids) 
      A_hestmannoy["8L19506", "8L19505"] else "Not in valid_phenotyped_ids", "\n")
cat("9. Sample sibling (if available, e.g., 8L19507, 8L19508):",
    if ("8L19507" %in% valid_phenotyped_ids & "8L19508" %in% valid_phenotyped_ids) 
      A_hestmannoy["8L19507", "8L19508"] else "Not in valid_phenotyped_ids", "\n")


##### =======
##### Run down to here regardless of number of traits analyzed (as long as they include the wing length).
##### =======



# 
# # ====================
# # C) Run the analysis (univariate)
# # ====================
# 
# # Extract real data for Hestmannøy
trait_subset <- traitData %>%
  filter(ringnr %in% valid_phenotyped_ids, sted_r == "hestmannøy", !is.na(ving_h)) %>%
  mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr))))

No <- nrow(trait_subset)
Na <- length(valid_phenotyped_ids)
animal_map <- setNames(seq_len(Na), valid_phenotyped_ids)
animal <- animal_map[trait_subset$ringnr]
# Y <- trait_subset$ving_h

Y <- matrix(0, No, 2)
Y[,1] <- trait_subset$ving_h
Y[,2] <- trait_subset$scriptsex

# Bundle the data
dataset <- list(
  No = No,
  Y = Y,
  animal = animal,
  Na = Na,
  A = as.matrix(A_hestmannoy)  # Dense for Stan
)

## Analyse the data with the simple animal model
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# MCMC settings
nc <- 4
nw <- 1000
ni <- 2000
nt <- 1

# parameters to monitor
tomonitor <- c('mu', 'var_A', 'var_E', 'var_R', 'var_P', 'evolvability', 'heritability')

# Call Stan from R
out_1 <- stan(file = 'C:\\Users\\Elias Ovesen\\OneDrive\\Skrivebord\\Prosjektoppgave\\animal_model_univariate.stan',
            data = dataset,
            pars = tomonitor,
            chains = nc, iter = ni, warmup = nw, thin = nt,
            open_progress = FALSE,
            seed = 123) # Change seed if needed

# Save the results
save(out_1, file = 'Output_Thesis_Animal_Model_Second_Run.Rdata')

# load the results
load('Output_Thesis_Animal_Model_Second_Run.Rdata')

# Take a look at the result summary
out_1

# Take a closer look at the MCMC chains, model diagnostics, and outputs:
library(shinystan)
launch_shinystan(out_1)
























# ====================

## Bivariate animal model

###### (Added later to just include the individuals with both phenotypes in the A matrix)

# Step 1: Identify phenotyped individuals meeting criteria
phenotyped_ids_bivariate <- traitData %>%
  filter(sted_r == "hestmannøy",
         !is.na(ving_h),
         !is.na(nebb_l)
         )  %>%
  mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr)))) %>%  # Match pedigree preprocessing
  pull(ringnr) %>%
  unique()


# Step 2: Ensure phenotyped_ids are in pedigree_prepped
valid_phenotyped_ids_bivariate <- phenotyped_ids_bivariate[phenotyped_ids_bivariate %in% pedigree_prepped$ringnr]
# if (length(valid_phenotyped_ids) < length(phenotyped_ids)) {
#   cat("Warning: Some phenotyped IDs not in pedigree:", setdiff(phenotyped_ids, pedigree_prepped$ringnr), "\n")
# }


A_hestmannoy_bivariate <- A_small_nadiv[valid_phenotyped_ids_bivariate, valid_phenotyped_ids_bivariate]


######






# Extract real data for Hestmannøy (two traits: ving_h and another trait)
trait_subset_bivariate <- traitData %>%
  filter(ringnr %in% valid_phenotyped_ids_bivariate, sted_r == "hestmannøy", !is.na(ving_h), !is.na(nebb_l)) %>%  # Filter for non-missing both traits
  mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr))))

No_2 <- nrow(trait_subset_bivariate)
Na_2 <- length(valid_phenotyped_ids_bivariate)
animal_map_2 <- setNames(seq_len(Na_2), valid_phenotyped_ids_bivariate)
animal_2 <- animal_map_2[trait_subset_bivariate$ringnr]

# Response matrix (No x 3 for bivariate with sex as fixed effect)
Y_2 <- matrix(0, No_2, 3)
Y_2[,1] <- trait_subset_bivariate$ving_h
Y_2[,2] <- trait_subset_bivariate$nebb_l
Y_2[,3] <- trait_subset_bivariate$scriptsex

# Bundle the data
dataset_bivariate <- list(
  No = No_2,
  Y = Y_2,
  animal = animal_2,
  Na = Na_2,
  A = as.matrix(A_hestmannoy_bivariate)  # Dense for Stan
)

## Analyse the data with the bivariate animal model
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# MCMC settings (increased for bivariate convergence)
nc_2 <- 4
nw_2 <- 3000
ni_2 <- 6000
nt_2 <- 1

# Parameters to monitor
tomonitor_bivariate <- c('mu', 'Sigma_A', 'Sigma_E', 'Sigma_R', 'heritability1', 'heritability2', 'cor_A', 'cor_E', 'cor_R')

# Call Stan from R
out_bivariate <- stan(file = 'Bivariate_Animal_Model.stan',
            data = dataset_bivariate, 
            pars = tomonitor_bivariate,
            chains = nc_2, iter = ni_2, warmup = nw_2, thin = nt_2,
            open_progress = FALSE,
            refresh = 10,
            seed = 123) # Change seed if needed

# Save the results
saveRDS(out_bivariate, 'Output_Bivariate_Animal_Model.rds')


## --------------------------------------------------------------
## 2. Write a human-readable summary to CSV
summ_biv <- summary(out_bivariate)$summary
write.csv(
  as.data.frame(summ_biv),
  file = "Output_Bivariate_Summary.csv",
  row.names = TRUE
)

cat("Summary written to: Output_Bivariate_Summary.csv\n")

## --------------------------------------------------------------
## 3. Export a *static* shinystan HTML report (no GUI needed)
library(shinystan)

# Create the shinystan object
shy_biv <- as.shinystan(out_bivariate)

# Export to a folder called "shinystan_report_bivariate"
#   - `overwrite = TRUE`  → replace if it already exists
#   - `open = FALSE`      → **do NOT launch the browser**
# export_shinystan(
#   shy_biv,
#   dir = "shinystan_report_bivariate", # NOT SURE THIS WORKS 
#   overwrite = TRUE,
#   open = FALSE
# )

cat("Shinystan HTML report written to folder: shinystan_report_bivariate\n")
cat("   → Download the whole folder and open index.html in any browser.\n")

# # Take a look at the result summary

#  out_bivariate
# 
# # Take a closer look at the MCMC chains, model diagnostics, and outputs:
# library(shinystan)
# launch_shinystan(out_bivariate)







# library(shinystan)
# 
# out_22 <- readRDS('C:/Users/Elias Ovesen/OneDrive/Skrivebord/Prosjektoppgave/output_markov/Output_Bivariate_Animal_Model.rds')
# launch_shinystan(out_22)


























# OLD VERSION OF CODE 





# Test of univarite animal model
# 
# library(MCMCglmm)
# library(lme4)
# library(brms)
# library(AGHmatrix)
# library(pedtricks)
# library(MASS)  # for mvrnorm
# library(tidyverse)
# library(nadiv)
# 
# 
# 
# traitData <- read.csv("C:\\Users\\Elias Ovesen\\OneDrive\\Skrivebord\\Prosjektoppgave\\data\\morphology.txt",
#                       sep = ";",
#                       header = TRUE,
#                       na.strings = c("NA", ""),
#                       fileEncoding = "Windows-1252",
#                       stringsAsFactors = FALSE,
#                       colClasses = c(ringnr = "character"))
# 
# 
# pedData <- read.csv("C:\\Users\\Elias Ovesen\\OneDrive\\Skrivebord\\Prosjektoppgave\\data\\pedigree_data.txt",
#                       sep = " ",
#                       header = TRUE,
#                       na.strings = "NA")
# 
# 
# # pedigree_full <- pedData %>%
# #   dplyr::select(ringnr, dam, sire) %>% 
# #   mutate(
# #     dam = if_else(is.na(dam), "UNKNOWN", dam),
# #     sire = if_else(is.na(sire), "UNKNOWN", sire)
# # )
# #   #%>% mutate(across(c(dam, sire), ~replace_na(., "0")))
# 
# 
# 
# # # Replace NA with explicit missing code
# # pedigree_clean <- pedData %>%
# #   dplyr::select(ringnr, dam, sire) %>%
# #   dplyr::mutate(
# #     dam = if_else(is.na(dam), "UNKNOWN", dam),
# #     sire = if_else(is.na(sire), "UNKNOWN", sire)
# #   )
# 
# # Integrate this before other cleaning...
# 
# # pedigree_clean <- pedData %>%
# #   dplyr::select(ringnr, dam, sire) %>%
# #   dplyr::mutate(
# #     # Strip suffixes from parents (adjust regex if pattern varies, e.g., if some have no "_")
# #     dam = gsub("_.*$", "", dam),  # Remove _ and after
# #     sire = gsub("_.*$", "", sire),
# #     # Optional: Ensure ringnr is clean too, if inconsistent
# #     ringnr = gsub("_.*$", "", ringnr)
# #   ) %>%
# #   # Then proceed with NA replacement as in Issue 1
# #   dplyr::mutate(
# #     dam = if_else(is.na(dam) | dam == "", "UNKNOWN", dam),
# #     sire = if_else(is.na(sire) | sire == "", "UNKNOWN", sire)
# #   )
# # 
# # # Now apply external handling from Issue 2 and compute A
# # 
# # 
# # 
# # 
# # # After cleaning NA as above...
# # 
# # # Extract unique parents (after any suffix stripping—see Issue 3)
# # all_parents <- unique(c(pedigree_clean$dam, pedigree_clean$sire))
# # external_parents <- setdiff(all_parents, c(pedigree_clean$ringnr, "UNKNOWN"))  # Exclude real IDs and missing
# # 
# # # Create founder rows if any externals
# # if (length(external_parents) > 0) {
# #   founder_df <- data.frame(
# #     ringnr = external_parents,
# #     dam = "UNKNOWN",
# #     sire = "UNKNOWN"
# #   )
# #   pedigree_full <- bind_rows(founder_df, pedigree_clean) %>%
# #     distinct(ringnr, .keep_all = TRUE)  # Dedup if needed
# # } else {
# #   pedigree_full <- pedigree_clean
# # }
# # 
# # # Now compute—verify should pass if all are now included
# # A_full <- AGHmatrix::Amatrix(pedigree_full, missing = "UNKNOWN", verify = FALSE)
# 
# 
# 
# 
# # # Your preprocessing
# # pedigree_clean <- pedData %>%
# #   dplyr::select(ringnr, dam, sire) %>%
# #   mutate(
# #     dam = gsub("_.*$", "", dam),
# #     sire = gsub("_.*$", "", sire),
# #     ringnr = gsub("_.*$", "", ringnr),
# #     dam = if_else(is.na(dam) | dam == "", "UNKNOWN", dam),
# #     sire = if_else(is.na(sire) | sire == "", "UNKNOWN", sire)
# #   )
# # 
# # # Add external parents as founders
# # all_parents <- unique(c(pedigree_clean$dam, pedigree_clean$sire))
# # external_parents <- setdiff(all_parents, c(pedigree_clean$ringnr, "UNKNOWN"))
# # 
# # if (length(external_parents) > 0) {
# #   founder_df <- data.frame(
# #     ringnr = external_parents,
# #     dam = "UNKNOWN",
# #     sire = "UNKNOWN"
# #   )
# #   pedigree_full <- bind_rows(founder_df, pedigree_clean) %>%
# #     distinct(ringnr, .keep_all = TRUE)
# # } else {
# #   pedigree_full <- pedigree_clean
# # }
# # 
# # # Diagnostic checks
# # cat("Rows in pedigree:", nrow(pedigree_full), "\n")
# # cat("Unique ringnr:", length(unique(pedigree_full$ringnr)), "\n")
# # cat("Any duplicate ringnr:", any(duplicated(pedigree_full$ringnr)), "\n")
# # cat("Column types:\n")
# # str(pedigree_full)
# # cat("Sample of data:\n")
# # head(pedigree_full, 10)
# # cat("Unique dams:", length(unique(pedigree_full$dam)), "\n")
# # cat("Unique sires:", length(unique(pedigree_full$sire)), "\n")
# # 
# # # Try computing A matrix
# # tryCatch({
# #   A_full <- Amatrix(pedigree_full, missing = "UNKNOWN", verify = FALSE)
# #   cat("A matrix computed successfully. Dimensions:", dim(A_full), "\n")
# # }, error = function(e) {
# #   cat("Error:", e$message, "\n")
# # })
# 
# 
# 
# 
# 
# # 
# # pedigree_clean <- pedData %>%
# #   dplyr::select(ringnr, dam, sire) %>%
# #   dplyr::mutate(
# #     # Ensure character, trim whitespace, remove suffixes
# #     across(c(ringnr, dam, sire), ~ as.character(trimws(gsub("_.*$", "", .)))),
# #     # Replace NA, empty, or invalid
# #     dam = if_else(is.na(dam) | dam == "" | dam == "NA", "UNKNOWN", dam),
# #     sire = if_else(is.na(sire) | sire == "" | sire == "NA", "UNKNOWN", sire)
# #   )
# # 
# # # Add external parents
# # all_parents <- unique(c(pedigree_clean$dam, pedigree_clean$sire))
# # external_parents <- setdiff(all_parents, c(pedigree_clean$ringnr, "UNKNOWN"))
# # 
# # if (length(external_parents) > 0) {
# #   founder_df <- data.frame(
# #     ringnr = external_parents,
# #     dam = "UNKNOWN",
# #     sire = "UNKNOWN"
# #   )
# #   pedigree_full <- bind_rows(founder_df, pedigree_clean) %>%
# #     distinct(ringnr, .keep_all = TRUE)
# # } else {
# #   pedigree_full <- pedigree_clean
# # }
# 
# # # Additional diagnostics
# # cat("Any non-alphanumeric ringnr:", any(grepl("[^A-Za-z0-9]", pedigree_full$ringnr)), "\n")
# # cat("Any empty ringnr:", any(pedigree_full$ringnr == ""), "\n")
# # cat("Rows with both parents UNKNOWN:", sum(pedigree_full$dam == "UNKNOWN" & pedigree_full$sire == "UNKNOWN"), "\n")
# # 
# # # Try computing
# # tryCatch({
# #   A_full <- Amatrix(pedigree_full, missing = "UNKNOWN", verify = FALSE, sparseform = TRUE)  # Sparse for large pedigree
# #   cat("A matrix computed successfully. Dimensions:", dim(A_full), "\n")
# # }, error = function(e) {
# #   cat("Error:", e$message, "\n")
# # })
# # 
# # 
# # library(nadiv)
# # 
# # # Your existing preprocessing
# # pedigree_clean <- pedData %>%
# #   dplyr::select(ringnr, dam, sire) %>%
# #   mutate(
# #     across(c(ringnr, dam, sire), ~ as.character(trimws(gsub("_.*$", "", .)))),
# #     dam = if_else(is.na(dam) | dam == "" | dam == "NA", "UNKNOWN", dam),
# #     sire = if_else(is.na(sire) | sire == "" | sire == "NA", "UNKNOWN", sire)
# #   )
# # 
# # # Add external parents (as you did)
# # all_parents <- unique(c(pedigree_clean$dam, pedigree_clean$sire))
# # external_parents <- setdiff(all_parents, c(pedigree_clean$ringnr, "UNKNOWN"))
# # 
# # if (length(external_parents) > 0) {
# #   founder_df <- data.frame(
# #     ringnr = external_parents,
# #     dam = "UNKNOWN",
# #     sire = "UNKNOWN"
# #   )
# #   pedigree_full <- bind_rows(founder_df, pedigree_clean) %>%
# #     distinct(ringnr, .keep_all = TRUE)
# # } else {
# #   pedigree_full <- pedigree_clean
# # }
# # 
# # # Use prepPed to fix missing parents and sort
# # pedigree_prepped <- prepPed(pedigree_full, gender = NULL, check = TRUE)
# 
# 
# 
# 
# ######
# 
# 
# 
# # # Use pedigree_prepped from your successful nadiv run
# # overlap_ids <- intersect(
# #   pedigree_prepped$dam[!pedigree_prepped$dam %in% c("UNKNOWN", NA)],
# #   pedigree_prepped$sire[!pedigree_prepped$sire %in% c("UNKNOWN", NA)]
# # )
# # 
# # cat("Number of IDs appearing as both dam and sire (excluding UNKNOWN/NA):", length(overlap_ids), "\n")
# # if (length(overlap_ids) > 0) {
# #   cat("Overlapping IDs:\n")
# #   print(overlap_ids)
# #   overlap_rows <- pedigree_prepped %>%
# #     filter(dam %in% overlap_ids | sire %in% overlap_ids) %>%
# #     arrange(dam, sire)
# #   cat("Rows with overlapping IDs:\n")
# #   print(head(overlap_rows, 20))
# # } else {
# #   cat("No overlapping IDs found (excluding UNKNOWN/NA).\n")
# # }
# # 
# # # Check specifically for 8L19505
# # cat("Rows where 8L19505 is dam or sire:\n")
# # print(pedigree_prepped %>% filter(dam == "8L19505" | sire == "8L19505"))
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
# # # If overlaps found, replace sire instances with UNKNOWN
# # if (length(overlap_ids) > 0) {
# #   pedigree_fixed <- pedigree_prepped %>%
# #     mutate(
# #       sire = if_else(sire %in% overlap_ids, "UNKNOWN", sire)
# #     )
# # } else {
# #   pedigree_fixed <- pedigree_prepped
# # }
# # 
# # # Re-run prepPed and makeA
# # library(nadiv)
# # pedigree_fixed_prepped <- prepPed(pedigree_fixed, gender = NULL, check = TRUE)
# # tryCatch({
# #   A_nadiv <- makeA(pedigree_fixed_prepped)
# #   A_nadiv_dense <- as.matrix(A_nadiv)  # Optional
# #   cat("Fixed A matrix computed. Dimensions:", dim(A_nadiv_dense), "\n")
# #   cat("Diagonal summary:", summary(diag(A_nadiv_dense)), "\n")
# #   cat("Parent-offspring check (8L19506, 8L19505):", A_nadiv_dense["8L19506", "8L19505"], "\n")
# #   cat("Sibling check (8L19507, 8L19508):", A_nadiv_dense["8L19507", "8L19508"], "\n")
# # }, error = function(e) {
# #   cat("nadiv Error:", e$message, "\n")
# # })
# 
# 
# 
# 
# 
# 
# 
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
# library(nadiv)
# 
# # Run prepPed with NA as missing
# pedigree_prepped <- prepPed(pedigree_full, gender = NULL, check = TRUE)
# 
# # # Check for UNKNOWN in ringnr
# # cat("Any ringnr = UNKNOWN:", any(pedigree_prepped$ringnr == "UNKNOWN"), "\n")
# 
# # Re-run makeA
# tryCatch({
#   A_nadiv <- makeA(pedigree_prepped)
#   A_nadiv_dense <- as.matrix(A_nadiv)
#   cat("New A matrix computed. Dimensions:", dim(A_nadiv_dense), "\n")
#   cat("Diagonal summary:", summary(diag(A_nadiv_dense)), "\n")
#   cat("Parent-offspring check (8L19506, 8L19505):", A_nadiv_dense["8L19506", "8L19505"], "\n")
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
# # # Only those from hestmannøy
# 
# # # Step 1: Identify phenotyped individuals meeting criteria
# # phenotyped_ids <- traitData %>%
# #   filter(sted_r == "hestmannøy" & !is.na(ving_h)) %>%
# #   pull(ringnr) %>%
# #   unique()  # Ensure unique IDs
# # 
# # # Step 2: Prune pedigree to include only phenotyped individuals and their ancestors
# # # (Requires full pedigree; assumes pedigree_prepped has all individuals with ringnr, dam, sire)
# # # Use nadiv::prunePed() or manual recursion; here, a manual function for completeness
# # 
# # # Function to get all ancestors of a set of individuals
# # get_ancestors <- function(ids, pedigree) {
# #   ancestors <- ids
# #   to_process <- ids
# #   while (length(to_process) > 0) {
# #     parents <- unique(c(pedigree$dam[pedigree$ringnr %in% to_process], pedigree$sire[pedigree$ringnr %in% to_process]))
# #     parents <- parents[!parents %in% ancestors & !is.na(parents)]
# #     ancestors <- c(ancestors, parents)
# #     to_process <- parents
# #   }
# #   return(ancestors)
# # }
# # 
# # # Get all relevant IDs (phenotyped + ancestors)
# # all_relevant_ids <- get_ancestors(phenotyped_ids, pedigree_prepped)
# # 
# # # Step 3: Create smaller pedigree
# # smaller_pedigree <- pedigree_prepped %>%
# #   filter(ringnr %in% all_relevant_ids) %>%
# #   arrange(match(ringnr, all_relevant_ids))  # Optional: Order by appearance for efficiency
# # 
# # # Step 4: Compute smaller A matrix (subset to phenotyped)
# # library(nadiv)
# # smaller_ped_prepped <- prepPed(smaller_pedigree, gender = NULL, check = TRUE)
# # A_small_nadiv <- makeA(smaller_ped_prepped)
# # 
# # # Subset to phenotyped individuals only (final smaller A)
# # # A_small <- A_small_nadiv[phenotyped_ids, phenotyped_ids]
# # 
# # # Output summary
# # cat("Smaller pedigree rows:", nrow(smaller_pedigree), "\n")
# # cat("Smaller A dimensions:", dim(A_small_nadiv), "\n")  # Should be length(phenotyped_ids) x length(phenotyped_ids)
# # 
# # 
# 
# 
# 
# # GET ONLY THE BIRDS FROM HESTMANNØY
# 
# library(dplyr)
# library(nadiv)
# library(Matrix)
# 
# # Step 1: Identify phenotyped individuals meeting criteria
# phenotyped_ids <- traitData %>%
#   filter(sted_r == "hestmannøy" & !is.na(ving_h)) %>%
#   mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr)))) %>%  # Match pedigree preprocessing
#   pull(ringnr) %>%
#   unique()
# 
# # Step 2: Ensure phenotyped_ids are in pedigree_prepped
# valid_phenotyped_ids <- phenotyped_ids[phenotyped_ids %in% pedigree_prepped$ringnr]
# if (length(valid_phenotyped_ids) < length(phenotyped_ids)) {
#   cat("Warning: Some phenotyped IDs not in pedigree:", setdiff(phenotyped_ids, pedigree_prepped$ringnr), "\n")
# }
# 
# # Step 3: Get all ancestors for valid phenotyped IDs
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
# # Step 4: Create smaller pedigree
# smaller_pedigree <- pedigree_prepped %>%
#   filter(ringnr %in% all_relevant_ids) %>%
#   arrange(match(ringnr, all_relevant_ids))
# 
# # Step 5: Compute smaller A matrix
# smaller_ped_prepped <- prepPed(smaller_pedigree, gender = NULL, check = TRUE)
# A_small_nadiv <- makeA(smaller_ped_prepped)
# 
# # Step 6: Subset to valid phenotyped individuals
# A_hestmannoy <- A_small_nadiv[valid_phenotyped_ids, valid_phenotyped_ids]
# 
# # Step 7: Diagnostics
# cat("1. Phenotyped individuals (after filtering):", length(valid_phenotyped_ids), "\n")
# cat("2. Smaller pedigree rows:", nrow(smaller_pedigree), "\n")
# cat("3. A_small dimensions:", dim(A_hestmannoy), "\n")
# cat("4. A_small symmetric:", isSymmetric(A_hestmannoy), "\n")
# cat("5. Diagonal summary (inbreeding):", summary(diag(as.matrix(A_hestmannoy))), "\n")
# cat("6. All phenotyped IDs in A_small:", all(valid_phenotyped_ids %in% rownames(A_hestmannoy)), "\n")
# cat("7. Positive semi-definite:", all(eigen(as.matrix(A_hestmannoy), only.values = TRUE)$values >= 0), "\n")
# cat("8. Sample parent-offspring (if available, e.g., 8L19506, 8L19505):",
#     if ("8L19506" %in% valid_phenotyped_ids & "8L19505" %in% valid_phenotyped_ids) 
#       A_hestmannoy["8L19506", "8L19505"] else "Not in valid_phenotyped_ids", "\n")
# cat("9. Sample sibling (if available, e.g., 8L19507, 8L19508):",
#     if ("8L19507" %in% valid_phenotyped_ids & "8L19508" %in% valid_phenotyped_ids) 
#       A_hestmannoy["8L19507", "8L19508"] else "Not in valid_phenotyped_ids", "\n")
# cat("10. No UNKNOWN in ringnr:", !any(smaller_pedigree$ringnr == "UNKNOWN"), "\n")
# 
# 
# 
# 
# 
# 
# 
# 
# ## ====================
# ## C) Run the analysis
# ## ====================
# 
# # Bundle the data
# dataset <- list(
#   No = No,
#   Y = Y,
#   animal = animal,
#   Na = Na,
#   A = as.matrix(A_small)
# )
# 
# ## Analyse the data with the simple animal model
# library(rstan)
# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())
# 
# # MCMC settings
# nc <- 4
# nw <- 1000
# ni <- 2000
# nt <- 1
# 
# # parameters to monitor
# tomonitor <- c('mu', 'var_A', 'var_E', 'var_R', 'var_P', 'evolvability', 'heritability')
# 
# # Call Stan from R
# # It takes 2-3 mins to compile and start I think, 
# # And then it takes ~15 mins to run on my computer.
# out <- stan(file = 'C:\\Users\\Elias Ovesen\\OneDrive\\Skrivebord\\Prosjektoppgave\\Avansert biologi\\Session 3\\Session 3 Animal Model.stan',,
#             data = dataset, 
#             pars = tomonitor,
#             chains = nc, iter = ni, warmup = nw, thin = nt,
#             open_progress = F,
#             seed = seed) # change the seed if you don't want to replicate the same sampling
# 
# # Save the results
# save(out, file = 'Output Session 3 Animal Model.Rdata')
# 
# # load the results
# load('Output Session 3 Animal Model.Rdata')
# 
# # Take a look at the result summary
# out
# 
# # Take a closer look at the MCMC chains, model diagnostics, and outputs:
# library(shinystan)
# launch_shinystan(out)
# 
# 
# # 
# # library(dplyr)
# # 
# # # Use the second pedigree_prepped (with NA as missing)
# # tryCatch({
# #   # Recompute for clarity (should match your last run)
# #   A_nadiv <- makeA(pedigree_prepped)
# #   A_nadiv_dense <- as.matrix(A_nadiv)
# #   cat("A matrix computed. Dimensions:", dim(A_nadiv_dense), "\n")
# #   cat("Diagonal summary:", summary(diag(A_nadiv_dense)), "\n")
# #   cat("Parent-offspring check (8L19506, 8L19505):", A_nadiv_dense["8L19506", "8L19505"], "\n")
# #   cat("Sibling check (8L19507, 8L19508):", A_nadiv_dense["8L19507", "8L19508"], "\n")  # Full siblings via dam 8L30947
# #   cat("8L19505 self-relationship:", A_nadiv_dense["8L19505", "8L19505"], "\n")  # Should be 1 (founder)
# # }, error = function(e) {
# #   cat("nadiv Error:", e$message, "\n")
# # })
# 
# 
# # # Check for potential selfing (individuals as their own parents)
# # selfing_rows <- pedigree_prepped %>%
# #   filter(dam == ringnr | sire == ringnr)
# # cat("Rows with selfing (individual as own parent):\n")
# # print(selfing_rows)
# # 
# # # Check for parents sharing IDs (dam = sire for an individual)
# # same_parent_rows <- pedigree_prepped %>%
# #   filter(dam == sire & dam != "UNKNOWN" & !is.na(dam))
# # cat("Rows where dam = sire (non-missing):\n")
# # print(same_parent_rows)
# 
# 
# # pedigree_fixed <- pedigree_prepped %>%
# #   mutate(
# #     sire = if_else(dam == sire & dam != "UNKNOWN" & !is.na(dam), NA_character_, sire)
#   # )
# # pedigree_fixed_prepped <- prepPed(pedigree_fixed, gender = NULL, check = TRUE)
# # A_nadiv <- makeA(pedigree_fixed_prepped)
# 
# 
# 
# 
# # library(AGHmatrix)
# # tryCatch({
# #   A_full <- Amatrix(pedigree_prepped, missing = NA, verify = FALSE, sparseform = TRUE)
# #   cat("AGHmatrix A matrix computed. Dimensions:", dim(A_full), "\n")
# #   cat("Parent-offspring check (8L19506, 8L19505):", as.matrix(A_full)["8L19506", "8L19505"], "\n")
# # }, error = function(e) {
# #   cat("AGHmatrix Error:", e$message, "\n")
# # })
# 
# ######
# 
# 
# # # Compute A matrix
# # tryCatch({
# #   A_nadiv <- makeA(pedigree_prepped)
# #   # Convert to dense for validation (optional, sparse is memory-efficient)
# #   A_nadiv_dense <- as.matrix(A_nadiv)
# #   cat("nadiv A matrix computed. Dimensions:", dim(A_nadiv_dense), "\n")
# #   # Validate
# #   cat("Diagonal summary (inbreeding):", summary(diag(A_nadiv_dense)), "\n")
# #   cat("Parent-offspring check (8L19506, 8L19505):", A_nadiv_dense["8L19506", "8L19505"], "\n")
# # }, error = function(e) {
# #   cat("nadiv Error:", e$message, "\n")
# # })
# 
# 
# 
# 
# 
# ####
# 
# # Testing the custom A matrix function from the biology course
# # pedigree_custom <- pedigree_prepped %>% rename(animal = ringnr)
# # A_custom <- Amatrix(pedigree_custom)  # Your handcrafted function
# 
# # Add checks from Grok
# 
# ####
# 
# 
# 
# 
# #A_full[1:5, 1:5]
# 
# 
# # Make an A matrix that only contains the individuals that are in both the pedigree and morphology data
# pedigree_relevant <- pedData %>% 
#   dplyr::select(ringnr, dam, sire) %>%
#   filter(pedData$ringnr %in% traitData$ringnr) %>%
#   mutate(across(c(dam, sire), ~replace_na(., "0")))
# 
# A_relevant  <- AGHmatrix::Amatrix(pedigree_relevant)
# 
# 
# # Make an A matrix with just birds from Hestmannøy:
# 
# 
# 
# 
# # Some statistics on the data:
# 
# # There are 11048 birds in the pedigree, and 40813 birds in the morphology data, of which 22268 are unique individuals.
# 
# # How many birds are in both pedigree and morphology data (6181)
# length(intersect(traitData$ringnr, pedData$ringnr))
# 
# # Character list of the mutual ring numbers in the pedigree and the morphology data 
# # (i.e. the birds that are both in the pedigree and have measurements)
# mutualRingnr <- intersect(traitData$ringnr, pedData$ringnr) 
# 
# # How many measurements of the birds we know the pedigree to (14707)
# sum(traitData$ringnr %in% mutualRingnr) 
# 
# 
# 
# # Checking that Hestmannøy is the most frequent island in the morphology data (7507 occurrences).
# tab <- table(traitData$sted_r)
# most_frequent <- names(tab)[which.max(tab)]
# count <- max(tab)
# 
# most_frequent
# count





