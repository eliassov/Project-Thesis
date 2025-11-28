now

# Set up personal library
personal_lib <- "~/R/library"
if (!dir.exists(personal_lib)) {
  dir.create(personal_lib, recursive = TRUE)
}
.libPaths(c(personal_lib, .libPaths()))

# Install missing packages to personal library
required_packages <- c("rstan", "tidyverse", "nadiv") 
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, lib = personal_lib, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}


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


# Run prepPed with NA as missing
pedigree_prepped <- prepPed(pedigree_full, gender = NULL, check = TRUE)



# Re-run makeA
tryCatch({
  A_nadiv <- makeA(pedigree_prepped)
  A_nadiv_dense <- as.matrix(A_nadiv)
  cat("New A matrix computed. Dimensions:", dim(A_nadiv_dense), "\n")
  cat("Diagonal summary:", summary(diag(A_nadiv_dense)), "\n")
}, error = function(e) {
  cat("nadiv Error:", e$message, "\n")
})



# GET ONLY THE BIRDS FROM HESTMANNØY (that are phenotyped)

# Identify phenotyped individuals meeting criteria
phenotyped_ids <- traitData %>%
  filter(sted_r == "hestmannøy" & !is.na(ving_h)) %>%
  mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr)))) %>%  # Match pedigree preprocessing
  pull(ringnr) %>%
  unique()

# Ensure phenotyped_ids are in pedigree_prepped
valid_phenotyped_ids <- phenotyped_ids[phenotyped_ids %in% pedigree_prepped$ringnr]


# Get all ancestors for valid phenotyped IDs
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

# Create smaller pedigree
smaller_pedigree <- pedigree_prepped %>%
  filter(ringnr %in% all_relevant_ids) %>%
  arrange(match(ringnr, all_relevant_ids))

# Compute smaller A matrix
smaller_ped_prepped <- prepPed(smaller_pedigree, gender = NULL, check = TRUE)
A_small_nadiv <- makeA(smaller_ped_prepped)

# Subset to valid phenotyped individuals
A_hestmannoy <- A_small_nadiv[valid_phenotyped_ids, valid_phenotyped_ids]


# Diagnostics
# cat("1. Phenotyped individuals (after filtering):", length(valid_phenotyped_ids), "\n")
# cat("2. Smaller pedigree rows:", nrow(smaller_pedigree), "\n")
# cat("3. A_small dimensions:", dim(A_hestmannoy), "\n")
# cat("4. A_small symmetric:", isSymmetric(A_hestmannoy), "\n")
# cat("5. Diagonal summary (inbreeding):", summary(diag(as.matrix(A_hestmannoy))), "\n")
# cat("6. All phenotyped IDs in A_small:", all(valid_phenotyped_ids %in% rownames(A_hestmannoy)), "\n")
# cat("7. Positive semi-definite:", all(eigen(as.matrix(A_hestmannoy), only.values = TRUE)$values >= 0), "\n")







# ====================

## Bivariate animal model


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



A_hestmannoy_bivariate <- A_small_nadiv[valid_phenotyped_ids_bivariate, valid_phenotyped_ids_bivariate]





# Extract real data for Hestmannøy (two traits: ving_h and nebb_l)
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
                      seed = 123)

# Save the results
saveRDS(out_bivariate, 'Output_Bivariate_Animal_Model.rds')


summ_biv <- summary(out_bivariate)$summary

write.csv(
  as.data.frame(summ_biv),
  file = "Output_Bivariate_Summary.csv",
  row.names = TRUE
)

cat("Summary written to: Output_Bivariate_Summary.csv\n")






