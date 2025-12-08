
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


pedigree_clean <- pedData %>%
  dplyr::select(ringnr, dam, sire) %>%
  mutate(
    across(c(ringnr, dam, sire), ~ as.character(trimws(gsub("_.*$", "", .)))),
    dam = if_else(is.na(dam) | dam == "" | dam == "NA", NA_character_, dam),
    sire = if_else(is.na(sire) | sire == "" | sire == "NA", NA_character_, sire)
  )

# Add parents
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


pedigree_prepped <- prepPed(pedigree_full, gender = NULL, check = TRUE)



tryCatch({
  A_nadiv <- makeA(pedigree_prepped)
  A_nadiv_dense <- as.matrix(A_nadiv)
  cat("New A matrix computed. Dimensions:", dim(A_nadiv_dense), "\n")
  cat("Diagonal summary:", summary(diag(A_nadiv_dense)), "\n")
}, error = function(e) {
  cat("nadiv Error:", e$message, "\n")
})








# CREATE CLEAN DATA SUBSET
trait_subset_temp <- traitData %>%
  # Island and traits
  filter(sted_r == "hestmannøy", 
      !is.na(ving_h),         
      !is.na(nebb_l)) %>%  # Added this filter so that it is the same exact A matrix for comparison with bivariate and latent variable model
  
  # Ring number standardized
  mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr)))) %>%
  
  # Sex logic
  mutate(
    raw_sex = if_else(!is.na(scriptsex) & scriptsex != "u", scriptsex, fieldsex)  # priority: scriptsex > fieldsex
  ) %>%
  mutate(
    sex = case_when(
      raw_sex %in% c("m", "pm") ~ "m",
      raw_sex %in% c("f", "pf") ~ "f",
      TRUE ~ NA_character_ # Maps 'u', NA, and typos to NA
    )
  ) %>%
  
  # Drop rows where sex is missing
  filter(!is.na(sex)) 


# DEFINE VALID IDs BASED ON CLEAN DATA

# We only keep those birds that are also in the pedigree.
valid_phenotyped_ids <- intersect(trait_subset_temp$ringnr, pedigree_prepped$ringnr)

# Create the final data object using only these valid IDs
trait_subset <- trait_subset_temp %>%
  filter(ringnr %in% valid_phenotyped_ids) 





# BUILD A MATRIX using only the valid IDs


# Function to get all ancestors recursively
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


# Get ancestry for only the valid individuals
all_relevant_ids <- get_ancestors(valid_phenotyped_ids, pedigree_prepped)

# Create smaller pedigree
smaller_pedigree <- pedigree_prepped %>%
  filter(ringnr %in% all_relevant_ids) %>%
  arrange(match(ringnr, all_relevant_ids))

# Compute smaller A matrix
smaller_ped_prepped <- prepPed(smaller_pedigree, gender = NULL, check = TRUE)
A_small_nadiv <- makeA(smaller_ped_prepped)

# Subset A matrix to match trait_subset
A_hestmannoy <- A_small_nadiv[valid_phenotyped_ids, valid_phenotyped_ids]




# PREPARE STAN DATA

No <- nrow(trait_subset)
Na <- length(valid_phenotyped_ids)


# Every animal in the data must exist in the matrix
if (!all(unique(trait_subset$ringnr) %in% rownames(A_hestmannoy))) {
  stop("Error: Data contains individuals not in the A-matrix")
}

# Every animal in the matrix should have data 
if (!all(rownames(A_hestmannoy) %in% trait_subset$ringnr)) {
  warning("Warning: A matrix contains individuals with no observations")
}

# Map animals to integers 1:Na
matrix_ids <- rownames(A_hestmannoy)
animal_map <- setNames(seq_len(length(matrix_ids)), matrix_ids)
animal_idx <- animal_map[trait_subset$ringnr]

# Prepare Y (response)
Y <- as.numeric(trait_subset$ving_h) 

# Prepare fixed effects matrix (X) with an intercept column (all 1's) and a sex column (0/1)
X_mat <- model.matrix(~ sex, data = trait_subset)

# Get number of fixed effect predictors (K)
K <- ncol(X_mat) 

dataset <- list(
  No = nrow(trait_subset),       # Number of observations
  Y = Y,
  X = X_mat,                     # Fixed effects design matrix
  K = K,                         # Number of fixed effects
  animal = animal_idx, # Map IDs to integers
  Na = length(valid_phenotyped_ids),        # Number of unique individuals
  A = as.matrix(A_hestmannoy)    
)




tomonitor <- c('beta', 'var_A', 'var_E', 'var_R', 'var_P', 'heritability')

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# MCMC settings
nc <- 4    # Number of Chains
nw <- 3000 # Warmup iterations
ni <- 6000 # Total iterations
nt <- 1    # Thinning



out_1 <- stan(file = 'animal_model_univariate.stan',
              data = dataset,
              # pars = tomonitor,  # Optional, but it doesn't take too much longer to monitor all parameters, and then you get breeding values and environmental random effects
              chains = nc, iter = ni, warmup = nw, thin = nt,
              open_progress = FALSE,
              seed = 123) 


saveRDS(out_1, file = "Output_Univariate_Animal_Model_Correct_Sorting.rds")
summ_uni <- summary(out_1)$summary
write.csv(as.data.frame(summ_uni), "Output_Univariate_Summary_Correct_Sorting.csv", row.names = TRUE)


cat("Summary written to: Output_Univariate_Summary_Correct_Sorting.csv\n")











