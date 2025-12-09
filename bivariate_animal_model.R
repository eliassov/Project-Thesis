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


# Load and preprocess pedigree

pedData <- read.csv("pedigree_data.txt", sep = " ", header = TRUE, na.strings = "NA")

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
  founder_df <- data.frame(ringnr = external_parents, dam = NA_character_, sire = NA_character_)
  pedigree_full <- bind_rows(founder_df, pedigree_clean) %>% distinct(ringnr, .keep_all = TRUE)
} else {
  pedigree_full <- pedigree_clean
}

# Prepare pedigree for A matrix calculation
pedigree_prepped <- prepPed(pedigree_full, gender = NULL, check = TRUE)


# PREPARE TRAIT DATA (BIVARIATE)
traitData <- read.csv("morphology.txt", sep = ";", header = TRUE, 
                      na.strings = c("NA", ""), fileEncoding = "Windows-1252", 
                      stringsAsFactors = FALSE, colClasses = c(ringnr = "character"))

# Filter and clean data: filter for island + both traits + sex
trait_subset_temp <- traitData %>%
  filter(sted_r == "hestmannøy", 
         !is.na(ving_h), 
         !is.na(nebb_l)) %>%  
  
  mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr)))) %>%
  
  # Logic to clean sex column (priority: scriptsex > fieldsex)
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


# ALIGN DATA WITH PEDIGREE

# Identify valid IDs (intersection of clean data and pedigree)
valid_ids <- intersect(trait_subset_temp$ringnr, pedigree_prepped$ringnr)

# Create final data subset
trait_subset <- trait_subset_temp %>%
  filter(ringnr %in% valid_ids) 



#  BUILD subsetted A-matrix 

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

smaller_ped_prepped <- prepPed(smaller_pedigree, gender = NULL, check = TRUE)
A_small_nadiv <- makeA(smaller_ped_prepped)

# Subset A to valid phenotyped individuals
A_hestmannoy <- A_small_nadiv[valid_ids, valid_ids]


cat("Data and Matrix: No =", nrow(trait_subset), "Na =", length(valid_ids), "\n")


# PREPARE STAN DATA

# Map animals to integers

matrix_ids <- rownames(A_hestmannoy)
animal_map <- setNames(seq_len(length(matrix_ids)), matrix_ids)
animal_idx <- animal_map[trait_subset$ringnr]

# Prepare response Y matrix (No x 2)
Y_mat <- as.matrix(trait_subset[, c("ving_h", "nebb_l")])

# Prepare fixed effects matrix X (intercept + sex)
X_mat <- model.matrix(~ sex, data = trait_subset)
K <- ncol(X_mat)

dataset_bivariate <- list(
  No = nrow(trait_subset),
  Y = Y_mat,
  X = X_mat,       
  K = K,          
  animal = animal_idx,
  Na = length(valid_ids),
  A = as.matrix(A_hestmannoy)
)


# RUNNING STAN

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# MCMC settings
nc <- 4
nw <- 3000
ni <- 6000
nt <- 1

# Parameters to monitor 
tomonitor <- c('beta', 'Sigma_A', 'Sigma_E', 'Sigma_R', 
               'heritability1', 'heritability2', 
               'cor_A', 'cor_E', 'cor_R', 'sd_A', 'sd_E', 'sd_R')

out_bivariate <- stan(file = 'Bivariate_Animal_Model.stan',
                      data = dataset_bivariate, 
                      pars = tomonitor, # Optional, but it doesn't take too much longer to monitor all parameters, and then you get breeding values and environmental random effects
                      chains = nc, iter = ni, warmup = nw, thin = nt,
                      open_progress = FALSE,
                      seed = 123)

# Save outputs
saveRDS(out_bivariate, 'Output_Bivariate_Final_Relevant.rds')

summ_biv <- summary(out_bivariate)$summary
write.csv(as.data.frame(summ_biv), file = "Output_Bivariate_Summary_Final_Relevant.csv", row.names = TRUE)

cat("Bivariate analysis complete. Summary saved.\n")
