

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
# 2. LOAD & CLEAN DATA (OPTIMIZED ORDER)
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

# A. LOAD & LIGHT FORMATTING (Crucial for matching)
trait_raw <- traitData %>%
  filter(sted_r == "hestmannøy", 
         !is.na(ving_h), !is.na(nebb_l), !is.na(tars_h), !is.na(nebb_h), !is.na(vekt)) %>%
  mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr))))

# B. ALIGN WITH PEDIGREE (The "Gatekeeper" Step)
valid_ids <- intersect(trait_raw$ringnr, pedigree_prepped$ringnr)
trait_subset_temp <- trait_raw %>% filter(ringnr %in% valid_ids)





# ==============================================================================
# C. DEEP CLEANING (Sex, Age)
# ==============================================================================

trait_subset_temp <- trait_subset_temp %>%
  mutate(raw_sex = if_else(!is.na(scriptsex) & scriptsex != "u", scriptsex, fieldsex)) %>%
  mutate(sex = case_when(
    raw_sex %in% c("m", "pm") ~ "m", 
    raw_sex %in% c("f", "pf") ~ "f", 
    TRUE ~ NA_character_
  ))%>%
  filter(!is.na(sex))%>%
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
# ==============================================================================
# D. ROW-WISE STATISTICAL OUTLIER REMOVAL
# ==============================================================================
traits_to_check <- c("ving_h", "tars_h", "nebb_l", "nebb_h", "vekt")
bad_rows_global <- c() # We now collect row INDICES, not bird IDs

for(t in traits_to_check) {
  f <- as.formula(paste(t, "~ sex + age_class"))
  fit <- lm(f, data = trait_subset_temp, na.action = na.exclude)
  
  resids <- rstudent(fit)
  
  # Keeping your biologically sound 4 vs 6 SD rule
  cutoff <- if(t == "vekt") 6.0 else 4.0 
  
  # Find the exact row numbers that violated the cutoff
  bad_rows_trait <- which(abs(resids) > cutoff)
  
  if(length(bad_rows_trait) > 0) {
    bad_rows_global <- c(bad_rows_global, bad_rows_trait)
    cat("Trait:", t, "- Removing", length(bad_rows_trait), "specific measurements (rows)\n")
  }
}

# ==============================================================================
# E. FINAL DATASET
# ==============================================================================
# Remove only the specific rows with typos, keeping the rest of the bird's history
if(length(bad_rows_global) > 0) {
  trait_subset_lv <- trait_subset_temp[-unique(bad_rows_global), ]
} else {
  trait_subset_lv <- trait_subset_temp
}

cat("Final Clean N Rows: ", nrow(trait_subset_lv), "\n")




# NEW: Log-transform FIRST, then scale
trait_subset_lv <- trait_subset_lv %>%
  mutate(
    ving_h_std = as.numeric(scale(log(ving_h))),
    nebb_l_std = as.numeric(scale(log(nebb_l))),
    tars_h_std = as.numeric(scale(log(tars_h))),
    nebb_h_std = as.numeric(scale(log(nebb_h))),
    vekt_std   = as.numeric(scale(log(vekt)))
  )


cat("Data Loaded. N Rows:", nrow(trait_subset_lv), "\n")




# 1. Convert to Date (R handles YYYY-MM-DD automatically)
trait_subset_lv$date_obj <- as.Date(trait_subset_lv$dato)

# 2. Extract month
trait_subset_lv$month <- as.numeric(format(trait_subset_lv$date_obj, "%m"))



# 1. Extract the Year as a factor
trait_subset_lv <- trait_subset_lv %>%
  mutate(
    Year = factor(format(date_obj, "%Y")),
    # Create the 3-level Season we discussed
    Season_3 = factor(case_when(
      month %in% c(5, 6, 7) ~ "Breeding",
      month %in% c(8, 9)    ~ "Molt",
      TRUE                  ~ "Winter"
    ))
  )




# Grouping low-volume measurers using case_when
trait_subset_lv <- trait_subset_lv %>%
  group_by(init) %>%
  mutate(init_clean = case_when(
    n() < 10 ~ "OTHER",
    TRUE     ~ as.character(init)
  )) %>%
  ungroup() %>%
  mutate(init_clean = factor(init_clean))




# ==============================================================================
# 2.5 LOAD & CLEAN FITNESS DATA
# ==============================================================================

fitness_data <- read.csv("Fitness.csv", header = TRUE, stringsAsFactors = FALSE, na.strings = "NA", fileEncoding = "Windows-1252")

fitness_data_clean <- fitness_data %>%
  dplyr::select(Location, Year_ID, Year, Sex, Own_survival, N_Recruits, Least_age) %>%
  filter(grepl("^Hestmann", Location)) %>%
  mutate(
    # Extract the ring number
    ringnr = sub("^.*_", "", Year_ID),
    
    # NEW: Extract the actual calendar year and overwrite the 1-22 index
    Year = sub("_.*$", "", Year_ID), 
    
    age_class = case_when(
      Least_age <= 1 ~ "juvenile",
      Least_age > 1 ~ "adult",
      TRUE ~ NA_character_
    ),
    
    Own_survival = as.numeric(Own_survival)
  ) %>%
  # Gatekeeper: Must be in the prepared pedigree
  filter(ringnr %in% pedigree_prepped$ringnr)

# Split into two clean datasets (Stan cannot handle NA in the response variable)
data_repro <- fitness_data_clean %>% drop_na(N_Recruits)
data_surv  <- fitness_data_clean %>% drop_na(Own_survival)

cat("Morphology records: ", nrow(trait_subset_lv), "\n")
cat("Reproduction records: ", nrow(data_repro), "\n")
cat("Survival records: ", nrow(data_surv), "\n")




# ==============================================================================
# 3. PEDIGREE PROCESSING
# ==============================================================================
pedigree_working <- pedigree_prepped
colnames(pedigree_working) <- c("id", "dam", "sire") 
pedigree_working$id   <- as.character(pedigree_working$id)
pedigree_working$dam  <- as.character(pedigree_working$dam)
pedigree_working$sire <- as.character(pedigree_working$sire)

# NEW: Keep all birds that appear in ANY of the three datasets
birds_to_keep <- unique(c(trait_subset_lv$ringnr, data_repro$ringnr, data_surv$ringnr))

ped_pruned <- MCMCglmm::prunePed(pedigree_working, keep = birds_to_keep)

A_inv_obj <- MCMCglmm::inverseA(ped_pruned, nodes = "ALL", scale = FALSE)
ordered_ids <- rownames(A_inv_obj$Ainv)

ped_ordered <- data.frame(
  id = ordered_ids,
  dam = as.character(A_inv_obj$pedigree[, 2]), 
  sire = as.character(A_inv_obj$pedigree[, 3]),
  stringsAsFactors = FALSE
)

dii_vector <- A_inv_obj$dii
Na <- nrow(ped_ordered)

dam_idx_mapped  <- match(ped_ordered$dam, ped_ordered$id)
sire_idx_mapped <- match(ped_ordered$sire, ped_ordered$id)

dam_idx_mapped[is.na(dam_idx_mapped)]   <- Na + 1
sire_idx_mapped[is.na(sire_idx_mapped)] <- Na + 1

id_map <- setNames(1:Na, ped_ordered$id)





# ==============================================================================
# 4. DATA MAPPING & UNIVERSAL FACTORS
# ==============================================================================
# 1. Map IDs
animal_morph <- unname(id_map[trait_subset_lv$ringnr])
animal_repro <- unname(id_map[data_repro$ringnr])
animal_surv  <- unname(id_map[data_surv$ringnr])

# 2. Universal Year Map (Crucial so Year 1 is the same everywhere)
all_years <- sort(unique(c(as.character(trait_subset_lv$Year), 
                           as.character(data_repro$Year), 
                           as.character(data_surv$Year))))
year_map <- setNames(1:length(all_years), all_years)

year_morph <- unname(year_map[as.character(trait_subset_lv$Year)])
year_repro <- unname(year_map[as.character(data_repro$Year)])
year_surv  <- unname(year_map[as.character(data_surv$Year)])

# 3. Measurer Map
init_factor <- as.integer(factor(trait_subset_lv$init_clean))

X_morph <- model.matrix(~ 1, data = trait_subset_lv)
X_surv <- model.matrix(~ 1, data = data_surv)
X_repro <- model.matrix(~ 1, data = data_repro)

# 4. Trait Matrix
Y_mat <- as.matrix(trait_subset_lv[, c("tars_h_std", "ving_h_std", "nebb_l_std", "nebb_h_std", "vekt_std")])





# ==============================================================================
# 5. STAN EXECUTION FOR JOINT MODEL
# ==============================================================================
dataset_joint <- list(
  # Dimensions
  No_morph = nrow(trait_subset_lv),
  No_repro = nrow(data_repro),
  No_surv  = nrow(data_surv),
  Nt = 5, 
  Nlv = 3,
  Na = Na,
  
  # Responses
  Y_morph = Y_mat,
  Y_repro = as.array(as.integer(data_repro$N_Recruits)),
  Y_surv  = as.array(as.integer(data_surv$Own_survival)),
  
  # Animal Indices
  animal_morph = as.array(as.integer(animal_morph)),
  animal_repro = as.array(as.integer(animal_repro)),
  animal_surv  = as.array(as.integer(animal_surv)),
  
  # Pedigree
  dam  = as.array(as.integer(dam_idx_mapped)), 
  sire = as.array(as.integer(sire_idx_mapped)),
  dii  = as.array(dii_vector),
  
  # Random Effects Indices
  N_year = length(all_years),
  year_morph = as.array(as.integer(year_morph)),
  year_repro = as.array(as.integer(year_repro)),
  year_surv  = as.array(as.integer(year_surv)),
  
  N_init = max(init_factor),
  init_morph = as.array(as.integer(init_factor)),
  
  K_morph = ncol(X_morph),
  K_repro = ncol(X_repro),
  K_surv  = ncol(X_surv),
  X_morph = X_morph,
  X_repro = X_repro,
  X_surv  = X_surv
)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())



init_fn_joint <- function() {
  
  # 1. Initialize RAW Loadings (Instead of final Lambda)
  # Start Factor 1 with a positive push, and later factors near zero
  L_raw_init <- matrix(rnorm(dataset_joint$Nt * dataset_joint$Nlv, 0, 0.1), 
                       nrow = dataset_joint$Nt, ncol = dataset_joint$Nlv)
  L_raw_init[, 1] <- rnorm(dataset_joint$Nt, 0.5, 0.1) # Factor 1 (Size)
  
  list(
    # --- Raw Loadings and Gradients ---
    # lambda_raw = L_raw_init,
    lambda_raw_unanchored = L_raw_init,
    
    gamma_repro_raw = as.array(rnorm(dataset_joint$Nlv, 0, 0.1)),
    gamma_surv_raw  = as.array(rnorm(dataset_joint$Nlv, 0, 0.1)),
    
    # --- MGP Shrinkage Parameters (Must be > 0) ---
    # a_1 = runif(1, 1.5, 2.5),
    # a_2 = runif(1, 1.5, 2.5),
    delta = as.array(runif(dataset_joint$Nlv, 1.5, 2.5)),
    # phi   = matrix(runif(dataset_joint$Nt * dataset_joint$Nlv, 0.5, 1.5), 
    #                nrow = dataset_joint$Nt, ncol = dataset_joint$Nlv),
    # phi_r = as.array(runif(dataset_joint$Nlv, 0.5, 1.5)),
    # phi_s = as.array(runif(dataset_joint$Nlv, 0.5, 1.5)),
    
    # --- Base Biology & Variances ---
    h2_lv = as.array(runif(dataset_joint$Nlv, 0.4, 0.6)),
    sd_R  = as.array(runif(dataset_joint$Nt, 0.5, 1.0)),
    
    # --- Fixed Effects ---
    beta_morph = matrix(rnorm(dataset_joint$K_morph * dataset_joint$Nt, 0, 0.1), 
                        nrow = dataset_joint$K_morph, ncol = dataset_joint$Nt),
    beta_repro = as.array(rnorm(dataset_joint$K_repro, 0, 0.1)),
    beta_surv  = as.array(rnorm(dataset_joint$K_surv, 0, 0.1)),
    
    # --- Latent Genetic and Environmental Noise ---
    w_a  = matrix(rnorm(dataset_joint$Nlv * dataset_joint$Na, 0, 0.1), 
                  nrow = dataset_joint$Nlv, ncol = dataset_joint$Na),
    w_pe = matrix(rnorm(dataset_joint$Nlv * dataset_joint$Na, 0, 0.1), 
                  nrow = dataset_joint$Nlv, ncol = dataset_joint$Na),
    
    # --- Random Effect Variances ---
    sd_year_morph = as.array(runif(dataset_joint$Nt, 0.05, 0.2)),
    sd_init_morph = as.array(runif(dataset_joint$Nt, 0.05, 0.2)),
    
    # SCALARS 
    sd_year_repro = runif(1, 0.05, 0.2),
    sd_year_surv  = runif(1, 0.05, 0.2),
    
    # --- Random Effect Z-Scores ---
    z_year_morph = matrix(rnorm(dataset_joint$N_year * dataset_joint$Nt, 0, 0.1), 
                          nrow = dataset_joint$N_year, ncol = dataset_joint$Nt),
    z_year_repro = as.array(rnorm(dataset_joint$N_year, 0, 0.1)),
    z_year_surv  = as.array(rnorm(dataset_joint$N_year, 0, 0.1)),
    z_init_morph = matrix(rnorm(dataset_joint$N_init * dataset_joint$Nt, 0, 0.1), 
                          nrow = dataset_joint$N_init, ncol = dataset_joint$Nt)
  )
}






# Run the Joint Model
out_joint <- stan(
  file = "gamma_no_phi.stan",
  data = dataset_joint,
  init = init_fn_joint,        
  chains = 4, 
  pars = c(
             "beta_morph", "beta_repro", "beta_surv",
           "Lambda", "gamma_repro", "gamma_surv", 
           "sd_R", "h2_lv",
           
             "sd_year_morph", "sd_year_repro", "sd_year_surv", "sd_init_morph",
           
             "morph_prop_var_explained", 
           "morph_var_explained", "repro_var_explained", "surv_var_explained",
           
             "lambda_raw", "gamma_repro_raw", "gamma_surv_raw",
           
             # "a_1", "a_2", 
           "delta", "tau"
          # , "phi", "phi_r", "phi_s"
           , "lp__"
  ),
  control = list(adapt_delta = 0.95, max_treedepth = 15), 
  iter = 2000,          
  warmup = 1000
)

saveRDS(out_joint, 'Output_gamma_shrink_no_phi_or_a.rds')





# ==============================================================================
# 6. RESULTS & DIAGNOSTICS
# ==============================================================================
# out_lv <- readRDS('Output_gamma_shrink.rds')

# # Print the parameters specific to the Joint Starter model
# print(out_lv,
#       pars = c( "beta_morph", "beta_repro", "beta_surv",
#                 "Lambda", "gamma_repro", "gamma_surv", 
#                 "sd_R", "h2_lv",
#                 
#                 "sd_year_morph", "sd_year_repro", "sd_year_surv", "sd_init_morph",
#                 
#                 "morph_prop_var_explained", 
#                 "morph_var_explained", "repro_var_explained", "surv_var_explained",
#                 
#                 "lambda_raw", "gamma_repro_raw", "gamma_surv_raw",
#                 
#                 "a_1", "a_2", "delta", "tau", "phi", "phi_r", "phi_s", "lp__"),            
#       probs = c(0.025, 0.975),
#       digits_summary = 3)




