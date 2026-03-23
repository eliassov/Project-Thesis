# ==============================================================================
# SCRIPT 1: DATA PREPARATION
# ==============================================================================
# Run this locally to clean data and build the Stan dataset list.

library(tidyverse)
library(nadiv)
library(MCMCglmm)

# --- 1. PEDIGREE ---
pedData <- read.csv("data/pedigree_data.txt", sep = " ", header = TRUE, na.strings = "NA")
pedigree_clean <- pedData %>%
  dplyr::select(ringnr, dam, sire) %>%
  mutate(across(c(ringnr, dam, sire), ~ as.character(trimws(gsub("_.*$", "", .)))),
         dam = if_else(is.na(dam) | dam == "" | dam == "NA", NA_character_, dam),
         sire = if_else(is.na(sire) | sire == "" | sire == "NA", NA_character_, sire))

all_parents <- unique(c(pedigree_clean$dam, pedigree_clean$sire))
external_parents <- setdiff(all_parents, c(pedigree_clean$ringnr, NA_character_))
if (length(external_parents) > 0) {
  founder_df <- data.frame(ringnr = external_parents, dam = NA_character_, sire = NA_character_)
  pedigree_full <- bind_rows(founder_df, pedigree_clean) %>% distinct(ringnr, .keep_all = TRUE)
} else {
  pedigree_full <- pedigree_clean
}
pedigree_prepped <- prepPed(pedigree_full, gender = NULL, check = TRUE)

# --- 2. MORPHOLOGY ---
traitData <- read.csv("data/morphology.txt", sep = ";", header = TRUE, 
                      na.strings = c("NA", ""), fileEncoding = "Windows-1252", 
                      stringsAsFactors = FALSE, colClasses = c(ringnr = "character"))

trait_raw <- traitData %>%
  filter(sted_r == "hestmannøy", 
         !is.na(ving_h), !is.na(nebb_l), !is.na(tars_h), !is.na(nebb_h), !is.na(vekt)) %>%
  mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr))))

valid_ids <- intersect(trait_raw$ringnr, pedigree_prepped$ringnr)
trait_subset_temp <- trait_raw %>% filter(ringnr %in% valid_ids)

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

# Outliers
traits_to_check <- c("ving_h", "tars_h", "nebb_l", "nebb_h", "vekt")
bad_rows_global <- c() 
for(t in traits_to_check) {
  f <- as.formula(paste(t, "~ sex + age_class"))
  fit <- lm(f, data = trait_subset_temp, na.action = na.exclude)
  resids <- rstudent(fit)
  cutoff <- if(t == "vekt") 6.0 else 4.0 
  bad_rows_trait <- which(abs(resids) > cutoff)
  if(length(bad_rows_trait) > 0) {
    bad_rows_global <- c(bad_rows_global, bad_rows_trait)
  }
}

if(length(bad_rows_global) > 0) {
  trait_subset_lv <- trait_subset_temp[-unique(bad_rows_global), ]
} else {
  trait_subset_lv <- trait_subset_temp
}

# Transformations & Grouping
trait_subset_lv <- trait_subset_lv %>%
  mutate(
    ving_h_std = as.numeric(scale(log(ving_h))),
    nebb_l_std = as.numeric(scale(log(nebb_l))),
    tars_h_std = as.numeric(scale(log(tars_h))),
    nebb_h_std = as.numeric(scale(log(nebb_h))),
    vekt_std   = as.numeric(scale(log(vekt)))
  ) %>%
  mutate(
    date_obj = as.Date(dato),
    month = as.numeric(format(date_obj, "%m")),
    Year = factor(format(date_obj, "%Y")),
    Season_3 = factor(case_when(
      month %in% c(5, 6, 7) ~ "Breeding",
      month %in% c(8, 9)    ~ "Molt",
      TRUE                  ~ "Winter"
    ))
  ) %>%
  group_by(init) %>%
  mutate(init_clean = case_when(
    n() < 10 ~ "OTHER",
    TRUE     ~ as.character(init)
  )) %>%
  ungroup() %>%
  mutate(init_clean = factor(init_clean))

# --- 3. FITNESS ---
fitness_data <- read.csv("data/Fitness.csv", header = TRUE, stringsAsFactors = FALSE, na.strings = "NA", fileEncoding = "Windows-1252")
fitness_data_clean <- fitness_data %>%
  dplyr::select(Location, Year_ID, Year, Sex, Own_survival, N_Recruits, Least_age) %>%
  filter(grepl("^Hestmann", Location)) %>%
  mutate(
    ringnr = sub("^.*_", "", Year_ID),
    Year = sub("_.*$", "", Year_ID), 
    age_class = case_when(
      Least_age <= 1 ~ "juvenile",
      Least_age > 1 ~ "adult",
      TRUE ~ NA_character_
    ),
    Own_survival = as.numeric(Own_survival)
  ) %>%
  filter(ringnr %in% pedigree_prepped$ringnr)

data_repro <- fitness_data_clean %>% drop_na(N_Recruits)
data_surv  <- fitness_data_clean %>% drop_na(Own_survival)

# --- 4. PEDIGREE PRUNING & INVERSE A ---
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
Na <- nrow(ped_ordered)
dam_idx_mapped  <- match(ped_ordered$dam, ped_ordered$id)
sire_idx_mapped <- match(ped_ordered$sire, ped_ordered$id)
dam_idx_mapped[is.na(dam_idx_mapped)]   <- Na + 1
sire_idx_mapped[is.na(sire_idx_mapped)] <- Na + 1
id_map <- setNames(1:Na, ped_ordered$id)

# --- 5. UNIVERSAL MAPS & MATRICES ---
animal_morph <- unname(id_map[trait_subset_lv$ringnr])
animal_repro <- unname(id_map[data_repro$ringnr])
animal_surv  <- unname(id_map[data_surv$ringnr])

all_years <- sort(unique(c(as.character(trait_subset_lv$Year), 
                           as.character(data_repro$Year), 
                           as.character(data_surv$Year))))
year_map <- setNames(1:length(all_years), all_years)

year_morph <- unname(year_map[as.character(trait_subset_lv$Year)])
year_repro <- unname(year_map[as.character(data_repro$Year)])
year_surv  <- unname(year_map[as.character(data_surv$Year)])

init_factor <- as.integer(factor(trait_subset_lv$init_clean))

X_morph <- model.matrix(~ 1, data = trait_subset_lv)
X_surv <- model.matrix(~ 1, data = data_surv)
X_repro <- model.matrix(~ 1, data = data_repro)
Y_mat <- as.matrix(trait_subset_lv[, c("tars_h_std", "ving_h_std", "nebb_l_std", "nebb_h_std", "vekt_std")])

# --- 6. COMPILE LIST AND EXPORT ---
dataset_joint <- list(
  No_morph = nrow(trait_subset_lv),
  No_repro = nrow(data_repro),
  No_surv  = nrow(data_surv),
  Nt = 5, 
  Nlv = 2,
  Na = Na,
  Y_morph = Y_mat,
  Y_repro = as.array(as.integer(data_repro$N_Recruits)),
  Y_surv  = as.array(as.integer(data_surv$Own_survival)),
  animal_morph = as.array(as.integer(animal_morph)),
  animal_repro = as.array(as.integer(animal_repro)),
  animal_surv  = as.array(as.integer(animal_surv)),
  dam  = as.array(as.integer(dam_idx_mapped)), 
  sire = as.array(as.integer(sire_idx_mapped)),
  dii  = as.array(A_inv_obj$dii),
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

# Save this so the HPC script can just load it!
saveRDS(dataset_joint, file = "data/clean_stan_data.rds")


cat("Data successfully prepped and saved to data/clean_stan_data.rds\n")