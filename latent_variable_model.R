# Set up personal library
personal_lib <- "~/R/library"
if (!dir.exists(personal_lib)) {
  dir.create(personal_lib, recursive = TRUE)
}
.libPaths(c(personal_lib, .libPaths()))

# Install missing packages
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



# Load data
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

# Preprocess pedigree
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

# Prepare pedigree and compute full A matrix
pedigree_prepped <- prepPed(pedigree_full, gender = NULL, check = TRUE)
A_nadiv <- makeA(pedigree_prepped)

# Filter for Hestmannøy birds with traits
phenotyped_ids <- traitData %>%
  filter(sted_r == "hestmannøy" & !is.na(ving_h) & !is.na(nebb_l)) %>%
  mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr)))) %>%
  pull(ringnr) %>%
  unique()

valid_phenotyped_ids <- phenotyped_ids[phenotyped_ids %in% pedigree_prepped$ringnr]

# Get ancestors
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

# Subset A to phenotyped individuals
A_hestmannoy <- A_small_nadiv[valid_phenotyped_ids, valid_phenotyped_ids]

# Prepare data for Stan (two traits: ving_h and nebb_l)
trait_subset_lv <- traitData %>%
  filter(ringnr %in% valid_phenotyped_ids, sted_r == "hestmannøy", !is.na(ving_h), !is.na(nebb_l)) %>%
  mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr)))) %>%
  mutate(
    ving_h_std = (ving_h - mean(ving_h)) / sd(ving_h),
    nebb_l_std = (nebb_l - mean(nebb_l)) / sd(nebb_l)
  )

No_lv <- nrow(trait_subset_lv)
Nt_lv <- 2
Na_lv <- length(valid_phenotyped_ids)
animal_map_lv <- setNames(seq_len(Na_lv), valid_phenotyped_ids)
animal_lv <- animal_map_lv[trait_subset_lv$ringnr]

Y_lv <- as.matrix(trait_subset_lv[, c("ving_h_std", "nebb_l_std")])

# Bundle data
dataset_lv <- list(
  No = No_lv,
  Nt = Nt_lv,
  Y = Y_lv,
  animal = animal_lv,
  Na = Na_lv,
  A = as.matrix(A_hestmannoy)
)

# Set Stan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# MCMC settings
nc_lv <- 4
nw_lv <- 3000
ni_lv <- 6000
nt_lv <- 1

# Parameters to monitor
tomonitor_lv <- c(
  "mu",                  # trait means (standardised scale)
  "lambda",              # loadings 
  "sd_psi_a",            # LV genetic SD
  "sd_psi_e",            # LV non-genetic individual SD
  "sd_R",                # trait-specific residual SDs
  "var_psi_a",           
  "var_psi_e",
  "h2_psi",
  "h2_traits",           # observed trait heritabilities
  "h2_traits_no_residual"   # individual-level heritabilities (without residual error)
)

# Call Stan
out_lv <- stan(file = 'latent_variable_stan.stan',
               data = dataset_lv,
               pars = tomonitor_lv, # Can be removed to keep all parameters
               chains = nc_lv, iter = ni_lv, warmup = nw_lv, thin = nt_lv,
               open_progress = FALSE,
               refresh = 10,
               seed = 123)


# Save results
saveRDS(out_lv, 'Output_LV_Animal_Model.rds')


# Write summary to CSV
summ_lv <- summary(out_lv)$summary
write.csv(as.data.frame(summ_lv), file = "Output_LV_Summary.csv", row.names = TRUE)
cat("Summary written to: Output_LV_Summary.csv\n")




# =============================================================================
# BACK-TRANSFORMATION 

fit <- out_lv
samples <- rstan::extract(fit)

flip <- samples$lambda[,1] < 0
samples$lambda[flip, ] <- -samples$lambda[flip, ]

# Original trait means and SDs (used when standardizing)
mean_wing <- mean(trait_subset_lv$ving_h, na.rm = TRUE)
sd_wing   <- sd(trait_subset_lv$ving_h, na.rm = TRUE)
mean_beak <- mean(trait_subset_lv$nebb_l, na.rm = TRUE)
sd_beak   <- sd(trait_subset_lv$nebb_l, na.rm = TRUE)

# 1. Trait means
mu_wing  <- samples$mu[,1] * sd_wing + mean_wing
mu_beak  <- samples$mu[,2] * sd_beak + mean_beak

# 2. Loadings (scale-free)
loading_wing <- samples$lambda[,1]
loading_beak <- samples$lambda[,2]

# 3. Variance components on original scale
VA_wing   <- samples$lambda[,1]^2 * samples$var_psi_a * sd_wing^2   # additive-genetic
VA_beak   <- samples$lambda[,2]^2 * samples$var_psi_a * sd_beak^2
VInd_wing <- samples$lambda[,1]^2 * samples$var_psi_e * sd_wing^2   # individual-level non-genetic
VInd_beak <- samples$lambda[,2]^2 * samples$var_psi_e * sd_beak^2
VR_wing   <- samples$sd_R[,1]^2 * sd_wing^2                        # residual
VR_beak   <- samples$sd_R[,2]^2 * sd_beak^2

# 4. Final summary table
posterior_summary <- data.frame(
  Parameter = c(
    "mu_wing", "mu_beak",
    "loading_wing", "loading_beak",
    "VA_wing", "VA_beak",
    "VInd_wing", "VInd_beak",
    "VR_wing", "VR_beak",
    "h2_wing", "h2_beak",
    "h2_wing_indiv", "h2_beak_indiv",
    "h2_psi"
  ),
  Mean = c(
    mean(mu_wing), mean(mu_beak),
    mean(loading_wing), mean(loading_beak),
    mean(VA_wing), mean(VA_beak),
    mean(VInd_wing), mean(VInd_beak),
    mean(VR_wing), mean(VR_beak),
    mean(samples$h2_traits[,1]), mean(samples$h2_traits[,2]),
    mean(samples$h2_traits_no_residual[,1]), mean(samples$h2_traits_no_residual[,2]),
    mean(samples$h2_psi)
  ),
  SD = c(
    sd(mu_wing), sd(mu_beak),
    sd(loading_wing), sd(loading_beak),
    sd(VA_wing), sd(VA_beak),
    sd(VInd_wing), sd(VInd_beak),
    sd(VR_wing), sd(VR_beak),
    sd(samples$h2_traits[,1]), sd(samples$h2_traits[,2]),
    sd(samples$h2_traits_no_residual[,1]), sd(samples$h2_traits_no_residual[,2]),
    sd(samples$h2_psi)
  ),
  `X2.5%` = c(
    quantile(mu_wing, 0.025), quantile(mu_beak, 0.025),
    quantile(loading_wing, 0.025), quantile(loading_beak, 0.025),
    quantile(VA_wing, 0.025), quantile(VA_beak, 0.025),
    quantile(VInd_wing, 0.025), quantile(VInd_beak, 0.025),
    quantile(VR_wing, 0.025), quantile(VR_beak, 0.025),
    quantile(samples$h2_traits[,1], 0.025), quantile(samples$h2_traits[,2], 0.025),
    quantile(samples$h2_traits_no_residual[,1], 0.025), quantile(samples$h2_traits_no_residual[,2], 0.025),
    quantile(samples$h2_psi, 0.025)
  ),
  `X97.5%` = c(
    quantile(mu_wing, 0.975), quantile(mu_beak, 0.975),
    quantile(loading_wing, 0.975), quantile(loading_beak, 0.975),
    quantile(VA_wing, 0.975), quantile(VA_beak, 0.975),
    quantile(VInd_wing, 0.975), quantile(VInd_beak, 0.975),
    quantile(VR_wing, 0.975), quantile(VR_beak, 0.975),
    quantile(samples$h2_traits[,1], 0.975), quantile(samples$h2_traits[,2], 0.975),
    quantile(samples$h2_traits_no_residual[,1], 0.975), quantile(samples$h2_traits_no_residual[,2], 0.975),
    quantile(samples$h2_psi, 0.975)
  )
)

# Round and save
posterior_summary[,2:5] <- lapply(posterior_summary[,2:5], round, digits = 4)
write.csv(posterior_summary, "Output_LV_Final_Results.csv", row.names = FALSE)

cat("\n=== LATENT VARIABLE MODEL RESULTS (original scale) ===\n")
print(posterior_summary, row.names = FALSE)
cat("\nResults saved to: Output_LV_Final_Results.csv\n")